#!/usr/bin/env python3
#SBATCH --partition=long       ### Partition (like a queue in PBS)
#SBATCH --job-name=RRPS7          ### Job Name
#SBATCH --time=1-20:01:00       ### Wall clock 0ime limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=28     ### Number of tasks to be launched per Node
#SBATCH --mail-user=rarichardson92@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

# Don't forget to load modules in bash, easybuild, prl, python/3.6.0 before running code!

import argparse
import re

def getarguments():
    parser=argparse.ArgumentParser(description = "Removes PCR duplicates from single-end SAM files (UMI/randomer, chrom/linkage group, base position). Requires sorted SAM files, optional file with line seperated UMIs. Outputs file with first occuring read by UMI and start position. Designed to be piped in from bash program for appropriate sorting/size determination.")
    parser.add_argument("-f", "--file", help = "Defines name and path of sam file to use in program. Required, must be a string. Designed as input from dedup, which sorts and saves file as file.sorted.sam. Rename file in this format if running without dedup.script.", required = True, type = str)
    parser.add_argument("-u", "--umi", help = "Defines name and path of UMI file to use in program. Optional, must be a string.", required = False, type = str)
    parser.add_argument("-sizef", help = "Defines size of forward window to use in program. Required, must be an integer. Set at 100 in dedup.script (Accounts for all soft clipping possible when sorted per UMI).", required = True, type = int)
    parser.add_argument("-sizer", help = "Defines of reverse window to use in program. Required, must be an integer. Determined by file contents in dedupt.script when used together.", required = True, type = int)
    parser.add_argument("-keep", help = "Boolean. Determines if duplicates are kept in an out file. Default is False.", required = False, default = False, type = bool)
    parser.add_argument("-p", "--paired", help = "Paired end file boolean. Currently not supported by this program.", required = False, default = False, type = bool)
    return parser.parse_args()


def fwdclip(CIGAR, start):
    '''For forward strands, adds soft clipping to start position'''
    if "S" in CIGAR:
        softclip = CIGAR.split("S")
        if str.isdigit(softclip[0]) == True:
            start = start - int(softclip[0])
    return(start)

def rvsclip(CIGAR, start):
    '''For reverse strands, adds all relevant cigar notation to start point'''
    Addlist = re.findall(r'[0-9]*[MNPD]', CIGAR)
    Cliplist = re.findall(r'[0-9]*S$', CIGAR)
    Totaladd=0
    for found in Addlist:
        Totaladd += int(found[:-1])
    if len(Cliplist) == 1:
        Totaladd += int(Cliplist[0][:-1])
    start = start + Totaladd - 1
    # -1, since original start position takes up a slot
    return(start)

def windowcheck(window, newline):
    '''Checks if line to add to window is a duplicate of anything in the window'''
    dup = False
    if newline in window:
            dup = True
            if newline[0] in dupcount:
                dupcount[newline[0]]+=1
            else:
                dupcount[newline[0]]=1
    return(dup)

def UMIcheck(UMIs, adderUMI):
    '''Determines if UMI of SAM file is in UMI dictionary is applicable'''
    continueflag = False
    if adderUMI in UMIs:
        continueflag = True
    return(continueflag)

args=getarguments()

file=str(args.file)
UMIfile=str(args.umi)
sizef=int(args.sizef)
sizer=int(args.sizer)
keep=bool(args.keep)
p=bool(args.paired)

if p == True:
    raise Exception('Paired-end input not supported.')

UMIs = set()
dupcount={}
goodcount = 0

if UMIfile != "None":
    with open(UMIfile, "r") as UMIf:
        for line in UMIf:
            UMIs.add(line.strip('\n'))

if keep == True:
    bad = open(file[:-4]+"_duplicate.sam", "w")

with open(file+".process.sam", "r") as fh, open(file[:-4]+"_deduped.sam", "w") as good:
    fwdwin = []
    rvswin = []
    #Sets up window listss
    breakflag = False
    #flag for loop break if UMIs are not good at EOF

    for line in fh:
        full=line.strip('\n')
        #grabs sam record with stripped newline
        adder=full.split('\t')
        #splits full into adder, a list from full by SAM fields

        if UMIfile != "None":
            checkflag = UMIcheck(UMIs, adder[0][-8:])
            #check if new umi is in dictionary
            while checkflag == False:
            #will reset above variables if new UMI isn't valid,
                try:
                    full=next(fh).strip('\n')
                except StopIteration:
                    breakflag = True
                    break
                    #will break loop at EOF, flag breaks out of second loop
                adder=full.split('\t')
                checkflag = UMIcheck(UMIs, adder[0][-8:])

        if breakflag == True:
            break
        #if flagged true, EOF. Discontinue loop, current line is invalid UMI.

        start=int(adder[3])
        #sets start value to where the pos field is located

        if int(adder[1]) & 16 != 16:
        #Adds soft clipping for forward strands
            start=fwdclip(adder[5], start)
            adder=[ adder[0][-8:], adder[2], start]

            duplicate = windowcheck(fwdwin, adder)
            if duplicate == False and len(fwdwin) < sizef: #replace 3 with desired window size
                good.write(full+"\n")
                fwdwin = fwdwin + [adder]
                goodcount+=1
            elif duplicate == False:
                good.write(full+"\n")
                fwdwin = fwdwin[1:] + [adder]
                goodcount+=1
            elif keep == True:
                bad.write(full+"\n")

        else:
            start = rvsclip(adder[5], start)
            adder=[ adder[0][-8:], adder[2], start ]
            duplicate = windowcheck(rvswin, adder)
            if duplicate == False and len(rvswin) < sizer: #replace 3 with desired window size
                good.write(full+"\n")
                rvswin = rvswin + [adder]
                goodcount+=1
            elif duplicate == False:
                good.write(full+"\n")
                rvswin = rvswin[1:] + [adder]
                goodcount+=1
            elif keep == True:
                bad.write(full+"\n")

if keep == True:
    bad.close()

print()
print("Number of duplicates per UMI/randomer:")
for key, value in dupcount.items():
    print(key, value, sep = '\t')
print()
print("Total duplicates: "+str(sum(dupcount.values())))
print("Remaining reads: "+str(goodcount))
print()
print("Program time elapsed:")
