#!/usr/bin/env python

import re
import sys

try:
    IDENTITY = float(sys.argv[1])
except:
    print("usage: calc_identity.py IDENTITY_SCORE", file=sys.stderr)
    print("where IDENTITY_SCORE is [0-100] identity percentage", file=sys.stderr)
    sys.exit(1)

# MAIN Loop
for line in sys.stdin:

    # USE MD Field if we can
    match = re.search(r'MD:Z:([A-Z0-9^]*)', line)
    if match: # MD FIELD found, use this
        MD = match.group(1)
        M = 0
        S = 0
        for m in re.findall(r'([0-9]+)', MD):
            M += int(m)
        for m in re.findall(r'([ATCG^])', MD):
            S += 1
        identity = round(100 * M/(M+S), 2)
        if identity > IDENTITY:
            #print(f"{identity} : MD:Z:{MD}")
            print(line.rstrip())
    else: # Use CIGAR if MD not found, only an approximation
        cigar = re.split(r'\s+', line)[5]
        M = 0   # alignment (match or missmatch)
        N = 0   # skipped from reference
        I = 0   # insertion to reference
        for m in re.findall(r'([0-9]+)M', cigar):
            M += int(m)
        for m in re.findall(r'([0-9]+)N', cigar):
            N += int(m)
        for m in re.findall(r'([0-9]+)I', cigar):
            I += int(m)
        identity = round(100 * M/(M+N+I), 2)
        if identity > IDENTITY:
            #print(f"{identity} : {cigar}")
            print(line.rstrip())



# Col   Field   Type    Regexp/Range                Brief description
# 1     QNAME   String  [!-?A-~]{1,254}             Query template NAME
# 2     FLAG    Int     [0, 2^16 − 1]               bitwise FLAG
# 3     RNAME   String  \*|[:rname:∧*=][:rname:]*   Reference sequence NAME
# 4     POS     Int     [0, 2^31 − 1]               1-based leftmost mapping POSition
# 5     MAPQ    Int     [0, 2^8 − 1]                MAPping Quality
# 6     CIGAR   String  \*|([0-9]+[MIDNSHPX=])+     CIGAR string
# 7     RNEXT   String  \*|=|[:rname:∧*=][:rname:]* Reference name of the mate/next read
# 8     PNEXT   Int     [0, 2^31 − 1]               Position of the mate/next read
# 9     TLEN    Int     [−2^31 + 1, 2^31 − 1]       observed Template LENgth
# 10    SEQ     String  \*|[A-Za-z=.]+              segment SEQuence
# 11    QUAL    String  [!-~]+                      ASCII of Phred-scaled base QUALity+33
# 12    MD                                          OPTIONAL 

# CIGAR FIELD
# Op    Description                                             Consumes    Consumes
#                                                               query       reference
#
# M     alignment match (can be a sequence match or mismatch)   yes         yes
# I     insertion to the reference                              yes         no
# D     deletion from the reference                             no          yes
# N     skipped region from the reference                       no          yes
# S     soft clipping (clipped sequences present in SEQ)        yes         no
# H     hard clipping (clipped sequences NOT present in SEQ)    no          no
# P     padding (silent deletion from padded reference)         no          no
# =     sequence match                                          yes         yes
# X     sequence mismatch                                       yes         yes

# CIGAR: 100S84M16S
# 
# MD: 27A5G2T6G4T1A6T2C0A7A1A0A0C0G1G1A5
# 68+16 = 84