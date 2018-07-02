"""
UniqSeq
Author: Jordan A. Berg, University of Utah
'A primer designer tool for identifying unique regions
in genomes of a given species for PCR amplification of field samples.'
"""

import random
import argparse
import pandas as pd
import datetime

# Input and Output files in command line script call
parser = argparse.ArgumentParser(description="This script identifies unique regions in genomes of a given species for PCR amplification of field samples..")
parser.add_argument("-i","--input_species",help="Input directory of like species that primers will amplify (copy and paste working directory address, make sure no spaces are used, must be fasta files)")
parser.add_argument("-o","--input_other",help="Input directory of other species that are similar but not the target (copy and paste working directory address, make sure no spaces are used, must be fasta files)")
parser.add_argument("-w","--window",type=int, help="Window size for region identification (default=1200)", default=1200)
parser.add_argument("-a","--anneal_temp",type=int, help="Desired Taq PCR annealing temp (default=65)" default=65)
parser.add_argument("-r","--output_number",type=int, help="Desired number of primer pairs to output (default=10)", default=10)
args = parser.parse_args()

#initiate
window = args.window #get window from user as args input
anneal = args.anneal_temp #get annealing temp from user in args input
k = args.output_number #get number of primer pairs to output as args input

seq1 = [] #same species, populate from user input
seq2 = [] #dissimilar, populate from user input
initseqs = set() #starter library of sequences
uniqseqs = set() #library of unique sequences

#sequences to compare
for subdir, dirs, files in os.walk(args.input_species): #walk through raw data files within a given directory (first argument passed into program)
    for x in files:
        if x.startswith('.'): #ignore hidden files
            pass
        else:
            stri = ""
            openfilex = open(os.path.join(args.input_species, x)).readlines()
            for line in openfilex:
                if line.startswith('>'):
                    pass
                else:
                    line = line.rstrip("\n")
                    line = line.rstrip("\t")
                    line = line.rstrip(" ")
                    stri = stri + line
            seq1.append(stri)

for subdir, dirs, files in os.walk(args.input_other): #walk through raw data files within a given directory (first argument passed into program)
    for x in files:
        if x.startswith('.'): #ignore hidden files
            pass
        else:
            stri = ""
            openfilex = open(os.path.join(args.input_other, x)).readlines()
            for line in openfilex:
                if line.startswith('>'):
                    pass
                else:
                    line = line.rstrip("\n")
                    line = line.rstrip("\t")
                    line = line.rstrip(" ")
                    stri = stri + line
            seq2.append(stri)

#populate library of unique features of given window size in species list
m = seq1[0]
end = len(m) + 1
y = window
x = window - window
while y != end:
    initseqs.add(m[x:y])
    y = y + 1
    x = x + 1

for n in seq1[1:]:
    end = len(n) + 1
    y = window
    x = window - window
    while y != end:
        if n[x:y] not in initseqs:
            pass
        else:
            if n[x:y] in uniqseqs:
                pass
            else:
                uniqseqs.add(n[x:y])
        y = y + 1
        x = x + 1

#remove any features that are shared by dissimilar species
for p in seq2:
    end = len(p) + 1
    y = window
    x = window - window
    while y != end:
        if p[x:y] not in uniqseqs:
            pass
        else:
            uniqseqs.remove(p[x:y])
        y = y + 1
        x = x + 1

#remove any seqs with illegal chars
bad_seqs = []
for uniqs in uniqseqs:
    for char in list(uniqs):
        if char == 'a' or char == 't' or char == 'c' or char == 'g' or char == 'A' or char == 'T' or char == 'C' or char == 'G':
            continue
        else:
            bad_seqs.append(uniqs)
            break

for se in bad_seqs:
    uniqseqs.remove(se)

#for the future
#make sure that primers aren't found in dissimilar species (implement blast to makes sure they are at least 75% different than anything in the genome of something else)

#report is the primer sequence is found elsewhere in the similar species
    
#make primer dictionary
primers = {} #for and rev primers, key being start_loc, end_loc, and sequence

def rc(sequence):
    new_seq = []
    for cha in sequence:
        if cha == 'A':
            new_seq.append('T')
        elif cha == 'a':
            new_seq.append('t')
        elif cha == 'T':
            new_seq.append('A')
        elif cha == 't':
            new_seq.append('a')
        elif cha == 'G':
            new_seq.append('C')
        elif cha == 'g':
            new_seq.append('c')
        elif cha == 'C':
            new_seq.append('G')
        elif cha == 'c':
            new_seq.append('g')
        else:
            pass
    return ''.join(new_seq)

for q in uniqseqs:
    char = list(q)
    loc_temp = 0
    for_seq = []
    for ch in char:
        for_seq.append(ch)
        for_temp = loc_temp
        if ch == 'A' or ch == 'a':
            loc_temp = loc_temp + 2
            if loc_temp > anneal:
                break
            else:
                continue
        elif ch == 'T' or ch == 't':
            loc_temp = loc_temp + 2
            if loc_temp > anneal:
                break
            else:
                continue
        elif ch == 'G' or ch == 'g':
            loc_temp = loc_temp + 3
            if loc_temp > anneal:
                break
            else:
                continue
        elif ch == 'C' or ch == 'c':
            loc_temp = loc_temp + 3
            if loc_temp > anneal:
                break
            else:
                continue
        else:
            pass

    loc_temp = -5
    rev_seq = []
    for ch in reversed(char):
        rev_seq.append(ch)
        rev_temp = loc_temp
        if ch == 'A' or ch == 'a':
            loc_temp = loc_temp + 2
            if loc_temp > anneal:
                break
            else:
                continue
        elif ch == 'T' or ch == 't':
            loc_temp = loc_temp + 2
            if loc_temp > anneal:
                break
            else:
                continue
        elif ch == 'G' or ch == 'g':
            loc_temp = loc_temp + 3
            if loc_temp > anneal:
                break
            else:
                continue
        elif ch == 'C' or ch == 'c':
            loc_temp = loc_temp + 3
            if loc_temp > anneal:
                break
            else:
                continue
        else:
            pass

    #add both to dictionary with q:[forward_seq, reverse_seq, final_temp (min)]
    if rev_temp < for_temp:
        primers[q] = (''.join(for_seq), rc(rev_seq), rev_temp) #reversed but not complemented yet
    else:
        primers[q] = (''.join(for_seq), rc(rev_seq), for_temp) #reversed but not complemented yet

#output random primer sets from dictionary for trial -- include conversions for strand specificity (take reverse complement of reverse primer)
print(primers) #make sure it doesn't output the same twice, output to pandas matrix and export or just print out, either is fine

matrix = pd.DataFrame(random.sample(primers.items(), k))
now = datetime.datetime.now()
matrix.to_csv(args.input_species + '/uniqseq_primers' + now.strftime("%Y-%m-%d") + '.csv', header=False, index=True, sep=',', mode='a', na_rep='NA')
