"""
UniqSeq
Author: Jordan A. Berg, University of Utah
'A primer designer tool for identifying unique regions 
in genomes of a given species for PCR amplification of field samples.'

"""

import random

#initiate
window = 40 #get window from user as args input
anneal = 65 #get annealing temp from user in args input
k = 10 #get number of primer pairs to output as args input
#sequences to compare
str1 = "abbaaababsa???aabaaa$$$gcccgcgcggcggagaggcagcggcgcgtttcggccggcgcggctatatazzbaaabbaaaazzbaaabbbaababaaaaabbbababaaaabbbabababbbababbabbaaabababababababbabba"
str2 = "bbbbbbaa$$$gcccgcgcggcggagaggcagcggcgcgtttcggccggcgcggctatatabaaaazzbaaabbaaaazzbaaabbbaababaaabbabaaaabasndkanjsdnkjansdnakjsndnakjsnda"
str3 = "bbbbbbaaabaaaazzbaaabbaaaazzbaaabbbaababaaabbabaaaababbabababaaaaaskdnjans$$$gcccgcgcggcggagaggcagcggcgcgtttcggccggcgcggctatatkdajsda"
str4 = "aabcbcbbccbaabaaaaabbcbbcbcbaabababcbcbcbbcbcbcbbabbbaababaaabbabaaaabbbaababaaabbabaaaaabbabababababcbcbcbcbcbcbcbcbcbbabasdfbnksndfkjnsdnfjsd"
str5 = "bcdbcbcbcbdbdbdbaababddbdbcbcbbbsbbsnsnsjsnaaazzbaaabbbaababaaafksndkfnjsndkjfnjsdnfksdjfnksdfsbbaababaaabbabaaaababaaabbbaababaaabbaba"

seq1 = [str1,str2,str3] #same species, populate from user input
seq2 = [str4,str5] #dissimilar, populate from user input
initseqs = [] #starter library of sequences
uniqseqs = [] #library of unique sequences

#populate library of unique features of given window size in species list
m = seq1[0]
end = len(m) + 1
y = window
x = window - window
while y != end:    
    initseqs.append(m[x:y])
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
                uniqseqs.append(n[x:y])
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
    
    loc_temp = 0
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
random.sample(primers.items(), k)
