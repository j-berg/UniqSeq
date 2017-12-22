"""
UniqSeq
Author: Jordan A. Berg, University of Utah
'A primer designer tool for identifying unique regions 
in genomes of a given species for PCR amplification of field samples.'

"""

#initiate
window = 20 #get window from user as args input
anneal = 65 #get annealing temp from user in args input
#sequences to compare
str1 = "abbaaababsa???aabaaaazzbaaabbaaaazzbaaabbbaababaaaaabbbababaaaabbbabababbbababbabbaaabababababababbabba"
str2 = "bbbbbbaaabaaaazzbaaabbaaaazzbaaabbbaababaaabbabaaaabasndkanjsdnkjansdnakjsndnakjsnda"
str3 = "bbbbbbaaabaaaazzbaaabbaaaazzbaaabbbaababaaabbabaaaababbabababaaaaaskdnjanskdajsda"
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
        
 #make primer dictionary
