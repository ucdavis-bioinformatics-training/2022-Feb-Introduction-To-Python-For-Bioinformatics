seq = "CGGTAGTCGAGCTGCGGATATAATATGCATATAGATCGCACGCTAGCTCATAAAAGCATGCATGCGGCTAGCTGCTGATCGTGTCG"

counts={}
for letter in list(seq):
    if (letter not in counts):
        counts[letter] = 1
    else:
        counts[letter] += 1

for letter in counts.keys():
    print(letter+" "+str(counts[letter]))