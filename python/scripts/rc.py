seq = "CGGTAGTCGAGCTGCGGATATAATATGCATATAGATCGCACGCTAGCTCATAAAAGCATGCATGCGGCTAGCTGCTGATCGTGTCG"

compdict = {"A":"T", "C":"G", "G":"C", "T":"A"}

revlist=[]
seqlist = list(seq)
seqlist.reverse()
for letter in seqlist:
    revlist.append(compdict[letter])

print(''.join(revlist))


# seq.translate(seq.maketrans("ACGT","TGCA"))[::-1]