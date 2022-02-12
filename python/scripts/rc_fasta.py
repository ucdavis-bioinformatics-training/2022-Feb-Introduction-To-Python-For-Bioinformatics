def rc(seq):
    compdict = {"A":"T", "C":"G", "G":"C", "T":"A"}
    revlist=[]
    seqlist = list(seq)
    seqlist.reverse()
    for letter in seqlist:
        revlist.append(compdict[letter])

    return(''.join(revlist))


currseq=""
outf = open("out.fa","w")
seqf = open("seq.fa","r")
for line in seqf:
    line = line.rstrip("\n")
    if (line.startswith(">")):
        if (currseq != ""):
            outf.write(rc(currseq)+"\n")
            currseq=""
        outf.write(line+"_rc\n")
    else:
        currseq += line
outf.write(rc(currseq)+"\n")

seqf.close()
outf.close()