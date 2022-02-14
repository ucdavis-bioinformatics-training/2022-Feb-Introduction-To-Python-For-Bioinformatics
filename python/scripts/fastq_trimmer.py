def trim(qual,threshold):
    pos=0
    for base in list(qual):
        if (ord(base)-33 < threshold):
            return(pos)
        else:
            pos += 1


out = open("trimmed.fastq","w")
fq = open("samp1.fastq","r")

while (True):

    header = fq.readline().rstrip("\n")
    if not header:
        break
    seq = fq.readline().rstrip("\n")
    header2 = fq.readline().rstrip("\n")
    qual = fq.readline().rstrip("\n")

    trimpos = trim(qual,15)
    seq = seq[0:trimpos]
    qual = qual[0:trimpos]

    out.write(header+"\n"+seq+"\n"+header2+"\n"+qual+"\n")

fq.close()
out.close()