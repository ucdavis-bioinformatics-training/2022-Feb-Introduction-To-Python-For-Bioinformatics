import argparse

def trim(qual,threshold):
    pos=0
    for base in list(qual):
        if (ord(base)-33 < threshold):
            return(pos)
        else:
            pos += 1



parser = argparse.ArgumentParser()
parser.add_argument("--infile", help="Input file name", required=True)
parser.add_argument("--outfile", help="Output file name", required=True)
parser.add_argument("--qual_threshold", help="Quality threshold for trimming", required=True, type=int, default=20)
parser.add_argument("--length_threshold", help="Length threshold for keeping read", required=True, default=1, type=int)
args = parser.parse_args()

out = open(args.outfile,"w")
fq = open(args.infile,"r")

while (True):

    header = fq.readline().rstrip("\n")
    if not header:
        break
    seq = fq.readline().rstrip("\n")
    header2 = fq.readline().rstrip("\n")
    qual = fq.readline().rstrip("\n")

    trimpos = trim(qual,args.qual_threshold)
    seq = seq[0:trimpos]
    qual = qual[0:trimpos]

    if (len(seq) >= args.length_threshold):
        out.write(header+"\n"+seq+"\n"+header2+"\n"+qual+"\n")

fq.close()
out.close()