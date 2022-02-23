from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse


def find_adapter(seq,adap,minmatch):
    testfull = seq.find(adap)
    if (testfull != -1):
        return(testfull)

    lenadap = len(adap)
    for i in range(lenadap):
        if lenadap-i >= minmatch and seq.endswith(adap[0:lenadap-i]):
            return(len(seq)-(lenadap-i))

    return(len(seq))


# adap = "CTGTCTCTTATACACATCTCCGAGCCCACGAGACAACATCGCGCATCTCGTATGCCGT"
# #seq = "ACTACAAGGACGACGATGATAAGAAGCTTCTGTCTCTTATACACATCTCCGAGCCCACGAGACAACATCGCGCATCTCGTATGCCGTCTTCTGCTTGAATAAATCGGAA"
# # seq = "ACTACAAGGACGACGATGATAAGAAGCTTCTGTCTCTTATACACATCTCCGAGCCCACGAGACAA"
# # seq = "ACTACAAGGACGACGATGATAAGAAGCTTCTGTCTC"
# seq = "CATAGCTACTTACGTAGCATCGTAGCGATCGTACTATCGTCAGTAGCTCGCGGCGCGCG"

# pos = find_adapter(seq,adap,10)
# print(pos)
# print(seq[0:pos])


parser = argparse.ArgumentParser()
parser.add_argument("--infile", help="Input file name", required=True)
parser.add_argument("--adapter", help="Adapter file name", required=True)
parser.add_argument("--outfile", help="Output file name", required=True)
parser.add_argument("--minmatch", help="Minimum length match for adapter", required=True, type=int, default=20)
parser.add_argument("--length_threshold", help="Length threshold for keeping read", required=True, default=1, type=int)
args = parser.parse_args()


for record in SeqIO.parse(args.adapter, "fasta"):
    adap = record.seq

fp = open(args.outfile,"w")
for record in SeqIO.parse(args.infile, "fastq"):
    pos = find_adapter(record.seq, adap, args.minmatch)

    if pos >= args.length_threshold:
        la = record.letter_annotations.copy()
        la['phred_quality'] = la['phred_quality'][0:pos]
        sr = SeqRecord(Seq(record.seq[0:pos]), id=record.id, letter_annotations=la)
        SeqIO.write(sr, fp, "fastq")
