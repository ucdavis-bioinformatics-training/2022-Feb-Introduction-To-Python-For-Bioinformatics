import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def trim_adapter(seq, adap, misthresh, minmatch):
    for i in range(len(seq)):
        mismatch = 0
        match = 0
        for j in range(len(adap)):

            #print(i,j)

            if i+j >= len(seq):
                break

            if seq[i+j] != adap[j]:
                mismatch += 1
            else:
                match += 1
            
            if mismatch > misthresh:
                break

        if (mismatch <= misthresh) and (match + mismatch == len(adap)):
            return(i)

        if match >= minmatch and mismatch <= misthresh and i+j == len(seq):
            return(i)

    return(len(seq))


#seq = "ACTACAAGGACGACGATGATAAGAAGCTTCTGTCTCTTATACACATCTCCGAGCCCACGAGACAACATCGCGCATCTCGTATGCCGTCTTCTGCTTGAATAAATCGGAA"
#seq = "ACTACAAGGACGACGATGATAAGAAGCTTCTGTCTCTTATACACATCTCCGAGCCCACGAGACAA"
seq = "ACTACAAGGACGACGATGATAAGAAGCTTCTGTCTC"

adap = "CTGTCTCTTATGCACATCTCCGAGCCCACGAGATAACATCGCGCATCTCGTATGCCGT"

pos = trim_adapter(seq,adap,5,10)
print(pos)
print(seq[0:pos])