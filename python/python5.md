# Biopython and bamnostic/pysam

Now we will explore some packages for interfacing with common bioinformatics file types, such as fasta, fastq, and bam. 

## Biopython

[Biopython](https://biopython.org/) is a set of python libraries for software that provide a robust interface to various file types used in bioinformatics. If you haven't installed Biopython, then do so now by opening a Terminal in VSCode and typing this:

	pip install biopython

Now, let's take a look at some of the features. Biopython is based on classes and objects. Simply put, an **Object** is data bundled with functions that works on that data and a **Class** is a blueprint for an object. The first class we'll look at is the "Seq" class, which is a class that holds sequence data.

```
# first we need to import the class from Biopython
from Bio.Seq import Seq

# create a sequence object
dna = Seq("AGTCGGACTGACTACTGATCGTACGCTATTATT")

# you can reverse complement a sequence
print(dna.reverse_complement())

# you can transcribe a sequence
rna = seq.transcribe()
print(rna)

# you can translate RNA
protein = rna.translate()
print(protein)
```

Another commonly used class in Biopython is the "SeqIO" class. This class is the standard input/output interface to Biopython. SeqIO can interface with many common Bioinformatics formats, but will only use sequences as [SeqRecord](https://biopython.org/docs/1.75/api/Bio.SeqRecord.html) objects.

```
# first import the class
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# The "parse" method is the main way to extract fasta/fastq and other types of records.
# This method returns an "iterator" that you can use in a loop.
# The iterator returns SeqRecord objects.
for record in SeqIO.parse("samp1.fastq", "fastq"):
    print(record.id)
    print(record.letter_annotations)

# you can do the same with a fasta file
for record in SeqIO.parse("seq.fa", "fasta"):
	print(record.id)
	print(record.seq)


# You can also write records (files) using the "write" method
# First create a list of SeqRecords, then write the whole list at once
seqlist = []
seqlist.append(SeqRecord(Seq("ACTGATCGAGCTAGCTCATACGCTGACTGACTGATCGTCAGATGTATATATGCTATGCTGTAGCTCGATCGTCA"), id="contig1", description=""))
seqlist.append(SeqRecord(Seq("CGTACGTACGATCGATGCTAGCATAGATACGGCGCGCGCGGCGCAGATCGATGACT"), id="contig2", description=""))
fp = open("ref.fasta","w")
SeqIO.write(seqlist, fp, "fasta")
fp.close()
```


## Adapter Trimmer

Now, you should have enough knowledge to be able to write a simple adapter trimmer for fastq files. The adapter will be specified in a fasta file (with the one sequence) and the fastq file will be single-end (to keep things simple). The algorithm is to match all or part of the adapter sequence to the end of a read. If it matches, remove the adapter part of the sequence. Then write the new record to the output file. You should use argparse for parsing options.

Hints:
1. The trimming function should take in a read sequence and adapter sequence. You will loop through the read starting from the first base and try to match the adapter string to the sequence string, however, you will also allow a certain number of mismatches. E.g., the following would constitute a match if the mismatch threshold was 2 or more:

<pre>
Adapter:                               CTGTCTCTTAT<span style="color: red">G</span>CACATCTCCGAGCCCACGAGA<span style="color: red">T</span>AACATCGCGCATCTCGTATGCCGT
Sequence: ACTACAAGGACGACGATGATAAGAAGCTTCTGTCTCTTATACACATCTCCGAGCCCACGAGACAACATCGCGCATCTCGTATGCCGTCTTCTGCTTGAATAAATCGGAA
</pre>

If you find a full match (within the mismatch threshold), then you cut the portion of the read from the adapter to the end. I.e., for the above example the output sequence would be:

<pre>
Output Sequence: ACTACAAGGACGACGATGATAAGAAGCTT
</pre>

Also, for a partial match (within the mismatch threshold) that goes to the end of the sequence and is greater than some threshold for matches (say 10), you would trim at the adapter:

<pre>
Adapter:                               <span style="color: green">CTGTCTCTTA</span>TACACATCTCCGAGCCCACGAGACAACATCGCGCATCTCGTATGCCGT
Sequence: ACTACAAGGACGACGATGATAAGAAGCTTCTGTCTCTTATACACATCTCCGAGCCCACGAGACAA
Output:   ACTACAAGGACGACGATGATAAGAAGCTT
</pre>

The following scenario **would not** be trimmed:

<pre>
Adapter:                               <span style="color: green">CTGTCTC</span>TTATACACATCTCCGAGCCCACGAGACAACATCGCGCATCTCGTATGCCGT
Sequence: ACTACAAGGACGACGATGATAAGAAGCTTCTGTCTC
Output:   ACTACAAGGACGACGATGATAAGAAGCTTCTGTCTC
</pre>