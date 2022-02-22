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

Now, you should have enough knowledge to be able to write a simple adapter trimmer for fastq files. The adapter will be specified in a fasta file (with the one sequence) and the fastq file will be single-end (to keep things simple). The algorithm is to match all or part of the adapter sequence to the end of a read. If it matches, remove the adapter part of the sequence. Then write the new record to the output file. You should use argparse for parsing options. The trimming function should take in a read sequence and adapter sequence. You will loop through the read starting from the first base and try to match the adapter string to the sequence string, however, you will also allow a certain number of mismatches. E.g., the following would constitute a match if the mismatch threshold was 2 or more:

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

Hints:
1. Put your imports at the top of your code file.
1. Write a trimming function that takes in a sequence, an adapter, a mismatch threshold, and a minimum matching threshold. You will need to loop through each position of the sequence, and for each position you check the adapter against the sequence starting from that position. You will need to keep track of the matches and mismatches. If the mismatches exceed the threshold, then you do not trim. If the end of the sequence is reached and the number of matches does not exceed the minimum, then you do not trim. Otherwise, you trim from the adapter position to the end of the read. The function returns the position where to trim.
2. Test out your function with some example cases to make sure it works.
2. In the main part of your code, use argparse to create options for the input fastq file, the input fasta adapter file, the output file name, the mismatch threshold, the minimum matching threshold, and the minimum length threshold after trimming.
3. Open the adapter file and, using Biopython, read the adapter sequence into an object.
4. Open the input file and output file.
5. Using Biopython, read the fastq file in one record at at time. Use your function to get the trimming position. Trim both the sequence and the qualities using that position.
6. Create a new SeqRecord and write it to the output file if the trimmed sequence length is greater than or equal to the minimum length threshold.
7. Do this for all the records.
8. Close your files.


## bamnostic and pysam

For some reason, Biopython does not do sequence alignment files (different from multiple alignment files), so typically people use a package called "pysam" for those files. However, pysam (for some reason) does not seem to install on Windows, so we are going to use a package called "bamnostic" instead, which is very similar to pysam, but will run on Windows. 