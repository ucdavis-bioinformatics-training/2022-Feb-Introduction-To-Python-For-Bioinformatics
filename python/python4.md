# Putting it all together

## Exercise 1 - Counting bases

Now that you've learned the basics of programming, it's time to write some beginning code. Let's write a program that will print the count for each nucleotide in a DNA sequence. You will need to use the "list" function which, when applied to a string, will return the string as a list of characters. You will also need to use loops and a dictionary. Create the program in a new file called "counts.py". You can use this sequence to begin with:

**CGGTAGTCGAGCTGCGGATATAATATGCATATAGATCGCACGCTAGCTCATAAAAGCATGCATGCGGCTAGCTGCTGATCGTGTCG**

Hints: 
1. First put the sequence in a variable.
2. Use the "list" function to break apart the sequence into a list.
3. Iterate over the list using a loop.
4. Add the counts of each letter to a dictionary where the keys are the letters and the counts are the values.
5. Print out the dictionary using a loop.


## Exercise 2 - Get the nucleotide count from a fasta file.

Now we will take our previous code and add to it. Create a new file called "count_file.py". Download [this fasta file](https://github.com/ucdavis-bioinformatics-training/2022-Feb-Introduction-To-Python-For-Bioinformatics/raw/master/python/data/seq.fa). Take a look at it... notice it is a multiline fasta file. We are going to do the same count as previously, but we are now doing it for an entire file. You will need to skip the header lines and only count the sequence lines and instead of using the sequence, you will open the file and read sequences from there. You will also reuse the code you wrote in the previous example.

Hints:
1. Open the file and read it in line by line in a loop.
2. Check if the line begins with a ">". You will probably need to look at a [list of python string methods](https://www.w3schools.com/python/python_ref_string.asp) to find some that you need.
3. If the line does not begin with a ">", then you know it is a sequence line. Insert the count code from the previous example, but you'll have to tweak it to make it work here. You will also need to make sure to remove the newline character ("\n") at the end of each line using a string method, before counting.
4. Output the count.
5. As a bonus piece, keep count of the number of sequences and output the average bases per sequence, for each base.


## Exercise 3 - Find the reverse complement

Take the sequence from the first exercise and write code to find and print the reverse complement. Call it "rc.py".

Hints:
1. You will need to create a dictionary with the 4 nucleotides as keys, and their complements as values.
2. You will use the "list" function to break apart the sequence into a list and reverse the list using [built-in list methods](https://www.w3schools.com/python/python_ref_list.asp).
3. Loop through the bases and create a new list where each element of the new list is the complement of each from the old list. Use your dictionary to do this.
4. Print the new list using the "join" string method.
5. For a bonus, see if you can do all of steps 2-4 in just one line. You will need to look at [string slicing in python](https://www.geeksforgeeks.org/string-slicing-in-python/), and you will need to use some special string methods.


## Exercise 4 - Reverse Complement of each fasta sequence, output to a file

Create a new file called "rc_fasta.py". Use the reverse complementing code from the previous exercise as a function in this one. For the "seq.fa" file, output the reverse complement of each sequence to another file as a fasta file.

Hints:
1. First create a function that takes a sequence (string) as input and returns the reverse complement.
2. Open your output file for writing. Call it "out.fa".
3. Open the fasta file and read each line.
4. For each header line, append "\_rc" to the end.
5. For the sequence lines, you will have to gather them together until you have all of them for one fasta entry.
6. Output the reverse complemented fasta entry (using your function) to the output file. You will need to figure out the appropriate time to do this in the loop. Note that the "write" function does not output a newline, you will have to add a newline.
7. Do this for all the fasta records.
8. Close your files.


## Exercise 5 - Quality-based fastq trimmer

Okay, now let's try something more complex... a tool for trimming fastq sequences using a quality score cutoff. We will only do this for single-end Illumina reads just to make it easier. First, you need to understand the fastq format and how the quality scores are encoded. So let's take a look at some [bioinformatics file types](filetypes).

Let's [download a small fastq file](data/samp1.fastq) to use. Take a look at it. Open a new program file called "fastq_trimmer.py". What we want to do is trim the sequences (and the quality values) so that they are trimmed at the place where a base's quality value (reading from left to right) drops below a given threshold. Let's use a threshold of 15 to begin with. Then write the trimmed data to a new file.

Hints:
1. First write a function that takes two parameters, a quality string and a quality threshold, and returns the position where the sequence and quality lines will be trimmed. You will need to use the "ord" built-in function, which returns the decimal ASCII value of a character. You will need to use a loop to go through each quality value and check if it drops below the threshold. If it does, you will return the position for that value. If it doesn't then you should return the last position.
2. Open the input and output files. Loop through the input file, using a "while (True)" loop, reading 4 lines at a time into 4 separate variables (header, seq, header2, qual). Make sure to strip the newline characters for each line.
2. Right after reading in "header", check to see if it is False. If it is, this means that there is no more input, so put code to exit the loop.
3. Use your function to find the place to cut. 
4. Cut both the sequence and quality lines using the trimming position and [string slicing in python](https://www.geeksforgeeks.org/string-slicing-in-python/).
5. Write a new fastq record (4 lines) to the output file. Do this for all the records.
6. Close your files.


## Exercise 6 - Fixing logical errors and adding command-line options in the trimmer

If we take a look at the output from our previous exercise, we will see that there are some obvious problems we need to fix, and maybe some enhancements we could make. And it would be useful to be able to run our code on different files with different thresholds without having to change the code itself. The way to do that is to use **command-line arguments**. The command-line is a way of interacting with the computer through typing commands rather than using a mouse. You will have noticed that when you press the play button in VSCode, it actually runs the command in a terminal. We can add options to this command which will be sent to our code every time. There are multiple ways to do this, but we will be using a [module called "argparse"](https://docs.python.org/3/howto/argparse.html). A module is a library of functions (and other things) that you can import into your code to use. Modules should be imported at the very top of your code. argparse has many capabilities, but we will only be using one. The following code will set up command line arguments for you and put them into a variable called "args". To access the values, you simply use **args.infile**, **args.outfile**, and **args.qual_threshold**.

```
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--infile", help="Input file name", required=True)
parser.add_argument("--outfile", help="Output file name", required=True)
parser.add_argument("--qual_threshold", help="Quality threshold for trimming", required=True, type=int, default=20)
args = parser.parse_args()
```

Now, we also want to add a parameter that is a threshold length for keeping a read after trimming, which should be a default of 1. Call it "length_threshold". Then use that parameter to discard any reads (records) whose lengths are less than the threshold. Create a new file called "fastq_trimmer2.py" and copy over your previous code into this new program file and then add to this code. To test your program, you will need to run something like the following line in a Terminal in VSCode:

python3 fastq_trimmer2.py --infile samp1.fastq --outfile samp1.trimmed.fastq --qual_threshold 20 --length_threshold 5

Hints:
1. Add the import to the top of your code and add the parser lines before you start the main code, but after your function.
2. Add another argument to the parser for the length_threshold option.
2. Change the code so that it uses the **args.infile**, **args.outfile**, and **args.qual_threshold** variables instead of the input file, output file, and the quality threshold.
3. After trimming a read, check the length against the length threshold. Only write out records (4 lines) that are greater than or equal to the threshold.
4. Test the code in the terminal with various arguments and see the differences in the output.