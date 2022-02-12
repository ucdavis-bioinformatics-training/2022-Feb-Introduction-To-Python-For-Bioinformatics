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

Now we will take our previous code and add to it. Create a new file called "count_file.py". Download [this fasta file](https://github.com/ucdavis-bioinformatics-training/2022-Feb-Introduction-To-Python-For-Bioinformatics/raw/master/python/seq.fa). Take a look at it... notice it is a multiline fasta file. We are going to do the same count as previously, but we are now doing it for an entire file. You will need to skip the header lines and only count the sequence lines and instead of using the sequence, you will open the file and read sequences from there. You will also reuse the code you wrote in the previous example.

Hints:
1. Open the file and read it in line by line in a loop
2. Check if the line begins with a ">". You will probably need to look at a [list of python string methods](https://www.w3schools.com/python/python_ref_string.asp) to find some that you need.
3. If the line does not begin with a ">", then you know it is a sequence line. Insert the count code from the previous example, but you'll have to tweak it to make it work here. You will also need to make sure to remove the newline character ("\n") at the end of each line using a string method, before counting.
4. Output the count.
5. As a bonus piece, keep count of the number of sequences and output the average bases per sequence, for each base.