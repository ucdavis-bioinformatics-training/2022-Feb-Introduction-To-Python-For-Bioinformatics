# Part 2: Basics of Python, Variables, and Basic Data Types


A few notes before we get started:
- Learning all the nuances of python takes a long time! Our goal here is to introduce you to as many concepts as possible
but if you are serious about mastering python you will need to apply yourself beyond this introduction. 
- We bring up a lot of concepts to expose you to them but we encourage you to have a "scientific" mentality and highly
encourage you to continue testing the waters beyond these materials. 
- Python can be run interactively, however, we will be saving the code to files to be run in VSCode. We will be using the "print" function to output values to the terminal.

## Variables

What is a variable? A variable is named storage that can be used to store anything from numbers to characters to much more complex data. In other languages, one needs to specify the "type" of the variable before assigning values, however, in python that is not the case. In python, a variable's type is determined by the value it is assigned. Variable names can contain only alphanumeric characters (A-Z, a-z, 0-9) and underscore (\_) and no other special characters and cannot start with a number. In the following section, you should copy and paste (or type) the code into your "helloworld.py" file and run the code each time to see the output in the terminal.

# Basic Data Types: Integers, Floating-point numbers, booleans, strings.

### Integers
whole numbers, negative or positive

```
gene_count = 7
print(gene_count)
```

### Floats
basically Real numbers

```
gene_exp = 3.89
print(gene_exp)
```

### Booleans

"In computer science, the Boolean data type is a data type that has one of two possible values (usually denoted true 
and false) which is intended to represent the two truth values of logic and Boolean algebra. It is named after George 
Boole, who first defined an algebraic system of logic in the mid 19th century." [-wikipedia](https://en.wikipedia.org/wiki/Boolean_data_type)

```
control = False
treatment = True
print(control)
print(treatment)
```

### Strings
Any text (including all symbols and spaces) is called a "string". Strings are designated by quotes around the text.

```
gene_id = "BRCA2"
print(gene_id)
hw = "HELLO!!! WORLD!!!"
print(hw)
```