counts={}

fp = open("seq.fa","r")
for line in fp:
    if (not line.startswith(">")):
        line = line.rstrip("\n")

        for letter in list(line):
            if (letter not in counts):
                counts[letter] = 1
            else:
                counts[letter] += 1

for letter in counts.keys():
    print(letter+" "+str(counts[letter]))