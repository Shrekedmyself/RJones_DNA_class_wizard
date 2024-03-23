#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Introduction and dictionaries


# In[1]:


# Not sure if I removed all of the extraneous stuff but bear with me on it

# Inserting code from 'helpful_variables.txt':
standard_code = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGA": "*",
    "UGU": "C",
    "UGC": "C",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


aa_mol_weights = {
    "A": 89.09,
    "C": 121.15,
    "D": 133.1,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "K": 146.19,
    "L": 131.17,
    "M": 149.21,
    "N": 132.12,
    "P": 115.13,
    "Q": 146.15,
    "R": 174.2,
    "S": 105.09,
    "T": 119.12,
    "V": 117.15,
    "W": 204.23,
    "X": 0,
    "Y": 181.19,
}


# In[ ]:


# Class seq introduction

####Don't hate me if this isn't perfect


# In[2]:


# call the class seq
class seq:
    # instance attributes: attributes that will be different for each object made for the class
    # Instance attribute
    def __init__(self, name, sequence, organism, type):
        self.name = name
        self.organism = organism
        self.sequence = sequence
        self.type = type

    # define the info function
    def info(self):
        print(self.name)
        print(self.organism)
        print(self.type)
        print(self.sequence)

    # define the length function
    def length(self):
        x = len(self.sequence)
        print(x)

    # define the fasta_out function.
    # write to a file using the sequence name as part of the file name
    # I have written most of the code for the fasta_out function,
    # but you need to add something. What is missing?
    # this url may be helpful https://www.w3schools.com/python/python_file_write.asp
    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "\n"
            + self.sequence
        )
        f.close()


# In[ ]:


# Class protein introduction

##### Let's hope this works...


# In[84]:


# Write the new protein class here
class protein(seq):
    def __init__(self, name, sequence, organism, type, size):
        self.size = size

        super().__init__(name, sequence, organism, type)

    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "_"
            + self.size
            + "\n"
            + self.sequence
        )
        f.close()

    def mol_weight(self):
        results = 0
        for aa in self.sequence:
            if aa in aa_mol_weights:
                results += aa_mol_weights[aa]
            else:
                print("Unknown amino acid, ignoring: {aa}")
        print("The molecular weight is: {results}")
        return results


# In[ ]:


# Class nucleotide introduction


# In[6]:


# Write the new nucleotide class here
class nucleotide(seq):
    def __init__(self, name, sequence, organism, type):

        super().__init__(name, sequence, organism, type)

    def gc_content(self):
        total = len(self.sequence)
        gc = self.sequence.count("G") + self.sequence.count("C")
        gc_percent = 100 * gc / total
        print(gc_percent)


# In[ ]:


# DNA class introduction


# In[110]:


# Write the DNA class here
class DNA(nucleotide):
    def __init__(self, name, sequence, organism, type):

        super().__init__(name, sequence, organism, type)

    # Add a transcribe method
    # This site may be helpful https://www.geeksforgeeks.org/python-string-replace/
    def transcribe(self):
        transcript = self.sequence.replace("T", "U")
        print(transcript)

    def six_frames(self):
        six = (
            []
        )  # Creating an empty list because I tried it with a string and it did not work
        for i in range(0, len(self.sequence), 3):
            six.append(
                self.sequence[i : i + 3]
            )  # Extract a block of three characters from the sequence and append it to our list
            rev_six = six[
                ::-1
            ]  # I mean I can't think of a better way to do this lol... Simplicity is beauty??
        print(six)
        print(rev_six)

    def reverse_complement(self):
        reverse = ""  # Yeah we doing this manually because I'm awful at python
        for base in self.sequence:
            if base == "G":
                reverse += "C"
            elif base == "C":
                reverse += "G"
            elif base == "T":
                reverse += "A"
            elif base == "A":
                reverse += "T"
        print(reverse[::-1])
        return reverse[::-1]


# In[ ]:


# RNA class introduction


# In[13]:


# Write the RNA class here
class RNA(nucleotide):
    def __init__(self, name, sequence, organism, type):

        super().__init__(name, sequence, organism, type)

    # Add a start method
    # This site may be https://www.geeksforgeeks.org/python-string-find/
    def start(self):
        print(self.sequence.find("AUG"))

    def translate(self):
        start_codon_index = self.sequence.find("AUG")  # Finding the start codon
        if (
            start_codon_index == -1
        ):  # What to do if AUG is not found (as .find will return -1)
            return "No start codon (AUG) found in the sequence."  # Using return as it will end
        codons = [
            self.sequence[i : i + 3]
            for i in range(start_codon_index, len(self.sequence), 3)
        ]  # Creating a new variable, codons, that will
        # Use the same `self.sequence[i:i+3]` indexing as the DNA class and then implement a for loop and then check the range with using our
        # Previously called start_codon_index variable with the length of the input sequence and going by every three characters
        protein_sequence = ""  # Empty variable (for now)
        for (
            codon
        ) in (
            codons
        ):  # A for loop within a for loop (within a for loop, paraphrased from Frank Herbert)
            if codon in standard_code:
                protein_sequence += standard_code[
                    codon
                ]  # Using an if statement to check if our codon is in our dictionary 'standard_code' and then adding it to our protein sequence
            else:
                protein_sequence += "?"  # Placeholder for unknown codons (shouldn't happen hopefully...)
        print(protein_sequence)
        return protein_sequence


# In[ ]:


# Workshop introduction

####Yeah this isn't gonna completely work but oh well!!!!


# In[113]:


# 8. Starting at the end of the notebook, assign the sequence: "CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA" to a variable of class DNA called uidA.
# This variable attributes should be name uidA, organism bacteria, and type DNA.
uidA = DNA(
    name="uidA",
    organism="bacteria",
    sequence="CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA",
    type="DNA",
)

# 9. Use the fasta_out function to write the sequence and information for uidA DNA to a fasta file
uidA.fasta_out()

# 10. Use the six_frames function and reverse_complement functions to output the six coding frames and reverse complement of the uidA DNA sequence
uidA.six_frames()
uidA.reverse_complement()

# 11. Using the transcribe function for the DNA class, transcribe the uidA DNA sequence to an RNA sequence.
uidA.transcribe()
transcribed = uidA.transcribe()
print(uidA.transcribe())
print(transcribed)

# 12. Save this RNA sequence as a RNA class object called uidA_RNA with the same other attributes except the name should be uidA_RNA and the type should be RNA
uidA_RNA = RNA(
    name="uidA_RNA", organism="bacteria", type="RNA", sequence="{uidA.transcribe}"
)

# 13. Use the fasta_out() function to write the RNA sequence and information for uidA to a fasta file
uidA_RNA.fasta_out()

# 14. Use the translate method on the RNA object uidA_RNA to translate the RNA sequence
translated = uidA_RNA.translate()
print(uidA_RNA.translate())
print(translated)

# 15. Save this amino acid sequence as a protein class object called uidA_protein. Set the name as uidA_protein and the type as protein.
uidA_protein = protein(
    name="uidA_protein",
    organism="bacteria",
    sequence="{uidA_RNA.translated}",
    type="protein",
    size="55 BURGERS 55 FRIES 55 TACOS 55 PIES 55 COKES 100 TATER TOTS 100 PIZZA 100 TENDERS 100 MEATBALLS 100 COFFEES 55 WINGS 55 SHAKES 55 PANCAKES 55 PASTAS 55 PASTAS AND 155 TATERS ",
)
# You can set the size attribute as any value. Use the fasta_out() function to write this protein sequence and information to a new fasta file.
uidA_protein.fasta_out()

# 16. Use the method mol_weight to output the molecular weight of the amino acid sequence uidA_protein
uidA_protein.mol_weight()
