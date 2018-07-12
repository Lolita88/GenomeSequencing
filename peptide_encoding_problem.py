"""
Peptide Encoding Problem: Find substrings of a genome encoding a given amino acid sequence.
Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
Output: All substrings of Text encoding Peptide (if any such substrings exist).
"""

def peptide_encoding_problem(dna, peptide, rna_codon_table_array):
    #looking to output all possible encoded peptides
    #codon_table = rna_codon_table_array
    encoded_peptides = []

    #strings that will have to be searched for MA once encoded to peptide
    dna_5_3 = dna #coding strand, like mRNA strand
    dna_3_5 = get_reverse_complement(dna_5_3)
    #print("dna_5_3 " + str(dna_5_3))
    #print("dna_3_5 " + str(dna_3_5))
    rna_5_3 = dna_transcribed_to_rna(dna_5_3) #look in this strand for MA
    rna_3_5 = dna_transcribed_to_rna(dna_3_5)
    print("rna_5_3 " + str(rna_5_3))
    print("rna_3_5 " + str(rna_3_5))
    pattern_codons = get_codons(peptide, rna_codon_table_array)
    #print("pattern_codons " + str(pattern_codons))

    #try to find codon pattern in rna 5'to 3'
    search_strand_for_pattern(rna_5_3, peptide, pattern_codons)

def get_codons(peptide, rna_codon_table_array):
    encoded_peptides_dict = {}
    rna_dict = {}
    for each in rna_codon_table_array:
        try: #need to handle exceptions when there are "stop codons"
            #that do not translate into amino acids and serve to halt translation
            k,v = each.split(" ")
        except ValueError:
            v = ""
            #print("stop")
        if v != "":
            rna_dict[k] = v
    #print("rna_dict " + str(rna_dict)) #whole codon table array minus stops

    this_peptide = []
    for each in peptide: #peptide is peptide passed in
        this_peptide.append(each)
    #print(this_peptide) # M A

    for each in this_peptide:
        #print("peptide " + str(each))
        temp_list = []
        for key in rna_dict:
            #print(rna_dict[key])
            #print(key)
            if each == rna_dict[key]: #match!
                # add to list to add to dict key
                temp_list.append(key) # for each key, one or more codons
        #print(temp_list)
        encoded_peptides_dict[each] = temp_list
    #print("encoded_peptides_dict " + str(encoded_peptides_dict))
    return encoded_peptides_dict

def get_reverse_strand_dna(dna):
    return dna[::-1]

def get_reverse_complement(string):
    rna = False
    if "U" in string:
        rna = True
    #can take dna or rna
    reverse_string = ""
    for nucleotide in string:
        if nucleotide == "A":
            if rna == True:
                reverse_string += "U"
            else:
                reverse_string += "T"
        elif nucleotide == "T" or nucleotide == "U":
            reverse_string += "A"
        elif nucleotide == "C":
            reverse_string += "G"
        elif nucleotide == "G":
            reverse_string += "C"
        else:
            return "non nucleotide found" #should throw error
    #return reverse_string
    return reverse_string[::-1]
            
def dna_transcribed_to_rna(dna):
    rna = dna.replace("T", "U")
    return rna

def search_strand_for_pattern(strand, peptide, pattern_codons):#will pattern_codons always be in same order, since dict?
    #go through all reading frames, look for first codon, find, then look for second
    #reading frame
    #codon_list[]

    for i in range(len(strand)-2):
        temp_codon_location = 0
        #print(i)
        #print(strand[i:i+3])
        #for each in pattern_codons:
        for key, val in pattern_codons.items():
            #print(pattern_codons[each])
            #print(len(val))
            if(len(val) == 1):
                #print(str(val).strip('[]').replace("'",""))
                temp_val = str(val).strip('[]').replace("'","")
                if temp_val == strand[i:i+3]:
                    print("match")
            elif(len(val) > 1): #need to check all
                for each in val:
                    #print(each)
            #print(strand[i:i+3])
            """
            if temp_val == strand[i:i+3]:
                print("match")
                #print(temp_val)
                temp_codon_location = i
                #don't save anything yet, but search for second codon
                #maybe
                # do need to mark, so can break out of and get to next codon
                #break
            """


dna = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
peptide = "MA"
rna_codon_table_array = [line.strip() for line in open("files/rna_codon_table_array.txt")]
#print(rna_codon_table_array[0])

print(peptide_encoding_problem(dna, peptide, rna_codon_table_array))

"""
#Sample output:
ATGGCC
GGCCAT
ATGGCC
"""


