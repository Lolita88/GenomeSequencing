"""
Peptide Encoding Problem: Find substrings of a genome encoding a given amino acid sequence.
Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
Output: All substrings of Text encoding Peptide (if any such substrings exist).
"""

def peptide_encoding_problem(dna, peptide, rna_codon_table_array):
    #looking to output all possible encoded peptides
    #codon_table = rna_codon_table_array
    decoded_peptides_found = {}

    #strings that will have to be searched for MA once encoded to peptide
    dna_5_3 = dna #coding strand, like mRNA strand
    dna_3_5 = get_reverse_complement(dna_5_3)
    #print("dna_5_3 " + str(dna_5_3))
    #print("dna_3_5 " + str(dna_3_5))
    rna_5_3 = dna_transcribed_to_rna(dna_5_3) #look in this strand for MA
    rna_3_5 = dna_transcribed_to_rna(dna_3_5)
    print("rna_5_3 " + str(rna_5_3))
    #print("rna_3_5 " + str(rna_3_5))
    pattern_codons = get_codons(peptide, rna_codon_table_array)
    #print("pattern_codons " + str(pattern_codons))

    #try to find codon pattern in rna 5'to 3'
    decoded_peptides_found = search_strand_for_pattern(rna_5_3, peptide, pattern_codons)
    #print(decoded_peptides_found)
    get_rev_complement_codon_seqs(decoded_peptides_found)


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
    temp_peptide_dict = {}
    for i in range(len(strand)-2):
        temp_codon_location = 0
        
        #print(i)
        #print(strand[i:i+3])
        #for each in pattern_codons:
        for key, val in pattern_codons.items():
            #print(pattern_codons[each])
            #print(len(val))
            if(len(val) == 1):
                #if single codon for first position, check
                temp_val = str(val).strip('[]').replace("'","")
                if(temp_val == strand[i:i+3]):
                    #if match, move down 3, check next
                    #print("match")
                    temp_peptide_dict[i] = []
                    temp_peptide_dict[i].append(temp_val)
                    #move to next
                else:
                    #no match, move to next reading frame, need to make sure this is happening
                    break
            elif(len(val) > 1): #need to check all
                for each in val:
                    #print(each)
                    #print(strand[i+3:i+6])
                    if(each == strand[i+3:i+6]):
                        #print(strand[i+3:i+6])
                        #print("match2")
                        temp_peptide_dict[i].append(each)
            else:
                print("oh hello")
    #print(temp_peptide_dict)
    return temp_peptide_dict          

def get_peptide(codons, rna_codon_table_array):
    encoded_peptides_dict = {}
    for each in rna_codon_table_array:
        try: #need to handle exceptions when there are "stop codons"
            #that do not translate into amino acids and serve to halt translation
            k,v = each.split(" ")
        except ValueError:
            v = ""
            #print("stop")
        if v != "":
            rna_dict[k] = v
    return codon_str[::-1]
    this_peptide = []
    for each in peptide: #peptide is peptide passed in
        this_peptide.append(each)
    #print(this_peptide) # M A

    for each in codons:
        temp_list = []
        for val in rna_dict:
            print(val)
            #print(rna_dict[key])
            #print(key)
            """if each == rna_dict[val]: #match!
                # add to list to add to dict key
                temp_list.append(key) # for each key, one or more codons"""
        #print(temp_list)
        encoded_peptides_dict[each] = temp_list
    #print("encoded_peptides_dict " + str(encoded_peptides_dict))
    return encoded_peptides_dict

def get_rev_complement_codon_seqs(decoded_peptides_dict, rna_codon_table_array):
    complement_codon_dict = {}
    for key, val in decoded_peptides_dict.items():
        """for each in val:
            print(each)
            get_peptide(each, rna_codon_table_array)"""
        print(val)
        temp_val = str(val).strip('[]').replace("'","").replace(",","").replace(" ","")
        print(temp_val)
    #return codon_str[::-1]

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


