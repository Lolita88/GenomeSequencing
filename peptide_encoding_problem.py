"""
Peptide Encoding Problem: Find substrings of a genome encoding a given amino acid sequence.
Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
Output: All substrings of Text encoding Peptide (if any such substrings exist).
"""
"""
TAQTREAM

T 	ACC, ACG, ACA, ACU
A	GCC, CCG, ,GCA, GCU
Q	CAA, CAG

T 	ACC, ACG, ACA, ACU
R	CGC, CGG, CGA, CGU
E	GAA, GAG
A	GCC, CCG, ,GCA, GCU
M	AUG
"""

pattern_codons = [] # needs to be a list of dicts

def peptide_encoding_problem(dna, peptide, rna_codon_table_array):
    #looking to output all possible encoded peptides
    #codon_table = rna_codon_table_array
    #decoded_peptides_found = {}

    #strings that will have to be searched for MA once encoded to peptide
    # assuming if one string, it's coming in 5' to 3'
    dna_5_3 = dna #coding strand, like mRNA strand
    dna_3_5 = get_reverse_complement(dna_5_3)
    #print("dna_5_3 " + str(dna_5_3))
    #print("dna_3_5 " + str(dna_3_5))
    rna_5_3 = dna_transcribed_to_rna(dna_5_3) #look in this strand for MA
    rna_3_5 = dna_transcribed_to_rna(dna_3_5) #then look here for MA but going to have to rev translate back
    #print("rna_5_3 " + str(rna_5_3))
    #print("rna_3_5 " + str(rna_3_5))
    #rna_5_3 AUGGCCAUGGCCCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA
    #rna_3_5 UCACCCGUUAAUACGGGUACUAUUGAUCUCAGUUCUGGGGGCCAUGGCCAU
    # fix get_codons to handle list of dicts not just dicts

    #pattern_codons = get_codons(peptide, rna_codon_table_array)
    get_codons(peptide, rna_codon_table_array)
    #attrobj = attrdict.wrap(obj)
    #pattern_codons = attrdict.wrap(get_codons(peptide, rna_codon_table_array))
    #print("pattern_codons " + str(pattern_codons))
    # {'M': ['AUG'], 'A': ['GCA', 'GCC', 'GCG', 'GCU']} # wrong now
    # [[{'T': 'ACA'}, {'T': 'ACC'}, {'T': 'ACG'}, {'T': 'ACU'}], [{'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}], [{'Q': 'CAA'}, {'Q': 'CAG'}], [{'T': 'ACA'}, {'T': 'ACC'}, {'T': 'ACG'}, {'T': 'ACU'}], [{'R': 'AGA'}, {'R': 'AGG'}, {'R': 'CGA'}, {'R': 'CGC'}, {'R': 'CGG'}, {'R': 'CGU'}], [{'E': 'GAA'}, {'E': 'GAG'}], [{'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}], [{'M': 'AUG'}]]

    #try to find codon pattern in rna 5'to 3'
    decoded_peptides_found_rna_5_3 = search_strand_for_pattern(rna_5_3, peptide, True)
    #decoded_peptides_found_rna_3_5 = search_strand_for_pattern(rna_3_5, peptide, False)
    #{0: ['AUG', 'GCC'], 6: ['AUG', 'GCC']}
    #['GGC', 'CAU']

    #list_1 = convert_found_codons_for_output(decoded_peptides_found_rna_5_3)
    #list_2 = convert_found_codons_for_output(decoded_peptides_found_rna_3_5)
    #return(join_peptides_from_lists(list_1, list_2))

    #print(convert_found_codons_for_output(decoded_peptides_found_rna_5_3))
    #print(convert_found_codons_for_output(decoded_peptides_found_rna_3_5))
    #convert_found_codons_for_output(decoded_peptides_found_rna_5_3)
    #convert_found_codons_for_output(decoded_peptides_found_rna_3_5)
def join_peptides_from_lists(a,b):
    peptide = []
    peptide.append(str('\n'.join(a)))
    peptide.append(str('\n'.join(b)))
    return(str('\n'.join(peptide)))
    
def convert_found_codons_for_output(codon_dict):
    #print(codon_dict) 
    # {0: ['AUG', 'GCC'], 6: ['AUG', 'GCC']}
    # {43: ['GGC', 'CAU']}
    output_rna = []
    for key,val in codon_dict.items():
        #print("key " + str(key)) #key 0, key 6
        string_rna = ""
    
        for each in val:
            #print("val " + str(val)) #val ['AUG', 'GCC']
            #print("each " + str(each)) # this is right!!!!, not dna yet
            # AUG, then GCC before switches over to key 6
            string_rna = string_rna + each
            #print("string_rna " + str(string_rna)) #working
            #print("current key " + str(key))
            #print("current val " + str(val))
            #print("string_rna before " + str(string_rna)) #working
            #output_rna.append("hello")
        #print("string_rna before " + str(string_rna)) #working
        #print("before append output_rna " + str(output_rna))
        output_rna.append(rna_transcribed_to_dna(string_rna))
        #print("after append output_rna " + str(output_rna)) 
    #print("string_rna after " + str(string_rna)) #working
    #print("this output_rna " + str(output_rna))
    return output_rna

def get_codons(peptide, rna_codon_table_array):
    #print("length " + str(len(peptide)))
    rna_dict = {}
    # turn rna_codon_table_array into a dict without stops
    for each in rna_codon_table_array:
        try: #need to handle exceptions when there are "stop codons"
            k,v = each.split(" ")
        except ValueError:
            v = ""
            #print("stop")
        if v != "":
            rna_dict[k] = v
    # go through passed in peptide string and add to list in order of peptides
    for i in range(len(peptide)): 
        #print(i)
        list_of = []
        for key,val in rna_dict.items():
            #print("peptide[i] " + str(peptide[i]))
            if peptide[i] == val: # then match of codon to peptide
                peptide_dict = {val:key}
                #print("peptide_dict " + str(peptide_dict))
                list_of.append(peptide_dict)
                #print("pattern_codons " + str(pattern_codons)) # all dict objs
        pattern_codons.append(list_of)
    print("pattern_codons " + str(pattern_codons))
    # [[{'M': 'AUG'}], [{'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}]]
    #return pattern_codons #list of dicts k = peptide, v = codon, all possible for passed in peptide

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

def get_complement(string):
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
    return reverse_string
    #return reverse_string[::-1]
            
def dna_transcribed_to_rna(dna):
    rna = dna.replace("T", "U")
    return rna

def rna_transcribed_to_dna(rna):
    dna = rna.replace("U", "T")
    return dna

def search_strand_for_pattern(strand, peptide, coding_strand):
    #search_strand_for_pattern(rna_5_3, peptide, True)
    #will pattern_codons always be in same order, since dict?
    #go through all reading frames, look for first codon, find, then look for second
    #pattern_codons is now a list of dicts, not just a dict
    print("mother fuck!")
    possible_peptides = [] #list of codons, don't think I need a dict
    print(peptide)
    #print("pattern_codons " + str(pattern_codons))
    """
    TAQTREAM
    pattern_codons [[{'T': 'ACA'}, {'T': 'ACC'}, {'T': 'ACG'}, {'T': 'ACU'}], [{'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}], [{'Q': 'CAA'}, {'Q': 'CAG'}], [{'T': 'ACA'}, {'T': 'ACC'}, {'T': 'ACG'}, {'T': 'ACU'}], [{'R': 'AGA'}, {'R': 'AGG'}, {'R': 'CGA'}, {'R': 'CGC'}, {'R': 'CGG'}, {'R': 'CGU'}], [{'E': 'GAA'}, {'E': 'GAG'}], [{'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}], [{'M': 'AUG'}]]
    """
    #print("strand " + str(strand))
    #print(len(strand)-2) #8140
    #for each codon in rna/dna strand, reading frame - one at a time
    for i in range(len(strand)-2):
    #for i in range(50): # cut down to shorter data from beginning of strand
        # GCTACCCCCCTGGTGGAGTA
        #print(i) #number
        #print(strand[i:i+3]) # codons, first codon is GCU but doesn't show in editor because too much info - fucking python, grrrrr
        temp_strand_codon = strand[i:i+3]
        temp_list = []
        #print("temp_strand_codon " + str(temp_strand_codon))
        # look for match of "val" from pattern_codon in each reading frame
        for j in range(len(pattern_codons)): 
            #print("j " + str(j))
            #print("pattern_codes[j] " + str(pattern_codons[j].list_of[0]))
            #for m in range(len(pattern_codons[j].list_of)):
            for dict_list in pattern_codons[j]:
                #temp_list = []
                #print(dict_list) #loops through all lists of dicts in order
                for key,val in dict_list.items():
                    #print("key " + str(key)) #key R
                    #print("val " + str(val)) #val AGA
                    # T    ACC, ACG, ACA, ACU
                    if temp_strand_codon == val:
                        # match and add to list of dicts and keep going in same reading frame "i"
                        print("match at " + str(i))
                        #print(key)
                        #print(temp_strand_codon)
                        temp_list.append(val)
                        #print("temp_list " + str(temp_list))
                        # don't need an else or a break here because breaks only work on for and while's...
                        # go where from here???
                        # write a function and call from here
                break
                
            break
        

    #print("possible_peptides " + str(possible_peptides))

def reverse(s): 
    if len(s) == 0: 
        return s 
    else: 
        return reverse(s[1:]) + s[0] 

def get_peptide(codons, rna_codon_table_array):
    rna_dict = {}
    codon_str = ""
    this_peptide = []
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
    
def get_rev_complement_codon_seqs(decoded_peptides_dict, rna_codon_table_array):
    complement_codon_dict = {}
    for key, val in decoded_peptides_dict.items():
        """for each in val:
            print(each)
            get_peptide(each, rna_codon_table_array)"""
        #print(val)
        temp_val = str(val).strip('[]').replace("'","").replace(",","").replace(" ","")
        #print(temp_val)
    #return codon_str[::-1]

dna = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
peptide = "MA"

#gets the codon table
rna_codon_table_array = [line.strip() for line in open("files/rna_codon_table_array.txt")]

"""
data = [line.strip() for line in open("files/peptide_encoding_problem.txt")]
dna = data[0]
peptide = '\n'.join(data[1:])
print(peptide)
#print(dna)
"""

print(peptide_encoding_problem(dna, peptide, rna_codon_table_array))

"""
#Sample output:
ATGGCC
GGCCAT
ATGGCC
"""
#rna_5_3 AUGGCCAUGGCCCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA
#rna_3_5 UCACCCGUUAAUACGGGUACUAUUGAUCUCAGUUCUGGGGGCCAUGGCCAU

