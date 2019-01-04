"""
Peptide Encoding Problem: Find substrings of a genome encoding a given amino acid sequence.
Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
Output: All substrings of Text encoding Peptide (if any such substrings exist).
"""

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
    pattern_codons = get_codons(peptide, rna_codon_table_array)
    #print("pattern_codons " + str(pattern_codons))
    # {'M': ['AUG'], 'A': ['GCA', 'GCC', 'GCG', 'GCU']}
    # [{'T': 'ACA'}, {'T': 'ACC'}, {'T': 'ACG'}, {'T': 'ACU'}, {'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}, {'Q': 'CAA'}, {'Q': 'CAG'}, {'T': 'ACA'}, {'T': 'ACC'}, {'T': 'ACG'}, {'T': 'ACU'}, {'R': 'AGA'}, {'R': 'AGG'}, {'R': 'CGA'}, {'R': 'CGC'}, {'R': 'CGG'}, {'R': 'CGU'}, {'E': 'GAA'}, {'E': 'GAG'}, {'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}, {'M': 'AUG'}]

    #try to find codon pattern in rna 5'to 3'
    decoded_peptides_found_rna_5_3 = search_strand_for_pattern(rna_5_3, peptide, pattern_codons, True)
    decoded_peptides_found_rna_3_5 = search_strand_for_pattern(rna_3_5, peptide, pattern_codons, False)
    #{0: ['AUG', 'GCC'], 6: ['AUG', 'GCC']}
    #['GGC', 'CAU']

    list_1 = convert_found_codons_for_output(decoded_peptides_found_rna_5_3)
    list_2 = convert_found_codons_for_output(decoded_peptides_found_rna_3_5)
    return(join_peptides_from_lists(list_1, list_2))
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
    #print("peptide " + str(peptide))
    #print(type(peptide))
    #print("length " + str(len(peptide)))
    encoded_peptides_dict = {}
    encoded_peptides = [] # needs to be a list of dicts
    rna_dict = {}

    # turn rna_codon_table_array into a dict without stops
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

    # go through passed in peptide string and add to list in order of peptides
    # for T, ex. first in list, create a dict for the obj
    
    for i in range(len(peptide)): #TAQTREAM
        print(i) # correct 0-7
        #peptide_dict = {}
        #encoded_peptides[i].append(peptide_dict) # make dict here
        #print("peptide_dict " + peptide_dict)
        for key,val in rna_dict.items():
            #print("v " + str(v))
            #print("peptide[i] " + str(peptide[i]))
            if peptide[i] == val: # then match of codon to peptide
                peptide_dict = {}
                #print("match")
                #print(key)
                #print(val)
                peptide_dict = {val:key}
                #print("peptide_dict " + str(peptide_dict))
                encoded_peptides.append(peptide_dict) #seems to be breaking...
                print("encoded_peptides " + str(encoded_peptides)) # all dict objs
    return encoded_peptides #list of dicts k = peptide, v = codon, all possible for passed in peptide

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

def search_strand_for_pattern(strand, peptide, pattern_codons, coding_strand):#will pattern_codons always be in same order, since dict?
    #go through all reading frames, look for first codon, find, then look for second
    #reading frame
    #codon_list[]
    #pattern_codons is now a list of dicts, not just a dict
    temp_peptide_dict = {}
    # go through strand
    for i in range(len(strand)-2):
        #temp_codon_location = 0
        #print(peptide)
        #print(i) #number
        #print(strand[i:i+3]) # codons
        #for each in pattern_codons:
        for key, val in pattern_codons.items():
            # dict_items([('M', ['AUG']), ('A', ['GCA', 'GCC', 'GCG', 'GCU'])])
            #print(len(val))
            #print(key) # M or A
            #print(val) # ['AUG'] or ['GCA', 'GCC', 'GCG', 'GCU']
            if(len(val) == 1):
                #if single codon for first position, check
                temp_val = str(val).strip('[]').replace("'","")
                #print(temp_val) # AUG
                #print(strand[i:i+3])
                if(temp_val == strand[i:i+3]): #looping through strand ONE reading frame at a time
                    #if match, add to temp_peptide_dict[], move down 3, check next
                    #print("match at " + str(i))
                    #print(strand[i:i+3])
                    temp_peptide_dict[i] = [] # wtf? it's a dict above
                    #print(type(temp_peptide_dict))
                    temp_peptide_dict[i].append(temp_val)
                    #move to next
                else:
                    #no match, move to next reading frame, need to make sure this is happening
                    break
            elif(len(val) > 1): #need to check all, if in this case: ('A', ['GCA', 'GCC', 'GCG', 'GCU'])]
                for each in val:
                    #print(each) 
                    #print(strand[i+3:i+6])
                    if(each == strand[i+3:i+6]):
                        #print(strand[i+3:i+6])
                        #print("match2")
                        #print(i)
                        temp_peptide_dict[i].append(each)
                        #print(temp_peptide_dict)
            #else:
                # do nothing
                #print("oh hello")
    if coding_strand:
        #print("coding strand")
        #print(temp_peptide_dict)
        return temp_peptide_dict
    else:
        #print("non coding strand" + str(temp_peptide_dict))
        new_temp_peptide_dict = {} 
        for key,val in temp_peptide_dict.items():
            #rev codon letters
            #print(val)
            new_temp_peptide_dict[key] = []
            temp_list = []
            for each in val:
                #print("each " + str(each))
                new_val = each[::-1]
                #print("new_val " + str(new_val))
                #rev codon order 
                temp_list.append(new_val)
                #print("temp_iist working " + str(temp_list)) #working
                #temp_list[::-1]
            temp_list.reverse()
            #print("temp_iist rev" + str(temp_list)) #working
            new_temp_peptide_dict[key] = []
            new_temp_peptide_dict[key] = temp_list
            #print("new_temp_peptide_dict[key] " + str(new_temp_peptide_dict[key])) #need rev compliment GGC, GGC CAT
            for key,val in new_temp_peptide_dict.items():
                rev_each = []
                for each in val:
                    temp_each = get_complement(each)
                    #print(temp_each)
                    rev_each.append(temp_each)
                    #print(each)
                    #print("rev_each " + str(rev_each))
                new_temp_peptide_dict[key] = rev_each
                #print("new_temp_peptide_dict[key].value() " + str(new_temp_peptide_dict[key].val))
                #new_temp_peptide_dict[key] += get_reverse_complement(temp_list[i])
                #print("new_temp_peptide_dict " + str(new_temp_peptide_dict))
    return new_temp_peptide_dict  

def reverse(s): 
    if len(s) == 0: 
        return s 
    else: 
        return reverse(s[1:]) + s[0] 

def get_peptide(codons, rna_codon_table_array):
    rna_dict = {}
    encoded_peptides_dict = {}
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
    
    #for each in peptide: #peptide is peptide passed in
        #this_peptide.append(each)
    #print(this_peptide) # M A
    """
    for each in codons:
        temp_list = []
        for val in rna_dict:
            print(val)
            #print(rna_dict[key])
            #print(key)
            if each == rna_dict[val]: #match!
                # add to list to add to dict key
                temp_list.append(key) # for each key, one or more codons
        #print(temp_list)
        encoded_peptides_dict[each] = temp_list
    print("encoded_peptides_dict " + str(encoded_peptides_dict))
    return encoded_peptides_dict"""

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
"""
dna = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
peptide = "MA"
"""
rna_codon_table_array = [line.strip() for line in open("files/rna_codon_table_array.txt")]

data = [line.strip() for line in open("files/peptide_encoding_problem.txt")]
dna = data[0]
peptide = '\n'.join(data[1:])
#print(peptide)
#print(dna)

print(peptide_encoding_problem(dna, peptide, rna_codon_table_array))

"""
#Sample output:
ATGGCC
GGCCAT
ATGGCC
"""
#rna_5_3 AUGGCCAUGGCCCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA
#rna_3_5 UCACCCGUUAAUACGGGUACUAUUGAUCUCAGUUCUGGGGGCCAUGGCCAU

