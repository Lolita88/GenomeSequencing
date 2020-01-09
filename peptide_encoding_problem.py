"""
Peptide Encoding Problem: Find substrings of a genome encoding a given amino acid sequence.
Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
Output: All substrings of Text encoding Peptide (if any such substrings exist).
"""

from pympler import tracker
import copy
tr = tracker.SummaryTracker()

#global vars
cur_frame = 0 #this is for incrementing codon window on strand
reading_frame = 0
peptide_frame = 0
strand_len = 0
end_reading_frame = 0
peptide_codon_len = 0

pattern_codons = [] #a list which hold dicts (pairs)
#pattern_codons = [[{'M': 'AUG'}], [{'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}]]
#pattern_codons = [[{'T': 'ACA'}, {'T': 'ACC'}, {'T': 'ACG'}, {'T': 'ACU'}], [{'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}], [{'Q': 'CAA'}, {'Q': 'CAG'}], [{'T': 'ACA'}, {'T': 'ACC'}, {'T': 'ACG'}, {'T': 'ACU'}], [{'R': 'AGA'}, {'R': 'AGG'}, {'R': 'CGA'}, {'R': 'CGC'}, {'R': 'CGG'}, {'R': 'CGU'}], [{'E': 'GAA'}, {'E': 'GAG'}], [{'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}], [{'M': 'AUG'}]]
temp_list = []
perm_peptide_list = []

def peptide_encoding_problem(dna, peptide, rna_codon_table_array):
    
    print("dna len " + str(len(dna))) #8140 WTF!!!!
    print("peptide " + str(peptide))
    #print("rna_codon_table_array " + str(rna_codon_table_array))
    global strand_len
    global end_reading_frame
    end_reading_frame = len(dna) - len(peptide)*3
    print("end_reading_frame " + str(end_reading_frame))
    strand_len = len(dna)

    # assuming if one string, it's coming in 5' to 3'
    dna_5_3 = dna #coding strand, like mRNA strand
    dna_3_5 = get_reverse_complement(dna_5_3) #rev and opp bp
    #print("dna_5_3 " + str(dna_5_3))
    #print("dna_3_5 " + str(dna_3_5))

    # transcribe to rna
    rna_5_3 = dna_transcribed_to_rna(dna_5_3)
    rna_3_5 = dna_transcribed_to_rna(dna_3_5) #then look here for MA but going to have to rev translate back
    print("rna_5_3 " + str(rna_5_3))
    print("rna_3_5 " + str(rna_3_5))
    
    #uncomment after debugging - DON'T FORGET # not hard coding the codon dict and won't read in for debugging
    get_codons(peptide, rna_codon_table_array) #saves in global var pattern_codons for now
    
    search_reading_frame_for_pattern(rna_5_3, peptide, True)
    #decoded_peptides_found_rna_3_5 = search_strand_for_pattern(rna_3_5, peptide, False)
    
    return(perm_peptide_list)

def search_reading_frame_for_pattern(strand, peptide, is_coding_strand):
    #rename manage_reading_frame?
    global reading_frame # starts at 0
    global peptide_frame # starts at 0
    global end_reading_frame # len of strand minus peptide len x3

    # if not at last reading frame
    if(reading_frame <= end_reading_frame):
        reading_frame += 1 #always increment 1, unless close to end
        peptide_frame = 0 #always start at beginning here
        search_next_peptide(strand, peptide, is_coding_strand)
    else:
        print("at end of strand")

def search_next_peptide(strand, peptide, is_coding_strand):
    global cur_frame
    global reading_frame
    global peptide_codon_len
    global strand_len
    global peptide_frame
    temp_strand_codon = strand[cur_frame:cur_frame+3]
    next_frame = cur_frame+3
    if(reading_frame <= (strand_len - peptide_codon_len)):
        #need to clear out the key, value pair values
        for key, value in pattern_codons[peptide_frame].items():
            print("pattern_codons[peptide_frame] " + str(pattern_codons[peptide_frame]))
            if temp_strand_codon == key:
                print("match at cur_frame " + str(cur_frame))
                #what do I need here? So far the rna peptide code, will need to change back to dna
                temp_list.append(key)
                # increment peptide_frame and cur_frame of strand
                peptide_frame += 1
                cur_frame += 3
                #print("peptide_frame " + str(peptide_frame))
                #print("cur_frame " + str(cur_frame))

                #move codon window down one
                temp_strand_codon = strand[cur_frame:cur_frame+3] 
                #if done, DO SOMETHING, EXIT LOOP TO SOMEWHERE, PROBLEM MAY BE HERE?
                if(len(temp_list) == len(peptide)):
                    #peptide matched, add to perm list, increment strand reading frame
                    perm_peptide_list = copy.deepcopy(temp_list)
                    temp_list.clear()
    
            #pattern_codons [[{'M': 'AUG'}], [{'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}]]
            
            #search next frame, reset peptide frame
       

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
    global pattern_codons
    rna_dict = {}
    # turn rna_codon_table_array into a dict without stops
    for each in rna_codon_table_array:
        try: #need to handle exceptions when there are "stop codons"
            k,v = each.split(" ")
        except ValueError:
            v = ""
        if v != "":
            rna_dict[k] = v
    #print("rna_dict " + str(rna_dict))
    # go through passed in peptide string and add to list in order of peptides
    for i in range(len(peptide)): #peptide is MA
        #print(i)
        # for each letter, create a dict list_of
        list_of = {}

        #loop through codon dict and look for the peptide like MA
        for key,val in rna_dict.items():
            #print("peptide[i] " + str(peptide[i]))
            if peptide[i] == val: # then match of codon to peptide
                #3 letter codon must be the key bc it's unique
                #1 letter abbreviation is not unique, so is val

                list_of[key] = val
                #print("list_of " + str(list_of)) # all dict objs
        pattern_codons.append(list_of)    
    #print("pattern_codons " + str(pattern_codons))
    
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
    return reverse_string[::-1] # reverses the string

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

#gets the codon table
#commented out for debugging - UNDO!!!!
rna_codon_table_array = [line.strip() for line in open("files/rna_codon_table_array.txt")]
#rna_codon_table_array = []
"""
data = [line.strip() for line in open("files/peptide_encoding_problem.txt")]
dna = data[0]
peptide = '\n'.join(data[1:])
"""

dna = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
peptide = "MA"
"""
dna = "GCTACCCCCCTGGTGGAGTATCATAAACTGCGCAGACCCGCGAAGCAATGTGTTCTACCTATACGCGAGCGGAGCGCGAATGAGGACAAAGTTGCAAGATGCCCTTCGAACAATATGTAAGACGGTAAAAAACCACTGATTCAGTAGCGCGGCGGTGCCTTAAAAGATGATCAACCATCCAGGCATGCTTCATAGCCTCCCTAGTTTGTGCCGTCCTAGTCGAACAGCTCAAACTCGGGAAGCGATGTTCATGGCCTCCCTTGTCTGAGCTGTAGTATATTCTCGCGTTAAGACCACGGGCGAGTTATCTGTGCCACTGTTATTTCGTAGGGTGTACGACTATTCTCGACAGCGACTTACATGTGCCGTGTCTCGTATACAGATCTTGCGAAAGCTAGGAGCGACCATTGAATCAGCAATTGCGAAATCCCATGAGACGGTGCTATTTAAGATCTAAAAGACGTAGATAGAGGTAAGAGAGCCGTTTGCCTTTCATAGTATACGGAATTGGGGCGGAAGGCGTCATTACCAAGGCATAGTACTGCTGATCTATCGCTATGAGCTACTAGTCAAGGGTGGCGTCAAACGGTTAATTGGCATTGCTTCCCGAGTCTGAGCCGTTGGTCGGGGTTGCGGCGGCATGTTCACGACATATCCGAAAAAATGGAAACGAAGTCAATTTACACCGGCCGGTATCTTGATTCCTCACCATCAAGCTAATTGAGCAGTTCAGTCCGGTCGACCAGAAGAGCCTCCGCTGTCATGCGTTACTAGTACGACACCATTTAAGTGAGAGAATTTTGGCCTACATCGCTATTACCGACTGCTGCTACTCTCACCCGCCTCCGTTTCTCAGGACGACCGCGTTTGATAGACGTATCCGGCCGATGTAGGCCGTTTGTGATCCCGGGCGAGCTACTGCCGGGGTTTAAGGCATCTTATGCCTAGACCCCACAATAGGGGTGTCTGCAAACTCGCCCTCGTCTAGAGTGAGTTTCGCTGAATTTTACTCCCCATAAACCCGCCTGCCGCCGACCGATGTCCAATTACAGGGTTCGAGTCTCGAGATACCTGTTAAGGGCACGCGTTTCCGCTCGGATGCAACACGGAGAACGTCGTTTCCTCCTAGCGCGAGGAAACATCTCGACGCATGCGGAGGTACTGTCAACGATCTCACTTGGATTGCAGAAGCCTGTGTTTCAACATAAGTGGAACTCCGCCGGGAGACTCGGTGGTTCTAAAGATGTAGTTCCCTTCGTAGATCGTTTATGAGAATAACTTAGTCATTAGCCGCCCTTCGTCCAATATGCCTATGGATTTCGGTGGTCAAGCTGATTTGTTTCTTCGCTGACAGTTACGGCACACGCGGGTCCAGCCTGTCTCACCCTATTGTACTTCGAGGAGTGGGTGTAACTGCGCAAACGAGGGAAGCCATGATTCTTCTCAACCTTGCCATGGCTTCGCGTGTCTGGGCAGTTGGTTCAATATTTGTGGTCCCTTACATTGGAAACGGACCTCAGATTCCCGAAGCGGCGCAATACTTCTGCGTAACATGTGCACCATGGCCTCGCGGGTTTGTGCCGTTGATTTATCTTGCTGTGGGTACGAGAAAGAGGCAGAGTACCACAATAAAGATTGCTGTTCGAATACCCCCCATACTCCCTCTACGTGTGAGCAGGCACGGAAGGCCTCGTGTCTCGTCCTCAAGCTCTTCTTTAGCAACCTGGCCGTTAAACCAGCGACGGCTGTCTAGTTAAAGTTGAGAGCAATTGTCTTGCGCGAGTCAAGCGTGGATCATTTGCAATGCGCAGCGGCTAATGGAATATTCCAAGTCCGTTTCCTAGGTAAGTGTCCTAGTAGACTTTTCTAACAATATACCTCACCATACCACCACTTATGAAGTTCGTACGTAATGCGAATGGTGTATGGACTGTAAACAGCATCTAGCCGTCATTTCGTATTCAGCTTCCCCACCGGAGTCGTGATATGGGTCTTCGAGTGTGATTTGTTCGCGAATATTCGTTGGTGGGATAACGTGGCAACTTGGCTACGGTCTAGAATTTAGGGTTATGGGCCTTCCCACGGGTTAGTTGTGATCCTCCAACTAACCTTGCGTCCGCTATTGGAATTGGCGCGAAGTGGCATGGGCTAATAATGACGGCTCAGAAGATGGCCCAAATTCTTAAATATGTATACAGTTACCATACGTTGGACACATCATGTGATATAGTTACTTTCCACGGACCAGACTACAAACTGTATTGGTGATCCCGACTAGGACGCACGCCTGTGTGAACCATGCCAGTGCGCCAACATATCCTAATAGTGCCGGCTTCGATAGGAAGACGTTACTGCAGAGCATGGTCCGGACATTGTCGTCACTCCAGAACGGGAATTCAGTGAATAGCCGTTTATTGTGCGTTTGGTAGACGAGTGGCGAGAGAGAATAATGGCCCGCTCGCAGTGGGCGTAGCTGGTCATTTTCATGTCATTGGAATACCGACTGCGAACGGCGTTCGACACGGTACTCGCACCGCATATCTGTATGGATCTGCCCTCGACGTTTATGTTGGTGGTATTCCTGTGGCTCCCTGGTTCAATTTATTGCGCGTAACGCACGTATAGCAAACATAACCTCGCGCATAAAGGGTCACCGAGGCCTGGTTGTAGCATTAGTATATTACACTCCAACCTTTAATAATTTTCAGACCGTTACCCGCGCTCATTTCCTCTTTAGGACTGCACAGACTCGCGAAGCTATGAGCCGAACTTTTACACCCCCATACCATTTCGTGCATCGGAGCCTCTTGGCTCAACACTGTTTGTGCTGGTACCGGAGGGGTGCTCTCAAATGACAACTCGCTGGATGCAGCAGCCTTGATTCACCCTGCATGGCCTCGCGTGTCTGTGCTGTGCCGCGGCGGTAGAACACATCTATTAAACTTGGTCACAACGCCGGCGAGTCGCGTTTTTAAGCCTCGGACGATGTACCTTGTGAATGTCGTGTCCCCCTGGAGTTAACAGATACGGCGTTTCCGGACACAGTTCAGCGTTGCATAGAAAAATCGTCCCTGATGTCATTGCGTTCGAGCCTACTAGCTAGATTCCCGATTAATCTAGACGGCGTCGTAGTGGCAGGTTTATCGTTAGATGGGGTCGTGCTACCTAGCAACAGTTATCCTTCCGGGAGGGCACATGGTCCGGAGTCCTTCATTGCCTCACGAGTTTGCGCCGTAAGGCCGCATATCGACCCTTCCCCGCGTCTTCTGTGCCTAAGAATTAAGCAATACCAACTAACCGGGTCCCTGACTCGCTGGCCTAGTATATAATTAAGATATATCGGGTATAACCGCCATCATGGTGCGTTCTCTCCGTCGCAGGACTTCAGACTAGAGACTTGTGTTCCGAAGAGTCTTTACGCACGTCTCGCGATCATTATAAGACACCGTACTGGTAAGTGCGGTATAACGCTGGCCCTGGAATCGCCTTTCATTGTACGTCATTGACGCTCTCGAAGCTGGATACTTACTTAACCAAGCATTACCTAGAGCAACGCAATCCCAGTACAGCTGGTTTCCATTCCTGCGCGGGTTTAATGGTCTTGCAGCACTCGGCGGTCCGCAATTCGGTCCGCGTACACGACCCACTCTAGCAGAGGCGAAGATCTCCGACATACGACCCGGTCCCGTAAGTAGTAAAGAAAAAGTCAAGCCTCCCGAGCCTGCGAGAGAGTGCAGTCTCCGACGGAACGGCGGTGGACGTAGCAATATACGGCCTTTATACGGGCGATTTAGGAATCGCCTTGCGAGCTACGGTGGAACCTAAGGATATAGTCTCGCCTAGCGGAAGTGGTATCCTGTAGTCGTCTACGCCCCGCTCGGTGTACCTCTTGACAAATTTACGACGAAAGAACCTTACCAGCCTGTTATAACCCTTGCGCTCGCTCTAGGCAAACCGGCGCAGTGCAAGATTTTATTCCCTGTGGGTTGCCGGCAGGGAGCACACGCAAATATGGGTTCTGAGCGTACTCTAGTTCAACAATAGGCCACTGCAAAAGAACAGATCTCGTAGTACTTGATATATGGACGTGCCGAAAAGCTATATTATTTGAAGGCAGTCCATGAATAGCCAGTGCTGACCGAAACAATTAACCACACCCAGACGGGGCTTTTGAGATGCGGGCGGCTATGTTAATCTTGCCAAAGGGTAGTGAACGGCACAAACCCGCGAGGCAATGGCACAGACAAGAGAGGCCATGCGGGAAGCGTGATGGAAAAAAATGACTACTGGGGAATTCTATATCCCCGCATAACCCGAGACTAAACGCGCCGTGACGCCAACGCCGTGTTAATTCCGCAAAATAGCAGAGGACGCGTGGAACCGCTGCAATCACGACGACCGCCCAAACACGGGAGGCTATGTGCAAAAGTGCTGACTCGGAGCCTCAGCCGGAGTCCGCCTATGCGACCGCTAAACGTCGCCAGTTTTCAGGTGAGCCACACATAAGGGTAGCATCCGTGGAGTGGGCAACTTGGCGAAGTGGAGCTCCATGCGGATGGGCCACGATCGTTTTAGCATATCCTTGTACGTGTTATGTTAGGCAATTGCCCCAAATGCGCGGTACCCCGACATCCGGTGACGCACCTTTCGCGGGTGCGAGACGATCGTAGGCCCGTTCCTGTGGGGGCGGTGTGATGTGCGAACACCCGTCTGCGGACGAAGTCAGCTATTTACAAAAGATTGTATACCGTATACGCTGCAGTTCTCAGAACTCACTTTTGGGAGCTGCTACTTGGAAAGCGCTAATCAGCAGCAATATTGAACATAGCCTCGCGCGTTTGCGCTGTGTACTCCTGAGATGCCCCTTGGCGGTTGCCTATAGGTAGTAGCCTCGACCCAAAGCTACCTGGAAACCATAAGGAGACCCCAGGCCGAACGGGTCAGTAGGTCCGTTGTTAGACGTCCTAGAAATCGACGCATAAAGTGGCGATCGGCTAATGAGCGATCACACCTTGGTGACTCAAATGTTCAGAGTATCTAGCGTAAAGACAGCCCATCTAGGCATCATTTAAGTCCGACCTGCTCCACGTAATAAGCGTTTGACATGGAGATACGGCGGACGGTCAAAGTCACGTGATTTCTACATCAATGCATTGCCTCGCGGGTTTGTGCCGTGTTCTGTGCGGTACAAGGTCCTACCAGATTTCTGGGGTAAGGCAACTCCCTTTGCCCCTGGGACGTTGCAGAGCGAATGTCAATCATGCTTCGACCCAGTTCATACAAGGTATTATATTTGTATTGCGTAGTGTTCATCACCGTCCTCGGTTGAGCATCCGCGCCTCTGAAGAGATATTCTACAAATGATCATGGCTTCTCTAGTCTGGGCGGTGGGCGCTAACGTAGTTTGTAGGTGACAAATATGAACAACGAGTACCGGCATTACATGGCGTATCTTGCTGGCTCCCGGCCCCGCGGGCGGTATCTGTGACACGTACCTAGGTTCCTCGTAGAGTTAACAGTTATGAACACAACTTCGAAAACCTGATACCACTTCTACCCGGACGTTCCGTAAAGATATATTTGATAGGGGTCTATCGATTTCGCTGATTTTTCGAGCTAGACACAACCGGCATCTATGCCACAATCAAACGGAAATATCACATTCACATAGTAGTATCAATGGTGTCCTTTCCTCGAGTAGGGATATAAAGTTTGCAATCGAAATGAAACTCGTTAACGTAAGAGCGACTACTGTGCATGCACGTACATCGCCTCCCGAGTCTGTGCCGTTTCCCCAAGATCACAACGCGTCAACTCCTCGTGCGGGATCACATGCAGAAGAGAGGTGCAAGAACCCCGGTGTTCCGCTAGGGAGTATTACGGCAAACCATCAGGTAGGAAGGGGCGAAGCGATAGTGAGATGGTTGCAATTCGGCGTCCCGCCCCGGGGAGTAACGACTTATACTGTAAGCTCACTGGGATAACACTTCTTTACCAGCCAATCAACATATAAATTCTGGCGATCTAGTGGAGTTCACTGAAGTTGATGATCGTTGCCAGGATATAATTGAGCGGCCAACAACCCATCCCGTGTAGAAGGCTACGATAAGTGCAAGGGGTTCGGTGTTACAAGTCCTTTACAGGGAACCGGACAGCCCCATTTGGACATCTTTAACTCCTTGGGGCCGGAAGTTCGCCTTGGGCCTATCCCCGTGCTAAGTGGCATGTAGGTGTTGGATATAAAAGAAGCCCCAGTTTGAATATCCAGATCCGAGTGGGGATCAGCAGATCGTTATGCGTACTTATGAAGTGATATCCTCTGCGGACATTAAACCGTGATCTCGGTATGTGTCGTCCCCAGGCACGCCGGGCGTAGTGGCCGAACGGCCTAGTGTTCACGAGACATCCTGGGCAGTTGAGGTTGGCTTCTCGGAATTGCTTGGTGCATCAATTTACAAGGAATTCCCTCTAAGTTGTAGGACTTGGCGACTAACAGCAAGTTTCAGGATCTGAAAACTGTCGCCCTCCATGACTCAGGACCCCACACGCGTCGCCATCGCAATCATCAGAGATACTGTGTTGGTCGACATCCGAGAGTGAGATATCCAATGCTTTGCCCGATACCACGGAACCAAATGAGACGGCACAGACACGCGAGGCCATGAAGCGGGAAAAGTCTCGCAGAATAGAGTGAAAAATGATCTCAGACCCTAGGTGAAACTTTGGGGCTCTCATATCCCGGGGAACCCGGTCGCTGGAACAGGGAACATAGCTTCCCTGGTCTGTGCGGTTATCGTCTGGGTATCCTCCCTTCGCCTTCATTCAAGTTGAGCACTGAAGCGTAATGACTGGCAACCCCGGCATAGAAAAAAACTCCTCGGTGTTGTGTTTAGAGGATCGTATGAGCCCTTTCCACTGAATTTCATGCTCCGGGCTGCACCTTTGCCTGCGCGCTGGGGACATGACACTAGACATCCCCGTGCCCTACGGTATACCGAATGGATAGAGGCGACCTTCTCTGAAACTTTAATCGGAATCGGCTTACTAGCATTGTCGATCTGAAAAGTTGCATATCAGGGAAACGCTGATGATCGAGGGGGGGAGGCTGTAGGACCGAGAAAGTATTACTGGGTGTTCGGGGGGTTAAGATGCGTGGCCCCGCTTGTGCCGCCTAACGCGACTCCTTGCTTCTAGAGCCGGTCGGCCCATACTCCGAAAATAACTCGCTCGGCTGTACCAACCTGTACAGAAATGAATCGGATGTAGCCGGTCAGTGGACATCGTTGCTCTCCAGCGGATACATCTGTTAACCAACTGTAGGATCTCTGAGCCAATGAGCTGGTCGCCCGGCACCACGCGTCGCGTGGTCTGGTCATTTGGCTTCCAAAGGACAGCGACGTTATCATCCTTAGCTTGCCTATCACACTCTCTAAGGATCTCTTACGGCCAAATCTCCTGAGCGAGGTATAAGCTTTCTATAGTTCGTGATCTCGGAGGAGTATGACCGAAGTGATTAAGAGTTTACTGGATCAACGACGGCACTGCATGTAACCCCCTATTGAGGTCGACGATTTCATCACTCTATGAGTACATCGTCCGCACGGGCACCTGTGTCCTACGTACTCCGTTTCGACCCTTGACGTACGTGTTACCCACGCAATTCCGGTAAGCATCTTCCTTACACTCATCTCGAGAGGCCGTTAATAAGGGTGTTGTTCTAGTCTGCTAGAGCGGCTGCCAATTTGCTGACTGGACAGCAAGTTTACCGGATACCAATAACTGGGAGAATAGAATCGTAAAAAGCCTATTTAACAAGTTGAAGTTTTTGCTGGGAGGTCCGGATGAGGCGTTAGTTGCAGACCTAGATAAATACGGACCAATTTGGGACAATCGATGTATCGTAGTTTAGACTAGGCCCGTCAAGTTGCGGTGCTCGAGGAGAGAGTACGTGCAGCCCCTTTGCATACGGAGAAGATGCCCCTTTAACTCGACGCCTGCATGATACTCGTGTTCTGTGCCTATATCTGGCAATTCCCGCCTCAATTCACATATCTCGAGTGGTGCACACTCCATTTAAGCCAGATTAGCTCTAACCCTTTGGGAGGCTGCAGAAGAGGGCGTTATGTCTAGTAGGAGCATTTT"
peptide = "TAQTREAM"
"""
print(peptide_encoding_problem(dna, peptide, rna_codon_table_array))

