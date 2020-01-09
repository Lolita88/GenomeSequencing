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

from pympler import tracker
import copy
tr = tracker.SummaryTracker()

#pattern_codons = [] # needs to be a list of dicts
#pattern_codons = [[{'M': 'AUG'}], [{'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}]]
pattern_codons = [[{'T': 'ACA'}, {'T': 'ACC'}, {'T': 'ACG'}, {'T': 'ACU'}], [{'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}], [{'Q': 'CAA'}, {'Q': 'CAG'}], [{'T': 'ACA'}, {'T': 'ACC'}, {'T': 'ACG'}, {'T': 'ACU'}], [{'R': 'AGA'}, {'R': 'AGG'}, {'R': 'CGA'}, {'R': 'CGC'}, {'R': 'CGG'}, {'R': 'CGU'}], [{'E': 'GAA'}, {'E': 'GAG'}], [{'A': 'GCA'}, {'A': 'GCC'}, {'A': 'GCG'}, {'A': 'GCU'}], [{'M': 'AUG'}]]
frame = 0 #this is for incrementing codon window on strand
reading_frame = 0
peptide_frame = 0
strand_len = 0
end_reading_num = 0
done = False

temp_list = []
perm_peptide_list = []

def peptide_encoding_problem(dna, peptide, rna_codon_table_array):
    #looking to output all possible encoded peptides
    #codon_table = rna_codon_table_array
    #decoded_peptides_found = {}
    print("dna len " + str(len(dna))) #8140 WTF!!!!
    global strand_len
    global end_reading_num
    end_reading_num = len(dna) - len(peptide)*3
    strand_len = len(dna)

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
    
    #uncomment after debugging - DON'T FORGET
    #get_codons(peptide, rna_codon_table_array)
    
    #try to find codon pattern in rna 5'to 3'
    decoded_peptides_found_rna_5_3 = search_strand_for_pattern(rna_5_3, peptide, True)
    #decoded_peptides_found_rna_3_5 = search_strand_for_pattern(rna_3_5, peptide, False)
    #{0: ['AUG', 'GCC'], 6: ['AUG', 'GCC']}
    #['GGC', 'CAU']

    #print("decoded_peptides_found_rna_5_3 " + str(decoded_peptides_found_rna_5_3))
    return(perm_peptide_list)

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
    global pattern_codons
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
    #print("rna_dict " + str(rna_dict))
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
    #print("pattern_codons " + str(pattern_codons))
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
    print("search_strand_for_pattern - from strand")
    return search_next_peptide_frame(strand, peptide, coding_strand) 

def stop_searching():
    global perm_peptide_list
    print("stop searching")
    return(perm_peptide_list)

def search_next_reading_frame(strand, peptide, coding_strand):
    global frame
    global reading_frame
    global peptide_frame
    global strand_len
    global end_reading_num
    #end_reading_num = len(dna) - len(peptide)*3 #set in prior function call
    strand_len = len(dna)
    #print("search_next_reading_frame")
    reading_frame += 1
    #print("reading_frame " + str(reading_frame))
    frame = reading_frame # this has to start at reading_frame, initially 0
    #print("global reset of frame " + str(frame))
    
    peptide_frame = 0 # reset to 0 each time reading_frame moves
    #if(reading_frame + (len(peptide)*3) <= len(strand)):
    if(end_reading_num >= reading_frame):
        print("search_next_peptide_frame - from search reading frame")
        search_next_peptide_frame(strand, peptide, coding_strand)
    else:
        print("should be done for real")

def search_next_peptide_frame(strand, peptide, coding_strand):
    global frame
    global reading_frame
    global peptide_frame
    global strand_len
    global end_reading_num
    global done
    done = False
    peptide_codon_len = len(peptide) * 3 
    print("search_next_peptide_frame - main")
    # look for match of "val" from pattern_codon in each reading frame
    temp_strand_codon = strand[frame:frame+3] # starts at 0 (technically, at reading_frame)
    #print("len strand " + str(strand_len)) 
    #next_frame = frame+3
    # have to start with a for here so I can control the break later down in order to stop the recursion, ahhhh
    # if not at end, check
    # look in list of lists of dicts
    if(reading_frame > end_reading_num):
        done = True
    if(reading_frame <= end_reading_num):
        for dict_list in pattern_codons[peptide_frame]:
            if(not done):
                for key,val in dict_list.items():
                    if(not done):
                        #print("key " + str(key)) #key R
                        print("val " + str(val)) #val AGA
                        print("temp_strand_codon " + str(temp_strand_codon)) #AUG
                        print("frame " + str(frame))
                        if temp_strand_codon == val:
                            # if match, add to temp_list of dicts and keep going in same reading_frame
                            print("match at " + str(frame))
                            temp_list.append(val)
                            print("temp_list " + str(temp_list))
                            peptide_frame += 1
                            #print("peptide_frame " + str(peptide_frame))
                            frame += 3 
                            #print("frame " + str(frame))
                            temp_strand_codon = strand[frame:frame+3] # if at end, shows error
                            if(len(temp_list) == len(peptide)): # match, but could be more if longer strand
                                #peptide matched, add to perm list, increment strand reading frame
                                perm_peptide_list = copy.deepcopy(temp_list)
                                """
                                if coding_strand == False:
                                    perm_peptide_list.reverse()
                                    print("perm_peptide_list reverse" + str(perm_peptide_list))
                                """
                                temp_list.clear()
                                print("perm_peptide_list " + str(perm_peptide_list))
                                # set peptide back to start - don't need to, will catch in following reading frames
                                print("calling 2")
                                done = True
                                search_next_reading_frame(strand, peptide, coding_strand)
                                #break
                                return
                            else:
                                # look for next match
                                print("search_next_peptide_frame - inside main 1")   
                                search_next_peptide_frame(strand, peptide, coding_strand)
                                break
                    else:
                        #break
                        return
            else:
                #break
                return
    # this is running for each dict obj
    if(reading_frame <= end_reading_num):
        temp_list.clear()
        search_next_reading_frame(strand, peptide, coding_strand) 
    else:
        done = True
        return
    

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
"""
dna = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
peptide = "MA"
"""
#gets the codon table
#commented out for debugging - UNDO!!!!
#rna_codon_table_array = [line.strip() for line in open("files/rna_codon_table_array.txt")]
rna_codon_table_array = []

"""
data = [line.strip() for line in open("files/peptide_encoding_problem.txt")]
dna = data[0]
peptide = '\n'.join(data[1:])
"""

dna = "GCTACCCCCCTGGTGGAGTATCATAAACTGCGCAGACCCGCGAAGCAATG"
"""
TGTTCTACCTATACGCGAGCGGAGCGCGAATGAGGACAAAGTTGCAAGATGCCCTTCGAACAATATGTAAGACGGTAAAAAACCACTGATTCAGTAGCGCGGCGGTGCCTTAAAAGATGATCAACCATCCAGGCATGCTTCATAGCCTCCCTAGTTTGTGCCGTCCTAGTCGAACAGCTCAAACTCGGGAAGCGATGTTCATGGCCTCCCTTGTCTGAGCTGTAGTATATTCTCGCGTTAAGACCACGGGCGAGTTATCTGTGCCACTGTTATTTCGTAGGGTGTACGACTATTCTCGACAGCGACTTACATGTGCCGTGTCTCGTATACAGATCTTGCGAAAGCTAGGAGCGACCATTGAATCAGCAATTGCGAAATCCCATGAGACGGTGCTATTTAAGATCTAAAAGACGTAGATAGAGGTAAGAGAGCCGTTTGCCTTTCATAGTATACGGAATTGGGGCGGAAGGCGTCATTACCAAGGCATAGTACTGCTGATCTATCGCTATGAGCTACTAGTCAAGGGTGGCGTCAAACGGTTAATTGGCATTGCTTCCCGAGTCTGAGCCGTTGGTCGGGGTTGCGGCGGCATGTTCACGACATATCCGAAAAAATGGAAACGAAGTCAATTTACACCGGCCGGTATCTTGATTCCTCACCATCAAGCTAATTGAGCAGTTCAGTCCGGTCGACCAGAAGAGCCTCCGCTGTCATGCGTTACTAGTACGACACCATTTAAGTGAGAGAATTTTGGCCTACATCGCTATTACCGACTGCTGCTACTCTCACCCGCCTCCGTTTCTCAGGACGACCGCGTTTGATAGACGTATCCGGCCGATGTAGGCCGTTTGTGATCCCGGGCGAGCTACTGCCGGGGTTTAAGGCATCTTATGCCTAGACCCCACAATAGGGGTGTCTGCAAACTCGCCCTCGTCTAGAGTGAGTTTCGCTGAATTTTACTCCCCATAAACCCGCCTGCCGCCGACCGATGTCCAATTACAGGGTTCGAGTCTCGAGATACCTGTTAAGGGCACGCGTTTCCGCTCGGATGCAACACGGAGAACGTCGTTTCCTCCTAGCGCGAGGAAACATCTCGACGCATGCGGAGGTACTGTCAACGATCTCACTTGGATTGCAGAAGCCTGTGTTTCAACATAAGTGGAACTCCGCCGGGAGACTCGGTGGTTCTAAAGATGTAGTTCCCTTCGTAGATCGTTTATGAGAATAACTTAGTCATTAGCCGCCCTTCGTCCAATATGCCTATGGATTTCGGTGGTCAAGCTGATTTGTTTCTTCGCTGACAGTTACGGCACACGCGGGTCCAGCCTGTCTCACCCTATTGTACTTCGAGGAGTGGGTGTAACTGCGCAAACGAGGGAAGCCATGATTCTTCTCAACCTTGCCATGGCTTCGCGTGTCTGGGCAGTTGGTTCAATATTTGTGGTCCCTTACATTGGAAACGGACCTCAGATTCCCGAAGCGGCGCAATACTTCTGCGTAACATGTGCACCATGGCCTCGCGGGTTTGTGCCGTTGATTTATCTTGCTGTGGGTACGAGAAAGAGGCAGAGTACCACAATAAAGATTGCTGTTCGAATACCCCCCATACTCCCTCTACGTGTGAGCAGGCACGGAAGGCCTCGTGTCTCGTCCTCAAGCTCTTCTTTAGCAACCTGGCCGTTAAACCAGCGACGGCTGTCTAGTTAAAGTTGAGAGCAATTGTCTTGCGCGAGTCAAGCGTGGATCATTTGCAATGCGCAGCGGCTAATGGAATATTCCAAGTCCGTTTCCTAGGTAAGTGTCCTAGTAGACTTTTCTAACAATATACCTCACCATACCACCACTTATGAAGTTCGTACGTAATGCGAATGGTGTATGGACTGTAAACAGCATCTAGCCGTCATTTCGTATTCAGCTTCCCCACCGGAGTCGTGATATGGGTCTTCGAGTGTGATTTGTTCGCGAATATTCGTTGGTGGGATAACGTGGCAACTTGGCTACGGTCTAGAATTTAGGGTTATGGGCCTTCCCACGGGTTAGTTGTGATCCTCCAACTAACCTTGCGTCCGCTATTGGAATTGGCGCGAAGTGGCATGGGCTAATAATGACGGCTCAGAAGATGGCCCAAATTCTTAAATATGTATACAGTTACCATACGTTGGACACATCATGTGATATAGTTACTTTCCACGGACCAGACTACAAACTGTATTGGTGATCCCGACTAGGACGCACGCCTGTGTGAACCATGCCAGTGCGCCAACATATCCTAATAGTGCCGGCTTCGATAGGAAGACGTTACTGCAGAGCATGGTCCGGACATTGTCGTCACTCCAGAACGGGAATTCAGTGAATAGCCGTTTATTGTGCGTTTGGTAGACGAGTGGCGAGAGAGAATAATGGCCCGCTCGCAGTGGGCGTAGCTGGTCATTTTCATGTCATTGGAATACCGACTGCGAACGGCGTTCGACACGGTACTCGCACCGCATATCTGTATGGATCTGCCCTCGACGTTTATGTTGGTGGTATTCCTGTGGCTCCCTGGTTCAATTTATTGCGCGTAACGCACGTATAGCAAACATAACCTCGCGCATAAAGGGTCACCGAGGCCTGGTTGTAGCATTAGTATATTACACTCCAACCTTTAATAATTTTCAGACCGTTACCCGCGCTCATTTCCTCTTTAGGACTGCACAGACTCGCGAAGCTATGAGCCGAACTTTTACACCCCCATACCATTTCGTGCATCGGAGCCTCTTGGCTCAACACTGTTTGTGCTGGTACCGGAGGGGTGCTCTCAAATGACAACTCGCTGGATGCAGCAGCCTTGATTCACCCTGCATGGCCTCGCGTGTCTGTGCTGTGCCGCGGCGGTAGAACACATCTATTAAACTTGGTCACAACGCCGGCGAGTCGCGTTTTTAAGCCTCGGACGATGTACCTTGTGAATGTCGTGTCCCCCTGGAGTTAACAGATACGGCGTTTCCGGACACAGTTCAGCGTTGCATAGAAAAATCGTCCCTGATGTCATTGCGTTCGAGCCTACTAGCTAGATTCCCGATTAATCTAGACGGCGTCGTAGTGGCAGGTTTATCGTTAGATGGGGTCGTGCTACCTAGCAACAGTTATCCTTCCGGGAGGGCACATGGTCCGGAGTCCTTCATTGCCTCACGAGTTTGCGCCGTAAGGCCGCATATCGACCCTTCCCCGCGTCTTCTGTGCCTAAGAATTAAGCAATACCAACTAACCGGGTCCCTGACTCGCTGGCCTAGTATATAATTAAGATATATCGGGTATAACCGCCATCATGGTGCGTTCTCTCCGTCGCAGGACTTCAGACTAGAGACTTGTGTTCCGAAGAGTCTTTACGCACGTCTCGCGATCATTATAAGACACCGTACTGGTAAGTGCGGTATAACGCTGGCCCTGGAATCGCCTTTCATTGTACGTCATTGACGCTCTCGAAGCTGGATACTTACTTAACCAAGCATTACCTAGAGCAACGCAATCCCAGTACAGCTGGTTTCCATTCCTGCGCGGGTTTAATGGTCTTGCAGCACTCGGCGGTCCGCAATTCGGTCCGCGTACACGACCCACTCTAGCAGAGGCGAAGATCTCCGACATACGACCCGGTCCCGTAAGTAGTAAAGAAAAAGTCAAGCCTCCCGAGCCTGCGAGAGAGTGCAGTCTCCGACGGAACGGCGGTGGACGTAGCAATATACGGCCTTTATACGGGCGATTTAGGAATCGCCTTGCGAGCTACGGTGGAACCTAAGGATATAGTCTCGCCTAGCGGAAGTGGTATCCTGTAGTCGTCTACGCCCCGCTCGGTGTACCTCTTGACAAATTTACGACGAAAGAACCTTACCAGCCTGTTATAACCCTTGCGCTCGCTCTAGGCAAACCGGCGCAGTGCAAGATTTTATTCCCTGTGGGTTGCCGGCAGGGAGCACACGCAAATATGGGTTCTGAGCGTACTCTAGTTCAACAATAGGCCACTGCAAAAGAACAGATCTCGTAGTACTTGATATATGGACGTGCCGAAAAGCTATATTATTTGAAGGCAGTCCATGAATAGCCAGTGCTGACCGAAACAATTAACCACACCCAGACGGGGCTTTTGAGATGCGGGCGGCTATGTTAATCTTGCCAAAGGGTAGTGAACGGCACAAACCCGCGAGGCAATGGCACAGACAAGAGAGGCCATGCGGGAAGCGTGATGGAAAAAAATGACTACTGGGGAATTCTATATCCCCGCATAACCCGAGACTAAACGCGCCGTGACGCCAACGCCGTGTTAATTCCGCAAAATAGCAGAGGACGCGTGGAACCGCTGCAATCACGACGACCGCCCAAACACGGGAGGCTATGTGCAAAAGTGCTGACTCGGAGCCTCAGCCGGAGTCCGCCTATGCGACCGCTAAACGTCGCCAGTTTTCAGGTGAGCCACACATAAGGGTAGCATCCGTGGAGTGGGCAACTTGGCGAAGTGGAGCTCCATGCGGATGGGCCACGATCGTTTTAGCATATCCTTGTACGTGTTATGTTAGGCAATTGCCCCAAATGCGCGGTACCCCGACATCCGGTGACGCACCTTTCGCGGGTGCGAGACGATCGTAGGCCCGTTCCTGTGGGGGCGGTGTGATGTGCGAACACCCGTCTGCGGACGAAGTCAGCTATTTACAAAAGATTGTATACCGTATACGCTGCAGTTCTCAGAACTCACTTTTGGGAGCTGCTACTTGGAAAGCGCTAATCAGCAGCAATATTGAACATAGCCTCGCGCGTTTGCGCTGTGTACTCCTGAGATGCCCCTTGGCGGTTGCCTATAGGTAGTAGCCTCGACCCAAAGCTACCTGGAAACCATAAGGAGACCCCAGGCCGAACGGGTCAGTAGGTCCGTTGTTAGACGTCCTAGAAATCGACGCATAAAGTGGCGATCGGCTAATGAGCGATCACACCTTGGTGACTCAAATGTTCAGAGTATCTAGCGTAAAGACAGCCCATCTAGGCATCATTTAAGTCCGACCTGCTCCACGTAATAAGCGTTTGACATGGAGATACGGCGGACGGTCAAAGTCACGTGATTTCTACATCAATGCATTGCCTCGCGGGTTTGTGCCGTGTTCTGTGCGGTACAAGGTCCTACCAGATTTCTGGGGTAAGGCAACTCCCTTTGCCCCTGGGACGTTGCAGAGCGAATGTCAATCATGCTTCGACCCAGTTCATACAAGGTATTATATTTGTATTGCGTAGTGTTCATCACCGTCCTCGGTTGAGCATCCGCGCCTCTGAAGAGATATTCTACAAATGATCATGGCTTCTCTAGTCTGGGCGGTGGGCGCTAACGTAGTTTGTAGGTGACAAATATGAACAACGAGTACCGGCATTACATGGCGTATCTTGCTGGCTCCCGGCCCCGCGGGCGGTATCTGTGACACGTACCTAGGTTCCTCGTAGAGTTAACAGTTATGAACACAACTTCGAAAACCTGATACCACTTCTACCCGGACGTTCCGTAAAGATATATTTGATAGGGGTCTATCGATTTCGCTGATTTTTCGAGCTAGACACAACCGGCATCTATGCCACAATCAAACGGAAATATCACATTCACATAGTAGTATCAATGGTGTCCTTTCCTCGAGTAGGGATATAAAGTTTGCAATCGAAATGAAACTCGTTAACGTAAGAGCGACTACTGTGCATGCACGTACATCGCCTCCCGAGTCTGTGCCGTTTCCCCAAGATCACAACGCGTCAACTCCTCGTGCGGGATCACATGCAGAAGAGAGGTGCAAGAACCCCGGTGTTCCGCTAGGGAGTATTACGGCAAACCATCAGGTAGGAAGGGGCGAAGCGATAGTGAGATGGTTGCAATTCGGCGTCCCGCCCCGGGGAGTAACGACTTATACTGTAAGCTCACTGGGATAACACTTCTTTACCAGCCAATCAACATATAAATTCTGGCGATCTAGTGGAGTTCACTGAAGTTGATGATCGTTGCCAGGATATAATTGAGCGGCCAACAACCCATCCCGTGTAGAAGGCTACGATAAGTGCAAGGGGTTCGGTGTTACAAGTCCTTTACAGGGAACCGGACAGCCCCATTTGGACATCTTTAACTCCTTGGGGCCGGAAGTTCGCCTTGGGCCTATCCCCGTGCTAAGTGGCATGTAGGTGTTGGATATAAAAGAAGCCCCAGTTTGAATATCCAGATCCGAGTGGGGATCAGCAGATCGTTATGCGTACTTATGAAGTGATATCCTCTGCGGACATTAAACCGTGATCTCGGTATGTGTCGTCCCCAGGCACGCCGGGCGTAGTGGCCGAACGGCCTAGTGTTCACGAGACATCCTGGGCAGTTGAGGTTGGCTTCTCGGAATTGCTTGGTGCATCAATTTACAAGGAATTCCCTCTAAGTTGTAGGACTTGGCGACTAACAGCAAGTTTCAGGATCTGAAAACTGTCGCCCTCCATGACTCAGGACCCCACACGCGTCGCCATCGCAATCATCAGAGATACTGTGTTGGTCGACATCCGAGAGTGAGATATCCAATGCTTTGCCCGATACCACGGAACCAAATGAGACGGCACAGACACGCGAGGCCATGAAGCGGGAAAAGTCTCGCAGAATAGAGTGAAAAATGATCTCAGACCCTAGGTGAAACTTTGGGGCTCTCATATCCCGGGGAACCCGGTCGCTGGAACAGGGAACATAGCTTCCCTGGTCTGTGCGGTTATCGTCTGGGTATCCTCCCTTCGCCTTCATTCAAGTTGAGCACTGAAGCGTAATGACTGGCAACCCCGGCATAGAAAAAAACTCCTCGGTGTTGTGTTTAGAGGATCGTATGAGCCCTTTCCACTGAATTTCATGCTCCGGGCTGCACCTTTGCCTGCGCGCTGGGGACATGACACTAGACATCCCCGTGCCCTACGGTATACCGAATGGATAGAGGCGACCTTCTCTGAAACTTTAATCGGAATCGGCTTACTAGCATTGTCGATCTGAAAAGTTGCATATCAGGGAAACGCTGATGATCGAGGGGGGGAGGCTGTAGGACCGAGAAAGTATTACTGGGTGTTCGGGGGGTTAAGATGCGTGGCCCCGCTTGTGCCGCCTAACGCGACTCCTTGCTTCTAGAGCCGGTCGGCCCATACTCCGAAAATAACTCGCTCGGCTGTACCAACCTGTACAGAAATGAATCGGATGTAGCCGGTCAGTGGACATCGTTGCTCTCCAGCGGATACATCTGTTAACCAACTGTAGGATCTCTGAGCCAATGAGCTGGTCGCCCGGCACCACGCGTCGCGTGGTCTGGTCATTTGGCTTCCAAAGGACAGCGACGTTATCATCCTTAGCTTGCCTATCACACTCTCTAAGGATCTCTTACGGCCAAATCTCCTGAGCGAGGTATAAGCTTTCTATAGTTCGTGATCTCGGAGGAGTATGACCGAAGTGATTAAGAGTTTACTGGATCAACGACGGCACTGCATGTAACCCCCTATTGAGGTCGACGATTTCATCACTCTATGAGTACATCGTCCGCACGGGCACCTGTGTCCTACGTACTCCGTTTCGACCCTTGACGTACGTGTTACCCACGCAATTCCGGTAAGCATCTTCCTTACACTCATCTCGAGAGGCCGTTAATAAGGGTGTTGTTCTAGTCTGCTAGAGCGGCTGCCAATTTGCTGACTGGACAGCAAGTTTACCGGATACCAATAACTGGGAGAATAGAATCGTAAAAAGCCTATTTAACAAGTTGAAGTTTTTGCTGGGAGGTCCGGATGAGGCGTTAGTTGCAGACCTAGATAAATACGGACCAATTTGGGACAATCGATGTATCGTAGTTTAGACTAGGCCCGTCAAGTTGCGGTGCTCGAGGAGAGAGTACGTGCAGCCCCTTTGCATACGGAGAAGATGCCCCTTTAACTCGACGCCTGCATGATACTCGTGTTCTGTGCCTATATCTGGCAATTCCCGCCTCAATTCACATATCTCGAGTGGTGCACACTCCATTTAAGCCAGATTAGCTCTAACCCTTTGGGAGGCTGCAGAAGAGGGCGTTATGTCTAGTAGGAGCATTTT"
"""
peptide = "TAQTREAM"

print(peptide_encoding_problem(dna, peptide, rna_codon_table_array))

"""
#Sample output:
ATGGCC
GGCCAT
ATGGCC
"""
#rna_5_3 AUGGCCAUGGCCCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA
#rna_3_5 UCACCCGUUAAUACGGGUACUAUUGAUCUCAGUUCUGGGGGCCAUGGCCAU

