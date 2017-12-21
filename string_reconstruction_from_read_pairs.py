# Taking paired prefixes and paired suffixes. For example, say I have a read AAAT|TTTA. 
# This would be broken down into two nodes, (AAA|TTT) and (AAT|TTA). Then, the read-pair 
# would represent an edge between them, i.e., (AAA|TTT) -> (AAT|TTA)
# Then, we can find an Eulerian path in this resulting paired-end de Bruijn graph, 
# and then we can reconstruct a string from the "read1" path and another string from 
# the "read2" path, and we can overlay them because we know d

from copy import deepcopy
from random import randint
import random

def string_reconstruction_from_read_pairs(k, d, kmer_pairs):
    #print(kmer_pairs)
    de_bruijn = paired_de_bruijn_graph_from_kmers(k, d, kmer_pairs)
    #print(de_bruijn)
    eulerian = find_eulerian_cycle(de_bruijn)
    genome = genome_path_from_eulerian_path(k, d, eulerian)
    print('\n')
    print(genome)

def find_eulerian_cycle(my_list):
    #create adj my_list
    adj_list, circuit_max = create_adjacency_list(my_list)

    #reduced adj my_list to keep track of traveled edges
    red_adj_list = {}
    red_adj_list = deepcopy(adj_list) #exact copy of dict
   
    #Find start node (graph must be directed/ubalanced)
    red_adj_list, start_node = find_start_end_nodes(red_adj_list)
   
    start = start_node
    curr_vert = start_node
    path = [curr_vert]
    stack = []
    circuit = []
    while len(circuit) != circuit_max:
        if red_adj_list[curr_vert] != []: #if neighbor nodes exist
            stack.append(curr_vert)
            pick = randint(0,len(red_adj_list[curr_vert])-1) #what is pick?
            temp = deepcopy(curr_vert) #why use deepcopy
            curr_vert = red_adj_list[temp][pick]
            red_adj_list[temp].remove(curr_vert)
        else:
            circuit.append(curr_vert)
            curr_vert = stack[len(stack)-1]
            stack.pop()

    #formatting
    path = start + '->'
    circuit_end = len(circuit)
    count = 0
    path = start + '->'
    
    # uses a reverse for loop. i starts at the high end and counts backwards.
    for i in range(len(circuit) - 1, -1, -1):
        if i == 0:
            path += circuit[i]
        else:
            path += (circuit[i] + '->')
    #print("path " + str(path))
    return path
    #print(path.strip('->'))

def create_adjacency_list(my_list):
    #saving my_list into a dictionary
    adj_list = {}
    circuit_max = 0
    for line in my_list:
        node = line.strip('\n')
        node = node.replace(' -> ', ' ')
        node = node.split(' ')
        adj_list.setdefault(node[0],[]) #adj my_list gets the start of the node pair(ie. first num)
        for num in node[1].split(','): #num gets assigned the end node of the pair. Split on comma needed when multiple end nodes
            adj_list[node[0]].append(num)
            circuit_max += 1
    return adj_list, circuit_max

def find_start_end_nodes(red_adj_list):
    start = {} 
    # set the length of the start node neighbors in each value of the dict 
    for each in red_adj_list:
        start.setdefault(each, 0)
        start[each] += len(red_adj_list[each])
    end = {}
    # loop through and count the end nodes appearances, ex. 3 appears twice
    for key in red_adj_list:
        for value in red_adj_list[key]:
            end.setdefault(value, 0)
            end[value] += 1
    for key in end:
        try: #need to handle exceptions in case keys aren't found
            # where there is a mismatch in path numbers, that signifies start/end node
            if start[key] != end[key]:
                if start[key] > end[key]:
                    start_node = key # out degree is one greater than in degree
                if start[key] < end[key]:
                    end_node = key # out degree is one less than in degree
        except KeyError: 
            end_node = key
    for key in start:
        try:
            if end[key] != start[key]:
                if end[key] < start[key]:
                    start_node = key
                if end[key] > start[key]:
                    end_node = key
        except KeyError:
            start_node = key
    red_adj_list[end_node] = []
    return red_adj_list, start_node

def paired_de_bruijn_graph_from_kmers(k, d, kmer_pairs): # don't need d or k
    # kmer_pair ex: AAAT|TTTA
    # pre_nodes ex: AAA|TTT
    # suf_nodes ex: AAT|TTA
    adjacency_list = {}
    pre_nodes = []
    suf_nodes = []
   
    for i in range(len(kmer_pairs)):
        temp_pair1, temp_pair2 = kmer_pairs[i].split("|")
        # is a pre/suf always one less than the full kmer????
        print("temp_pair1 " + str(temp_pair1))
        print("temp_pair2 " + str(temp_pair2))
        pre_nodes = temp_pair1[0:-1] + "|" + temp_pair2[0:-1]
        suf_nodes = temp_pair1[1:] + "|" + temp_pair2[1:]
        
        print("pre_nodes " + str(pre_nodes))
        print("suf_nodes " + str(suf_nodes))
        
        #print(pre_nodes[1:k])
        print(str(pre_nodes[1:k-2])) + " " + str(suf_nodes[0:2]))
        #print(str(suf_nodes[0:k-d]) + " " + str(suf_nodes[-k+1:-1]))
        #print("\n")        
        
        #if((pre_nodes[1:k-1] == suf_pair[0:k-1]) and (pre_nodes[-d:] == suf_pair[-k+1:-1])):
        if((pre_nodes[1:k-d+1] == suf_nodes[0:k-d]) and (pre_nodes[-d:] == suf_nodes[-k+1:-1])):
            #add to dict for output
            if pre_nodes in adjacency_list.keys(): # if already there, append
                adjacency_list[pre_nodes] += ("," + suf_nodes)
            else: # write anew
                adjacency_list[pre_nodes] = pre_nodes + " -> " + suf_nodes
    print("adjacency_list.values " + str(adjacency_list.values()))
    return adjacency_list.values()

def genome_path_from_eulerian_path(k, d, eulerian_path):
    # takes in something like this:
    # GTG|GTG->TGG|TGA->GGT|GAG->GTC|AGA->TCG|GAT->CGT|ATG->GTG|TGT->TGA|GTT->GAG|TTG->AGA|TGA
    # returns a genome sequence like this GTGGTCGTGAGATGTTGA
    print(eulerian_path)
    kmers = eulerian_path.split("->")
    print(kmers)
    pre_string = ""
    suf_string = ""
    genome = ""
    for i in range(len(kmers)):
        pre, suf = kmers[i].split("|")
        #print("pre " + pre[0])
        #print("suf " + suf[0])
        pre_string = pre_string + pre[0]
        suf_string = suf_string + suf[2]
        """
        if(i == len(kmers)-1): # if at end of suf_string, need to bring down the last chars
            suf_string = suf_string + suf
        else:
            suf_string = suf_string + suf[2]"""
    print(pre_string)
    print(suf_string)
    genome = pre_string + suf_string[d:]
    print("genome is " + genome)
    print("should be GTGGTCGTGAGATGTTGA")
    return genome
"""
        FirstPatterns ? the sequence of initial k-mers from GappedPatterns
        SecondPatterns ? the sequence of terminal k-mers from GappedPatterns
        PrefixString ? StringSpelledByPatterns(FirstPatterns, k)
        SuffixString ? StringSpelledByPatterns(SecondPatterns, k)
        for i = k + d + 1 to |PrefixString|
            if the i-th symbol in PrefixString does not equal the (i - k - d)-th symbol in SuffixString
                return "there is no string spelled by the gapped patterns"
        return PrefixString concatenated with the last k + d symbols of SuffixString
"""

def hamming_distance(p, q):
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1
    return count

# sample input
k = 4
d = 2
k_pairs = [
    "GAGA|TTGA",
    "TCGT|GATG",
    "CGTG|ATGT",
    "TGGT|TGAG",
    "GTGA|TGTT",
    "GTGG|GTGA",
    "TGAG|GTTG",
    "GGTC|GAGA",
    "GTCG|AGAT"
]
# sample output: GTGGTCGTGAGATGTTGA

"""
raw_data = [line.strip() for line in open("files/kmer_pairs_dataset3.txt")]
if(raw_data[0] == "Input"):
    del raw_data[0]
k, d = map(int, raw_data[0].split(" "))
del raw_data[0] # deletes k,d before import of kmer pairs
correct_answer = ""
k_pairs = []
for line in raw_data:
    if(line == "Output"):
        break
    else:
        k_pairs.append(line)
"""
my_answer = string_reconstruction_from_read_pairs(k, d, k_pairs)
print(my_answer)
"""
#test code if you have the answer to compare
correct_answer = "GAAAGGTACAAATACTGGCGACCTCGCTGTTCGACACTTCATCACTGCTCCGGGGCGCTCAGGAGGGACGGTTCCCTGTACCATTGGAAGTCAATAGTCTAAGGTACAAAGAGAAGACCCGACCCGACAGAGGGGGTTCTGCGCCGGGTTTCGAGCTTGTAACCCCCCAGAGAATTAGATCCACCGTCTGTGTGGACAAAGTAGTAAAGCTAGCATACCAAATTGAAATTCGGAGTTTGACTACCAGATCCACGCATACGCTGCACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTAGAAATTCAGAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGGGGGTAATTCGTAGTTAGGTACAGAAAACTCCCGGACAGAACCGCATATAACCGATGAAGCAAGGGTTCTTCATTTAATACGACCCTAACCGGTATTGCTGCTAGCTTGATTTTCCTAGCAATCTAAACTCTATGTATGAGGCCACTCGGACGCCCGCTAGTGCCGGCAGCTAGCTACTGCCCTTCACCAGGAGCACGCACTATGCCTATCGGGCAATGCTGATCATACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGATCCCTCTGCAGAAAGCGGTGGCGGCGGGTCTAAGCAAGTCCAACGCAATACCAGGAAATCACCGTATCGTTAGCGACCAGTAGGTGATGGTTTGTAAGTTCGGACTACAGGCGGATGTGTCCCCGCCAGTTAAAAGTCGACTTTCTGTTACAACTGCTCCCTACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTCTAATGATCCCCACAGATCGTGTTTCAACGTTGAAGTCTGAATGGGTTCGTGAATATAGCCATCCAACGTGGACAATAAGATGAGCTTTATAGTTTCCGATCCTCATGGCGATCGAATAAGATCTATCCGCTTGTGTGTGTACGAGTCGCCGACTAACCGGTCTTGGGATATATACGTCACAGATTAAGTACTCGTCACGAGCTTGAATGGGAAGATAAGTAGACTCTTGTCGGGCACACACAGAGACTCCGACGCATCGAGATCGCAAACACTGCCTCCAGCCGGGGGATGCTAATCGTCGCGGTCGGTCCGAGCTTTATTCTACATCGTGGTGTTTCCGACCGAGCCATAAGAACAGTGTCCAAGTCACAAGAGGAGCACGCGGTGGAGGTCGTTCGCTATACAATATATTTGCAACTGTGTCTGGCATCACGCGCATTTCTCACACTTCCAAACGTGCTGCATTCTAAATGATTTCATGAATAGATTGTCTACTAGTTCACCCAAGGTATTACAGCACTGGTCATGTGCCGCTCTGGCACGGCTAGTATCAGGGCCGACTGTGTCCTAGGCGGCTGTTTTCGGGAGCCCAAGGGAAGCAATCAATGCGTTGACTGACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGATGACCCTAGGGCGGGCTGTTTGAGTGGGCTATCGGCGACCATTTCGCCCCGCAAGCCCCCGTCACGATACCGAGACCCTGAAGCTCATAAACGCCTATCTTTGTTGCATGAACAACGGGAGTAAGCGAGGCCAGGCCATACGGTTTCGGAGCCGCAGAATAGCCTTTACACGACCTCTACACAACCCAAAGTGAATATCCACGGGGTATTGTTTGTGCTGCTACCATGCGCAGTCAACATCGCCCACCGGCGATGTGTTTCAGATCTGCAGGCCCATCAACCGTTGTGACACCACCCCGGCTTTCAGAAGCAGTATGTCGGACACATTGACCTGTAGCGTTAGTTCTGTACAAGGGACCCTGCTCACTCGAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTTGGGAATCTAATGCGGTCTGCCATGGGACCCCTAACAACAAGGGACGCCTCCACCTGTCTAGGAGGAAAGACTTTACACACATCTTCTAGTTTCGAGAAGCACCGTAGCCAGTGGACCCTGAGTGGTTACAAACAAATCGCAGTTTAGCGCTTACCGACAAAGGCGGGAGCTTCGTTCACATTAGGATTGAAAGAACTTAAGAGTCTGTAGGCTCGGAGGTCTCTATATATACCATCTAGTCGTCCGGGGCATGATTAACTAAGAGTTGATCTGAGTCGGAACATAATGCCTGATCTGACCCCAATTCACTACGGTCGCACGTTCCGGGAACACCTACCGATCAGTCCGGAACTGTGACCTAAGAAGTCTCAAGCCTTTACGTCAATGTTCCCGGTGAAGGACTGTGTAACGGTCGCCTTCGCGCCCCCCCATAGGCCGGTCCTTCTCGTTGCAGGATAGCTAAGTCCCATATAGAGTTGTTGGTGTACCATTACGCTGATTTTACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTGGTGTGTTCAACTCAGCATACCCGGTTAGTCTGGAGCACTCCCCGTGCCTACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTGCCGAAATCCTTGTCAAGTAGTGGCATCTACTGCCGCGGGGCAGGACTGATGCTGACCCAAGACCACGCTCCTATCAGCCGGTGGCGCATCAGGGTGGGCTATAAGTTATATTCCTACTGTACGGCTGAAGTCAGCTGTAGTCAGGGAGCGGTTCCTGAGCCGGCTGATTCCGCTCGTAATGCGCTATGTAGAAGCATAGTTAGCCTCGCGCTCGTGTGTGGGCAGTCGTATAAGTAGTTTAGCTCCCGATGCGAAGGAGTTGCAAGTACCTACAAACTCGTAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGTCCGGTTCATATTAACCATGCACCAAGGTTTGACTAAATCAACTCGTGGGAATCCGACGTGACAAAATCCCCAGATATGCCGGGGGTGCACGTGAATACGTCGTAAGTTGAGCGTCCTATGACGGGAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGCCAAAGTAATACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGAATCCTTGACTGGCCTGCCGAGTGTTATCTCGCTGCTATCCCCCCCCAAGGTAGAAATGGAAGTGGGATCCAGCGCACCAAGGCACTTCACACAGGCATTACCCCAGCACCACGAATTAGCTTGCAGCTAAAGACAGGGTATTTTACGGAGTATATGATCTCTGTGAGGTACCGTATTCACACATCGTGGGATGTCTGCCATGAGCTTTTCCATTAGTATCCGGCGAGTTTTGATCCAAGTTACCAAACAAGGTTGTCCTCCAGGTCCTACGTGCTGAACGGCCAACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGCCGGATCGGCCCTGACTCAAGTAAGGTCTGGTTGCTTGTCACTACATAAAGCCACGGAAGGGTGCGCGGCCCCAAAGCTGCGTCCGGATTCGACTCCCGTTTGCCTGGCTTCTCGACGAATCTAACGTTCTCATTAACCGAAAACCCTGAGCGGCTTAACCTCATTCGTCCCAGAATCAAACCCATCGTGTATCACCGTTGGCCCAGCAGGGAAGACAAGACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGTGTGGCTCCCCAACCGAAGAATAAGATCCCTTCGCCGCCACAGAAGCAACACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGATTTAAAGTCTACCGTGGGGAGCCGGACGAGAACAAGACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGTAGGCGGCTACCTTCTGCCCATATCTCGGAATCCTTAGGGTCTTTAGCTGCTGGCAACACGGGGATGCCTTAGTCGCGGGGGAGCCGTTAGATCGGTTTCAGACTTGCCGACACCGTTCACCGTGGCCGGCCAACCGGCCGGTTCTGTGCTCCACTGAGTTTAGATAGGAATCCATACCTATAGTTTTACGTACACCAACTGGTTAACAAGCCGTCTCCGCGATGACAAACGGTGGGGGCACGAGCTCGAGGTAAGAGGTTGCGATCCGATTACAATGTATGACTACTTATAATGGTCTTACCTTAAATAGGGAACGGGTTTACAGATTACATGTCGGCGAAGACGTTACTGTATTTCGGCCAATGAGCAAATTCCCACCAGGCTCGCTCGGCTTAAATAGAATAGTCAGATTGGCTCTGAATGCTTTGGGCGTGCATTGAGAGCAGCCATATGTAATATTAAATGGCAGTACAAATCATACAGTTCAGAACTGCCGACAGCGCAGGAGTTTAAGGGTATCGAATATTGCGCTATCCGTGAGTGCTCTTAGCGATGGGGGGGCGGACCCTAAGTCTGACCCCCTCTCCTACCTTCTACGGATTACTATTATTGGCACTTGATGAGTAATCATTTCTAGCAAGAGTCTTATAAGGTAAACAATAACTTAGTAAGTGAGTCATGTAGTGTGCTTCCAGGACGAGTCGGCAAACTCTGTAGTCTTATGCTCATGTCTGACCTGCTGTGCCCAAAATTCTCTTCGTAAGGAGGGCTTTATAATGTTATGGGCACGACTTCGCATTGGGTCCACGCCCCAGGACTTCAGCATGTTATTTTGGGTTGCAGGATTTAAGAGAGCCTCATGCGTTGATAAGCCAAAGTGGGGGTATGGTGGGACCTCTCACCATGAGAGTTAAGTTAACTCACCGTGGCTCAAAAAAAGCTGGTTAGAATCTGCGAGTAATACGAGCGGGAAAATCTGGAATAACAGAAGCGACACCCTGACCTACAGTCGTTCAGTACTAGGTTACAAGTGAACCACTCGCGGATATAGTCAGGCGGGGATGTCCCGCGCGTTGATTAGACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGCGAGTTAATTGAAGTTTGCCTAGACACCGCTGAGGCTGGTTCGACATACCCTTAGGGAGGCCAAGCTATATAAAACCAAGATCATTGACCCCCTACGTGATACGTGATTTCAAACTTTACAATCATTAGGGTCGCCAGTGGAGAATCTATAGAATCTTTTCTACAGGCTACAGAGAAGCATTTTTCACAGGACCGCGTGGCGCAAACAATCCGATGGGGACCATCTGTGAACTCCCATACGTGACTATTCTGTGTCACATGAGGGGAGCTAGGGGGATTGAGTGCTCATGTCGGTTGGAGACCATTTTGAGTGCACAAGGGACCCTGCTCACTCGATTGGGAATCTAATGCGGTCTGCCATGGGGCGAAGCATACTTACCTTGATCAACGCAGTGATTATTCATCTGAAGAGGATTGGGATAATACTGCGAACATATTGGAAAATTAACTGATTTATCTTCTGATCGATTCCCACACTCCACGAATTGGGGTGCCATGCTCCCATAGTAGGCCCTAGAGATGCCGATCATTCCGCAGGTGTGCCTAAGTGGACAGTCACTTGGCACTTAGGCCAATAAGTACAACAAAGGGATCAGTGGGCAAATTATCAGCGTACAATTCCCAGATATATAGGCGGCGAGAAAAGCTTCAAAAGACTTAATTTACTAGCCTCCTACAAACTCTAGATGAGGATTGGCTCTTGATGCTAGCGTTTTCATTTTCCATTACAAGACATTAGGCTGATAATTGCAGAGATTGGCGGCGTAGACTGACAGTCGCGATCAATCTGCGTGTTA"
print("len of my ans: " + str(len(my_answer)))
print("len of correct ans: " + str(len(correct_answer)))
print ('hamming distance:', hamming_distance(my_answer, correct_answer))
"""
