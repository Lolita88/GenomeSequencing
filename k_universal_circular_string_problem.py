# Function finds a universal circular string in binary numbers
# example uses eight binary 3-mers (000, 001, 011, 111, 110, 101, 010, and 100) exactly once
# output is 00011101
from copy import deepcopy
from random import randint
import random
import itertools

def k_universal_circular_string_problem(k):
    # get length of a binary (k)mer
    # if k == 3, then length is 8, 8 binary 3-mers possible
    # if k == 4, then there would be 16 
    # The only thing we need to do is solve the k-Universal Circular String Problem 
    # is to find an Eulerian cycle in DeBruijn(BinaryStringsk). Note that the nodes 
    # of this graph represent all possible binary (k - 1)-mers. A directed edge 
    # connects (k - 1)-mer Pattern to (k - 1)-mer Pattern' in this graph if there 
    # exists a k-mer whose prefix is Pattern and whose suffix is Pattern'. 

    # turn k into binary version of kmers
    # itertools.product - cartesian product of input iterables., equivalent to a nested for-loop
    # product(range(2), repeat=3) or can use product([0,1], ...)
    lst = [list(i) for i in itertools.product([0, 1], repeat=k)]
    binary_list = []
    for each in lst:
        binary_list.append(''.join(map(str, each)))
    #print(binary_list)
    de_bruijn = de_bruijn_graph_from_kmers(binary_list)
    eulerian = find_eulerian_cycle(de_bruijn)
    genome = genome_path_from_eulerian_path(eulerian)
    return genome

def de_bruijn_graph_from_kmers(kmers):
    adjacency_list = {}
    kmers.sort()
    for i in range(len(kmers)):
        pre = kmers[i][0:-1]
        suf = kmers[i][1:]
        if(pre[1:] == suf[:-1]):
            #add to dict for output
            if pre in adjacency_list.keys(): # if already there, append
                adjacency_list[pre] += ("," + suf)
            else: # write anew
                adjacency_list[pre] = pre + " -> " + suf
    print(adjacency_list.values())
    return adjacency_list.values()

def create_adjacency_my_list(my_list):
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

def find_eulerian_cycle(my_list):
    #create adj my_list - key/value pairs
    # keys are the pre, values are the sufs they can point to
    adj_list, circuit_max = create_adjacency_my_list(my_list)

    #reduced adj my_list to keep track of traveled edges
    red_adj_list = {}
    red_adj_list = deepcopy(adj_list) #exact copy of dict
    #arbitrary starting point (if graph is directed/balanced)
    start = random.choice(list(my_list))
    start = start.split(" ")
    start = start[0]
    curr_vert = start
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
            #stack.pop()
    """print(len(circuit))
    circuit.pop(0)
    print(len(circuit))"""
    print(circuit)
    #formatting
    path = start + '->'
    path = ""
    for vert in circuit[:0:-1]: # changing this from [::-1] fixed the error of having 
        # the start node repeated at the end - due to cycle. [:-1] leaves off last element
        path += (vert + '->')
    #print(path.strip())

    return path.strip('->')
    #genome_path_from_eulerian_path(path.strip('->'))

def genome_path_from_eulerian_path(eulerian_path):
    # takes in something like this GGC->GCT->CTT->TTA->TAC->ACC->CCA
    # returns a genome sequence like this GGCTTACCA
    print(eulerian_path)
    kmers = eulerian_path.split("->")
    genome = ""
    for i in range(len(kmers)):
        genome = genome[:i] + kmers[i]
    return genome

k = 3
print(k_universal_circular_string_problem(k))
# sample output for binary 4mer: 0000110010111101
