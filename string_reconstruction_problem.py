# genomePath(EulerianPath(DeBruijnGraph(kmers)))
# This code takes k length and list of kmers, creates a DeBruijn graph from them (Define DeBruijn here).
# Then, with the DeBruijn data, find a Eulerian path and stitch that into a genome path. 
from copy import deepcopy
from random import randint
import random

def string_reconstruction_problem(k, kmers):
    #change this, it's bad format, duh
    return genome_path_from_eulerian_path(find_eulerian_cycle(de_bruijn_graph_from_kmers(kmers)))

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
    return adjacency_list.values()

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
    #print(adj_list)
    return adj_list, circuit_max

#Find start/end nodes
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
    for vert in circuit[::-1]:
        path += (vert + '->')
    return path.strip('->')

def genome_path_from_eulerian_path(eulerian_path):
    # takes in something like this GGC->GCT->CTT->TTA->TAC->ACC->CCA
    # returns a genome sequence like this GGCTTACCA
    kmers = eulerian_path.split("->")
    genome = ""
    for i in range(len(kmers)):
        genome = genome[:i] + kmers[i]
    return genome
"""    
k = 4
kmers = [
    "CTTA",
    "ACCA",
    "TACC",
    "GGCT",
    "GCTT",
    "TTAC"
]
"""
data = [line.strip() for line in open("files/string_reconstruction_problem.txt")]
k = data[0]
kmers = data[1:]

print(string_reconstruction_problem(k, kmers))