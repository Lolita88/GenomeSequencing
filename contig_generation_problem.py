def contig_generation_problem(kmers):
    """
    Contig Generation Problem: Generate the contigs from a collection of reads (with imperfect coverage).
     Input: A collection of k-mers Patterns.
     Output: All contigs in DeBruijn(Patterns).
    """
    de_bruijn = de_bruijn_graph_from_kmers(kmers)
    paths = maximal_non_branching_paths(de_bruijn)
    contigs = convert_paths_to_contigs(paths)
    return contigs

def convert_paths_to_contigs(paths):
    contigs = []
    for path in paths:
        nodes = path.strip('\n')
        nodes = nodes.replace(' -> ', ' ')
        nodes = nodes.split(' ')
        # go through, if not at end, compare overlap and add on
        text = ""
        for i in range(len(nodes)):
            text = text[0:i] + nodes[i] # start with prior node, then add on beginning at overlap each time
        contigs.append(text)
    contigs.sort()
    return contigs

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

def maximal_non_branching_paths(adjacency_list):
    """
    Input: The adjacency list of a adjacency_list whose nodes are integers.
    Output: The collection of all maximal nonbranching paths in this adjacency_list.    
    """
    """
    This code takes an adjacency list, formats to graph/dict, determines which nodes are 1-in-1-out.
    From that, creates paths that do not branch (non-maximal-branching-paths) aka paths consisting
    of single edges. 
    """
    
    #reformats data to remove -> and save to graph
    graph = {}
    for line in adjacency_list:
        node = line.strip('\n')
        node = node.replace(' -> ', ' ')
        node = node.split(' ')
        graph.setdefault(node[0],[]) #graph gets the start of the node pair(ie. first num)
        for num in node[1].split(','): #num gets assigned the end node of the pair. Split on comma needed when multiple end nodes
            graph[node[0]].append(num)
    
    degree = in_and_out_degree(graph)
    paths = []
    # for each node in gragh
    for v,e in graph.items():
        # if v is not a 1-in-1-out node
        if not one_in_one_out(v, degree):
            if degree[1].get(v,0) > 0:
                # for each outgoing edge (v, w) from v
                for w in graph[v]: # for v = 3 it will be 2 paths
                    # NonBranchingPath ← the path consisting of the single edge (v, w)
                    non_branching_path = str(v) + " -> " + str(w)
                    # while w is a 1-in-1-out node
                    while one_in_one_out(w, degree): #!!!! w is not a 1in out node!!!!
                        # extend NonBranchingPath by the outgoing edge (w, u) from w 
                        for u in graph[w]:
                            #extend NonBranchingPath by the outgoing edge (w, u) from w 
                            non_branching_path += " -> " + u
                            w = u # breaks infinite loop but overwriting w prevents that
                    paths.append(non_branching_path)
        else: # if isololated path (all nodes should be 1-in-1-out)
            if not visited(v, paths): 
                temp_cycle = isolated_cycle(v, degree, graph)
                if temp_cycle:
                    paths.append(temp_cycle)
    return paths

def one_in_one_out(vertex, degree):
    
    deg_in = degree[0].get(vertex,0)
    deg_out = degree[1].get(vertex,0)
    #print("deg_in " + str(deg_in))
    #print("deg_out " + str(deg_out))
    return (deg_in == 1) and (deg_out == 1)

def visited(vertex, paths):
    # checks for vertex in paths
    #print("paths " + str(paths))
    # this works for large data set
    for path in paths:
        # because this is a circular path, starting point does not matter - only order
        if vertex in path:
            return True
    return False

def isolated_cycle(vertex, degree, graph):
    cycle = [vertex]
    while one_in_one_out(cycle[-1], degree):
        cycle.append(graph[cycle[-1]][0])
        if cycle[0]==cycle[-1]: # if first == last, end of
            # format this
            cycle_path = ""
            for each in cycle:
                cycle_path += each + " -> "
            return cycle_path[:-4] # removes last ->
    return None
    
def in_and_out_degree(graph):
    '''
    return the in and out degree lists for a given adjacency_list's adjacency list
    '''
    ind = {}
    outd = {}
    for key, value in graph.items():
        outd[key] = len(value)
        for kk in value:
            ind[kk] = ind.get(kk,0)+1
    return (ind,outd)
"""
kmers = [
    "ATG",
    "ATG",
    "TGT",
    "TGG",
    "CAT",
    "GGA",
    "GAT",
    "AGA"
]
"""
kmers = [line.strip() for line in open("files/contig_generation_data.txt")]

print('\n'.join(contig_generation_problem(kmers)))
#print(contig_generation_problem(kmers))

"""
Sample output:
AGA ATG ATG CAT GAT TGGA TGT
"""