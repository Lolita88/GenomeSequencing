
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
                    #print("w " + str(w))
                    non_branching_path = str(v) + " -> " + str(w)
                    #print("non_branching_path " + str(non_branching_path))
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
            #print(cycle_path)
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
adjacency_list = [
    "1 -> 2",
    "2 -> 3",
    "3 -> 4,5",
    "6 -> 7",
    "7 -> 6"
]
"""
adjacency_list = [line.strip() for line in open("files/maximal_non_branching_paths_data.txt")]

print('\n'.join(maximal_non_branching_paths(adjacency_list)))

"""
sample output:
1 -> 2 -> 3
3 -> 4
3 -> 5
7 -> 6 -> 7
"""