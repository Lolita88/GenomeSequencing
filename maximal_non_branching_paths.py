"""
def maximal_non_branching_paths(graph):
    paths = {}
    for line in graph:
        node = line.strip('\n')
        node = node.replace(' -> ', ' ')
        node = node.split(' ')
        paths.setdefault(node[0],[]) #adj my_list gets the start of the node pair(ie. first num)
        for num in node[1].split(','): #num gets assigned the end node of the pair. Split on comma needed when multiple end nodes
            paths[node[0]].append(num)
    ind,outd = in_and_out_degree(paths)
    print(ind,outd)

def in_and_out_degree(adjacency_list):
    '''
    return the in and out degree lists for a given graph's adjacency list
    '''
    ind = {}
    outd = {}
    print(adjacency_list)
    for key, value in adjacency_list.items():
        #print("k " + str(key))
        #print("v " + str(value))
        outd[key] = len(value) #counts the out-degree of th3 keys(start nodes)
        for kk in value:
            ind[kk] = ind.get(kk,0)+1 #counts in-degree at index[value] 
            # get(kk,0) 0 is value to be returned if key is not found
            # since ind is empty at that index, the +1 is adding the single count for a degree
            #print("ind[kk] " + str(ind[kk]))
    return (ind,outd)
  
    MaximalNonBranchingPaths(Graph)
        Paths ← empty list
        for each node v in Graph
            if v is not a 1-in-1-out node
                if out(v) > 0
                    for each outgoing edge (v, w) from v
                        NonBranchingPath ← the path consisting of the single edge (v, w)
                        while w is a 1-in-1-out node
                            extend NonBranchingPath by the outgoing edge (w, u) from w 
                            w ← u
                        add NonBranchingPath to the set Paths
        for each isolated cycle Cycle in Graph
            add Cycle to Paths
        return Paths
    
"""
def maximal_non_branching_paths(adjacency_list):
    ''' 
    Implement MaximalNonBranchingPaths.
    Input: The adjacency list of a graph whose nodes are integers.
    Output: The collection of all maximal nonbranching paths in this graph.    
    '''
    #print(adjacency_list)
    #tempad = adjacency_list
    #why have dadj? same exact as adjacency list
    """dadj = {}
    for key, value in adjacency_list.items():
        dadj[key] = dadj.get(key,[]) + value[:]
    """
    degree = in_and_out_degree(adjacency_list)
    paths = []
    for v, e in adjacency_list.items():
    #        # cut-off optimization, skip v node if already in path list
    #        if visited(v) : 
    #            continue
        deg_in = degree[0].get(v,0)
        deg_out = degree[1].get(v,0)
        if (deg_in == 1 and deg_out == 1):
            # vertex v is 1-in-1-out node
            # could be part of a new isolated cycle, check this...
            if not visited(v, paths):
                cycle = isolated_cycle(v, degree, adjacency_list)
                if cycle:
                    paths.append(cycle)
        elif (deg_out>0):
            # explore vertex v outgoing branches
            for w in e:
                paths.append(non_branching_path([v,w], degree, adjacency_list))
    return paths

def is_one_in_one_out(vertex, degree):
    deg_in = degree[0].get(vertex,0)
    deg_out = degree[1].get(vertex,0)
    return (deg_in == 1) and (deg_out == 1)

def visited(vertex, paths):
    # linear time visited-vertex implementation
    # fixme : use dict to get constant time
    for path in paths:
        if vertex in path:
            return True
    return False

def isolated_cycle(vertex, degree, adjacency_list):
    '''
    return isolated cycle including vertex v if any
    '''
    cycle = [vertex]
    #print(cycle[-1])
    while is_one_in_one_out(cycle[-1], degree):
        cycle.append(adjacency_list[cycle[-1]][0])
        if cycle[0]==cycle[-1]:
            return cycle
    return None
def non_branching_path(edge, degree, adjacency_list):
    '''
    return the non-branching path starting with edge edge
    '''
    branch = edge[:]
    #print(branch[-1])
    while is_one_in_one_out(branch[-1], degree):     
        branch.append(adjacency_list[branch[-1]][0])
    return branch
    

def in_and_out_degree(adjacency_list):
    '''
    return the in and out degree lists for a given graph's adjacency list
    '''
    ind = {}
    outd = {}
    for key, value in adjacency_list.items():
        outd[key] = len(value)
        for kk in value:
            ind[kk] = ind.get(kk,0)+1
    return (ind,outd)

graph = [
    "1 -> 2",
    "2 -> 3",
    "3 -> 4,5",
    "6 -> 7",
    "7 -> 6"
]
paths = {}
for line in graph:
    node = line.strip('\n')
    node = node.replace(' -> ', ' ')
    node = node.split(' ')
    paths.setdefault(node[0],[]) #adj my_list gets the start of the node pair(ie. first num)
    for num in node[1].split(','): #num gets assigned the end node of the pair. Split on comma needed when multiple end nodes
        paths[node[0]].append(num)

print(maximal_non_branching_paths(paths))
"""
sample output:
1 -> 2 -> 3
3 -> 4
3 -> 5
7 -> 6 -> 7
"""