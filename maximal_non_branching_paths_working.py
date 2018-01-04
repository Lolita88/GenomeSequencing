
def maximal_non_branching_paths(adjacency_list):
    """
    MaximalNonBranchingPaths(adjacency_list)
        Paths ← empty list
        for each node v in adjacency_list
            if v is not a 1-in-1-out node
                if out(v) > 0
                    for each outgoing edge (v, w) from v
                        NonBranchingPath ← the path consisting of the single edge (v, w)
                        while w is a 1-in-1-out node
                            extend NonBranchingPath by the outgoing edge (w, u) from w 
                            w ← u
                        add NonBranchingPath to the set Paths
        for each isolated cycle Cycle in adjacency_list
            add Cycle to Paths
        return Paths
    Input: The adjacency list of a adjacency_list whose nodes are integers.
    Output: The collection of all maximal nonbranching paths in this adjacency_list.    
    """
    #reformat data from list with pointers "->" to dict
    graph = {}
    for line in adjacency_list:
        node = line.strip('\n')
        node = node.replace(' -> ', ' ')
        node = node.split(' ')
        graph.setdefault(node[0],[]) #adj my_list gets the start of the node pair(ie. first num)
        for num in node[1].split(','): #num gets assigned the end node of the pair. Split on comma needed when multiple end nodes
            graph[node[0]].append(num)

    #adjacency_list = graph

    degree = in_and_out_degree(graph)
    paths = []
    for v, e in graph.items():
        deg_in = degree[0].get(v,0)
        deg_out = degree[1].get(v,0)
        if (deg_in == 1 and deg_out == 1):
            # vertex v is 1-in-1-out node
            # could be part of a new isolated cycle, check this...
            if not visited(v, paths):
                cycle = isolated_cycle(v, degree, graph)
                if cycle:
                    paths.append(cycle)
        elif (deg_out>0):
            # explore vertex v outgoing branches
            for w in e:
                paths.append(non_branching_path([v,w], degree, graph))
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

def isolated_cycle(vertex, degree, graph):
    '''
    return isolated cycle including vertex v if any
    '''
    cycle = [vertex]
    while is_one_in_one_out(cycle[-1], degree):
        cycle.append(graph[cycle[-1]][0])
        if cycle[0]==cycle[-1]:
            return cycle
    return None

def non_branching_path(edge, degree, graph):
    '''
    return the non-branching path starting with edge edge
    '''
    branch = edge[:]
    #print(branch[-1])
    while is_one_in_one_out(branch[-1], degree):     
        branch.append(graph[branch[-1]][0])
    return branch
    
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

adjacency_list = [
    "1 -> 2",
    "2 -> 3",
    "3 -> 4,5",
    "6 -> 7",
    "7 -> 6"
]

print(maximal_non_branching_paths(adjacency_list))
"""
sample output:
1 -> 2 -> 3
3 -> 4
3 -> 5
7 -> 6 -> 7
"""