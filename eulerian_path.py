from copy import deepcopy
from random import randint

def create_adjacency_list(list):
    #saving list into a dictionary
    adj_list = {}
    circuit_max = 0
    for line in list:
        node = line.strip('\n')
        node = node.replace(' -> ', ' ')
        node = node.split(' ')
        adj_list.setdefault(node[0],[]) #adj list gets the start of the node pair(ie. first num)
        for num in node[1].split(','): #num gets assigned the end node of the pair. Split on comma needed when multiple end nodes
            adj_list[node[0]].append(num)
            circuit_max += 1
    #print(adj_list)
    #print(circuit_max)
    return adj_list, circuit_max

def find_eulerian_cycle(list):
    #create adj list
    adj_list, circuit_max = create_adjacency_list(list)

    #reduced adj list to keep track of traveled edges
    red_adj_list = {}
    red_adj_list = deepcopy(adj_list) #exact copy of dict

    #arbitrary starting point (if graph is directed/balanced)
    #pick numbers that have multiple verts if poss, imo
    start = '6'
    curr_vert = '6'
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

list = [
    "0 -> 3",
     "1 -> 0",
     "2 -> 1,6",
     "3 -> 2",
     "4 -> 2",
     "5 -> 4",
     "6 -> 5,8",
     "7 -> 9",
     "8 -> 7",
     "9 -> 6"]

#list = open('files/eulerian_cycle.txt', 'r')

print(find_eulerian_cycle(list))