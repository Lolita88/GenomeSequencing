"""
Input: The adjacency list of an Eulerian directed graph.
Output: An Eulerian cycle in this graph.
"""

from copy import deepcopy
from random import randint
import random, sys


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
            pick = randint(0,len(red_adj_list[curr_vert])-1) 
            temp = deepcopy(curr_vert) 
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
"""
#Sample Input:
my_list = [
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
"""
#my_list = open('files/eulerian_cycle.txt', 'r')

data = []
for line in sys.stdin:
    data.append(line)
my_list = data

print(find_eulerian_cycle(my_list))

"""
Sample Output:
2->1->0->3->2->6->8->7->9->6->5->4->2
"""