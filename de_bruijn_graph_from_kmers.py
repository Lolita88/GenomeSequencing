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
"""
kmers = [
    "GAGG",
    "CAGG",
    "GGGG",
    "GGGA",
    "CAGG",
    "AGGG",
    "GGAG"
]
"""
kmers = [line.strip() for line in open("files/dataset_200_8.txt")]
#print(kmers)

print('\n'.join(de_bruijn_graph_from_kmers(kmers)))

