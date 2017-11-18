# Function takes in a collection of kmers, compares prefixes/suffixes for overlapping
# matches. Where matches are found, full kmers are paired and output. Sort makes it more
# efficient and cleaner, lexicographically. Plus, order doesn't matter. 
def overlap_graph(kmers):
    kmers.sort()
    patterns = []
    prefix_patterns = {}
    suffix_patterns = {}
    for i in range(len(kmers)):
        prefix_patterns[kmers[i][0:-1]] = kmers[i]
        suffix_patterns[kmers[i][1:]] = kmers[i]
    for pre in prefix_patterns:
        for suf in suffix_patterns:
            if pre == suf:
                patterns.append(suffix_patterns[suf] + " -> " + prefix_patterns[pre])
    return patterns

"""
kmers = [
"ATGCG",
"GCATG",
"CATGC",
"AGGCA",
"GGCAT"
]
""" 
kmers = [line.strip() for line in open("files/dataset_198_10.txt")]
#print(overlap_graph(kmers))
print('\n'.join(overlap_graph(kmers)))