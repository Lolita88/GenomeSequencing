#Code takes 2 sets of kmer pairs and lines them up to match patters. 
#Does assume that they are already in order...
#Input: Integers k and d followed by a sequence of (k, d)-mers (a1|b1), … , (an|bn) such that Suffix(ai|bi) = Prefix(ai+1|bi+1) for 1 ≤ i ≤ n-1.
#Output: A string Text of length k + d + k + n - 1 such that the i-th (k, d)-mer in Text is equal to (ai|bi)  for 1 ≤ i ≤ n (if such a string exists).

"""
Pseudocode
first_patterns = the sequence of initial k-mers from GappedPatterns
        second_patterns = the sequence of terminal k-mers from GappedPatterns
        PrefixString ← StringSpelledByPatterns(FirstPatterns, k)
        SuffixString ← StringSpelledByPatterns(SecondPatterns, k)
        for i = k + d + 1 to |PrefixString|
            if the i-th symbol in PrefixString does not equal the (i - k - d)-th symbol in SuffixString
                return "there is no string spelled by the gapped patterns"
        return PrefixString concatenated with the last k + d symbols of SuffixString
"""

def string_spelled_by_gapped_patterns(k, d, kmer_pairs):
    first_patterns = []
    second_patterns = []
    for i in range(len(kmer_pairs)):
        x,y = kmer_pairs[i].split("|")
        first_patterns.append(x)
        second_patterns.append(y)
    prefix_string = string_spelled_by_patterns(first_patterns, k)
    suffix_string = string_spelled_by_patterns(second_patterns, k)
    
    #The solution actually is to generate De Brujin graph for paired edges first and find whether it is connected or not. If not connected then find start and finish nodes. If connected use random nodes to simulate start and finish. Then use these paired start/finish when constructing prefix and suffix strings.
    offset = k + d
    
    for i in range(offset,len(prefix_string)-1):
        
        if(prefix_string[i] != suffix_string[i-offset]):
            return "there is no string spelled by the gapped patterns"
    return prefix_string + suffix_string[-(k + d):]

def string_spelled_by_patterns(patterns, k):
    return string_reconstruction(patterns)
    
def string_reconstruction(kmers):
    text = ""
    for i in range(len(kmers)):
        text = text[:i] + kmers[i]
    return text

def hamming_distance(p, q):
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1
    return count  
"""
#simple test data
k = 4
d = 2
kmer_pairs = [
"GACC|GCGC",
"ACCG|CGCC",
"CCGA|GCCG",
"CGAG|CCGG",
"GAGC|CGGA",
]
print(string_spelled_by_gapped_patterns(k, d, kmer_pairs))
"""
#raw_data = sys.stdin.read().splitlines()
#raw_data = open('files/open('files/eulerian_cycle.txt', 'r')

raw_data = [line.strip() for line in open("files/kmer_pairs_dataset.txt")]
if(raw_data[0] == "Input"):
    del raw_data[0]
k, d = map(int, raw_data[0].split(" "))
del raw_data[0] # deletes k,d before import of kmer pairs

kmer_pairs = []
for line in raw_data:
    if(line == "Output"):
        break
    else:
        #kmer_pair = line.split("|")
        kmer_pairs.append(line)

#test code if you have the answer to compare
my_answer = string_spelled_by_gapped_patterns(k, d, kmer_pairs)

print(my_answer)
"""
correct_answer = "CCAATTGTTGGCAACAAAGAATCGCTTATGCTAGGGTGACGTGCCAATCGACTGATTTGACTGGCCGGGGGATCGGCTGCGTAAAACCGGTGTCAGAATAAATAGTCATGGCCGGCGTCGACAGGCGCCCCGAGGGATAGGTAACGGGCGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTCGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTGTCAACTACGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTGAGCCTGAGGCCCGTGAAGAAGCCGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTGGTTCTGGGTGCATAGCCGCGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTGACGCCACGAAGTGTCAACTAGTGTTGTCATGAGAGAGTTATTATAGCAGGCCTACTTGTAGGTAAATACACTCTAGGTTATTCGCTCTGCTCCCCTCCTGCGTAACCCCTACCGTGAAGAAGCGGTCGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTGTGTTACTACCCATAGCGTCGGCCTCGTGAAGAAGCGGTTCTGGGCGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTTGCATAGCCGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGACGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTAGTGTCAACTCGGACGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTCGCCACGAAGTGTCAACTACGTGGCAATCATCGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTGTACTAGTTTAGCTGTAGGGCTTGAGGCAATTCCACGATCAGCGGGAACAGCGATATAACCCTTACATATCTAAACGCTGGACTGCATAAAGTAAGCAAGGAAATTGACTGAGGCGCTTACCCCGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTCCAGTATCAAGCCGCAACCGGGCCCGTGACTCATCCTCCTGCATACCCGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTGAACGGGGCCTGGTCCCGTTTTCGAAGGGTGAGTTCTGCTTAGCGTTGTCTTTCATTCGCTCAAAAGTCCCGCGTAAGAGCATCCTGGATTGTTCGCCCTGTAAGCGGGACTACGCGTGCCGATGGTGGGCTTGCAATTATCATAGCGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTTCCTGTTCCGTCAATTCCTCTCTAAATACTATCTAACCTGGTCGCAGAACTCGAAGAACTACCGGCCGTCAGCAATTCTAGCTTAATACCTCGTCGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTTGAATAGTGCGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTGCCCCTCGGAACGGTATGTACTGCAAGCGTAGAAACCCTGATAGCTTGGATGACGAAACTGTTAGATGTACTGCCAACGGTTAGTCGCGCTGTCGGTTTCGTTAACGATGCATTAAGTCGAACTCGTACCTAGAAACGTGAAGAAGCGGTTCTGGGTGCATAGCCGGACGCCACGAAGTGTCAACTAGTGGGATATTGGTGAAGCAGAGGACGAATTGCGATATCCAAGATGAGAACTGTTTGTCAGTCGGGGAAGACCCAGCTGACTACGCTCAGAGCCCGGTCATGTGTCTGAATCAATCTAAAAACGTATAGTTTGGCTACTGGGGCGCTAGGTGC"
print("len of my ans: " + str(len(my_answer)))
print("len of correct ans: " + str(len(correct_answer)))
print ('hamming distance:', hamming_distance(my_answer, correct_answer))
"""