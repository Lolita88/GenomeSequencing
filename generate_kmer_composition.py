def generate_kmer_composition(k, d, genome):
    basepairs = []
    offset = 2 * k + d - 1
    for i in range(len(genome) - offset):
        kmer1 = genome[i:i+3]
        kmer2 = genome[i + 2 + k: i + 2 + 2*k]
        print(kmer1 + " " + kmer2)
        basepairs.append(kmer1 + "|" + kmer2)
    basepairs.sort()
    return basepairs

print (generate_kmer_composition(3, 2, "TAATGCCATGGGATGTT"))
