"""
Protein Translation Problem: Translate an RNA string into an amino acid string.
Input: An RNA string Pattern and the array GeneticCode.
Output: The translation of Pattern into an amino acid string Peptide.
"""

def protein_translation_problem(rna, rna_codon_table_array):
    rna_dict = {}
    pattern = ""
    for each in rna_codon_table_array:
        try: #need to handle exceptions when there are "stop codons"
        #  that do not translate into amino acids and serve to halt translation
            k,v = each.split(" ")
        except ValueError: 
            #if v is a stop codon, don't save it
            v = ""
            #print("stop") 
        if v != "":
            # v is an AA abbreciaton, so save
            rna_dict[k] = v
       
    for i in range(0,len(rna), 3):
        # loop through the rna string 3 at a time. Only working w one reading frame.
        #print(rna[i:i+3])
        curr_codon = rna[i:i+3]
        if curr_codon in rna_dict:
            # add codon's single letter version to string pattern from dict look up
            pattern = pattern + rna_dict[curr_codon]
    return pattern       

rna_codon_table_array = [line.strip() for line in open("files/rna_codon_table_array.txt")]
rna = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
#rna = "AUGGCCCAGUCAAGAAGGGCUGAGGUUGCGGCGAAGGCCUGGCAGAGAGAAAAGCCAACGCGAAGGGUGGGCGUGGCACUCUUGGAUUCAAUAAAACGGCAAGCCCUCCUCCCGGACGUUGAGGGUGUAGAACAAUUGGUAAAGGAAAGAAGUGGACGGGGUGUCGCACCGUCAAACAUCGCCGCGAUCUGUCGCACUGACGUGUCACGUGUGUCGACUCCGUGCCAAGUGGCAUGUCAACAUUAUGUAGAACCCGUGACGUUUCAUUUCUCGGUCAUCCGGGCGGCGGCGCCGGACGAUAUUCUGGUAGCGCCUUCGCGUGCAUUGCUCGCAAAGGUUUGGAGUGCGUGCAUGCAUAUGUUGGAACAACAGUGCUGCCACUGGCGUUUAAAGAGUGAAAAUACAGAUCCUUUCACGACCUUGCAAAGUCGUACAUGGCUACCCUGGAUAAGAGAGGACAAAGCUCUGAAGAUAAGUACUCGUCAUAUCAUUCGGCCUACUAAGAUAUGUAAGGGGAACGGUCCCCGUAUAUCUUGUCCAUACGACUGCGGAAACCCAACGAACAUUCGUGAUAUGUAUCGGUCCUCGGGCUGGGUCGAUGGAAACUCAGAUGUGAGCCCCCUGCUUCAGUCACAUAUCAGGACGCUCAAAAAAAGAACGGUUGAUCAAAGCGUAACUCAUAAAACACCGGCCGGAGACUAUCGCUCAAUGAGGACCCACCUGUUCUGCGUCAAGAUUUGCAACCAUCUGUUAGAACUACUCGUAAUGACUAGGUUCCCACAUGUAAAUGCAAUGGGCUCAAUUAGUGGCACUAUUGACGUCGGUGCGCUGCACCCUACCACCCUACUACGGACGAAGAUGGAUGUGUACGGUUACAUCCUGCCUACGCUUCGGGCCAGCUCAGCUCUGUGUCCGUGUGGCGACAGGAAAGAGCUUAAUGUUCAUUUCCUUAGGACUAAUCACGCCCGAUGUAACUCACCAAACCGCCAGGGACUGACUUCGACAGCUUGCGUGUCGCUCCCCCUUAUUAUAACCUCCAGCCGAAGCUCUUGUGGAGGAGGGCGAUACGUAGUUCCGGGAUUUGCGCGCAGCGCGAACACAUUCCGCCUGAGUUUUAAGAUAGAUCCUCGGGUGUACCUUGUAGUAGACUAUCGCGUAACUCUACGCGAGGGCGGACCAUGCGAGGGCUUCGACAACGACGACUUGGCUCGCCGACAUCUAGGGGCACUAGCGACUGCCUACUUUGACAGUGGAAACCCUCAUCAGUAUGCUCAGUCGCAUGAGUUUCUACUGGUAUUUCAAUUGAUAAUGCAUCAUCAUUUUAAGCAUCUCCACGGCUUAUGGUCCGCGCAGCAUAAAACUUCCCGUUUGGACCCCAUGUCUGUACAGUGCCCCACCAAGUUUCGGGCUAAGGUUGGAGGCGACAAUCGGGGCCCCCUAUUUCCUCUGUCUGAUUUGAAUAUGAACGAGACUGCGCCAAACAGUGUCCGGUGGCUGUUAUCGAACUGGGUGAAUGCAUGGAUGACCAGUAACGGCGAGCAUACAUCUGCAGUAGAACUCUUCGGACAGUGCUUAACAGGCCCGGGCCAAGUGCAUGACGAGGGACCUUUAAAUCGGAGUCGAACUUCCCCCUACCCGGGCACCCACUAUUGUGGAGAGGAGGAUUAUUUGGCAGGUAGAUAUAAUCAGAGCCAUAAUUCCUCUGUACGGCGCUGCGGGGCAACGGGCGAGAAGCAUUUGCUUGUUAACAUGACGUCGACUUAUGGCACACACGGUUGCCUCCAUCAACUGACUCCAACGGACAAUUCAGACCUGAUAUACAACGCGCCUGAAAAGGGAGAGGCCCCCGACAAGGCCACUUCAGCAGAUUUAAGAUAUAACAUGUCCGCGCCGCAACCAACGGUUUUCCACAUCAGCCCAAGCGGGCAUGUUUAUUCGUCGACGGGGGCUUAUCUAGAUGGACCUAGGCGCUAUUUCUCGCGCUUUAUCCGGGGAGACGCAAGAUGUAUGCACUGGAUCGAAGGCAAAUAUUAUAUGUAUGGUUUGUUCGCACGGAUAAACACAGAUAACUGUUGCUGCCGUUCCUACUGUCGUGCACUAGCGCAUAUUAAUAGCGUUCGCGGAGUUGGUAACGUACCCCACGUCACGAAAUGUCGUCUAGGCUGUCCAACGGCACGGGAGACGCACGCAGUAAUAGCGACGGACUUACUGGGUGGCAAUCUUUUUGCUGUGUGCAUAGAACUAAACAGUGAUAUGCCCACCCAGAGCUUCCACCUAAGCGAGGCUGCAGAUCGAUCUAACCUCUCACUCUUUUUAAAUAGCGUAUUGUAUUUGCAUGAGCACACACUUAGUGCUCAAGGCUAUGGACAUGUUACCGACAGACAUAGGCGAGUAUAUGGACGCUACGGGGCGAGCCCCUUGUGGCCUGGAGGGGAGUCGUUUCAUAACCGUGGUAGGCUAAUGUGUUCCGCUGGCGCUGAUACAAACGCCACGCGGGCAGCAGGCCUUCAAUACCCUGCGAUGCGACUUCGGGGAUAUUUUGGUCGUCAAUCUAAAUAUAUAUCGACUUGCUCCAGGGUAAAACAAUACUCAGGCGGCAUCGUAGAUCUAUUUAACGCUAACUGCUUCGGAAAAAGGGACUGUGAUAGGAAUCUCAGUUCAAUUGACUGCAAAAUAAGAACAAUUCGUGCCGCGUACGACUCGUGGAGAUUUAGACCCGGCGAGGAACCGGAUUCCCCUCGCCUCUCUGGACCACCCCGUAAACAUUCCUUGUCAUUGGACUUCAUGUCACCUCGUAUUCAAGUGCCGUACAGGAAUUCUUCCUGUCGAUGGCAUGAUAGAUCCCUUCGCACUCGACUUCGGACCAGGAGGAUCCUCGCAAACGAGGCAGAGAGCAGGUCACUCCACAGCCUAGACAAACACGGUAUAUGGGGCGAAAAGUCACAGGGAGACGUCCAGUCGUAUCAGAUGAUUUGUAGCUUCGUGGAACGGGGGAGCCUUUCGAAGAUCUACAGUUCCGUCCCUCGGUCUAUCCAAAGUUUACCUAAAAGAACAAGGAGUGGUGUGGACGGGCUCAUGGCACGUUGUCAUUCGUCGAAUAAGCCAACAUUGCUUAUCAGGACUGAACUAACGAGGGACGGAGAACGCAUCGAGAAGUAUAACGCACAUUGCGACCCUCGUGAUUUCAAUCGCAUCUUCCACGCUAUACGUGAUAGACUCAUCACGGCUACCGUUUCAAUAGGGCUGGGUGCGAAGAAGCCUGCGGGAGUGUAUGGUCACGGAUUAGGGGACGCUGGUGCGAUGCUGGACUUGCACAACUCGGGAAUUGACUUGAGAUUGCGAUUACAGGACCCAACGGGGGCAGAAUACCCGCCCGGCCAGGCAGGGGCAGGCCUAGUUCAUCUGACAACAUCUGUGCCGCUCUGCAUCCCCCGUCCAGCCUCGGCCUACAUGAAGUGGCGGCUUCCUAAGCUUCAUGCUCGUCUUGAAGCUGGCGCUCCAGUGCUUGAUCUAGCAAGAAGCCCUCCCAAACAUGCUUGGUCUUCCGGCCUGCACCUUCUAGUCAACUUACUUGCUGUCAGAAGGAGCUUGAACUUCACUACGUGGGGGGUUGUUCAGAUUUUACAAAGGCCAGAUGCUCCGGUCUUAAUCGAACUACCAGGCGCCGAUAGACCAACAACAACUUUAAACCGUAUCAGCAAGCUCAGUUCUAAGCCAUUCGAGUACCCCUGCUUCAUAUUGAUCCGGGUGGAGCAUGAGGCGCUGGGCUGCAGCCAUCUGUGCACGUUGCUUGUUCAUAUGCCUGGACCGCCCAUGCGAUCAAAACUACUAGUUAGGCGUCUAUAUAUAGACUCCGAGCUGAUGAACCAACGAAUAGUCGCGGGUAGUAUUUUUUCAUAUGGCAUCAGACUUCGUUACAGAAAUCCCGAUAUACAUCCCAGCGUUGCUGAUGUGGUACUCAGGCCGGACGACCCGGCAUUAAUGCCCCACCGGCAACCGGUCGGACAAAGGUGCUCGUCGGUUGACAGUCUGAAUACCCCGGUUUUGCGAUCGCGGCGCGUAGAGGACAAACUCGAAGGUGAUCCGGAUGAAGGGACGCCUUGCUCGACGUGGGAACGGGAACCAACGUUACAGACCGACGAGCCCCGGAACAGUAGAAGUGCGUUGAAACUUCGUGCACGUCCAAGUAACGUACCGAUGCGUGGGCAUUCUAAUCCAGUAUUAUUAACUCUUCUAGGCACAUUUUUUUCCGACCUCGUAUUGCGAGGCCUUCGAUCCUCCGCCCUAGCGUCCGUAAAGGAGGAUCAGCAUCCCUGUGACGUUAGUAUCGGAGUCACGUCGCUCUGUGGCCUUUCGACCCUUUGCUUUCUCACACGAACUGAAUUUUCCUAUUUCCUAUGUAGGGUUUUGGUCCCACCCCCACUAGGCUCAGUAUUCCUUACCGAGUCAACCUGCCUAGGAGAAUCUUACUUGUCAGUCAUUCGACCAUCUGCGCCAGACGCAUCAUUGGUACGCGAGCUAGAUUCGAUCCUAGUUGGAGACUCUGGCAAGCUAUCCCUUCCAACAGGCUCGCGCAAUUUACCAGAUAGUAAGGGAUGUGUACAGCAUAUCGUACUGCUCUGCGACUGCGUAUUGAGCGAAACUCCUGAUGGCAAGUUAAAUUCGAUGACUGAGGAGAGCGGACCGCCCACCUCGCGCAAAAGACAAUCCUGCCACUUAACAAGUUCCCAACUGGAAGACCGUCGGGUUAACGGUCGCGAUGGAAAGUACAACGCCCGAACGAGUUCCCACUACUAUAAGAGUCGCAUAAAGUAUCCCUCAUGGGGUUUACUACCCUGUACUACGGGCCGACGCCUCAUAAUCCUUGACGCUCGGAAGUCGGCCCUGCUCCGCGGCAUAUUGACAUGUGAAGUACUGACACCACCGUCCGUCCGGCCGAUGAACCCGUUAACUCGUCUCAGUCCAGACAUUCUGCGACACGGUCGGUCUAAGCGUAUGGGGCUCUACCUACUGAGACAGGCAUCCAUUAUUUUCGACGGUACUUCGAUCCCAUUAACGCUACUCCACUCGCGUUGUCUAAAAACGCAGCAUCUCGAAAAACCCGUAGCUGAGCAAUGUAGUUGUCCUGCUAGCUCUUUGGGAACCCUGAGCAUAUUUCCCUGCACUAUAGUAACACCUCGGCUGGACAAAAUAUGGGAGCUGCUAAGCUUAUCGUCAAGCUCCAUUCAAAUCAACCACAUGGAUACACUCAGUCCACGGCUACACGAGACAUUGGUCGGACCGAUUGUGGAUCGAAAUUGUCUUCUGAUAUCAGUACGCUACGAUCAAUGCCGGAAGUGCAAAAUUCCUCCCCACAAUGCCCGGAAGGGCCUACCUGUAGCGUCGCGUUUGGCUUCGACGAACUCGAAGCUGGCUAAAUGUGAAUUUACAAUGUCUCACAAGAUCGCAAAAUUCCACGCCAGUCAGCUAUUGCGGACGUGGUUUAGAACCCGUUGUCCAAUGAAUGAAAUUGGACAUAUCAUAAGCUCUAGAAGCUUUUCGCAUGCUCGUCCCCCCCCUAUUUACUGCAGAUGCGCCUUCUCACGUCGAUCGGGGCCCAAGCAUGUGCAAGGUGUACUACAAGACCACUCUCACGGCUGUCUGUGUCCACACCCCGAUUUCGACUUCGAUGGCACUCCGAAAAGCUUGGAUGUAAACUUUGCUCAAAAUCCUGGUAGAGACUGCGGAACAUGCCAAGACACUCCCUACCUGCGUCAUGCUAGGUGCAUAAGACAGUUAACCUCGAAACGAAGACCAGUUAUCUACCUGGUCGGAUGGUCUAGAUUUUGGCAUGCGAGGCGAACCCGAAGUGACCCUUUGCAUACUGACUACCAGGGAAAACUAAGGAAUCAAGAACGGCACGUGCGACCUCCACUAGCGAUACAAUCCCGGAGGGCCAAUACGAGCAGCAGCAACGAUUCUAACCUAAUAGCAUCCCUUCCCUUGUCCCACAGCUCGACCUUUCACCGAUCAUCGUUAGGCGGUAAUAUAGCGGAGUCGCCACAUUACCCGUAUACGCCAAUCAUGAGUACCGCCAUCAAGGAUACGUAUCUCACUGUGCUCAAUAAAAUCAGUCCCUGUCAUCUUUGGGCGGGAGGGGCCCACACUCGAGUGCUGCGUGACACGAAGUGGUCAACAGUUUACGAAUGGUCAAUACGCCCAUGGCGCCGAGAAAUUGGACUAAGGGCCUGCCAUAGUUUCCCAGUUUCAAGUUCGUGUUGCAAAUCGCGUAAUUCUCAUAAUGUCAGGUGUACCUGUUUUAAGCGACUUACCGUGUGGACGUUUCGUCAGGCAGGUCCGGUCGACACCUUGCUCAAUGAUGUCUACAAUUUAGUGCCUGCAGGGGACAUAGAAUCGUACAGGCCUCAUAUUGUCGUCAUCCGCGCUGCUGUGCAGUGGACCAGAGUCCUGUGGGGGUGUAGGUGGACCAUACAACCCAUCACCUCAGCAAUCCGCUCCUCGCGGAUCUCGAUAGCGAGUGUGAGGCCGACACUUUACAUUUUCCAGGCCCUGGUUCACGUAUCCGUGAUAUCACGCCGCCGUAUGUGCCUAAUUUUUAUAGGACAUGCUGCUUGGAGCGUAACUAUCAGGGUAUCGCUUAUCUCCAAACCGUUACCAUGCGCAACCGUGCCUCUACUUUGGCUGAUCCACCUCGUGCAGAUUCCUGAAAAGCCCAACACCUUAGCAAGAAACCGCGUCAUCGAUACUACUCACCGGAUCGGACUUUUAAUAGUCAAACCUGGACAUAUCGAGCGCGCGGUUGGAGAUGCAUACUUGAAAAGGCGGACGUGCGAAGAUUGUCUUAGAGCUUGUCGGCAAUCAGCGUGGCUAGCUGUUGCAAACGAGCGGAAGUUACUAUGCCUCGAAUGGGCAUUGGAUAGAGGCGUUGAAACCGUUGUCACCACUCUCCCGGUCGAGUUUGGCGAUCAAGCUGAUGUAUUGGUCGAUAGCACGGUGCUUAACACUGGAGGAGACGCUACGAGUGAAAGGACUCGAGCGGACGCACGCGGCUUUCACAGGUACUUAAAAACGGUGCGAAAAAGCAGUCCCCCGCAUAGACGCACGGUGUUGAUGAAACGUAUAUUGUAUAAGGCCACUUCCGCUGUUGUCCCUGAACCUUCACAUCGAAAAGCGGAUAGACGUGCAACAGAGGAGCCCUCGGGACCAGAAGCCUUAUGCUCCGCCCAUACUCUGAUUUUCAAGUCCAUACAACUUAGUAAAAACCGGCAGAUAUUGCGCGUAGAGGAUCGCUGCCGCUACCGCAUGUGGGAUUUGAGCCUAAGCUUGCUAGCUGUCGGCAGAUUUGGAACAUGGACUGAGGAACAUGGUGCCAGCAUAAGUGGGCCCUUGGAUCACUGUAGACUCUGUACUAUGAAUCACGGCGGUCAACGUCGGACGACGUCAAGCGGAAGAGUUCGUGUAUGGAAGCGUGAUAGGCCCGAAAAAUACACCGAUAGCGUUAAACAAGUAGUUGGGGCGCGCAUUUACAACAGUAUGUCGUUUUUCCGCCAAGCCCGAAGUAGUAGCAUUCGCGGUAUCUCUGAUGUUCAAUACCGUUGCCUUCUCAUGAAAACUCCUACAAAUUUUGCCUCGAGAGAAGGGGCCCCGUUGGUCCCGCUAGUAACGUUGUACAGACGAGGGCCUGCCUCUGAAGCCUCUGCUCCCGUACGAUUACUCAAACGCGACGCGCGACACCUGCGUCAUACAUUUGCGCAGGACGCAUGUACCCCUAGGUCUGACACCUCGUAUGCGGUACGAGCCCCUUUAAUCGCAAACAACACGCUUGGCCCGUGUGCAGCAAUUCGGGGGCCUGGUCUUAUGGACGUAUAUGUAGAACAAAGUUACAGAGAUGAACCGGUCAGACAUACCCCGGUCCCAAACGCCUUGAGUACCCUAAGGUCUAGCAAACCGUACUGCUUAAGGCUGACAACGUUCACCGCCACCGACCGAAGAAGCGCAGCAGCCGUACAUAUGCCCGUCCGACGACCUUGGAAACGACGAAGCAGAGCUGUCCAUAAGAUGGUUAGUAGUGGUCCGGGACUACGUCCCCGGGAGCCAUUUCAGAGGGAGAGACCCCUUGGCAAAUUACUUAAUUCGCUUGCACAGAUAAAGCAGUCAACAACUAGAGACCGUAGGAUAGUAAAUAGUGAGUGUGUCGAUUUGACGAAGAAGGGGAUCGCCAAGUUUUGGCUGCGACUAGUUGAGAGUUCAGAUGGAUCUGGGGAACCAGUGAAUAUGAUUCACUGCGCGGAAUUGCCGAUCCUGGAACAAACGGAGGCAAGAGUAUGGGAGAGCCGACUGAAGCUGUUAGUGCACGUCUUCACGCAUAUCUCGUGUGCUCCUGACUAUAUCCCCAUGCGAAGCCCUUACUUUCGAAAAGACAAGGUCGCGAGUAGUCUCGGCCUUCGCAUGAGUCAUGCACUGCACAGAUUAGCGUGUUUAUCUCGUGCGUCAUCUCCAUCAAACCUCAAUGACCGGUCGCGAUUUUUCUGUUUUUUCGCACUCUGUACUGUGAGGGUUCCGAGCUCUAAUGCGCGUUUACGUAAAUCGACGCUCCAGCGUCAGCCCUAUCGUGUUUGCUCCAGGGCUCUAAUUAGAUACCUCCUAUCUCCAGAUUUCAAGGGGCAUGAUAUCAACCGCACCGUACACUUUAGAAUGUGUCGCCGUGUUAACGUUACGACUCCGGAUGAACUUACACUCGCCUUACGUGGAGGCAUCGAUGGCAGCUCGCUAAGUUGUUCGUCCAUCUCCGAGAGGGUGCGGUCCGAACCCACGGUACUUGAAGCCUUAUACGUCGUAACUGACGAUAUGAUAUAG"
#print(rna_codon_table_array)
#Sample output: MAMAPRTEINSTRING

print(protein_translation_problem(rna, rna_codon_table_array))

