#import sys
def string_reconstruction(kmers):
    text = ""
    for i in range(len(kmers)):
        text = text[:i] + kmers[i]
    return text

"""
kmers = [
    "ACCGA",
    "CCGAA",
    "CGAAG",
    "GAAGC",
    "AAGCT"
]
"""

data = [line.strip() for line in open("files/dataset_198_3.txt")]

#print(string_reconstruction(kmers))
print(string_reconstruction(data))
        
        
