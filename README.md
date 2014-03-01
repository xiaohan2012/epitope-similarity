epitope-similarity
==================

Epitope similarity between two structures


#How to use

from epitope_similarity import similarity
from epitope_similarity.util.load_pdb import load_pdb_struct

p1 = load_pdb_struct ('path/to/struct_1.pdb') #load the first structure
p2 = load_pdb_struct ('path/to/struct_2.pdb') #load the second structure

score1, score2, score3 = similarity (p1, p2) #calculate the similarity score