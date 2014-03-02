epitope-similarity
==================

Epitope similarity calculation between two structures

#How to use

1. Without epitope specification: `python get_sim_score.py --query-pdb data/sample1.pdb  --against-pdb data/sample2.pdb`
2. With specification: `python get_sim_score.py --query-pdb data/sample1.pdb  --query-epitope 211,213,214,224,225,226,227,228,229 --against-pdb data/sample2.pdb --against-epitope 216,217,218,219,220,221`
3. With spinimage configuration: `python get_sim_score.py --query-pdb data/sample1.pdb  --against-pdb data/sample2.pdb  --spin-image-radius-range=20,40 --spin-image-radius-step=2 --spin-image-height-range=10,60 --spin-image-height-step=5 --sphere-radius-range=20,40 --sphere-radius-step=2`
