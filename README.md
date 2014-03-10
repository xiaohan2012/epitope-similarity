epitope-similarity
==================

Epitope similarity calculation between two structures

#Before you run

```
> cd epitope-similarity
> virtualenv venv
> source venv/bin/activate
> pip install -r requirements.txt
```

#Check if everything is ok

`> python test/similarity_test.py`

	
#How to use

1. Without epitope specification: `python get_sim_score.py --query-pdb data/sample1.pdb  --against-pdb data/sample2.pdb`
2. With specification: `python get_sim_score.py --query-pdb data/sample1.pdb  --query-epitope 211,213,214,224,225,226,227,228,229 --against-pdb data/sample2.pdb --against-epitope 216,217,218,219,220,221`
3. With spinimage configuration: `python get_sim_score.py --query-pdb test/data/sample1.pdb --against-pdb test/data/sample2.pdb --spin-image-radius-range=-10,10 --spin-image-radius-step=2 --spin-image-height-range=-20,20 --spin-image-height-step=5 --sphere-radius-range=-20,0 --sphere-radius-step=2`
