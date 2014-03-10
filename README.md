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

1. Without any specification: `python get_sim_score.py --query-pdb test/data/sample1.pdb  --against-pdb test/data/sample2.pdb`
2. With epitope specification: `python get_sim_score.py --query-pdb test/data/sample1.pdb  --query-epitope 211,213,214,224,225,226,227,228,229 --against-pdb test/data/sample2.pdb --against-epitope 216,217,218,219,220,221`
3. With spinimage parametrization: `python get_sim_score.py --query-pdb test/data/sample1.pdb --against-pdb test/data/sample2.pdb  --spin-image-radius-step=2 --spin-image-height-step=5 --sphere-radius-step=2`
4. With cutoff parametrization: `python get_sim_score.py --query-pdb test/data/sample1.pdb --against-pdb test/data/sample2.pdb  --cutoff=10`

You can specify the fingerprint path and skip the fp calcualtion step by using the option --query-fp and --against_fp, for example: 

`python get_sim_score.py --query-pdb test/data/sample1.pdb --against-pdb test/data/sample2.pdb --query-fp=fp/h2-v2/sample1.pdb.fp --against-fp=fp/h2-v2/sample2.pdb.fp`

