from __future__ import division

import os

from similarity import FPWithComplex, similarity_between

def help ():
    print   """
    Usage:

    Without any specification: python get_sim_score.py --query-pdb test/data/sample1.pdb  --against-pdb test/data/sample2.pdb 
    
    With specification: python get_sim_score.py --query-pdb test/data/sample1.pdb  --query-epitope 211,213,214,224,225,226,227,228,229 --against-pdb test/data/sample2.pdb --against-epitope 216,217,218,219,220,221

    With spinimage configuration: python get_sim_score.py --query-pdb test/data/sample1.pdb --against-pdb test/data/sample2.pdb  --spin-image-radius-step=2 --spin-image-height-step=5 --sphere-radius-step=2

    You can specify the fingerprint path and skip the fp calcualtion step by using the option --query-fp and --against_fp, for example: 
    
    python get_sim_score.py --query-pdb test/data/sample1.pdb --against-pdb test/data/sample2.pdb --query-fp=fp/h2-v2/sample1.pdb.fp --against-fp=fp/h2-v2/sample2.pdb.fp
    """

if __name__ == "__main__":
    import sys, getopt

    try:
        #parse the cmd argument
        optlist, args = getopt.getopt (sys.argv[1:], "", ['query-pdb=', 'against-pdb=', 'query-fp=', 'against-fp=', 'query-epitope=', 'against-epitope=', 'spin-image-radius-step=', 'spin-image-height-step=', 'sphere-radius-step=', 'cutoff='])
    except:
        help ()
        sys.exit (-1)

    #and make them into the right data type
    spin_image_radius_range = (0, 20)
    spin_image_height_range =  (-30, 10)
    sphere_radius_range = (0, 20)
    
    spin_image_height_step = 5
    spin_image_radius_step = 2
    sphere_radius_step = 2

    cutoff = 20.0
    
    query_epitope = []
    against_epitope = []

    query_fp_path, against_fp_path = None, None

    while len(optlist) > 0:
        opt, val = optlist.pop ()
        if opt == '--spin-image-radius-step':
            spin_image_radius_step = float (val)
        elif opt == '--spin-image-height-step':
            spin_image_height_step = float (val)
        elif opt == '--sphere-radius-step':
            sphere_radius_step = float (val)
        elif opt == '--cutoff':
            cutoff = float (val)
        elif opt == '--query-pdb':
            query_pdb_path = val
        elif opt == '--against-pdb':
            against_pdb_path = val
        elif opt == '--query-fp':
            query_fp_path = val
        elif opt == '--against-fp':
            against_fp_path = val
        elif opt == '--query-epitope':
            query_epitope = map (int, val.split (','))
        elif opt == '--against-epitope':
            against_epitope = map (int, val.split (','))
        else:
            raise Exception ("Invalid option")

    #calculate the finger print
    from get_fp import Complex
    from Bio.PDB.PDBParser import PDBParser
    p = PDBParser(PERMISSIVE=1)

    query_struct = p.get_structure(os.path.basename (query_pdb_path), query_pdb_path)
    against_struct = p.get_structure(os.path.basename (against_pdb_path), against_pdb_path)

    query_complex = Complex (query_struct, query_epitope)
    against_complex = Complex (against_struct, against_epitope)
    
    if query_fp_path is None or  against_fp_path is None:#if fp is not given
        query_complex.get_fp(spin_image_radius_step = spin_image_radius_step, spin_image_height_step = spin_image_height_step, sphere_radius_step = sphere_radius_step)
        against_complex.get_fp(spin_image_radius_step = spin_image_radius_step, spin_image_height_step = spin_image_height_step, sphere_radius_step = sphere_radius_step)
        
        query_fp_string = query_complex.fp2str ()
        against_fp_string = against_complex.fp2str ()
    else:
        #if fp is given, read them
        with open (query_fp_path, 'r') as f1, open(against_fp_path, 'r') as f2:
            query_fp_string = f1.read ()
            against_fp_string = f2.read ()
        
    query = FPWithComplex (query_complex, query_fp_string)
    against = FPWithComplex (against_complex, against_fp_string)
    
    score1, score2, score3 = similarity_between (query, against, cutoff = cutoff)
    #z1, z2, z3 = similarity_between (query, query, cutoff = cutoff) #the normalization constant
    #print score1, score2, score3
    print score1, score2, score3
