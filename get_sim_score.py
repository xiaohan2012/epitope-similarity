import os

from similarity import FPWithComplex, similarity_between

#from concat import already_concated, concat_files, concated_file_path

if __name__ == "__main__":
    """
    Usage:

    python get_sim_score.py --query-pdb data/sample1.pdb  --against-pdb data/sample2.pdb  --spin-image-radius-range=20,40 --spin-image-radius-step=2 --spin-image-height-range=10,60 --spin-image-height-step=5 --sphere-radius-range=20,40 --sphere-radius-step=2
    
    python get_sim_score.py --query-pdb data/sample1.pdb  --query-epitope 211,213,214,224,225,226,227,228,229 --against-pdb data/sample2.pdb --against-epitope 216,217,218,219,220,221

    python get_sim_score.py --query-pdb ../web/tests/data/sample3_A.pdb:../web/tests/data/sample3_B.pdb  --query-epitope 1,2,3,4,5,6,7,8,9,10 --against-pdb ../web/tests/data/sample3_A.pdb:../web/tests/data/sample3_B.pdb --against-epitope 1,2,3,4,5,6,7,8,9,10
    """
    import sys, getopt

    #parse the cmd argument
    optlist, args = getopt.getopt (sys.argv[1:], "", ['query-pdb=', 'against-pdb=', 'query-epitope=', 'against-epitope=', 'spin-image-radius-range=', 'spin-image-radius-step=', 'spin-image-height-range=', 'spin-image-height-step=', 'sphere-radius-range=', 'sphere-radius-step=',])
    
    #and make them into the right data type
    spin_image_radius_range = (0, 20)
    spin_image_radius_step = 2
    spin_image_height_range =  (-30, 10)
    spin_image_height_step = 5
    sphere_radius_range = (0, 20)
    sphere_radius_step = 2

    query_epitope = []
    against_epitope = []

    while len(optlist) > 0:
        opt, val = optlist.pop ()
        if opt == '--spin-image-radius-range':
            spin_image_radius_range = map (float,  val.split (','))
        elif opt == '--spin-image-radius-step':
            spin_image_radius_step = float (val)
        elif opt == '--spin-image-height-range':
            spin_image_height_range = map (float,  val.split (','))
        elif opt == '--spin-image-height-step':
            spin_image_height_step = float (val)
        elif opt == '--sphere-radius-range':
            sphere_radius_range = map (float,  val.split (','))
        elif opt == '--sphere-radius-step':
            sphere_radius_step = float (val)
        elif opt == '--query-pdb':
            paths = val.split (':')
            if len (paths) == 1:
                query_pdb_path = val
            elif len (paths) > 1:
                if not already_concated (paths):
                    concat_files (paths)
                query_pdb_path = concated_file_path (paths)
        elif opt == '--query-fp':
            query_fp_path = val
        elif opt == '--against-pdb':
            paths = val.split (':')
            if len (paths) == 1:
                against_pdb_path = val
            elif len (paths) > 1:
                if not already_concated (paths):
                    concat_files (paths)
                against_pdb_path = concated_file_path (paths)

        elif opt == '--against-fp':
            against_fp_path = val
        elif opt == '--query-epitope':
            query_epitope = map (int, val.split (','))
        elif opt == '--against-epitope':
            against_epitope = map (int, val.split (','))
        else:
            raise Exception ('Invalid option')

    #calculate the finger print
    from get_fp import Complex
    from PDB.PDBParser import PDBParser
    p = PDBParser(PERMISSIVE=1)

    query_struct = p.get_structure(os.path.basename (query_pdb_path), query_pdb_path)
    against_struct = p.get_structure(os.path.basename (against_pdb_path), against_pdb_path)

    query_complex = Complex (query_struct, query_epitope)
    against_complex = Complex (against_struct, against_epitope)
    
    query_complex.get_fp ()
    against_complex.get_fp ()
    
    query_fp_string = query_complex.fp2str ()
    against_fp_string = against_complex.fp2str ()

    query = FPWithComplex (query_complex, query_fp_string)
    against = FPWithComplex (against_complex, against_fp_string)
    
    score1, score2, score3 = similarity_between (query, against)
    
    result = {
        'status': 0,
        'score1': score1, 
        'score2': score2, 
        'score3': score3
    }

    from simplejson import dumps
    print dumps (result)    
