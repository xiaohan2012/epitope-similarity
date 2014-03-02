from setting import *

import unittest, os
from get_fp import Complex
from PDB.PDBParser import PDBParser
from similarity import FPWithComplex, similarity_between


class SimilarityTest (unittest.TestCase):
    def test_basic (self):
        path1 = DIRNAME + '/data/sample1.pdb'
        path2 = DIRNAME + '/data/sample2.pdb'
        
        p = PDBParser(PERMISSIVE=1)

        query_struct = p.get_structure(os.path.basename (path1), path1)
        against_struct = p.get_structure(os.path.basename (path2), path2)

        query_complex = Complex (query_struct)
        against_complex = Complex (against_struct)
    
        query_complex.get_fp ()
        against_complex.get_fp ()
    
        query_fp_string = query_complex.fp2str ()
        against_fp_string = against_complex.fp2str ()

        query = FPWithComplex (query_complex, query_fp_string)
        against = FPWithComplex (against_complex, against_fp_string)
    
        score1, score2, score3 = similarity_between (query, against)
        
        expected = {"score1": 118.00269647021572, "score2": -8.0, "score3": 20.0}
        actual = {"score1": score1, "score2": score2, "score3": score3}
        
        self.assertEqual (actual, expected)

    def test_basic_with_epitope (self):
        path1 = DIRNAME + '/data/sample1.pdb'
        path2 = DIRNAME + '/data/sample2.pdb'
        
        p = PDBParser(PERMISSIVE=1)

        query_struct = p.get_structure(os.path.basename (path1), path1)
        against_struct = p.get_structure(os.path.basename (path2), path2)

        query_complex = Complex (query_struct, epitope = [211,213,214,224,225,226,227,228,229])
        against_complex = Complex (against_struct, epitope = [216,217,218,219,220,221])
    
        query_complex.get_fp ()
        against_complex.get_fp ()
    
        query_fp_string = query_complex.fp2str ()
        against_fp_string = against_complex.fp2str ()

        query = FPWithComplex (query_complex, query_fp_string)
        against = FPWithComplex (against_complex, against_fp_string)
    
        score1, score2, score3 = similarity_between (query, against)
        
        expected = {'score1': 34.705754203703862, 'score2': 0, 'score3': 6}
        actual = {"score1": score1, "score2": score2, "score3": score3}
        
        self.assertEqual (actual, expected)

if __name__ == '__main__':
    unittest.main ()
