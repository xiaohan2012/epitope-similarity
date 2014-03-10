from setting import *

import unittest, os
from get_fp import Complex
from Bio.PDB.PDBParser import PDBParser
from similarity import FPWithComplex, similarity_between

class SimilarityTest (unittest.TestCase):
    def test_bitlength_110 (self):
        """
        nothing is specified, test whether the legnth is 110
        """
        path = DIRNAME + '/data/sample1.pdb'
        
        p = PDBParser(PERMISSIVE=1)

        struct = p.get_structure(os.path.basename (path), path)

        complex = Complex (struct)
    
        complex.get_fp ()

        actual_length = len(complex.fp2str ().split ('\n') [0].split (' '))
        expected = 110
        self.assertEqual (actual_length, expected)
    
        actual = complex.fp2str ()


    def test_bitlength_230 (self):
        """
        radius step is set to 2, test whether the legnth is 230
        """
        path = DIRNAME + '/data/sample1.pdb'
        
        p = PDBParser(PERMISSIVE=1)

        struct = p.get_structure(os.path.basename (path), path)

        complex = Complex (struct)
    
        complex.get_fp (spin_image_height_step = 2)
    
        actual_length = len(complex.fp2str ().split ('\n') [0].split (' '))
        expected = 230
        self.assertEqual (actual_length, expected)

if __name__ == '__main__':
    unittest.main ()
