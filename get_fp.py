"""
Finger print generation(110 bit version)
"""

import os,  sys, glob, warnings, math
import numpy  as np

from collections import defaultdict

from vec_geom import angle, length

from util.dist import res_ca_dist



class Residue(object):
    """
    The residue class used for fingerprint generation
    """
    abbrev_mapping = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "ASX": "B", "CYS": "C", "GLN": "Q", "GLU": "E", "GLX": "Z", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "UNK": "X"}
    
    hydro_dict = {'A':0.61,'C':1.07,'D':0.46,'E':0.47,'F':2.02,'G':0.07,'H':0.61,'I':2.22,'K':1.15,'L':1.53,'M':1.18,'N':0.06,'P':1.95,'Q':0.0,'R':0.6,'S':0.05,'T':0.05,'V':1.32,'W':2.65,'Y':1.88, 
                  'X': 0.9974999999999999, 'B': 0.26, 'Z': 0.235, 'J': 1.875}
    charged_dict = {'A':-0.01,'C':0.12,'D':0.15,'E':0.07,'F':0.03,'G':0.0,'H':0.08,'I':-0.01,'K':0.0,'L':-0.01,'M':0.04,'N':0.06,'P':0.0,'Q':0.05,'R':0.04,'S':0.11,'T':0.04,'V':0.01,'W':0.0,'Y':0.03, 
                    'X':0.04000000000000001, 'B': 0.105, 'Z': 0.06, 'J': -0.01}
    h_bond_dict = {'A':0,'C':0,'D':1,'E':1,'F':0,'G':0,'H':1,'I':0,'K':2,'L':0,'M':0,'N':2,'P':0,'Q':2,'R':4,'S':1,'T':1,'V':0,'W':1,'Y':1, 
                   'X':0.85, 'B': 1.5, 'Z': 1.5, 'J': 0.0}

    def __init__(self, res, comp):

        self.c = comp
        self.fp = None #to be assigned later
        self.ca = None
        self.body = res
        self.resnum = res.get_id()
        if isinstance(self.resnum, tuple):
            self.resnum = self.resnum[1]
        assert type(self.resnum) == int,  "resnum should be integer but is '%r' of value '%r'" %(type(self.resnum), self.resnum)

        for a in res: #iterate over the atoms in the residue
            if a.get_name () == 'CA':
                self.ca = a
                break
            
        if self.ca is None:
            sys.stderr.write("resnum %d has no CA atom\n" %self.resnum)
            
    def get_id (self):
        """ residue id"""
        return self.body.get_id() [1]

    def get_list (self):
        """
        return the atom list
        """
        return self.body.get_list ()

    def get_resname_abbrev (self):
        """
        get the abbreviation code of residue name
        """
        return Residue.abbrev_mapping[self.body.get_resname()]

    def set_fp_length (self, length):
        """
        set the length of the fingerprint
        """
        self.fp = [0] * length
        
    def is_valid (self):
        """
        if the residue is valid (contain CA atom or not)
        """
        return self.ca is not None
        
    def turn_on_bit(self,bit_num,count):
        """
        for the bit position, `bit_num`, set the value to `count`
        """
        self.fp[bit_num] = count

    def get_surrounding_res(self):
        """
        get surrouding residues
        """
        for res in self.c.get_residues ():
            if res.ca == self.ca:continue
            yield res

    def get_struct_fp(self, spin_image_radius_step, spin_image_height_step, 
                      spin_image_radius_ind_min, spin_image_radius_ind_max, 
                      spin_image_height_ind_min, spin_image_height_ind_max,
                      spin_image_radius_seg_cnt, spin_image_height_seg_cntm,
                      offset = 0):
        """
        get the structural based fingerprints, given a set of spinimage parameters and the fingerprint bit starting position
        """
        def get_region_number(res_ca):
            """
            according to the given residue CA atom, calculate which spin image region this residue should belong to
            """
            my_point = np.array(self.ca.get_coord ())
            v1 = my_point - np.array([0,0,0])#the original point
            v2 = np.array(res_ca.get_coord ()) - my_point
            
            ang = angle(v1, v2)

            v2_len = length(v2)
            
            x, y = v2_len * math.sin(ang), v2_len * math.cos(ang)
            
            if x < 0:
                raise ValueError("x cannot be negative")

            x_ind  = math.floor( x / spin_image_radius_step)
            y_ind  = math.floor( y / spin_image_height_step)
            
            if x_ind >= spin_image_radius_ind_min and x_ind <= spin_image_radius_ind_max and\
               y_ind >= spin_image_height_ind_min and y_ind <= spin_image_height_ind_max:
                #within the range
                return int((y_ind - spin_image_height_ind_min) * spin_image_radius_seg_cnt + x_ind)
            else:
                return None
            
        #type one fp
        d_ = defaultdict(int)
        for res in self.get_surrounding_res():
            num = get_region_number(res.ca)
            if num is not None:
                d_[num] += 1
        
        for bit_num,count in d_.iteritems():
            self.turn_on_bit(offset + bit_num, count)

        return self.fp           
    
    def get_surrounding_fp(self, dist_step, dist_ind_min, dist_ind_max, offset = 0):
        """
        get the surrouding-residue based fingerprints
        """
        #type two fp
        d_ = defaultdict(list)
        for other in self.get_surrounding_res():
            dist = res_ca_dist(other,self)
            dist_ind = math.floor( dist / dist_step )

            if dist_ind > dist_ind_min and dist_ind < dist_ind_max:
                #within range
                d_[dist_ind].append(other)
        
        slice_count = dist_ind_max - dist_ind_min #how many slices for the sphere
        
        for i in xrange(slice_count + 1):
            # for layer i
            h_bond, charged , hydro = 0 , 0 , 0
            for res in d_[i]:
                # for res in layer i
                code = res.get_resname_abbrev()
                hydro += Residue.hydro_dict[code]
                charged += Residue.charged_dict[code]
                h_bond += Residue.h_bond_dict[code]
            
            #fp for layer i,in the 3 aspects
            self.turn_on_bit(offset + i , hydro)
            self.turn_on_bit(offset + slice_count + i , charged)
            self.turn_on_bit(offset + slice_count*2 + i , h_bond)

    def __repr__(self):            
        return "Residue#%d" %(self.get_id())


class Complex(object):
    def __init__(self , pdb_st, epitope = []):
        """
        pdb_st: the pdb structure
        epitope: list of int, the epitope residue number list, also they are the residues to be considered
        """
        self.st = pdb_st

        all_epitope = (True if len(epitope) == 0 else False)
        
        self.residues = []
            
        for res in self.st.get_residues():
            res = Residue(res, self)
            # if res.is_valid () and (all_epitope or (len(epitope) != 0 and res.resnum in epitope)):
            if res.is_valid() and (all_epitope or res.get_id() in epitope): #if the res is valid and filter those residue in the epitope
                self.residues.append(res)

        assert len(self.residues) > 0, "%s should have at least one residue" %(pdb_st.id)

    def get_residues (self):
        return self.residues
        
    def get_fp(self, 
               spin_image_radius_range = (0, 20),
               spin_image_height_range =  (-30, 10),
               sphere_radius_range = (0, 20),
               spin_image_radius_step = 2,               
               spin_image_height_step = 5,
               sphere_radius_step = 2):
        """
        get the 110-bit fingerprint
        """
        #radius part
        spin_image_radius_min , spin_image_radius_max = spin_image_radius_range
        spin_image_radius_seg_cnt = ( spin_image_radius_max - spin_image_radius_min ) / spin_image_radius_step
        spin_image_radius_ind_min , spin_image_radius_ind_max = int(spin_image_radius_min / spin_image_radius_step), int(spin_image_radius_max / spin_image_radius_step - 1) #index min and max must be integer. Maybe some warning should be put here.

        #height part
        spin_image_height_min , spin_image_height_max = spin_image_height_range
        spin_image_height_seg_cnt = ( spin_image_height_max - spin_image_height_min ) / spin_image_height_step
        spin_image_height_ind_min , spin_image_height_ind_max = int(spin_image_height_min / spin_image_height_step), int(spin_image_height_max / spin_image_height_step - 1) #index min and max must be integer. Maybe some warning should be put here.

        cylinder_slice_count = (spin_image_height_ind_max - spin_image_height_ind_min + 1) * (spin_image_radius_ind_max - spin_image_radius_ind_min + 1)
        
        #sphere part
        dist_min, dist_max = sphere_radius_range
        dist_step = sphere_radius_step
        
        dist_ind_min, dist_ind_max =  int(dist_min / dist_step), int(dist_max / dist_step - 1)#index min and max must be integer. Maybe some warning should be put here.
        
        sphere_slice_count = dist_ind_max - dist_ind_min + 1
        
        fp_size = cylinder_slice_count + sphere_slice_count * 3
        self.res_list = []
        
        for i , res in enumerate(self.residues):
            
            res.set_fp_length (fp_size)#initialize the fingerprint
            
            res.get_surrounding_fp(dist_step, dist_ind_min, dist_ind_max, offset = cylinder_slice_count)
            
            res.get_struct_fp(spin_image_radius_step, spin_image_height_step, 
                              spin_image_radius_ind_min, spin_image_radius_ind_max, 
                              spin_image_height_ind_min, spin_image_height_ind_max,
                              spin_image_radius_seg_cnt, spin_image_height_seg_cnt)
            
            self.res_list.append(res)
        
        return self.res_list
    
    def fp2str (self):
        """
        print the fingerprint string representation
        """
        return "\n".join("%s\t%s" %(res.get_id()," ".join(map(str,res.fp))) for res in self.res_list)

    def write_fp_to_file(self,path):
        """
        write the fp string to file given path
        """
        if not os.path.exists( os.path.dirname(path) ):
            os.mkdir(os.path.dirname(path))
        with open(path,'w') as f:
            f.write(self.fp2str ())
            
if __name__ == "__main__":
    import sys
    
    from Bio.PDB.PDBParser import PDBParser
    p = PDBParser(PERMISSIVE=1)
    
    try:
        path = sys.argv[1]
        horizontal = int(sys.argv[2])
        vertical = int(sys.argv[3])
        
        epitope = [] if len(sys.argv) == 4 else map(int, sys.argv[-1].split (','))
    except:
        print "Expecting format like: \n\tpython get_fp.py path/to/pdb/file horizontal-step  vertical-step [a list of epitope residue numbers]"
        print 'For example: \n\tpython get_fp.py test/data/sample1.pdb  5 4 211,213,214,224,225,226,227,228,229'
        sys.exit (-1)
        
    struct = p.get_structure(os.path.basename (path), path)
    c = Complex (struct, epitope)
    
    try:
        c.get_fp (spin_image_radius_step = horizontal,
                  spin_image_height_step = vertical)
    except:
        "Error"
    
    print c.fp2str ()
