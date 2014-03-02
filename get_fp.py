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
    abbrev_mapping = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "ASX": "B", "CYS": "C", "GLN": "Q", "GLU": "E", "GLX": "Z", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
    
    def __init__(self , res , comp, 
                 spin_image_radius_range = (0, 20),
                 spin_image_radius_step = 2,
                 spin_image_height_range =  (-30, 10),
                 spin_image_height_step = 5,
                 sphere_radius_range = (0, 20),
                 sphere_radius_step = 2):
        self.hydro_dict = {'A':0.61,'C':1.07,'D':0.46,'E':0.47,'F':2.02,'G':0.07,'H':0.61,'I':2.22,'K':1.15,'L':1.53,'M':1.18,'N':0.06,'P':1.95,'Q':0.0,'R':0.6,'S':0.05,'T':0.05,'V':1.32,'W':2.65,'Y':1.88}
        self.charged_dict={'A':-0.01,'C':0.12,'D':0.15,'E':0.07,'F':0.03,'G':0.0,'H':0.08,'I':-0.01,'K':0.0,'L':-0.01,'M':0.04,'N':0.06,'P':0.0,'Q':0.05,'R':0.04,'S':0.11,'T':0.04,'V':0.01,'W':0.0,'Y':0.03}
        self.h_bond_dict={'A':0,'C':0,'D':1,'E':1,'F':0,'G':0,'H':1,'I':0,'K':2,'L':0,'M':0,'N':2,'P':0,'Q':2,'R':4,'S':1,'T':1,'V':0,'W':1,'Y':1}

        self.c = comp
        self.fp = [0] * 110
        self.ca = None
        self.body = res
        # self.resnum = res.resnum 
        self.resnum = res.get_id ()
        
        self.spin_image_radius_min , self.spin_image_radius_max = spin_image_radius_range
        self.spin_image_radius_step = spin_image_radius_step
        self.spin_image_radius_seg_cnt = ( self.spin_image_radius_max - self.spin_image_radius_min ) / self.spin_image_radius_step
        self.spin_image_radius_ind_min , self.spin_image_radius_ind_max = self.spin_image_radius_min / self.spin_image_radius_step, self.spin_image_radius_max / self.spin_image_radius_step - 1

        self.spin_image_height_min , self.spin_image_height_max = spin_image_height_range
        self.spin_image_height_step = spin_image_height_step
        self.spin_image_height_seg_cnt = ( self.spin_image_height_max - self.spin_image_height_min ) / self.spin_image_height_step
        self.spin_image_height_ind_min , self.spin_image_height_ind_max = self.spin_image_height_min / self.spin_image_height_step, self.spin_image_height_max / self.spin_image_height_step - 1

        self.orig_point = np.array([0,0,0])

        self.dist_min, self.dist_max = sphere_radius_range
        self.dist_step = sphere_radius_step
        self.dist_ind_min = self.dist_min / self.dist_step
        self.dist_ind_max = self.dist_max / self.dist_step - 1

        
        for a in res: #iterate over the atoms in the residue
            # if a.pdbname.strip() == "CA":
            if a.get_name () == 'CA':
                self.ca = a
                break
            
        if self.ca is None:
            warnings.warn ("resnum %d has no CA atom" %res.get_id ())
            
    def get_id (self):
        """ residue id"""
        return self.body.get_id () [1]

    def get_list (self):
        """
        return the atom list
        """
        return self.body.get_list ()

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

    def get_region_number(self,res_ca):
        """
        according to the given residue CA atom, calculate which spin image region this residue should belong to
        """
        my_point = np.array(self.ca.get_coord ())
        v1 = my_point - self.orig_point 
        v2 = np.array(res_ca.get_coord ()) - my_point
        
        ang = angle(v1, v2)

        v2_len = length(v2)
        
        x, y = v2_len * math.sin(ang), v2_len * math.cos(ang)
        
        if x < 0:
            raise ValueError("x cannot be negative")

        x_ind  = math.floor( x / self.spin_image_radius_step)
        y_ind  = math.floor( y / self.spin_image_height_step)
        
        if x_ind >= self.spin_image_radius_ind_min and x_ind <= self.spin_image_radius_ind_max and\
            y_ind >= self.spin_image_height_ind_min and y_ind <= self.spin_image_height_ind_max:
            #within the range
            return int((y_ind - self.spin_image_height_ind_min) * self.spin_image_radius_seg_cnt + x_ind)
        else:
            return None
            
    def get_resname_abbrev (self):
        """
        get the abbreviation code of residue name
        """
        return Residue.abbrev_mapping[self.body.get_resname()]

    def get_struct_fp(self):
        """
        get the structural based fingerprints
        """
        #type one fp
        d_ = defaultdict(int)
        for res in self.get_surrounding_res():
            num = self.get_region_number(res.ca)
            if num is not None:
                d_[num] += 1
        
        for bit_num,count in d_.iteritems():
            self.turn_on_bit(bit_num,count)

        return self.fp           
    
    def get_surrounding_fp(self):
        """
        get the surrouding-residue based fingerprints
        """
        #type two fp
        d_ = defaultdict(list)
        for other in self.get_surrounding_res():
            dist = res_ca_dist(other,self)
            dist_ind = math.floor( dist / self.dist_step )

            if dist_ind > self.dist_ind_min and dist_ind < self.dist_ind_max:
                #within range
                d_[dist_ind].append(other)

        for i in xrange(self.dist_ind_max - self.dist_ind_min + 1):
            # for layer i
            h_bond, charged , hydro = 0 , 0 , 0
            for res in d_[i]:
                # for res in layer i
                code = res.get_resname_abbrev()
                hydro += self.hydro_dict[code]
                charged += self.charged_dict[code]
                h_bond += self.h_bond_dict[code]
            
            #fp for layer i,in the 3 aspects
            self.turn_on_bit(80 + i , hydro)
            self.turn_on_bit(90 + i , charged)
            self.turn_on_bit(100 + i , h_bond)

    def last_30_bits(self):
        """
        get the type-2 fp, which is the last 30 bits
        """
        self.get_surrounding_fp()
        return self.fp[-30:]

    def __repr__(self):            
        return "ca atom index:%d" %(self.ca.index)


class Complex(object):
    def __init__(self , pdb_fp, epitope = []):
        """
        pdb_fp: the pdb structure
        epitope: the residue number list
        """
        self.st = pdb_fp

        all_epitope = (True if len(epitope) == 0 else False)

        self.residues = []
            
        for res in self.st.get_residues ():
            res = Residue(res, self)
            # if res.is_valid () and (all_epitope or (len(epitope) != 0 and res.resnum in epitope)):
            if res.is_valid () and (all_epitope or (len(epitope) != 0 and res.get_id () in epitope)): #if the res is valid and filter those residue in the epitope
                self.residues.append(res)

    def get_residues (self):
        return self.residues
        
    def get_fp(self):
        """
        get the 110-bit fingerprint
        """
        self.res_list = []
        for i , res in enumerate(self.residues):
            #print "residue %d" %i
            res.get_surrounding_fp()
            res.get_struct_fp()
            #print res.fp,len(res.fp)
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
    
    from PDB.PDBParser import PDBParser
    p = PDBParser(PERMISSIVE=1)
    
    path = sys.argv[1]
    epitope = [] if len(sys.argv) == 2 else map(int, sys.argv[2].split (','))

    struct = p.get_structure(os.path.basename (path), path)
    c = Complex (struct, epitope);
    
    c.get_fp ()
    print c.fp2str ()
