"""
define the `similarity` of two chains
"""

import re

from collections import OrderedDict,defaultdict
from UserDict import UserDict
import math
from numpy import corrcoef,array,vstack,zeros

class residue_fp(object):
    """
    The fingerprint for  residue
    """
    def __init__(self,fp_str,comp,residue_id):
        """
        Parameter:
        fp_str: the string representation of the fingerprint
        comp: the complex associated with the residue
        residue_id: the residue id
        """
        self.fp_str = fp_str
        self.complex = comp
        self.residue_id = residue_id

        self.res = filter (lambda r: r.get_id() == self.residue_id, self.complex.residues) [0] #search for the residue in the complex that has id = self.residue_id

    def get_edit_dist_to(self,residue_fp1):
        """
        Measure the disntace of self.fp_str to the other finger print using the correlation coefficient

        Paramter: 
        residue_fp1: the other finger print

        Return:
        The distance, float
        """
        return corrcoef(vstack((self.fp_str,residue_fp1.fp_str)))[0,1]

    def within_range_of(self, o_res, range_dist = 20):
        """
        determine if the residue is within distance `range_dist` of the `o_res`

        Parameter:
        o_res: the object residue
        range_dist: the distance, int

        Return: 
        boolean
        """
        def distance(xyz1,xyz2):
            x1,y1,z1 = tuple(xyz1)
            x2,y2,z2 = tuple(xyz2)
            return math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)

        outer_c = 0
        
        o_res_total = len(o_res.get_list())
        res_total = len(self.res.get_list())
        
        # o_res_total = len(o_res.atom)
        # res_total = len(self.res.atom)
        for a1 in self.res.get_list ():
            inner_c = 0
            for a2 in o_res.get_list ():
                if distance(a1.get_coord(), a2.get_coord()) < range_dist:
                    inner_c += 1
            if inner_c >= o_res_total / 2.:
                outer_c += 1
        if outer_c >= res_total / 2.:
            return True
        return False
    
    def get_string(self):
        """
        Get the fingerprint string representation
        """
        return self.fp_str

class residue_fp_list(UserDict):
    """
    The fingerprints of a list of residues
    """
    def __init__(self, fp_string = '', comp = None):
        """
        Parameter:
        fp_string: the total/big string of the fingerprints, in the format: res_id1\t1,0,1,1\nres_id2\t1,1,0,0
        comp: the associating complex
        """
        UserDict.__init__(self)

        self.comp = comp
        self.data = OrderedDict()#the place to store the fingerprints
        
        from StringIO import StringIO
        for line in StringIO (fp_string).readlines():
            s_line = re.split(r"[ ,\t]",line) 
            residue_id,fp = int(s_line[0]), map(float, s_line[1:])
            
            self.data[residue_id] = residue_fp(fp, self.comp, residue_id)


class dist_mat(UserDict):
    """
    The distance matrix that measures pairwise distance between residue fingerprints from complex A and residue fingerprints from complex B 
    """
    
    def __init__(self,fp1,fp2):
        """
        Parameter:
        fp1: the first fingerprint list 
        fp2: the second fingerprint list
        """
        self.data = defaultdict(dict)
        
        self.fp1 = fp1.fp
        self.fp2 = fp2.fp
        
        for res1, fp1 in self.fp1.items():
            for res2, fp2 in self.fp2.items():
                self.data[res1][res2] = fp1.get_edit_dist_to(fp2)

        self.clustered_fp1_res = set()
        self.clustered_fp2_res = set()
        self.clusters = []

    def find_closest_tuple(self):
        """
        find the closest residue pairs in terms of their fingerprint distance meanwhile ignoring those already in the clusters or cosidered not suitable
        
        Return: 
        the residue pair with the **maximum** fingerprint similarity with their similarity score
        """
        from itertools import chain
        
        all_pairs = list(chain.from_iterable(map(lambda (row, cols): map (lambda (col, val): (row, col, val), cols.items ()), self.data.items ())))#from defaultdict(<type 'dict'>, {res1: {res2: 3, res3: 3}, res2: {res4: 6}}) to [(res1, res2, 3), (res1, res3, 3), (res2, res4, 6)]
        
        tuples = filter(lambda tpl: 
                        tpl not in self.not_suitable_tuple and 
                        tpl [0] not in self.clustered_fp1_res and 
                        tpl [1] not in self.clustered_fp2_res,
                        all_pairs)  #filter out those in the not_suitable list
                
        if len (tuples) == 0:
            return None
        else:
            return max (tuples, key = lambda (_,__,num): num)

    def find_cluster_helper(self):
        """
        Perform one round of residue cluster discovery and add the discoverd cluster to the class variable
        
        steps: 
        1. find the closest pair as the cluster base
        2. use the base to find within-range pairs in a greedy-manner, always consider those pairs whose are more similar with each other
        """
        
        cluster = set()
        
        self.not_suitable_tuple = set()
        
        self.extending_tuple = self.find_closest_tuple() #we can the closest residue as the expansion point
        
        if self.extending_tuple is None: #if no more, it's done
            return

        cluster.add(self.extending_tuple); ext_res1 , ext_res2 , dist = self.extending_tuple #add the first tuple to the cluster and unpack it
        
        self.clustered_fp1_res.add(ext_res1); self.clustered_fp2_res.add(ext_res2)#add residues to corresponding sides
        
        center_res1, center_res2 = self.fp1[ext_res1].res, self.fp2[ext_res2].res #get the expansion residues
        
        while True:
            t = self.find_closest_tuple() #get the next closest(edit distance) tuple
            if t is None:#cannot find any appropriate tuple, quit
                break
            
            res1 , res2 , dist = t
            if self.fp1[res1].within_range_of(center_res1) and self.fp2[res2].within_range_of(center_res2):
                #check if it is within the range of extending tuple,if so, add it to the cluster
                cluster.add(t)

                self.clustered_fp1_res.add(res1)
                self.clustered_fp2_res.add(res2)
                #print "cluster size: %d,total residue number: f1 = %d, f2 = %d" %(len(cluster),len(self.fp1),len(self.fp2))
            else:
                #print "out of range"
                self.not_suitable_tuple.add(t)
                            
        self.clusters.append(cluster)
        self.extending_tuple = None
        
        return cluster

    def find_clusters(self):
        while self.find_cluster_helper():pass
        return self.clusters
        
res_sim_mat = {"AA" : 1,"AC" : 0,"AD" : 0,"AE" : 0,"AF" : 0,"AG" : 0,"AH" : 0,"AI" : 0,"AK" : 0,"AL" : 0,"AM" : 0,"AN" : -2,"AP" : -1,"AQ" : -1,"AR" : -1,"AS" : 1,"AT" : 0,"AV" : 0,"AW" : -3,"AY" : -2,"CC" : 9,"CD" : -3,"CE" : -4,"CF" : -2,"CG" : -3,"CH" : -3,"CI" : -1,"CK" : -3,"CL" : -1,"CM" : -1,"CN" : -3,"CP" : -3,"CQ" : -3,"CR" : -3,"CS" : -1,"CT" : -1,"CV" : -1,"CW" : -2,"CY" : -2,"DD" : 6,"DE" : 2,"DF" : -3,"DG" : -1,"DH" : -1,"DI" : -3,"DK" : -1,"DL" : -4,"DM" : -3,"DN" : 1,"DP" : -1,"DQ" : 0,"DR" : -2,"DS" : 0,"DT" : -1,"DV" : -3,"DW" : -4,"DY" : -3,"EE" : 5,"EF" : -3,"EG" : -2,"EH" : 0,"EI" : -3,"EK" : 1,"EL" : -3,"EM" : -2,"EN" : 0,"EP" : -1,"EQ" : 2,"ER" : 0,"ES" : 0,"ET" : -1,"EV" : -2,"EW" : -3,"EY" : -2,"FF" : 6,"FG" : -3,"FH" : -1,"FI" : 0,"FK" : -3,"FL" : 0,"FM" : 0,"FN" : -3,"FP" : -4,"FQ" : -3,"FR" : -3,"FS" : -2,"FT" : -2,"FV" : -1,"FW" : 1,"FY" : 3,"GG" : 6,"GH" : -2,"GI" : -4,"GK" : -2,"GL" : -4,"GM" : -3,"GN" : 0,"GP" : -2,"GQ" : -2,"GR" : -2,"GS" : 0,"GT" : -2,"GV" : -3,"GW" : -2,"GY" : -3,"HH" : 8,"HI" : -3,"HK" : -1,"HL" : -3,"HM" : -2,"HN" : 1,"HP" : -2,"HQ" : 0,"HR" : 0,"HS" : -1,"HT" : -2,"HV" : -3,"HW" : -2,"HY" : 2,"II" : 4,"IK" : -3,"IL" : 2,"IM" : 1,"IN" : -3,"IP" : -3,"IQ" : -3,"IR" : -3,"IS" : -2,"IT" : -1,"IV" : 3,"IW" : -3,"IY" : -1,"KK" : 5,"KL" : -2,"KM" : -1,"KN" : 0,"KP" : -1,"KQ" : 1,"KR" : 2,"KS" : 0,"KT" : -1,"KV" : -2,"KW" : -3,"KY" : -2,"LL" : 4,"LM" : 2,"LN" : -3,"LP" : -3,"LQ" : -2,"LR" : -2,"LS" : -2,"LT" : -1,"LV" : 1,"LW" : -2,"LY" : -1,"MM" : 5,"MN" : -2,"MP" : -2,"MQ" : 0,"MR" : -1,"MS" : -1,"MT" : -1,"MV" : 1,"MW" : -1,"MY" : -1,"NN" : 6,"NP" : -2,"NQ" : 0,"NR" : 0,"NS" : 1,"NT" : 0,"NV" : -3,"NW" : -4,"NY" : -2,"PP" : 7,"PQ" : -1,"PR" : -2,"PS" : -1,"PT" : -1,"PV" : -2,"PW" : -4,"PY" : -3,"QQ" : 5,"QR" : 1,"QS" : 0,"QT" : -1,"QV" : -2,"QW" : -2,"QY" : -1,"RR" : 5,"RS" : -1,"RT" : -1,"RV" : -3,"RW" : -3,"RY" : -2,"SS" : 4,"ST" : 1,"SV" : -2,"SW" : -3,"SY" : -2,"TT" : 5,"TV" : 0,"TW" : -2,"TY" : -2,"VV" : 4,"VW" : -3,"VY" : -1,"WW" : 11,"WY" : 2,"YY" : 7}

class FPWithComplex(object):
    """
    Fingerprint with complex
    """
    def __init__(self, complex, fp_string):
        """
        Parameter:
        complex: the complex structure
        fp_string: the complex's residues' fingerprint string
        """
        self.complex = complex #reads the structure
        self.fp = residue_fp_list(fp_string, complex) #instantiates the fp list
    
def similarity_between(c1,c2):
    """
    Parameter: 
    c1: complex fingerprint 1
    c2: complex fingerprint 2

    Return:
    Three part similiarity scores
    """
    #print "generating distance matrix"
    pair = dist_mat(c1,c2)
    
    #print "finding clusters"
    clusters = pair.find_clusters()
    
    
    #val1, val2 and val3 starts
    val1,val2,val3 = 0, 0, 0

    #print "calculating value1"
    for c in clusters:
        for t in c:
            val1 += t[2] * 10#the edit distance
    
    #print "calculating value2"
    for c in clusters:
        for t in c:
            res1,res2,dist = t
            res1_code = pair.fp1[res1].res.get_resname_abbrev()
            res2_code = pair.fp2[res2].res.get_resname_abbrev()
            try:
                val2 += res_sim_mat[res1_code + res2_code]
            except KeyError:
                val2 += res_sim_mat[res2_code + res1_code ]

    #print "calculating value3"
    for c in clusters:
        val3 += len(c)#pair count

    return val1 , val2 , val3
