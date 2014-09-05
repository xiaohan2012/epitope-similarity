# Similarity matrix utility functions
"""
List of util:
1. load a simmat in txt format 
2. print the similarity matrix to csv format
"""
from lmatrix import lmatrix

def load_cids(path):
    """() => set of str
    return a set of the complex ids
    >>> len(load_cids("fp_370_atg.txt"))
    236
    """
    return sorted(list(set([c.strip()
                            for l in open(path)
                            for c in l.split()[:2]])))


def load_simmat(path):
    """(str) => lmatrix

    load the lmatrix from file
    """
    


    #init matrixc
    m = lmatrix(load_cids(path))

    #set values
    for l in open(path).readlines():
        c1,c2,value = l.split()
        value = float(value)
        m[c1,c2] = value
        m[c2,c1] = value

    return m

def print_simmat_in_csv(path):
    mat = load_simmat(path)
    #mat = mat / mat.diagonal()
    #mat = mat 
    #print mat[:10,:10]
    print mat.to_csv_str()

if __name__ == "__main__":
    import sys
    path = sys.argv[1]
    print_simmat_in_csv(path)
