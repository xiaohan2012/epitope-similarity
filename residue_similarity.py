
"""
Defines similarity between different residues
"""
from copy import deepcopy
from csv import parse_csv_table

def _fill(matrix, extra_ordinary = [], extra_special = []):
    """
    Fill in missing fields 
    
    `extra_ordinary` specifies the residue id, whose value are averaged over ALL values
    `extra_special`  specifies the residue id, whose value are averaged over CERTAIN values

        >>> from StringIO import StringIO
        >>> data = parse_csv_table(StringIO(",A,B\\nA,1,2\\nB,3,4"))
        >>> newdata = _fill(data, ["X"], [("C", ("A", "B"))])
        >>> newdata["X"]["A"]
        2.0
        >>> newdata["X"]["B"]
        3.0
        >>> newdata["X"]["X"]
        2.5
        >>> newdata["C"]["A"]
        2.0
        >>> newdata["C"]["B"]
        3.0
        >>> newdata["C"]["C"]
        2.5
        >>> newdata["C"]["X"]
        2.5
        >>> newdata["X"]["C"]
        2.5
    """
    matrix_copy = deepcopy(matrix)
    original_keys = matrix_copy.keys()
    
    def block_avg(mat, row_names = None, col_names = None):
        """Calculate the average value of the submatrix specified by `row_names` and `col_names`
        Default to be the whole row or column.
        """
        if not col_names:
            col_names = mat.keys()
        if not row_names:
            row_names = mat.keys()

        return float(sum([mat[row_name][col_name] 
                    for col_name in col_names
                    for row_name in row_names])) / len(col_names) / len(row_names)
        
    for extra_key in extra_ordinary:
        for key in  original_keys:
            matrix[extra_key][key] = block_avg(matrix_copy, col_names = [key])
            matrix[key][extra_key] = block_avg(matrix_copy, row_names = [key])

        matrix[extra_key][extra_key] = block_avg(matrix_copy)
        
    for extra_key, fields in extra_special:
        for key in  original_keys:
            matrix[extra_key][key] = block_avg(matrix_copy, 
                                               row_names = fields, col_names = [key])
            matrix[key][extra_key] = block_avg(matrix_copy, 
                                               row_names = [key], col_names = fields)

        matrix[extra_key][extra_key] = block_avg(matrix_copy, fields, fields)

    # extra ordinary key VS extra special key
    for ordinary_key in extra_ordinary:
        for special_key, fields in extra_special:
            matrix[ordinary_key][special_key] = block_avg(matrix_copy, col_names = fields)
            matrix[special_key][ordinary_key] = block_avg(matrix_copy, row_names = fields)
            
    return matrix
    
    
def load_similarity_maitrx(path):
    """
    Load the similarity matrix with missing fields filled in
    
    >>> m = load_similarity_maitrx("ref/blosum.csv")
    >>> m["X"]["X"]
    -1.045
    >>> m["X"]["A"]
    -0.85
    >>> m["B"]["A"]
    -2.0
    >>> m["A"]["B"]
    -1.5
    >>> m["X"]["B"]
    -1.0
    >>> m["B"]["X"]
    -1.1
    """
    
    matrix = parse_csv_table(open(path))
    return _fill(matrix, ["X"], [("B", ("N", "D"))])

if __name__ == "__main__":
    import doctest
    doctest.testmod()

    
