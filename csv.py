from collections import defaultdict

def parse_csv_table(readable, separator = ","):
    """
    Parse csv table into dictionary
    
        >>> from StringIO import StringIO
        >>> dict(parse_csv_table(StringIO(",A,B\\nA,1,2\\nB,3,4")))
        {'A': {'A': 1.0, 'B': 2.0}, 'B': {'A': 3.0, 'B': 4.0}}
    """
    data = defaultdict(dict)
    
    # First line is the columns
    columns = filter(None, readable.readline().strip().split(separator))

    for line in readable:
        row, rest = line.strip().split(separator, 1)
        for value, col in zip(rest.split(separator), columns):
            data[row][col] = float(value)

    return data
    


if __name__ == "__main__":
    import doctest
    doctest.testmod()            
