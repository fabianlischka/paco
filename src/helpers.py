import logging
import numpy as np

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s. %(levelname)s in %(funcName)s: %(message)s')

# usage:
# with open('..\generate\u.txt', 'r') as f:
#    x = readVectFF(f)

# with open('..\generate\stiff.txt', 'r') as f:
#    s = readSparseFF(f)

def readVectFF(f):
    """Reads a vector from a file as saved by FreeFem++, returns ndarray

    First line contains the number of entries: nbcoef.
    Next lines contain the entries, with 5 entries per line,
    separated by tab.
    Entries are read left to right, then line by line.
    """

    nbcoef = int(f.next().strip())
    logging.info("attempt to read %s values", nbcoef)
    arr = np.loadtxt(f)
    arr = arr.flatten(order='C')
    if len(arr) != nbcoef:
        errMsg = "readVectFF: Inconsistent dimension, expected %d, got %d" %(nbcoef, len(arr))
        logging.error(errMsg)
        raise RuntimeError(errMsg)
    logging.info("successfully read %s values", nbcoef)
    return arr

def readSparseFF(f):
    """Reads a sparse matrix (Morse) from a file saved by FreeFem++, returns ndarray

    First line contains: n m (is symmetic) nbcoef
    Next lines contain, for each nonzero coefficient:   i j a_ij where (i,j) \in  {1,...,n}x{1,...,m}
    """

    # strip comment header
    row = f.next().strip()
    while row[0] == '#':
        row = f.next().strip()
    # read shape
    n, m, issym, nbcoef = map(int, row.split())
    logging.info("sparse (%s x %s), sym? %s, attempt to read %s values", n, m, issym, nbcoef)
    # read nzeros
    ijv = np.loadtxt(f)
    nrows, ncols = ijv.shape
    # check dims
    if nrows != nbcoef:
        errMsg = "Inconsistent dimensions, expected %d, got %d" %(nbcoef, nrows)
        logging.error(errMsg)
        raise RuntimeError(errMsg)
    if ncols != 3:
        errMsg = "Inconsistent dimensions, expected %d, got %d" %(3, ncols)
        logging.error(errMsg)
        raise RuntimeError(errMsg)
    # populate matrix
    res = np.zeros((n,m))
    for ifloat,jfloat,v in ijv:
        i, j = int(ifloat) - 1, int(jfloat) - 1  # python is 0-based
        res[i,j] = v
        if issym != 0:
            res[j,i] = v

    logging.info("successfully read %s values", nbcoef)
    return res
