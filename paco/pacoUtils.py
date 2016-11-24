from distutils import spawn
import logging
import numpy as np
import os
import os.path
import string


# Config
_extraPaths = ['/u1/local/ff++/bin', '/usr/local/bin']
_searchPath = os.environ['PATH']
_pathFF = None
for extraPath in _extraPaths:
    _searchPath += (os.pathsep + extraPath)

def getPathFF():
    global _extraPaths, _searchPath, _pathFF
    if _pathFF is None:
        logging.debug("attempt to find FreeFem++-nw binary...")
        p =  spawn.find_executable("FreeFem++-nw", _searchPath)
        if p is None:
            errMsg = "Cannot find 'FreeFem++-nw', looked in %s" % _searchPath
            logging.error(errMsg)
            raise RuntimeError(errMsg)

        _pathFF = p
        logging.debug("successfully found FreeFem++-nw binary at %s.", _pathFF)
    return _pathFF


def readVectFF(f):
    """Read a vector in the form output by FreeFem++. Return ndarray.

    First line contains the number of entries: nbcoef.
    Next lines contain the entries, with 5 entries per line,
    separated by tab.
    Entries are read left to right, then line by line.
    Usage:
    with open('..\generate\u.txt', 'r') as f:
        x = readVectFF(f)
    """
    logging.debug("attempt to read FF vect...")
    nbcoef = int(f.next().strip())
    logging.debug("attempt to read %s values from FF vect", nbcoef)
    # note: can't use loadtxt or genfromtxt, as last row might have fewer vals
    # if number of coeffs is not multiple of 5...
    arr = np.empty(nbcoef)
    oldc = 0
    for row in f:
        vals = map(float, row.strip().split())
        c = oldc + len(vals)
        if c > nbcoef:
            break
        arr[oldc:c] = vals
        oldc = c

    if c != nbcoef:
        errMsg = "readVectFF: Inconsistent dimension, expected %d, got %d" % (
            nbcoef, c)
        logging.error(errMsg)
        raise RuntimeError(errMsg)
    logging.info("successfully read %s values from FF vect.", nbcoef)
    return arr


def createDenseFromIJX(iijjxx, M=-1, N=-1, issym=False, ij_onebased=True):
    """Create a dense matrix S from (nnz x 3) array with i, j, x. Return ndarray.

    If M, N are provided, they specify size of S as m-by-n.
    Note: Elements in v that have duplicate coordinates in i,j are *added*.
    If issym, for given (i,j), transposed (j,i) value is also set.
    If ij_onebased, an offset of 1 is subtracted from i, j before placing into
    the (zero-based) ndarray.
    """
    # populate matrix
    if M == -1:
        M = max(ii)
    if N == -1:
        N = max(jj)
    if ij_onebased:
        offset = 1
    else:
        offset = 0
    res = np.zeros((M, N))
    for ifloat, jfloat, x in iijjxx:
        i, j = int(ifloat) - offset, int(jfloat) - offset  # python is 0-based
        res[i, j] += x
        if issym:
            res[j, i] += x
    return res


def readSparseFF(f):
    """Read a sparse matrix (Morse) from file saved by FreeFem++. Return ndarray.

    First line contains: n m (is symmetric) nbcoef
    Next lines contain, for each nonzero coefficient: i j a_ij
                        where (i,j) \in  {1,...,n}x{1,...,m}
    Returns the matrix as a dense ndarray.
    Usage:
    with open('..\generate\stiff.txt', 'r') as f:
        s = readSparseFF(f)
    """

    # strip comment header
    logging.debug("attempt to read FF sparse matrix...")
    row = f.readline().strip()
    while row[0] == '#':
        row = f.readline().strip()
    # read shape
    n, m, issym, nbcoef = map(int, row.split())
    logging.debug(
        "read sparse FF (%s x %s), sym? %s, %s values...", n, m, issym, nbcoef)
    # read nzeros
    ijv = np.loadtxt(f)
    nrows, ncols = ijv.shape
    # check dims
    if nrows != nbcoef:
        errMsg = "Inconsistent dimensions, expected %d, got %d" % (
            nbcoef, nrows)
        logging.error(errMsg)
        raise RuntimeError(errMsg)
    if ncols != 3:
        errMsg = "Inconsistent dimensions, expected %d, got %d" % (3, ncols)
        logging.error(errMsg)
        raise RuntimeError(errMsg)
    # populate matrix
    res = createDenseFromIJX(ijv, M=n, N=m,  # note weird FF convention
                             issym=(issym != 0), ij_onebased=True)
    # res = np.zeros((n,m))
    # for ifloat,jfloat,v in ijv:
    #     i, j = int(ifloat) - 1, int(jfloat) - 1  # python is 0-based
    #     res[i,j] = v
    #     if issym != 0:
    #         res[j,i] = v

    logging.info("successfully read %s values from FF sparse matrix.", nbcoef)
    return res


def readMeshFF(f):
    """Read a mesh from file saved by FreeFem++. Return (vv, tt, ee).

    Fileformat:
    First line contains: v t e
        v = number of vertices, t = number of triangles, e = number of edges
    Next v lines contain: x y l_b
        x,y = coord of vertex (both floats), l_b = boundary label (int)
    Next t lines contain: v1 v2 v3 l_r
        v1,v2,v3 = indices of vertices constituting triangle, l_r = region label (int)
    Next e lines contain: v1 v2 l_b
        v1,v2 = indices of vertices constituting boundary edge, l_b = boundary label (int)

    Return (vv, tt, ee), where
    vv is a v x 3 ndarray of floats: x, y, l_b
    tt is a t x 4 ndarray of ints:   v1, v2, v3, l_r
    ee is a e x 3 ndarray of ints:   v1, v2, l_b
    """
    logging.debug("attempt to read FF mesh...")
    # strip comment header
    row = f.readline().strip()
    while row[0] == '#':
        row = f.readline().strip()
    # read shape
    v, t, e = map(int, row.split())
    logging.debug(
        "read FF mesh with %s vertices, %s triangles, %s edges...", v, t, e)
    # read vertices
    vv = np.empty((v, 3))
    for vi in range(v):
        row = f.next().strip()
        x, y, l_b = map(float, row.split())
        vv[vi] = x, y, l_b
    # read triangles or faces
    tt = np.empty((t, 4), dtype=int)
    for ti in range(t):
        row = f.next().strip()
        v1, v2, v3, l_r = map(int, row.split())
        tt[ti] = v1, v2, v3, l_r
    # read edges
    ee = np.empty((e, 3), dtype=int)
    for ei in range(e):
        row = f.next().strip()
        v1, v2, l_b = map(int, row.split())
        ee[ei] = v1, v2, l_b

    logging.info(
        "successfully read mesh from FF with %s vertices, %s triangles, %s edges.",
        v, t, e)

    return (vv, tt, ee)


def with_file(filename, callback):
    """Call `callback` with the opened file `filename`. Return result of call.

    See http://softwareengineering.stackexchange.com/questions/262346/should-i-pass-in-filenames-to-be-opened-or-open-files
    """
    with open(filename, 'r') as f:
        return callback(f)


def write_csr_matrix(A, f):
    """Write dense A in compressed sparse row format (zero-based) to file f.

    Text file format:
    - one number per line
    - First number is m = number of rows
    - Second number is n = number of columns
    - Next m+1 numbers are the row pointers
    - Next nnz numbers are the column indices
    - Next nnz numbers are the matrix values
    """
    m, n = A.shape
    logging.debug("write (%s x %s) matrix in CSR format...", m, n)
    f.write('%d\n' % m)
    f.write('%d\n' % n)
    counter = 0
    for i in range(m):
        f.write('%d\n' % counter)
    #    counter = counter + nnz(A(i,:));
        counter += len(A[i].nonzero()[0])
    f.write('%d\n' % counter)

    counter = 0
    for i in range(m):
        for k in A[i].nonzero()[0]:
            #        fprintf(fid,'%d\n',k-1)
            f.write('%d\n' % k)

    for i in range(m):
        for k in A[i].nonzero()[0]:
            #         fprintf(fid,'%.16e\n',full(A(i,k)));
            f.write('%r\n' % A[i, k])
            # note: use %r = repr() to reduce rounding errors
            #  %s or str() outputs fewer digits
            #  in particular, float(str(x)) might not be == x
            #  while float(repr(x)) == x, AFAIK
    logging.info("successfully wrote (%s x %s) matrix in CSR format.", m, n)


def read_csr_matrix(f):
    """Read A in compressed sparse row format from file f. Return dense ndarray.

    Text file format:
    - one number per line
    - First number is m = number of rows
    - Second number is n = number of columns
    - Next m+1 numbers are the row pointers (int, zero-based)
    - Next nnz numbers are the column indices (int, zero-based)
    - Next nnz numbers are the matrix values (float)
    """
    logging.debug("attempting to read matrix in CSR format...")
    # strip comment header
    row = f.next().strip()
    while row[0] == '#':
        row = f.next().strip()
    # read shape
    m = int(row)
    row = f.next().strip()
    n = int(row)
    logging.debug("attempting to read (%s x %s) matrix in CSR format...", m, n)

    # read m+1 row pointers
    counter = 0
    rowPointers = np.empty(m + 1, dtype=int)  # we store nnz in the last place
    for i in range(m + 1):
        rowPointers[i] = int(f.next().strip())

    nnz = rowPointers[m]

    # read nnz colum indices
    colIndices = np.empty(nnz, dtype=int)
    for i in range(nnz):
        colIndices[i] = int(f.next().strip())
        if colIndices[i] >= n:
            errMsg = "Inconsistent dims, col idx %d > %d" % (colIndices[i], n)
            logging.error(errMsg)
            raise RuntimeError(errMsg)

    # read nnz values
    values = np.empty(nnz)
    for i in range(nnz):
        values[i] = float(f.next().strip())

    # populate matrix
    res = np.zeros((m, n))
    for i in range(m):
        for nzi in range(rowPointers[i], rowPointers[i + 1]):
            res[i, colIndices[nzi]] = values[nzi]

    logging.info(
        "successfully read (%s x %s) matrix in CSR format with nnz %s.",
        m, n, nnz)
    return res


def writeABCD(A, B, C, D, destDir=".", prefix=""):
    ensurePath(destDir)
    fpre = prefix + "s2_"
    with open(os.path.join(destDir, fpre+"A.txt"), 'w') as f:
        write_csr_matrix(A, f)
    with open(os.path.join(destDir, fpre+"B.txt"), 'w') as f:
        write_csr_matrix(B, f)
    with open(os.path.join(destDir, fpre+"C.txt"), 'w') as f:
        write_csr_matrix(C, f)
    with open(os.path.join(destDir, fpre+"D.txt"), 'w') as f:
        write_csr_matrix(D, f)


def readABCD(sourceDir=".", prefix=""):
    fpre = prefix + "s2_"
    A = with_file(os.path.join(sourceDir, fpre+'A.txt'), read_csr_matrix)
    B = with_file(os.path.join(sourceDir, fpre+'B.txt'), read_csr_matrix)
    C = with_file(os.path.join(sourceDir, fpre+'C.txt'), read_csr_matrix)
    D = with_file(os.path.join(sourceDir, fpre+'D.txt'), read_csr_matrix)
    return A, B, C, D


def writeParamsFile(destDir=".", prefix="", params):
    fpre = prefix + "s3_"
    with open('template_params.txt', 'r') as f:
        template = string.Template(f.read())
    outfile = template.substitute(params)
    with open(os.path.join(destDir, fpre+"params.txt"), 'w') as f:
        f.write(outfile)


def ensurePath(path):
    """Create path if it does not exist."""
    if len(path) == 0:
        return
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
