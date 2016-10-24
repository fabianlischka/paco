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
    Returns the vector as an ndarray
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

    First line contains: n m (is symmetric) nbcoef
    Next lines contain, for each nonzero coefficient:   i j a_ij where (i,j) \in  {1,...,n}x{1,...,m}
    Returns the matrix as a dense ndarray.
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

def readMeshFF(f):
    """Reads a mesh from a file saved by FreeFem++.

    First line contains: v t e
    where v = number of vertices, t = number of triangles, e = number of edges
    Next v lines contain: x y l_b
    where x,y = x,y coord of vertex (both floats), l_b = boundary label (int)
    Next t lines contain: v1 v2 v3 l_r
    where v1,v2,v3 = indices of vertices constituting triangle, l_r = region label (int)
    Next e lines contain: v1 v2 l_b
    where v1,v2 = indices of vertices constituting boundary edge, l_b = boundary label (int)
    """
    # strip comment header
    row = f.next().strip()
    while row[0] == '#':
        row = f.next().strip()
    # read shape
    v, t, e = map(int, row.split())
    logging.info("mesh with %s vertices, %s triangles, %s edges: attempting to read...", v, t, e)
    # read vertices
    vv = np.empty((v,3))
    for vi in range(v):
        row = f.next().strip()
        x, y, l_b = map(float, row.split())
        vv[vi] = x, y, l_b
    # read triangles or faces
    tt = np.empty((t,4), dtype = int)
    for ti in range(t):
        row = f.next().strip()
        v1, v2, v3, l_r = map(int, row.split())
        tt[ti] = v1, v2, v3, l_r
    # read edges
    ee = np.empty((e,3), dtype = int)
    for ei in range(e):
        row = f.next().strip()
        v1, v2, l_b = map(int, row.split())
        ee[ei] = v1, v2, l_b

    logging.info("successfully read %s lines", v+t+e)
    return (vv, tt, ee)


def fvm(vv,q,u,where): # finite volume method
    ii = np.zeros(12, dtype = int)
    jj = np.zeros(12, dtype = int)
    ss = np.zeros(12)

    count = 0
    for i in range(3): # python 0-based unlike matlab...
        ip  =  (i+1) % 3
        ipp = (ip+1) % 3
        unL = -( (q[ip,1] + q[i,1] - 2*q[ipp,1]) * u[0]
                -(q[ip,0] + q[i,0] - 2*q[ipp,0]) * u[1] )/6
        ii[count],   ss[count]   = vv[i],   unL
        ii[count+1], ss[count+1] = vv[ip], -unL
        if unL > 0:
            jj[count]   = vv[i]
            jj[count+1] = vv[i]
        else:
            jj[count]   = vv[ip]
            jj[count+1] = vv[ip]
        count += 2
        if where[i] != 0 and where[ip] != 0:
            unL = +( (q[ip,1] - q[i,1] ) * u[0]
                    -(q[ip,0] - q[i,0] ) * u[1] )/2
            if unL > 0:
                ii[count],   ss[count]   = vv[i], unL
                ii[count+1], ss[count+1] = vv[i], unL
                jj[count]   = vv[i]
                jj[count+1] = vv[i]
                count += 2
    return ii,jj,ss
