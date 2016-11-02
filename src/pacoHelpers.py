from distutils import spawn
import logging
import numpy as np
import os
import os.path
import shutil
import subprocess


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
        logging.info("successfully found FreeFem++-nw binary at %s.", _pathFF)
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


def fvm(vv, q, u, where):  # finite volume method
    """Finite volume method.

    Inputs:
    vv: int,       3. Vertices defining this face: v1, v2, v3. Zero-based (!)
    q : float, 3 x 2. For each of those vertices, x and y coord
    u : float,     2. Flow. x and y coord
    where: float   3. Vertice label (should be int)
    """
    ii = np.zeros(12, dtype=int)
    jj = np.zeros(12, dtype=int)
    ss = np.zeros(12)

    count = 0
    for i in range(3):  # python 0-based unlike matlab...
        ip = (i + 1) % 3
        ipp = (ip + 1) % 3
        unL = -( (q[ip, 1] + q[i, 1] - 2 * q[ipp, 1]) * u[0]
                -(q[ip, 0] + q[i, 0] - 2 * q[ipp, 0]) * u[1]) / 6
        ii[count],     ss[count]     = vv[i],   unL
        ii[count + 1], ss[count + 1] = vv[ip], -unL
        if unL > 0:
            jj[count]     = vv[i]
            jj[count + 1] = vv[i]
        else:
            jj[count]     = vv[ip]
            jj[count + 1] = vv[ip]
        count += 2
        if where[i] != 0 and where[ip] != 0:
            unL = +( (q[ip, 1] - q[i, 1]) * u[0]
                    -(q[ip, 0] - q[i, 0]) * u[1]) / 2
            if unL > 0:
                ii[count],     ss[count]     = vv[i], unL
                ii[count + 1], ss[count + 1] = vv[i], unL
                jj[count]     = vv[i]
                jj[count + 1] = vv[i]
                count += 2
    return ii, jj, ss


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


def getAB(mm, kk, bb, u, v, vv, ff, nu):
    """Mangle stuff, using the finite volume method. Return A, B.

    Inputs:
    mm, float, N_v x N_v. Diagonal. Mass.
    kk, float, N_v x N_v. Stiffness
    bb, float, N_b x N_v. Rih
    u,  float, N_f x   1. Flow (x-coord)
    v,  float, N_f x   1. Flow (y-coord)
    vv, float, N_v x   3. Vertices: x, y, label
    ff, int,   N_f x   4. Faces: v1, v2, v3, label. NOTE: 1-based
    nu, float,   1 x   1.

    Return A, B as dense ndarray matrices.
    """
    # check that non-zero values in bb are all 1
    if np.any(bb[bb.nonzero()] != 1):
        errMsg = "Some non-zero values in Rih/bb are not one. Fixing."
        logging.warn(errMsg)
        # raise RuntimeError(errMsg)
        # just fix it for now
        bb[bb.nonzero()] = 1.0

    N_f = ff.shape[0]  # number of faces/triangles
    N_v = vv.shape[0]  # number of vertices
    N_b = bb.shape[0]  # number of ?
    ii = np.empty(12 * N_f, dtype=int)
    jj = np.empty(12 * N_f, dtype=int)
    ss = np.empty(12 * N_f)

    for k in range(N_f):  # for each face
        # [ik,jk,sk] = fvm(ff(k,:),vv(ff(k,:),1:2),[u(k), v(k)], vv(ff(k,:),3));
        ff_idx = ff[k, 0:3] - 1  # python is zero-based, -1 to index into vv
        ik, jk, sk = fvm(ff_idx, vv[ff_idx, 0:2], [u[k], v[k]], vv[ff_idx, 2])
        ii[12 * k: 12 * k + 12] = ik
        jj[12 * k: 12 * k + 12] = jk
        ss[12 * k: 12 * k + 12] = sk

    # orig: A = sparse(ii,jj,ss,size(vv,1),size(vv,1),size(ii,1));  %    Nv x Nv
    # Note: we call it _ns here, not scaled by M yet
    A_ns = createDenseFromIJX(zip(ii, jj, ss), M=N_v, N=N_v, ij_onebased=False)

    # orig: M = sparse(mi,mj,ms,size(vv,1),size(vv,1));             %    Nv x Nv
    # orig: M2 = spdiags(sqrt(diag(M)),0,size(vv,1),size(vv,1));    %    Nv x Nv
    # M2 = np.diagflat(np.sqrt(np.diag(M)))                    # not needed here
    M2di = 1 / np.sqrt(np.diag(mm))  # vector of diagonal, sqrt, inverted

    # orig: K = sparse(ki,kj,ks,size(vv,1),size(vv,1));             %    Nv x Nv
    K = kk

    # orig: B = sparse(bj,bi,ones(length(bi),1),size(vv,1),max(bi)); %   Nv x Nb
    # create the normal way, then transpose.
    # Note: checked that all nz == 1 above.
    B_ns = bb.transpose()

    # orig: AA = A + nu * K;
    AA = A_ns + nu * K
    # Am = M2\AA/M2;  % this is A
    A = M2di * AA * M2di.reshape(-1, 1)

    # Bm = M2\B;      % this is B
    B = B_ns * M2di.reshape(-1, 1)  # multiply the rows

    return A, B


def writeABCD(A, B, C, D, destPath="", prefix = ""):
    ensurePath(destPath)
    fpre = prefix + "s2_"
    with open(os.path.join(destPath, fpre+"A.txt"), 'w') as f:
        write_csr_matrix(A, f)
    with open(os.path.join(destPath, fpre+"B.txt"), 'w') as f:
        write_csr_matrix(B, f)
    with open(os.path.join(destPath, fpre+"C.txt"), 'w') as f:
        write_csr_matrix(C, f)
    with open(os.path.join(destPath, fpre+"D.txt"), 'w') as f:
        write_csr_matrix(D, f)


def ensurePath(path):
    """Create path if it does not exist."""
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise


def runS1(sourcePath=None, prefix=None, pathFF=None):
    if sourcePath is None:
        sourcePath = "stokes.edp"
        logging.warn("stage 1: no sourcePath provided, assuming %s", sourcePath)
    if pathFF is None:
        pathFF = getPathFF()

    # logic here:
    # if sourcePath is a directory, not a file: add stokes.edp
    if os.path.isdir(sourcePath):
        sourcePath = os.path.join(sourcePath, "stokes.edp")
        logging.warn("stage 1: sourcePath only directory, assuming %s", sourcePath)

    # now it's a file
    # if the file is not there, abort
    if not os.path.isfile(sourcePath):
        errMsg = "stage 1: sourcePath not found: %s" % sourcePath
        logging.error(errMsg)
        raise RuntimeError(errMsg)

    # if the file does not end with edp, warn
    fileExtension = os.path.splitext(sourcePath)[1]
    if fileExtension != ".edp":
        logging.warn("stage 1: sourcePath does not end with .edp: %s", sourcePath)

    # if prefix given and stokes file does not start with prefix+"s0_", copy to there
    # then use it
    if prefix is not None:
        fpre = prefix + "s0_"
        dirname, basename = os.path.split(sourcePath)
        if basename.startswith(fpre):
            pass # nothing to do
        else:
            newBasename = fpre + basename
            oldSourcePath = sourcePath
            sourcePath = os.path.join(dirname, newBasename)
            # copy over
            logging.info("stage 1: copying from %s to %s", oldSourcePath, sourcePath)
            shutil.copyfile(oldSourcePath, sourcePath)

    # execute FreeFem++-nw sp
    logging.info("Executing FreeFem++ at %s with file %s", pathFF, sourcePath)
    logLevelFF = 0
    if logging.getLogger().getEffectiveLevel() < 20:
        logLevelFF = 2

    subprocess.check_call([pathFF, "-v", str(logLevelFF), sourcePath])

    # maybe check that expected files are there?
    pass


def runS2(sourcePath="..\\run1\\", destPath=None, prefix="", nu=0.01):
    if destPath is None:
        destPath = sourcePath # os.path.join(sourcePath, "data")
    fpre = prefix + "s1_" # full prefix

    # read intermediate output from FF
    mm = with_file(os.path.join(sourcePath, fpre+'mass.txt'), readSparseFF)
    kk = with_file(os.path.join(sourcePath, fpre+'stiff.txt'), readSparseFF)
    bb = with_file(os.path.join(sourcePath, fpre+'Rih.txt'), readSparseFF)
    u  = with_file(os.path.join(sourcePath, fpre+'u.txt'), readVectFF)
    v  = with_file(os.path.join(sourcePath, fpre+'v.txt'), readVectFF)
    vv, ff, ee = with_file(os.path.join(sourcePath, fpre+'stokes.msh'), readMeshFF)

    # compute A, B using fvm
    A, B = getAB(mm, kk, bb, u, v, vv, ff, nu)
    # NOTE hardcoded parameters: C, D
    C = np.zeros((1, 736))
    C[0, 699] = 2.4325069738165400e+01
    D = C

    # write next intermediate files
    writeABCD(A, B, C, D, destPath=destPath, prefix=prefix)
