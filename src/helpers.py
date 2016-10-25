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
    """Reads a vector in the form output by FreeFem++.

    First line contains the number of entries: nbcoef.
    Next lines contain the entries, with 5 entries per line,
    separated by tab.
    Entries are read left to right, then line by line.
    """

    nbcoef = int(f.next().strip())
    logging.info("attempt to read %s values", nbcoef)
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
        errMsg = "readVectFF: Inconsistent dimension, expected %d, got %d" %(nbcoef, c)
        logging.error(errMsg)
        raise RuntimeError(errMsg)
    logging.info("successfully read %s values", nbcoef)
    return arr

def createDenseFromIJX(iijjxx, M = -1, N = -1, issym = False, ij_onebased = True):
    """Creates a dense matrix from (nnz x 3) array with i, j, x
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
    res = np.zeros((M,N))
    for ifloat, jfloat, x in iijjxx:
        i, j = int(ifloat) - offset, int(jfloat) - offset  # python is 0-based
        res[i,j] = x
        if issym:
            res[j,i] = x
    return res

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
    res = createDenseFromIJX( ijv, M = n, N = m,  # note weird FF convention
                              issym = (issym!=0), ij_onebased = True)
    # res = np.zeros((n,m))
    # for ifloat,jfloat,v in ijv:
    #     i, j = int(ifloat) - 1, int(jfloat) - 1  # python is 0-based
    #     res[i,j] = v
    #     if issym != 0:
    #         res[j,i] = v

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

def with_file(filename, callback):
    """Call `callback` with the opened file `filename`. Return result.

    See http://softwareengineering.stackexchange.com/questions/262346/should-i-pass-in-filenames-to-be-opened-or-open-files
    """
    with open(filename, 'r') as f:
        return callback(f)

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

def write_csr_matrix(A, f):
    """
    Write dense A in compressed sparse row format (zero-based) to file f

    Text file format:
    - one number per line
    - First number is m = number of rows
    - Second number is n = number of columns
    - Next m+1 numbers are the row pointers
    - Next nnz numbers are the column indices
    - Next nnz numbers are the matrix values
    """
    m, n = A.shape
    f.write('{}\n'.format(m))
    f.write('{}\n'.format(n))

    counter = 0
    for i in range(m):
        f.write('{}\n'.format( counter ))
    #    counter = counter + nnz(A(i,:));
        counter += len( A[i].nonzero()[0] )
    f.write('{}\n'.format( counter ))

    counter = 0
    for i in range(m):
        for k in A[i].nonzero()[0]:
    #        fprintf(fid,'%d\n',k-1)
            f.write('{}\n'.format( k ))

    for i in range(m):
        for k in A[i].nonzero()[0]:
    #         fprintf(fid,'%.16e\n',full(A(i,k)));
            f.write('{}\n'.format( A[i,k] ))



def getA():  # FIXFIXFIX
    base = "..\\run1\\"

    mm = with_file(base+'mass.txt', readSparseFF)
    kk = with_file(base+'stiff.txt', readSparseFF)
    bb = with_file(base+'Rih.txt', readSparseFF)
    u  = with_file(base+'u.txt', readVectFF)
    v  = with_file(base+'v.txt', readVectFF)
    vv, ff, ee = with_file(base+'stokes.msh', readMeshFF)

    # check that non-zero values in bb are all 1
    if np.any(bb[bb.nonzero()] != 1):
        errMsg = "Some non-zero values in Rih/bb are not one. Fixing."
        logging.warn(errMsg)
        # raise RuntimeError(errMsg)
        # just fix it for now
        bb[bb.nonzero()] = 1.0

    nu = 0.01

    N_f = ff.shape[0]  # number of faces/triangles
    N_v = vv.shape[0]  # number of vertices
    N_b = bb.shape[0]  # number of ?
    ii = np.empty(12*N_f, dtype = int)
    jj = np.empty(12*N_f, dtype = int)
    ss = np.empty(12*N_f)

    for k in range(N_f):
        # [ik,jk,sk] = fvm(ff(k,:),vv(ff(k,:),1:2),[u(k), v(k)], vv(ff(k,:),3));
        ff_idx = ff[k,0:3] - 1  # python is zero-based
        ik, jk, sk = fvm( ff_idx, vv[ ff_idx, 0:2], [u[k], v[k]], vv[ ff_idx, 2] )
        ii[12*k : 12*k+12] = ik
        jj[12*k : 12*k+12] = jk
        ss[12*k : 12*k+12] = sk

    # A = sparse(ii,jj,ss,size(vv,1),size(vv,1),size(ii,1));   % A  Nv x Nv
    A = createDenseFromIJX(zip(ii,jj,ss), M = N_v, N = N_v, ij_onebased = True)
    # M = sparse(mi,mj,ms,size(vv,1),size(vv,1));              %    Nv x Nv
    M = mm
    # M2 = spdiags(sqrt(diag(M)),0,size(vv,1),size(vv,1));     %    Nv x Nv
    M2 = np.diagflat(np.sqrt(np.diag(M)))
    M2di = 1/np.sqrt(np.diag(M)) # vector
    # K = sparse(ki,kj,ks,size(vv,1),size(vv,1));              %    Nv x Nv
    K = kk
    # B = sparse(bj,bi,ones(length(bi),1),size(vv,1),max(bi)); % B  Nv x Nb
    # create the normal way, then transpose. Note: checked that all nz == 1 above.
    Bt = bb
    B = Bt.transpose()


    # AA = A + nu * K;
    AA = A + nu*K
    # Am = M2\AA/M2;  % this is A
    Am = M2di * AA * M2di.reshape(-1,1)

    # Bm = M2\B;      % this is B
    Bm = B * M2di.reshape(-1,1)  # multiply the rows
