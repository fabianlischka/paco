import logging
import numpy as np
import os
import os.path
import pacoUtils as pu
import subprocess
import yaml



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
    A_ns = pu.createDenseFromIJX(zip(ii, jj, ss), M=N_v, N=N_v, ij_onebased=False)

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


def runStage1(params):
    """Run given file through FreeFem++ to obtain s1 files.

    Prefixes the given .edp file with given prefix + "s0_",
        unless it's of that form already.
    Then copies to destDir.
    Then, runs FreeFem++ on it.
    """

    pathFF = pu.getPathFF(params['Config']['DirsBin'])
    destDir = params['Config']['DirData']
    prefix  = params['Config']['Prefix']
    # find the script
    scriptName = params['runStage1']['Script']
    scriptPath = None
    scriptDirs = [params['Config']['DirConfig'], params['Config']['DirCurrent'],
                 params['Config']['DirData'] ]
    for scriptDir in scriptDirs:
        scriptPath = os.path.join(scriptDir, scriptName)
        if os.path.isfile(scriptPath):
            break
    if scriptPath is None:
        msg = "Stage 1: could not find script %s in directories %s." % (
            scriptName, ", ".join(scriptDirs) )
        print(msg)
        logging.critical(msg)
        raise RuntimeError(msg)

    # if the file does not end with edp, warn
    fileExtension = os.path.splitext(scriptPath)[1]
    if fileExtension != ".edp":
        logging.warn("stage 1: scriptPath does not end with .edp: %s", scriptPath)

    # determine destination directory, copy source file over if necessary
    newBasename = pu.copyIfNecessary(scriptPath, destDir, prefix + "s0_")

    logLevelFF = 0
    if logging.getLogger().getEffectiveLevel() < 20:
        logLevelFF = 2

    # execute FreeFem++-nw sp
    logging.info("Executing FreeFem++ (at %s) in dir %s on file %s",
                    pathFF, destDir, newBasename)
    subprocess.check_call([pathFF, "-v", str(logLevelFF), newBasename],
                          cwd=destDir)

    # maybe check that expected files are there?
    pass


def runStage2(params):
    sourceDir = params['Config']['DirData']
    prefix = params['Config']['Prefix']
    destDir = sourceDir
    fpre = prefix + "s1_" # full prefix

    # read intermediate output from FF
    mm = pu.with_file(os.path.join(sourceDir, fpre+'mass.txt'), pu.readSparseFF)
    kk = pu.with_file(os.path.join(sourceDir, fpre+'stiff.txt'), pu.readSparseFF)
    bb = pu.with_file(os.path.join(sourceDir, fpre+'Rih.txt'), pu.readSparseFF)
    u  = pu.with_file(os.path.join(sourceDir, fpre+'u.txt'), pu.readVectFF)
    v  = pu.with_file(os.path.join(sourceDir, fpre+'v.txt'), pu.readVectFF)
    vv, ff, ee = pu.with_file(os.path.join(sourceDir, fpre+'stokes.msh'), pu.readMeshFF)

    nu = params['runStage2']['nu']
    # compute A, B using fvm
    A, B = getAB(mm, kk, bb, u, v, vv, ff, nu)
    C = pu.matrixFromConfig(params['runStage2']['C'], A.shape[1])
    D = pu.matrixFromConfig(params['runStage2']['D'], A.shape[1])

    # write next intermediate files
    pu.writeABCD(A, B, C, D, destDir=destDir, prefix=prefix)

def runStage3(params):
    destDir = params['Config']['DirData']
    binDir = params['Config']['DirPacoBin']
    prefix = params['Config']['Prefix']

    # Generate params.txt
    pu.writeParamsFile(params['runStage3']['Params'],
                       params['Config']['DirPacoRoot'],
                       destDir=destDir, prefix=prefix)

    # Generate y_hat if necessary
    if 'yhat' in params['runStage3']['Params']:
        yhat = np.array(params['runStage3']['Params']['yhat'])
        fn = os.path.join(destDir, prefix+"s2_yhat.txt")
        yhat.tofile(fn, sep="\n", format="%r")

    # execute targets
    for subStage in params['runStage3']['SubStages']:
        funcToCall = subStageMethods[subStage['Method']]
        binaryFullPath = os.path.join(binDir, subStage['Target'])
        logging.info("Executing %s in dir %s with prefix %s.",
                    binaryFullPath, destDir, prefix)
        funcToCall(params, binaryFullPath, destDir, prefix)

def subStageDirect(params, binaryFullPath, destDir, prefix):
    subprocess.check_call([binaryFullPath, destDir, prefix],
                          cwd=destDir)

def subStageBatch(params, binaryFullPath, destDir, prefix):
    # Substitute:
    # Email: ${Email}
    # BinaryFullPath: ${BinaryFullPath}
    # DirData: ${DirData}
    # Prefix: ${Prefix}

    substitutions = params['Config'] # contains Email, DirData, Prefix
    substitutions['BinaryFullPath'] = binaryFullPath
    outPath = pu.writeBatchFile( substitutions,
                       params['Config']['DirPacoRoot'],
                       destDir=destDir, prefix=prefix )
    subprocess.check_call(['qsub', outPath])

subStageMethods = { "direct" : subStageDirect,
                    "batch"  : subStageBatch }
