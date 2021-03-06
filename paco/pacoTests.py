# import logging
# import numpy as np
import logging
import numpy as np
import os.path
import pacoUtils as pu

def checkS0():
    pass

def checkS1(sourceDir=".", prefix):
    pass

def checkS2(sourceDir=".", prefix):
    A, B, C, D = pu.readABC(sourceDir=sourceDir, prefix=prefix)
    return checkConsistency(A, B, C, D)

def checkConsistency(A, B, C, D):
    msg = "" # no problems
    if A.ndim != 2:
        msg = "A must have 2 dimensions."
    else if A.shape[0] != A.shape[1]:
        msg = "A must be square."
    else if A.shape[0] != B.shape[1]:
        msg = "A and B must have same number of rows."
    else if A.shape[1] != C.shape[1] or A.shape[1] != D.shape[1]:
        msg = "A, C and D must have the same number of columns."

    if len(msg) > 0:
        logging.error(msg)
        return False
    else
        logging.debug("Matrices consistent.")
        return True


def compareMatrices(A, B):
    if A.shape != B.shape:
        return (1000000, "Different shape.")
    if np.array_equal(A, B):
        return (0, "Arrays equal.")
    diff = np.linalg.norm(A-B)
    if np.allclose(A,B):
        return(diff, "Arrays different, but close.")
    return(diff, "Arrays different.")


def compareS0(sourceFullPath1, sourceFullPath2):
    # compare stokes.edp
    pass

def compareS1(sourceDir1, prefix1, sourceDir2, prefix2):
    # compare u, v, mass, stiff, Rih, stokes.msh
    if sourceDir2 is None:
        sourceDir2 = sourceDir1
    basenames = ['mass', 'stiff', 'Rih', 'u', 'v']
    readFuncs = [readSparseFF, readSparseFF, readSparseFF, readVectFF, readVectFF]
    fpre1 = prefix1 + "s1_"
    fpre2 = prefix2 + "s1_"

    for basename, readFunc in zip(basenames, readFuncs):
        obj1 = with_file(os.path.join(sourceDir1, fpre1+basename), readFunc)
        obj2 = with_file(os.path.join(sourceDir2, fpre2+basename), readFunc)
        diff, res = compareMatrices(obj1, obj2)
        d += diff
        msg = "For %s: %s. Diff = %s" % (basename, res, diff)
        print(msg)
        logging.info(msg)

    vv1, ff1, ee1 = with_file(os.path.join(sourceDir1, fpre1+'stokes.msh'), readMesFF)
    vv2, ff2, ee2 = with_file(os.path.join(sourceDir2, fpre2+'stokes.msh'), readMesFF)
    for basename, obj1, obj2 in [   ("vv",vv1,vv2),
                                    ("ff",ff1,ff2),
                                    ("ee", ee1, ee2)]:
        diff, res = compareMatrices(obj1, obj2)
        d += diff
        msg = "For %s: %s. Diff = %s" % (basename, res, diff)
        print(msg)

    msg = "Total diff: %s." % ( d )
    print(msg)
    logging.info(msg)
    return d

def compareS2(sourceDir1, prefix1, sourceDir2, prefix2):
    # compare A, B, C, D
    if sourceDir2 is None:
        sourceDir2 = sourceDir1
    basenames = ["A", "B", "C", "D"]
    ext = ".txt"
    readFunc = pu.read_csr_matrix
    fpre1 = prefix1 + "s2_"
    fpre2 = prefix2 + "s2_"

    msg = "Compare Stage 2: %s and %s." % ( os.path.join(sourceDir1, fpre1),
                                            os.path.join(sourceDir2, fpre2))
    print(msg)
    logging.info(msg)

    d = 0.0
    for basename in basenames:
        obj1 = pu.with_file(os.path.join(sourceDir1, fpre1+basename+ext), readFunc)
        obj2 = pu.with_file(os.path.join(sourceDir2, fpre2+basename+ext), readFunc)
        diff, res = compareMatrices(obj1, obj2)
        d += diff
        msg = "For %s: %s. Diff = %s" % (basename, res, diff)
        print(msg)
        logging.info(msg)
    msg = "Total diff: %s." % ( d )
    print(msg)
    logging.info(msg)
    return d

def compareS3(prefix1, prefix2, sourceDir1=".", sourceDir2=None):
    if sourceDir2 is None:
        sourceDir2 = sourceDir1
