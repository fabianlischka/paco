# import logging
# import numpy as np
import os.path




def compareS0(sourceFullPath1, sourceFullPath2):
    # compare stokes.edp
    pass

def compareS1(prefix1, prefix2, sourceDir1=".", sourceDir2=None):
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
        if obj1.shape == obj2.shape:
            # look at np.array_equal(A,B), np.allclose(A,B)
            pass

    vv1, ff1, ee1 = with_file(os.path.join(sourceDir1, fpre1+'stokes.msh'), readMesFF)
    vv2, ff2, ee2 = with_file(os.path.join(sourceDir2, fpre2+'stokes.msh'), readMesFF)

def compareS2(prefix1, prefix2, sourceDir1=".", sourceDir2=None):
    # compare A, B, C, D
    if sourceDir2 is None:
        sourceDir2 = sourceDir1


def compareS3(prefix1, prefix2, sourceDir1=".", sourceDir2=None):
    if sourceDir2 is None:
        sourceDir2 = sourceDir1
