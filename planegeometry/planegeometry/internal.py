import numpy as np

def distance_(A, B):
    return np.linalg.norm(B-A)

def vlength_(v):
    return np.linalg.norm(v)

def dot_(u, v=None):
    if v is None:
        v = u
    return np.vdot(u, v)
