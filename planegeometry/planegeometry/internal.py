import numpy as np

def distance_(A, B):
    return np.linalg.norm(B-A)

def vlength_(v):
    return np.linalg.norm(v)

def dot_(u, v=None):
    if v is None:
        v = u
    return np.vdot(u, v)

def det2x2_(row1, row2):
    return row1[0]*row2[1] - row2[0]*row1[1]

def line_line_intersection_(P1, P2, Q1, Q2):
    dx1 = P1[0] - P2[0]
    dx2 = Q1[0] - Q2[0]
    dy1 = P1[1] - P2[1]
    dy2 = Q1[1] - Q2[1]
    D = det2x2_((dx1, dy1), (dx2, dy2))
    if D == 0:
        return None
    D1 = det2x2_(P1, P2)
    D2 = det2x2_(Q1, Q2)
    return np.array([
      det2x2_((D1, dx1), (D2, dx2)),
      det2x2_((D1, dy1), (D2, dy2))
    ]) / D
