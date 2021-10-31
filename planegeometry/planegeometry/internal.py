from math import cos, sin
import numpy as np

epsilon_ = np.finfo(float).eps

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

def unit_vector_(beta):
    return np.array([cos(beta), sin(beta)])

def ellipse_points_(t, O, a, b, alpha):
    x = a * np.cos(t)
    y = b * np.sin(t)
    cosalpha = cos(alpha)
    sinalpha = sin(alpha)
    return np.column_stack(
        (
         O[0] + cosalpha*x - sinalpha*y,
         O[1] + sinalpha*x + cosalpha*y
        )
    )

def circle_points_(t, O, r):
    return np.column_stack(
        (
         O[0] + r * np.cos(t),
         O[1] + r * np.sin(t)
        )
    )
  
def collinear_(A, B, C, tol = 0):
    notdistinct = np.allclose(A, B) or np.allclose(A, C) or np.allclose(B, C)
    if notdistinct:
        return True
    AB = B - A
    AC = C - A
    z1 = complex(*AB)
    z2 = complex(*AC)
    z = z1.conjugate() * z2
    re = z.real
    re2 = re*re
    im = z.imag
    return re2 / (re2 + im*im) >= 1 - tol

def circle_as_ellipse_(C):
    r = C.radius
    return Ellipse(C.center, r, r, 0)