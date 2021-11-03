from math import cos, sin
import numpy as np
from fractions import Fraction as Fr
import numbers

def is_point_(P):
    return (isinstance(P, tuple) or isinstance(P, list) or isinstance(P, np.ndarray)) and len(P) == 2

def error_if_not_point_(**kwargs):
    key = list(kwargs.keys())[0]
    P = kwargs[key]
    if not is_point_(P):
        # frame = currentframe()
        # P = getargvalues(frame)#.locals["P"]
        raise ValueError("`%s` is not a point." % key)
    return

def is_number_(x):
    return isinstance(x, numbers.Number)

def error_if_not_number_(**kwargs):
    key = list(kwargs.keys())[0]
    x = kwargs[key]
    if not is_number_(x):
        raise ValueError("`%s` is not a number." % key)
    return

def error_if_not_positive_(**kwargs):
    key = list(kwargs.keys())[0]
    x = kwargs[key]
    if x <= 0:
        raise ValueError("`%s` is not positive." % key)
    return

def farey_(n):
    return [Fr(0, 1)] + sorted(
        {Fr(m, k) for k in range(1, n+1) for m in range(1, k+1)}
    )

def farey_stack_(n):
    fractions = farey_(n)
    lists = [[f.denominator, f.numerator] for f in fractions]
    return [lists[i] + lists[i+1] for i in range(len(lists)-1)]

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

def det2x2_mat_(M):
    return det2x2_(M[0, :], M[1, :])

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

def from_complex_(z):
    return np.array([z.real, z.imag])

def mod2_(z):
    re = z.real
    im = z.imag
    return re*re + im*im