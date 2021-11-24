# -*- coding: utf-8 -*-
from planegeometry.geometry import Mobius, Circle, cross_ratio
from planegeometry.internal import approx_equal_
import numpy as np

def test_cross_ratio_real():
    circ = Circle((2, 3), 4)
    A = circ.point_from_angle(10, degrees=True)
    B = circ.point_from_angle(90, degrees=True)
    C = circ.point_from_angle(180, degrees=True)
    D = circ.point_from_angle(270, degrees=True)
    cr = cross_ratio(A, B, C, D)
    assert approx_equal_(cr.imag, 0)
    
def test_cross_ratio_Mobius_invariant():
    np.random.seed(666)
    Pts = np.random.rand(4, 2)
    cr = cross_ratio(*Pts)
    A, B, C, D = [*Pts]
    mob = Mobius([[2 + 3j, 1j], [-1j, 2]])
    cr2 = cross_ratio(
        mob.transform(A),    
        mob.transform(B),    
        mob.transform(C),    
        mob.transform(D)
    )
    assert approx_equal_(cr.real, cr2.real) and approx_equal_(cr.imag, cr2.imag)
