# -*- coding: utf-8 -*-
from planegeometry.geometry import Triangle, intersection_circle_circle
from planegeometry.internal import is_point_
from numpy import random, allclose

def test_contains():
    t = Triangle((0,0), (1,5), (3,3))
    random.seed(666)
    pts = t.random_points(3, "in")
    assert t.contains(pts[0, :])
    assert t.contains(pts[1, :])
    assert t.contains(pts[2, :])

def test_malfatti_tangent():
    t = Triangle((0,0), (1,5), (3,3))
    Mcircles = t.malfatti_circles()
    C1, C2, C3 = Mcircles["circles"].values()
    I1 = intersection_circle_circle(C1, C2)
    assert is_point_(I1)
    I2 = intersection_circle_circle(C1, C3)
    assert is_point_(I2)
    I3 = intersection_circle_circle(C2, C3)
    assert is_point_(I3)
    TA, TB, TC = Mcircles["tangency_points"].values()
    assert allclose(TA, I3)
    assert allclose(TB, I2)
    assert allclose(TC, I1)
    