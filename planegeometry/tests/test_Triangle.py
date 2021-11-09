# -*- coding: utf-8 -*-
from planegeometry.geometry import Triangle, intersection_circle_circle, Line
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
    
def test_excircles():
    t = Triangle((0,0), (1,5), (3,3))
    AB = Line(t.A, t.B)
    AC = Line(t.A, t.C)
    BC = Line(t.B, t.C)
    excircles = t.excircles()
    IAB = excircles["A"].intersection_with_line(AB)
    assert is_point_(IAB)
    IAC = excircles["A"].intersection_with_line(AC)
    assert is_point_(IAC)
    IBA = excircles["B"].intersection_with_line(AB)
    assert is_point_(IBA)
    IBC = excircles["B"].intersection_with_line(BC)
    assert is_point_(IBC)
    ICA = excircles["C"].intersection_with_line(AC)
    assert is_point_(ICA)
    ICB = excircles["C"].intersection_with_line(BC)
    assert is_point_(ICB)