# -*- coding: utf-8 -*-
from planegeometry.geometry import Mobius, Circle, Triangle, Line
from numpy import allclose, array_equal

def test_image_circle():
    circ = Circle((2, 3), 4)
    # a case with c = 0
    Mob = Mobius([[1+1j, 2], [0, 3-2j]])
    circ1 = Mob.transform_circle(circ)
    P = Mob.transform(circ.point_from_angle(0))
    Q = Mob.transform(circ.point_from_angle(60))
    R = Mob.transform(circ.point_from_angle(120))
    circ2 = Triangle(P, Q, R).circumcircle()
    assert circ1.is_equal(circ2)
    # a case with c != 0
    Mob = Mobius([[1+1j, 2], [2-3j, 3-2j]])
    circ1 = Mob.transform_circle(circ)
    P = Mob.transform(circ.point_from_angle(0))
    Q = Mob.transform(circ.point_from_angle(60))
    R = Mob.transform(circ.point_from_angle(120))
    circ2 = Triangle(P, Q, R).circumcircle()
    assert circ1.is_equal(circ2)
    # a case with -d/c on the circle
    Mob = Mobius([[1+1j, 2], [1, -6-3j]])
    line1 = Mob.transform_circle(circ)
    P = Mob.transform(circ.point_from_angle(30))
    Q = Mob.transform(circ.point_from_angle(60))
    line2 = Line(P, Q)
    assert line1.is_equal(line2)

def test_image_line():
    line = Line((2,3), (1,5))
    # in case c=0, the image is a line
    Mob = Mobius([[1+1j, 2], [0, 3-2j]])
    assert isinstance(Mob.transform_line(line), Line)
    # case of a circle image
    Mob = Mobius([[1+1j, 2], [2-3j, 3-2j]])
    circ1 = Mob.transform_line(line)
    P = Mob.transform(line.A)
    Q = Mob.transform(line.B)
    R = Mob.transform((line.A + line.B)/2)
    circ2 = Triangle(P, Q, R).circumcircle()
    assert circ1.is_equal(circ2)
    
def test_mapping_three_points():
    P1 = (1, 2)
    P2 = (3, 4)
    P3 = (7, 7)
    Q1 = (2, 1)
    Q2 = (4, 4)
    Q3 = (-7, 7)
    Mob = Mobius.from_mapping_three_points(P1, P2, P3, Q1, Q2, Q3)
    R1 = Mob.transform(P1)
    R2 = Mob.transform(P2)
    R3 = Mob.transform(P3)
    assert allclose(Q1, R1)
    assert allclose(Q2, R2)
    assert allclose(Q3, R3)


def test_generalized_power():
    # case 1 : diag(c(lambda,lambda))
    Mat = [[1+2j, 0], [0, 1+2j]]
    Mob = Mobius(Mat)
    pow2 = Mob.power(2).M
    gpow2 = Mob.gpower(2).M
    assert allclose(pow2, gpow2)
    powminus2 = Mob.power(-2).M
    gpowminus2 = Mob.gpower(-2).M
    f = powminus2[0, 0] / gpowminus2[0, 0]
    assert allclose(powminus2, f * gpowminus2)
    # case 2: tr(Mat)^2 = 4*det(Mat)
    Mat = [[1, 3j], [0, 1]]
    Mob = Mobius(Mat)
    pow2 = Mob.power(2).M
    gpow2 = Mob.gpower(2).M
    assert allclose(pow2, gpow2)
    powminus2 = Mob.power(-2).M
    gpowminus2 = Mob.gpower(-2).M
    assert allclose(powminus2, gpowminus2)
    # case 3: two distinct eigenvalues
    Mat = [[1, 3j], [1+2j, 1]]
    Mob = Mobius(Mat)
    pow2 = Mob.power(2).M
    gpow2 = Mob.gpower(2).M
    assert allclose(pow2, gpow2)
    powminus2 = Mob.power(-2).M
    gpowminus2 = Mob.gpower(-2).M
    f = powminus2[0, 0] / gpowminus2[0, 0]
    assert allclose(powminus2, f * gpowminus2)

