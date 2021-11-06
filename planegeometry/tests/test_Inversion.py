from planegeometry.geometry import Inversion, Circle
import numpy as np

def test_invert_circle():
    circ0 = Circle((0,2), 3)
    iota = Inversion((5,5), 6)
    circ1 = iota.invert_circle(circ0)
    A = iota.invert((3,2))
    B = iota.invert((0,5))
    C = iota.invert((-3,2))
    assert circ1.includes(A)
    assert circ1.includes(B)
    assert circ1.includes(C)

def test_swapping_two_circles():
    ## non-intersecting circles, external
    circ1 = Circle((0,2), 3)
    circ2 = Circle((4,5), 1)
    #
    iota = Inversion.from_swapping_two_circles(circ1, circ2, True)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    #
    iota = Inversion.from_swapping_two_circles(circ1, circ2, False)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    #
    iota = Inversion.from_swapping_two_circles(circ2, circ1, True)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    #
    iota = Inversion.from_swapping_two_circles(circ2, circ1, False)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    ## non-intersecting circles, internal
    circ1 = Circle((0,0), 3)
    circ2 = Circle((0.7,1), 0.7)
    #
    iota = Inversion.from_swapping_two_circles(circ1, circ2, True)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    #
    iota = Inversion.from_swapping_two_circles(circ1, circ2, False)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    #
    iota = Inversion.from_swapping_two_circles(circ2, circ1, True)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    #
    iota = Inversion.from_swapping_two_circles(circ2, circ1, False)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    ## intersecting circles
    circ1 = Circle((5,4), 2)
    circ2 = Circle((8,5), 3)
    #
    iota = Inversion.from_swapping_two_circles(circ1, circ2, True)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    #
    iota = Inversion.from_swapping_two_circles(circ1, circ2, False)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    #
    iota = Inversion.from_swapping_two_circles(circ2, circ1, True)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    #
    iota = Inversion.from_swapping_two_circles(circ2, circ1, False)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    ## tangent circles
    circ1 = Circle((5,4), 2)
    circ2 = Circle((8,4), 1)
    #
    iota = Inversion.from_swapping_two_circles(circ1, circ2, True)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    #
    iota = Inversion.from_swapping_two_circles(circ1, circ2, False)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    #
    iota = Inversion.from_swapping_two_circles(circ2, circ1, True)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))
    #
    iota = Inversion.from_swapping_two_circles(circ2, circ1, False)
    assert circ1.is_equal(iota.invert_circle(circ2))
    assert circ2.is_equal(iota.invert_circle(circ1))

def test_fixing_two_circles():
    circ1 = Circle((1, 0), 3)
    circ2 = Circle((2, 2), 1)
    iota = Inversion.from_fixing_two_circles(circ1, circ2)
    newcirc1 = iota.invert_circle(circ1)
    newcirc2 = iota.invert_circle(circ2)
    assert circ1.is_equal(newcirc1)
    assert circ2.is_equal(newcirc2)

def test_fixing_three_circles():
    circ1 = Circle((1, 0), 3)
    circ2 = Circle((2, 2), 1)
    circ3 = Circle((3, 5), 2)
    iota = Inversion.from_fixing_three_circles(circ1, circ2, circ3)
    newcirc1 = iota.invert_circle(circ1)
    newcirc2 = iota.invert_circle(circ2)
    newcirc3 = iota.invert_circle(circ3)
    assert circ1.is_equal(newcirc1)
    assert circ2.is_equal(newcirc2)
    assert circ3.is_equal(newcirc3)
    
def test_compose_inversions():
    iota1 = Inversion((1, 1), 2)
    iota2 = Inversion((3, 2), 4)
    M = (4, 5) # np.array([4, 5], dtype=float)
    P = iota1.invert(iota2.invert(M))
    Mob = iota2.compose(iota1)
    Q = Mob.transform(M)
    assert np.allclose(P, Q)
    # with a negative power
    iota2 = Inversion((3, 2), -4)
    M = (4, 5) # np.array([4, 5], dtype=float)
    P = iota1.invert(iota2.invert(M))
    Mob = iota2.compose(iota1)
    Q = Mob.transform(M)
    assert np.allclose(P, Q)
    
