from planegeometry.geometry import Inversion, Circle, Line, mid_circles, Reflection

def test_mid_circles():
    # two congruent circles which do not intersect each other
    circ1 = Circle(center = (1,2), radius = 3)
    circ2 = Circle(center = (14,15), radius = 3)
    mc = mid_circles(circ1, circ2)
    assert isinstance(mc, Line)
    R = Reflection(mc)
    assert R.reflect_circle(circ1).is_equal(circ2)
    assert R.reflect_circle(circ2).is_equal(circ1)
    # two congruent circles which intersect each other
    circ1 = Circle(center = (1,2), radius = 3)
    circ2 = Circle(center = (2,2), radius = 3)
    mc = mid_circles(circ1, circ2)
    assert isinstance(mc, list)
    R = Reflection(mc[1])
    assert R.reflect_circle(circ1).is_equal(circ2)
    assert R.reflect_circle(circ2).is_equal(circ1)
    iota = Inversion.on_circle(mc[0])
    assert iota.invert_circle(circ1).is_equal(circ2)
    assert iota.invert_circle(circ2).is_equal(circ1)
    # two tangent congruent circles
    circ1 = Circle(center = (1,2), radius = 3)
    circ2 = Circle(center = (7,2), radius = 3)
    mc = mid_circles(circ1, circ2)
    assert isinstance(mc, Line)
    R = Reflection(mc)
    assert R.reflect_circle(circ1).is_equal(circ2)
    assert R.reflect_circle(circ2).is_equal(circ1)
    # two non-intersecting circles outside each other
    circ1 = Circle((0,2), 3)
    circ2 = Circle((4,5), 1)
    mc = mid_circles(circ1, circ2)
    assert isinstance(mc, Circle)
    iota = Inversion.on_circle(mc)
    assert iota.invert_circle(circ1).is_equal(circ2)
    assert iota.invert_circle(circ2).is_equal(circ1)
    # two non-intersecting circles, one inside the other
    circ1 = Circle((0,0), 3)
    circ2 = Circle((0.7,1), 0.7)
    mc = mid_circles(circ1, circ2)
    assert isinstance(mc, Circle)
    iota = Inversion.on_circle(mc)
    assert iota.invert_circle(circ1).is_equal(circ2)
    assert iota.invert_circle(circ2).is_equal(circ1)
    # two intersecting circles
    circ1 = Circle((5,4), 2)
    circ2 = Circle((8,5), 3)
    mc = mid_circles(circ1, circ2)
    assert isinstance(mc, list)
    iota = Inversion.on_circle(mc[0])
    assert iota.invert_circle(circ1).is_equal(circ2)
    assert iota.invert_circle(circ2).is_equal(circ1)
    iota = Inversion.on_circle(mc[1])
    assert iota.invert_circle(circ1).is_equal(circ2)
    assert iota.invert_circle(circ2).is_equal(circ1)
    # two tangent circles, externally
    circ1 = Circle((5,4), 2)
    circ2 = Circle((8,4), 1)
    mc = mid_circles(circ1, circ2)
    assert isinstance(mc, Circle)
    iota = Inversion.on_circle(mc)
    assert iota.invert_circle(circ1).is_equal(circ2)
    assert iota.invert_circle(circ2).is_equal(circ1)
    # two tangent circles, internally
    circ1 = Circle((5,4), 2)
    circ2 = Circle((6,4), 1)
    mc = mid_circles(circ1, circ2)
    assert isinstance(mc, Circle)
    iota = Inversion.on_circle(mc)
    assert iota.invert_circle(circ1).is_equal(circ2)
    assert iota.invert_circle(circ2).is_equal(circ1)
