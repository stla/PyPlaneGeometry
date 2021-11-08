from planegeometry.geometry import Circle
from math import isclose, pi

def test_orthogonal():
    circ1 = Circle((0, 0), 3)
    circ2 = circ1.orthogonalThroughTwoPointsOnCircle(0.5, 1, arc = False)
    assert circ1.is_orthogonal(circ2)
    assert isclose(circ1.angle(circ2), pi/2)