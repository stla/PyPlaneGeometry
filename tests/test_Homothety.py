from planegeometry.geometry import Homothety, Circle

def test_homothety_transform_circle():
    H = Homothety((1, 1), -3)
    circ = Circle((2, 0), 4)
    M = circ.point_from_angle(20)
    imcirc = H.transform_circle(circ)
    P = H.transform(M)
    assert imcirc.includes(P)
    
