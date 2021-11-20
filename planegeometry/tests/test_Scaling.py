from planegeometry.geometry import Scaling, Circle

def test_scale_circle():
    S = Scaling((1,1), (2,3), 2)
    f = S.as_affine()
    circ = Circle((-1, 2), 4)
    ell1 = S.scale_circle(circ)
    ell2 = f.transform_ellipse(circ)
    assert ell1.is_equal(ell2)
