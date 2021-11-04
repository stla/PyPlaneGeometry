from planegeometry.geometry import Ellipse

def test_ellipse_equation():
    ell = Ellipse((4,3), 5, 1, 100)
    A, B, C, D, E, F = ell.equation().values()
    ell2 = Ellipse.from_equation(A, B, C, D, E, F)
    assert ell.is_equal(ell2)
