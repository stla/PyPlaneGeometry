from planegeometry.geometry import Ellipse, Affine
import numpy as np

def test_ellipse_equation():
    ell = Ellipse((4,3), 5, 1, 100)
    A, B, C, D, E, F = ell.equation().values()
    ell2 = Ellipse.from_equation(A, B, C, D, E, F)
    assert ell.is_equal(ell2)
    
def test_affine_transforms_ellipse():
    ell0 = Ellipse((1,1), 5, 2, 30)
    f = Affine([[3.5, 2], [0, 4]], [-1, 1.25])
    ell1 = f.transform_ellipse(ell0)
    path0 = ell0.path(3)
    Q = f.transform(path0[0, :])
    assert ell1.includes(Q)
    Q = f.transform(path0[1, :])
    assert ell1.includes(Q)
    Q = f.transform(path0[2, :])
    assert ell1.includes(Q)

def test_ellipse_random_on():
    ell = Ellipse((4, 3), 5, 1, 100)
    rpoints = ell.random_points(3, "on")
    assert ell.includes(rpoints[0, :])
    assert ell.includes(rpoints[1, :])
    assert ell.includes(rpoints[2, :])

def test_ellipse_random_in():
    ell = Ellipse((4, 3), 5, 1, 100)
    rpoints = ell.random_points(3, "in")
    assert ell.contains(rpoints[0, :])
    assert ell.contains(rpoints[1, :])
    assert ell.contains(rpoints[2, :])
    
def test_ellipse_from_five_points():
    ell = Ellipse((2, 3), 5, 4, 50)
    np.random.seed(666)
    pts = ell.random_points(5, "on")
    ell2 = Ellipse.from_five_points(*pts)
    assert ell.is_equal(ell2)

def test_ellipse_matrix():
    ell = Ellipse((2, 3), 5, 4, 50)
    S = ell.shape_matrix()
    ell2 = Ellipse.from_center_and_matrix((2, 3), S)
    assert ell.is_equal(ell2)


