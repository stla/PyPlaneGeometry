from math import pi
from planegeometry.geometry import Rotation, Ellipse
import numpy as np

def test_rotation_transform_matrix():
    R = Rotation((1,3), 30)
    Ms = np.random.rand(3, 2)
    Ps = R.transform(Ms)
    assert np.array_equal(Ps[0, :], R.transform(Ms[0, :]))
    assert np.array_equal(Ps[1, :], R.transform(Ms[1, :]))
    assert np.array_equal(Ps[2, :], R.transform(Ms[2, :]))
    
def test_rotate_ellipse():
    ell = Ellipse((5,4), 3, 2, pi/6, False)
    R = Rotation((1,2), 40)
    f = R.as_affine()
    ell1 = R.rotate_ellipse(ell)
    ell2 = f.transform_ellipse(ell)
    assert ell1.is_equal(ell2)
    #
    ell = Ellipse((5,4), 3, 2, 30)
    R = Rotation((1,2), pi/5, False)
    f = R.as_affine()
    ell1 = R.rotate_ellipse(ell)
    ell2 = f.transform_ellipse(ell)
    assert ell1.is_equal(ell2)

