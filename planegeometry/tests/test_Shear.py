from planegeometry.geometry import Shear
import numpy as np

def test_shear_transform_matrix():
    S = Shear((1,1), (1,3), 0.5, 30)
    Ms = np.random.rand(3, 2)
    Ps = S.transform(Ms)
    assert np.array_equal(Ps[0, :], S.transform(Ms[0, :]))
    assert np.array_equal(Ps[1, :], S.transform(Ms[1, :]))
    assert np.array_equal(Ps[2, :], S.transform(Ms[2, :]))
    
