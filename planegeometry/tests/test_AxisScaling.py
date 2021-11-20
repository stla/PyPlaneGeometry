from planegeometry.geometry import ScalingXY
import numpy as np

def test_axisscaling_transform_matrix():
    S = ScalingXY((1,1), 4, 2)
    Ms = np.random.rand(3, 2)
    Ps = S.transform(Ms)
    assert np.array_equal(Ps[0, :], S.transform(Ms[0, :]))
    assert np.array_equal(Ps[1, :], S.transform(Ms[1, :]))
    assert np.array_equal(Ps[2, :], S.transform(Ms[2, :]))
    
