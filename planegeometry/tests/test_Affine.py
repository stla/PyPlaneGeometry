from planegeometry.geometry import Affine
import numpy as np

def test_affine_transform_matrix():
    f = Affine([[3.5, 2], [0, 4]], [-1, 1.25])
    Ms = np.random.rand(3, 2)
    Ps = f.transform(Ms)
    assert np.array_equal(Ps[0, :], f.transform(Ms[0, :]))
    assert np.array_equal(Ps[1, :], f.transform(Ms[1, :]))
    assert np.array_equal(Ps[2, :], f.transform(Ms[2, :]))
    
def test_affine_from_mapping_three_points():
    P1 = (1, 2)
    P2 = (3, 4)
    P3 = (7, 7)
    Q1 = (2, 1)
    Q2 = (4, 4)
    Q3 = (-7, 7)
    f = Affine.from_mapping_three_points(P1, P2, P3, Q1, Q2, Q3)
    assert np.array_equal(Q1, f.transform(P1))
    assert np.array_equal(Q2, f.transform(P2))
    assert np.array_equal(Q3, f.transform(P3))
