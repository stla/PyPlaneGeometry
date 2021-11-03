import planegeometry.geometry as g
import numpy as np

A = (0, 0)
B = (1, 0)
C = (2, 1)
t = g.Triangle(A, B, C)
print("flatness", t.flatness)
t.show()
print(t.equal_detour_point())
print(t.incircle())
circ = t.circumcircle()

print("******************************")
ell1 = g.Ellipse((1,1), 5, 1, 30)
ell2 = g.Ellipse((4,-1), 3, 2, 50) 
f = g.Affine.from_ellipse_to_ellipse(ell1, ell2)
print(f.transform_ellipse(ell1)) # should be ell2

print("******************************")
A = np.array([
        [complex(1,1), complex(0,2)],
        [complex(2,1), complex(1,0)]
    ])
M = g.Mobius(A)
print(M.power(2))
print(M.gpower(2))

print("******************************")
circ = g.Circle((0,0), 3)
circ2 = M.transform_circle(circ)
A = circ.point_from_angle(1)
B = M.transform(A)
print(circ2.includes(B))

print("******************************")
L = g.Line((0,0), (2,2))
P = (-1, 5)
M = L.projection(P)
print(L.includes(M))

print("******************************")
L = g.Line((0,0), (2,2))
P = (-1, 5)
Q = L.reflection(P)
R = L.reflection(Q)
print(P-R)

circ = g.Circle((2,3), 2)
iota = g.Inversion((1,2), 3)
