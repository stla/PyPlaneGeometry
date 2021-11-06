from math import sqrt
import planegeometry.geometry as g
import numpy as np
import matplotlib.pyplot as plt

R = g.Rotation((1,1), 30)
M = (2,1)
R.rotate(M)

M = [[2,1], [2,1], [2,1]]
R.rotate(M)

print("##############################")
A = (0, 0)
B = (1, 0)
C = (2, 1)
t = g.Triangle(A, B, C)
print("flatness", t.flatness)
t.show()
print(t.equal_detour_point())
print(t.incircle())
circ = t.circumcircle()
print(circ.includes(A))
        
steinerEllipse = t.steiner_ellipse()
steinerInellipse = t.steiner_inellipse()

figure, axes = plt.subplots(figsize=(10, 10))
axes.set_aspect(1)
plt.plot([A[0], B[0], C[0], A[0]], [A[1], B[1], C[1], A[1]])
axes.add_artist(
    plt.Polygon(
        steinerEllipse.path(), closed=True, fill=False, 
        edgecolor="red", linewidth=2
    )
)
axes.add_artist(
    plt.Polygon(
        steinerInellipse.path(), closed=True, fill=False, 
        edgecolor="green", linewidth=2
    )
)
plt.xlim(-1, 3)
plt.ylim(-1, 3)
plt.show()

print("--------------------------------")
circ = g.Circle((1,1), 5)
line = g.Line((2,-2), (0,4))
Is = g.intersection_circle_line(circ, line)
print(circ.includes(Is[0]) and circ.includes(Is[1]))

print("--------------------------------")
ell = g.Ellipse((1,1), 5, 1, 30)
line = g.Line((2,-2), (0,4))
Is = g.intersection_ellipse_line(ell, line)
print(ell.includes(Is[0]) and ell.includes(Is[1]))

figure, axes = plt.subplots(figsize=(10, 10))
axes.set_aspect(1)
#axes.axline(line.A, line.B, linewidth=2, color="red")
axes.add_artist(
    plt.Polygon(
        ell.path(), closed=True, fill=False, 
        edgecolor="green", linewidth=2
    )
)
plt.plot([Is[0][0], Is[1][0]], [Is[0][1], Is[1][1]])
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.show()

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
