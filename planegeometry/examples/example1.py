from math import sqrt
import planegeometry.geometry as g
import numpy as np
import matplotlib.pyplot as plt

pts = [
 [0.320245806478104, 0.847395475730481],
 [-2.72896368534611, 1.50534613738611],
 [1.69205843484972, -0.170790753814987],
 [0.420225244618854, 0.82579609214895],
 [-0.589010821409943, -1.76884017226363],
 [-0.417306283096679, 0.374124193982946],
 [0.481812291426394, -0.0593464493617124],
 [0.77950047792719, -0.25745159458116],
 [-2.84646134512524, -0.354800896200286],
 [2.56857349981754, 1.47786749664524],
 [2.32825781223248, 0.686633419065656],
 [-0.518519102008833, -0.296478138003052],
 [0.732161106464013, -0.282239632083836],
 [2.06995626950519, 1.55030888578695],
 [-0.885630303841598, 0.748142680098466],
 [1.12034686089333, 1.03949215112694],
 [-0.869227436389726, 0.490566063720416],
 [-1.12432300759563, -0.681038095223974],
 [2.62795848265206, -0.844981559109443],
 [-0.146512330304941, -0.295855688164278],
 [-1.05735878109256, -1.86414829634254],
 [-0.0563971519881219, -1.0075352613717],
 [-0.869424038317555, 0.327981441656761],
 [-1.34729475620851, -0.997518894240321],
 [2.10160701737814, -1.31216251071656],
 [-0.2569011202426, 0.694065859250904],
 [-1.99573307078148, -0.123594195830325],
 [2.81005363680948, -0.568438457376384],
 [2.07695896209, -0.783720711661815],
 [-0.666999989780655, 0.353084824771068]
]

pts = np.asarray(pts)

ell = g.Ellipse.LownerJohnEllipse(pts)

ell = g.Ellipse((1,1), 5, 1, 30)
points = ell.random_points(4, "on")

ell2 = g.Ellipse.from_boundary4(points)
S = ell2[0]
c = ell2[1]

print("::::::::::::::::::::::::::::::::::")
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
