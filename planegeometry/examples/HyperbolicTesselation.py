import planegeometry.geometry as g
from math import cos, sin, atan2
import matplotlib.pyplot as plt
import seaborn as sb

Thetas0 = [0, 2, 3.8]

circ = g.Circle((0,0), 3)

depth = 5
colors = sb.color_palette(palette="Set1", n_colors=depth)

arcs = [circ.orthogonalThroughTwoPointsOnCircle(
    Thetas0[i], Thetas0[(i+1) % 3], arc = True
) for i in range(3)]
  
inversions = [g.Inversion(arc.center, arc.radius**2) for arc in arcs]
  
Ms = [None]*depth
  
Ms[0] = [(cos(theta), sin(theta)) for theta in Thetas0]
  
Ms[1] = [None]*3
for i in range(3):
     im1 = (i-1) % 3
     M = inversions[i].invert(Ms[0][im1])
     Ms[1][i] = (M, i)

for d in range(depth)[2:]:
    n1 = len(Ms[d-1])
    n2 = 2*n1
    Ms[d] = [None]*n2
    k = 0
    while k < n2:
        for j in range(n1):
            M, iota = Ms[d-1][j]
            for i in range(3):
                if i != iota:
                    newM = inversions[i].invert(M)
                    Ms[d][k] = (newM, i)
                    k += 1

figure, axes = plt.subplots(facecolor="black", figsize=(10, 10))
axes.set_aspect(1)

def draw_arc(arc, color):
    axes.add_artist(
        plt.Polygon(arc.path(), fill=False, color=color, closed=False)
    )

[draw_arc(arc, colors[0]) for arc in arcs]

for i in range(3):
    x, y = Ms[0][i]
    Ms[0][i] = atan2(y,x)
for j in range(depth)[1:]:
    for k in range(len(Ms[j])):
        a, iota = Ms[j][k]
        x, y = a
        Ms[j][k] = atan2(y,x)
        
thetas = Ms[0]
for d in range(depth)[1:]:
    thetas = thetas + Ms[d]
    thetas.sort()
    for i in range(len(thetas)):
        ip1 = (i+1) % len(thetas)
        arc = circ.orthogonalThroughTwoPointsOnCircle(
            thetas[i], thetas[ip1], arc = True
        )
        draw_arc(arc, colors[d])



plt.xlim(-4, 4)
plt.ylim(-4, 4)
plt.axis("off")
plt.show()
