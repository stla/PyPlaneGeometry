# -*- coding: utf-8 -*-
from math import atan2
from planegeometry.geometry import Circle, Inversion, Arc
from planegeometry.internal import unit_vector_
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb

def tessellation2(depth, Thetas0):
    circ = Circle((0,0), 3)
    arcs = [
        circ.orthogonal_through_two_points_on_circle(
            Thetas0[i], Thetas0[(i+1) % 3], arc=True
        ) for i in [0, 1, 2]
    ]
    inversions = [Inversion(arc.center, arc.radius**2) for arc in arcs]
    Ms = [None]*depth
    Ms[0] = [(unit_vector_(theta), None) for theta in Thetas0]
    Ms[1] = [None]*3
    for i in range(3):
        im1 = 2 if i==0 else i-1
        M = inversions[i].invert(Ms[0][im1][0])
        Ms[1][i] = (M, i)
    for d in range(2, depth):
        n1 = len(Ms[d-1])
        n2 = 2 * n1
        Ms[d] = [None]*n2
        k = 0
        while k < n2:
            for j in range(n1):
                M = Ms[d-1][j]
                for i in range(3):
                    if i != M[1]:
                        newM = inversions[i].invert(M[0])
                        Ms[d][k] = (newM, i)
                        k += 1
    path = [arc.path() for arc in arcs]
    path.reverse()
    path = np.vstack(tuple(path))
    Paths = [(path, 0)]
    for arc in arcs:
        path1 = arc.path()
        xb, yb = arc.starting_point()
        xa, ya = arc.ending_point()
        alpha1 = atan2(ya, xa)
        alpha2 = atan2(yb, xb)
        path2 = Arc((0,0), 3, alpha1, alpha2, False).path()
        Paths.append((np.vstack((path1, path2)), 1))
    Thetas = [None]*len(Ms)
    for i, M in enumerate(Ms):
        Thetas[i] = []
        for pt in M:
            Thetas[i].append(atan2(pt[0][1], pt[0][0]))
    for d in range(2, depth):
        thetas = np.sort(np.concatenate(Thetas[:d]))
        n = len(thetas)
        for i in range(n):
            arc = circ.orthogonal_through_two_points_on_circle(
                thetas[i], thetas[(i+1) % n], arc = True
            )
            path1 = arc.path()
            xb, yb = arc.starting_point()
            xa, ya = arc.ending_point()
            alpha1 = atan2(ya, xa)
            alpha2 = atan2(yb, xb)
            path2 = Arc((0,0), 3, alpha1, alpha2, False).path()
            Paths.append((np.vstack((path1, path2)), d))
    return Paths


depth = 6
Paths = tessellation2(depth, [0, 2, 3.8])

colors = sb.color_palette(palette="bright", n_colors=depth)
  
figure, axes = plt.subplots(facecolor="black", figsize=(10, 10))
axes.set_aspect(1)
for path in Paths:
    axes.add_artist(
        plt.Polygon(
            path[0].tolist(), closed=True, fill=True, 
            facecolor=colors[path[1]], edgecolor="black", linewidth=2
        )
    )
plt.xlim(-4, 4)
plt.ylim(-4, 4)
plt.axis("off")
plt.show()

