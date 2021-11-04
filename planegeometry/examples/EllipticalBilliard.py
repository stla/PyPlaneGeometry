# -*- coding: utf-8 -*-
from math import sqrt, atan2, cos, sin, pi
import planegeometry.geometry as geom
import numpy as np
import matplotlib.pyplot as plt

def draw_segment(line):
    x1, y1 = line.A
    x2, y2 = line.B
    plt.plot([x1, x2], [y1, y2], color="blue")

def reflect(incidentDir, normalVec):
    return incidentDir - 2*np.vdot(normalVec, incidentDir) * normalVec


# n: number of segments; P0: initial point; v0: initial direction
def trajectory(n, P0, v0):
    out = [None]*n
    L = geom.Line(P0, P0+v0)
    inters = geom.intersection_ellipse_line(ell, L)
    Q0 = inters[1]
    out[0] = geom.Line(inters[0], inters[1])
    for i in range(1, n):
        theta = atan2(Q0[1], Q0[0])
        t = ell.theta2t(theta, degrees = False)
        nrmlVec = ell.normal(t)
        v = reflect(Q0-P0, nrmlVec)
        inters = geom.intersection_ellipse_line(ell, geom.Line(Q0, Q0+v))
        out[i] = geom.Line(inters[0], inters[1])
        P0 = Q0
        Q0 = inters[1] if np.allclose(Q0, inters[0]) else inters[0]
    return out


ell = geom.Ellipse((0,0), 6, 3, 0)

P0 = ell.point_from_angle(60)
v0 = np.array([cos(pi+0.8), sin(pi+0.8)])
traj = trajectory(150, P0, v0)

figure, axes = plt.subplots(figsize=(10, 10))
axes.set_aspect(1)
axes.add_artist(
    plt.Polygon(
        ell.path(), closed=True, fill=False,  
        edgecolor="black", linewidth=2
    )
)
for line in traj:
    draw_segment(line)
plt.title("Elliptical billiard", fontdict = {"fontsize": 40})
plt.xlim(-6.1, 6.1)
plt.ylim(-6.1, 6.1)
plt.axis("off")
plt.show()

