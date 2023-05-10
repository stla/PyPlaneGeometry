# -*- coding: utf-8 -*-
from math import sqrt, cos, pi
from planegeometry.geometry import Triangle, Circle, Inversion
import pyvista as pv
import numpy as np


# outer Soddy circle ----------------------------------------------------------
def soddyOuterCircle(ABC):
    a = ABC.a 
    b = ABC.b
    c = ABC.c
    p = a + b + c
    s = p / 2
    area = sqrt(s*(s-a)*(s-b)*(s-c))
    R = a*b*c / (4*sqrt(s*(a+b-s)*(a+c-s)*(b+c-s)))
    r = area / s
    radius = area / (4*R + r - p)
    p1 = a - area/(s-a)
    p2 = b - area/(s-b)
    p3 = c - area/(s-c)
    center = (p1*ABC.A + p2*ABC.B + p3*ABC.C) / (p1 + p2 + p3)
    return Circle(center, abs(radius))

def iteration(circles0, inversions):
    out = []
    for mycircle in circles0:
       circle = mycircle[0]
       ji = mycircle[1]
       for i in range(4):
         if i != ji:
          out.append((inversions[i].invert_circle(circle), i))
    return out

def gasket(pltr, circles0, inversions, depth, colors):
    if depth > 0:
       circles = iteration(circles0, inversions)
       for circle in circles:
            circ = circle[0]
            cx, cy = circ.center
            sphere = pv.Sphere(circ.radius, center = (cx, cy, 0))
            pltr.add_mesh(sphere, color=colors[0], specular=20, smooth_shading=True)
       del colors[0]
       gasket(pltr, circles, inversions, depth-1, colors)

def drawGasket(pltr, t, depth, colors):
    Mcircles = list(t.malfatti_circles()["circles"].values())
    triangle = Triangle(Mcircles[0].center, Mcircles[1].center, Mcircles[2].center)
    soddyO = soddyOuterCircle(triangle);
    Mcircles.append(soddyO);
    for i in range(4):
        circ = Mcircles[i]
        cx, cy = circ.center
        circle = pv.Circle(circ.radius)
        circle.translate((cx,cy,0))
        pltr.add_mesh(circle, color="#363940", show_edges=True, opacity=1, line_width=5, edge_color="black")
    inversions = [None]*4
    circles0 = [None]*4
    inversions[0] = Inversion.from_fixing_three_circles(soddyO, Mcircles[1], Mcircles[2])
    inversions[1] = Inversion.from_fixing_three_circles(soddyO, Mcircles[0], Mcircles[2])
    inversions[2] = Inversion.from_fixing_three_circles(soddyO, Mcircles[0], Mcircles[1])
    inversions[3] = Inversion.from_fixing_three_circles(Mcircles[0], Mcircles[1], Mcircles[2])
    for i in range(4):
        circ = (inversions[i].invert_circle(Mcircles[i]), i)
        cx, cy = circ[0].center
        sphere = pv.Sphere(circ[0].radius, center = (cx, cy, 0))
        pltr.add_mesh(sphere, color=colors[0], specular=20, smooth_shading=True)
        circles0[i] = circ
    del colors[0]
    gasket(pltr, circles0, inversions, depth, colors)
    
depth = 3
A = np.array([0,0], dtype=float)
B = np.array([1,0], dtype=float) # do not change
C = np.array([cos(2*pi*0/180), sqrt(3)/2], dtype=float)
t = Triangle(A, B, C)
Mcircles = list(t.malfatti_circles()["circles"].values())
triangle = Triangle(Mcircles[0].center, Mcircles[1].center, Mcircles[2].center)
soddyO = soddyOuterCircle(triangle)
center = soddyO.center
radius = soddyO.radius
A1 = (A-center)/radius
B1 = (B-center)/radius
C1 = (C-center)/radius
t1 = Triangle(A1, B1, C1)
pltr = pv.Plotter()
pltr.set_background("#363940")
drawGasket(pltr, t1, depth, ["yellow", "orange", "magenta", "forestgreen"])
pltr.show()
