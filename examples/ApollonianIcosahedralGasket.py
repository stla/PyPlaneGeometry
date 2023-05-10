# -*- coding: utf-8 -*-
from math import sqrt, inf, copysign, acos, cos, sin
from functools import reduce
import cmath
import planegeometry.geometry as g
import numpy as np
import pyvista as pv
import seaborn as sb

def vlength(v):
    return np.linalg.norm(v) 

def transfoMatrix(A, B, C):
    Bprime = B - A
    Cprime = C - A
    n = np.cross(Bprime, Cprime)
    n = n / vlength(n)
    q = np.cross(n, Bprime)
    q = q / vlength(q)
    p = Bprime / vlength(Bprime)
    R = np.column_stack((p, q, n))
    M = np.hstack(
        (
            np.vstack(
                (R, np.zeros((1,3), dtype=float))
            ), 
            np.concatenate((A, [1])).reshape((4, 1))
        )
    )
    tR = np.transpose(R)
    return {
     "Mat": M,
     "P": tR.dot(Bprime)[:2],
     "Q": tR.dot(Cprime)[:2]
    }

A = np.array([0, 0, 1])
B = np.array([10, 2, 0])
C = np.array([5, 10, 3])

def distance(A, B):
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    return np.linalg.norm(B - A)

def pair_to_triplet(A):
    return tuple(A) + (0,)

def innerSoddyRadius(r1, r2, r3):
    return 1/(1/r1 + 1/r2 + 1/r3 + 2*sqrt(1/r1/r2 + 1/r2/r3 + 1/r3/r1))

def pretty_sphere(radius, center):
    sphere = pv.Sphere(
        radius, center=pair_to_triplet(center), 
        theta_resolution=60, phi_resolution=60
    )
    sphere.scale([1.0, 1.0, 0.3])
    return sphere

def innerSoddyCircle(plotter, c1, c2, c3, transfo, **kwargs):
    r1 = c1.radius
    r2 = c2.radius
    r3 = c3.radius 
    radius = innerSoddyRadius(r1, r2, r3)
    O1 = c1.center
    O2 = c2.center
    O3 = c3.center
    center = g.Triangle(O1, O2, O3).equal_detour_point()[0]
    c123 = g.Circle(center, radius)
    sphere = pretty_sphere(
        radius, center
    ).transform(transfo)
    plotter.add_mesh(sphere, smooth_shading=True, specular=10, **kwargs)
    return (
      ("ccc", c123, c1, c2),
      ("ccc", c123, c2, c3),
      ("ccc", c123, c1, c3)
    )

def side_circle_circle(plotter, A, B, cA, cB, transfo, **kwargs):
    if A[1] > B[1]:
        return side_circle_circle(plotter, B, A, cB, cA, transfo, **kwargs)
    rA = cA.radius
    rB = cB.radius
    oA = cA.center
    oB = cB.center
    zoA = complex(*oA)
    zoB = complex(*oB)
    zB = complex(*A)
    alpha = acos((B[0]-A[0]) / distance(A, B))
    zX1 = cmath.exp(-alpha*1j) * (zoA-zB)
    zX2 = cmath.exp(-alpha*1j) * (zoB-zB)
    #stopifnot(sign(Im(zX1))==sign(Im(zX2)))
    soddyR = innerSoddyRadius(rA, rB, inf)
    if zX1.real < zX2.real: 
        Y = complex(
            2*rA*sqrt(rB)/(sqrt(rA)+sqrt(rB)) + zX1.real,
            copysign(soddyR, zX1.imag)
        )
    else:
        Y = complex(
            2*rB*sqrt(rA)/(sqrt(rA)+sqrt(rB)) + zX2.real,
            copysign(soddyR, zX1.imag)
        )
    Z = cmath.exp(1j*alpha)*Y + zB
    c_x = Z.real
    c_y = Z.imag
    center = (c_x, c_y)
    cAB = g.Circle(center, soddyR)
    sphere = pretty_sphere(
        soddyR, center
    ).transform(transfo)
    plotter.add_mesh(sphere, smooth_shading=True, specular=10, **kwargs)
    return (
      ("ccc", cAB, cA, cB),
      ("ccl", cA, cAB, A, B),
      ("ccl", cAB, cB, A, B)
    )

def side_side_circle(plotter, A, B, C, circle, transfo, **kwargs):
    zA = complex(*A)
    zO = complex(*(circle.center))
    vec = zA - zO
    zP = zO + circle.radius * vec/abs(vec)
    P = np.asarray((zP.real, zP.imag))
    OP = P - circle.center
    tangent = g.Line(P, P + np.array([-OP[1], OP[0]]))
    lineAC = g.Line(A, C)
    lineAB = g.Line(A, B)
    P1 = g.intersection_line_line(tangent, lineAC)
    P2 = g.intersection_line_line(tangent, lineAB)
    incircle = g.Triangle(A, P1, P2).incircle()
    sphere = pretty_sphere(
        incircle.radius, incircle.center
    ).transform(transfo)
    plotter.add_mesh(sphere, smooth_shading=True, specular=10, **kwargs)
    return (
        ("cll", A, B, C, incircle),
        ("ccl", circle, incircle, A, B),
        ("ccl", circle, incircle, A, C)        
    )

def new_holes(plotter, holes, transfo, color):
    newholes = ()
    for i in range(3):
        hole = holes[i]
        htype = hole[0]
        if htype == "ccc":
            newholes = newholes + (
                innerSoddyCircle(
                    plotter, hole[1], hole[2], hole[3], transfo, color=color
                ),
            )
        elif htype == "ccl":
            newholes = newholes + (
                side_circle_circle(
                    plotter, hole[3], hole[4], hole[1], hole[2], transfo, 
                    color=color
                ),
            )
        elif htype == "cll":
            newholes = newholes + (
                side_side_circle(
                    plotter, hole[1], hole[2], hole[3], hole[4], transfo, 
                    color=color
                ),
            )
    return newholes

def c_tuples(tup1, tup2):
    return tup1 + tup2

def my_triangle(A, B, C):
    points = np.array([
        list(pair_to_triplet(A)), 
        list(pair_to_triplet(B)), 
        list(pair_to_triplet(C))    
    ])
    faces = [3, 0, 1, 2]
    return pv.PolyData(points, faces = faces)

def drawFace(plotter, A, B, C, colors, depth):
    
    tm = transfoMatrix(A, B, C)
    A = (0.0, 0.0)
    B = tm["P"]
    C = tm["Q"]
    transfo = tm["Mat"]
  
    ABC = g.Triangle(A, B, C) 
    mcircles = ABC.malfatti_circles()["circles"]
    C1 = mcircles["cA"]
    C2 = mcircles["cB"]
    C3 = mcircles["cC"]
    plotter.add_mesh(
        pretty_sphere(
            C1.radius, C1.center
        ).transform(transfo), 
        specular=10, color=colors[0], smooth_shading=True
    )    
    plotter.add_mesh(
        pretty_sphere(
            C2.radius, C2.center
        ).transform(transfo), 
        specular=10, color=colors[0], smooth_shading=True
    )    
    plotter.add_mesh(
        pretty_sphere(
            C3.radius, C3.center
        ).transform(transfo), 
        specular=10, color=colors[0], smooth_shading=True
    )    
    holes = (
      side_circle_circle(plotter, A, B, C1, C2, transfo, color = colors[1]),
      side_circle_circle(plotter, B, C, C2, C3, transfo, color = colors[1]),
      side_circle_circle(plotter, C, A, C3, C1, transfo, color = colors[1]),
      innerSoddyCircle(plotter, C1, C2, C3, transfo, color = colors[1]),
      side_side_circle(plotter, A, B, C, C1, transfo, color = colors[1]),
      side_side_circle(plotter, B, A, C, C2, transfo, color = colors[1]),
      side_side_circle(plotter, C, A, B, C3, transfo, color = colors[1])
    )
  
    for d in range(depth):
        n_holes = len(holes)
        Holes = ()
        for i in range(n_holes):
            Holes += (new_holes(plotter, holes[i], transfo, colors[d+2]), )
        holes = reduce(c_tuples, Holes)

###############################################################################

vertices = np.array([
    [0, 0.618033988749895, 1],
    [0, 0.618033988749895, -1],
    [0, -0.618033988749895, 1],
    [0, -0.618033988749895, -1],
    [0.618033988749895, 1, 0],
    [0.618033988749895, -1, 0],
    [-0.618033988749895, 1, 0],
    [-0.618033988749895, -1, 0],
    [1, 0, 0.618033988749895],
    [-1, 0, 0.618033988749895],
    [1, 0, -0.618033988749895],
    [-1, 0, -0.618033988749895]
])

faces = [
    [0, 2, 8],
    [0, 8, 4],
    [0, 4, 6],
    [0, 6, 9],
    [0, 9, 2],
    [3, 11, 1],
    [3, 1, 10],
    [3, 10, 5],
    [3, 5, 7],
    [3, 7, 11],
    [8, 2, 5],
    [4, 8, 10],
    [6, 4, 1],
    [9, 6, 11],
    [2, 9, 7],
    [1, 11, 6],
    [10, 1, 4],
    [5, 10, 8],
    [7, 5, 2],
    [11, 7, 9]
]

depth = 1
colors = sb.color_palette(palette="Set1")

pltr = pv.Plotter(window_size = [512, 512])
pltr.set_background("#363940")
for face in faces:
    ABC = vertices[face]
    triangle = pv.PolyData(ABC, faces = [3, 0, 1, 2])
    pltr.add_mesh(
        triangle, show_edges=True, edge_color="yellow", 
        line_width=1, color="yellow", render_lines_as_tubes=True
    )
    A, B, C = ABC
    drawFace(pltr, A, B, C, colors, depth)
pltr.camera.zoom(1.2)
pltr.show()
