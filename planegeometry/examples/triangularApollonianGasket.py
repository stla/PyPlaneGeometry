from math import sqrt, inf, copysign
from functools import reduce
import cmath
import planegeometry.geometry as g
import numpy as np
import pyvista as pv

def distance(A, B):
    return np.linalg.norm(B - A)

def innerSoddyRadius(r1, r2, r3):
    return 1/(1/r1 + 1/r2 + 1/r3 + 2*sqrt(1/r1/r2 + 1/r2/r3 + 1/r3/r1))

def innerSoddyCircle(plotter, c1, c2, c3, **kwargs):
    r1 = c1.radius
    r2 = c2.radius
    r3 = c3.radius 
    radius = innerSoddyRadius(r1, r2, r3)
    O1 = c1.center
    O2 = c2.center
    O3 = c3.center
    center = g.Triangle(O1, O2, O3).equalDetourPoint()[0]
    c123 = g.Circle(center, radius)
    sphere = pv.Sphere(radius, center = tuple(center) + (0,))
    plotter.add_mesh(sphere, smooth_shading=True, specular=10, **kwargs)
    return (
      ("ccc", c123, c1, c2),
      ("ccc", c123, c2, c3),
      ("ccc", c123, c1, c3)
    )

def side_circle_circle(plotter, A, B, cA, cB, **kwargs):
    if A[1] > B[1]:
        return side_circle_circle(B, A, cB, cA, **kwargs)
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
    if zX1.real < zX2.real: #use math.copysign
        Y = complex(
            2*rA*sqrt(rB)/(sqrt(rA)+sqrt(rB)) + zX1.real,
            copysign(soddyR, Im(zX1))
        )
    else:
        Y = complex(
            2*rB*sqrt(rA)/(sqrt(rA)+sqrt(rB)) + zX2.real,
            copysign(soddyR, Im(zX1))
        )
    Z = exp(1j*alpha)*Y + zB
    c_x = Z.real
    c_y = Z.imag
    center = (c_x, c_y)
    cAB = Circle(center, soddyR)
    sphere = pv.Sphere(soddyR, center = (c_x, c_y, 0))
    plotter.add_mesh(sphere, smooth_shading=True, specular=10, ...)
    return (
      ("ccc", cAB, cA, cB),
      ("ccl", cA, cAB, A, B),
      ("ccl", cAB, cB, A, B)
    )

def side_side_circle(plotter, A, B, C, circle, **kwargs):
    zA = complex(*A)
    zO = complex(*(circle.center))
    vec = zA - zO
    zP = zO + circle.radius * vec/abs(vec)
    P = np.asarray((zP.real, zP.imag))
    OP = P - circle.center
    tangent = Line(P, P + np.array([-OP[1], OP[0]]))
    lineAC = Line(A, C)
    lineAB = Line(A, B)
    P1 = g.intersection_line_line(tangent, lineAC)
    P2 = g.intersection_line_line(tangent, lineAB)
    incircle = g.Triangle(A, P1, P2).incircle()
    sphere = pv.Sphere(incircle.radius, center = tuple(incircle.center) + (0,))
    plotter.add_mesh(sphere, smoothing=True, specular=10, **kwargs)
    return (
        ("cll", A, B, C, incircle),
        ("ccl", circle, incircle, A, B),
        ("ccl", circle, incircle, A, C)        
    )

def new_holes(holes, color):
    newholes = []
    for i in range(2):
        hole = holes[i]
        type = hole[0]
        if type == "ccc":
            newholes.append(
                innerSoddyCircle(hole[1], hole[2], hole[3], color = color)
            )
        elif type == "ccl":
            newholes.append(
                side_circle_circle(
                    hole[3], hole[4], hole[1], hole[2], color = color
                )    
            )
        elif type == "cll":
            newholes.append(
                side_side_circle(
                    hole[1], hole[2], hole[3], color = color
                )    
            )
    return newholes

def c_list(list1, list2):
    return list1 + list2

def pair_to_triplet(A):
    return tuple(A) + (0,)

def drawTriangularGasket(plotter, A, B, C, colors, depth):
    ABC = g.Triangle(A, B, C) 
    mcircles = ABC.malfatti_circles()["circles"]
    C1 = mcircles[0]
    C2 = mcircles[1]
    C3 = mcircles[2]
    triangle = pv.Triangle(
        [pair_to_triplet(A), pair_to_triplet(B), pair_to_triplet(C)]
    )
    plotter.add_mesh(triangle, show_edges=True, line_width=5)
    plotter.add_mesh(
        pv.Sphere(
            C1.radius, center=C1.center+(0,), specular=10, color=colors[0]
        )    
    ) 
    plotter.add_mesh(
        pv.Sphere(
            C2.radius, center=C2.center+(0,), specular=10, color=colors[0]
        )    
    ) 
    plotter.add_mesh(
        pv.Sphere(
            C3.radius, center=C3.center+(0,), specular=10, color=colors[0]
        )    
    ) 
    holes = (
      side_circle_circle(A, B, C1, C2, color = colors[1]),
      side_circle_circle(B, C, C2, C3, color = colors[1]),
      side_circle_circle(C, A, C3, C1, color = colors[1]),
      innerSoddyCircle(C1, C2, C3, color = colors[1]),
      side_side_circle(A, B, C, C1, color = colors[1]),
      side_side_circle(B, A, C, C2, color = colors[1]),
      side_side_circle(C, A, B, C3, color = colors[1])
    )
  
  for d in range(depth):
      n_holes = len(holes)
      Holes = tuple([None]*n_holes)
      for i in range(n_holes):
          Holes[i] = new_holes(holes[i], colors[d+2])
      holes = reduce(c_list, Holes)


