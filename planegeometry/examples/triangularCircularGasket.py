from math import sqrt, inf, copysign, acos, cos, sin, pi
from functools import reduce
import cmath
import planegeometry.geometry as g
import numpy as np
import pyvista as pv
import seaborn as sb

def distance(A, B):
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    return np.linalg.norm(B - A)

def pair_to_triplet(A):
    return tuple(A) + (0,)

def innerSoddyRadius(r1, r2, r3):
    return 1/(1/r1 + 1/r2 + 1/r3 + 2*sqrt(1/r1/r2 + 1/r2/r3 + 1/r3/r1))

def pretty_sphere(radius, center):
    return pv.Sphere(
        radius, center=pair_to_triplet(center), 
        theta_resolution=90, phi_resolution=90
    )

def innerSoddyCircle(plotter, c1, c2, c3, **kwargs):
    r1 = c1.radius
    r2 = c2.radius
    r3 = c3.radius 
    radius = innerSoddyRadius(r1, r2, r3)
    O1 = c1.center
    O2 = c2.center
    O3 = c3.center
    center = g.Triangle(O1, O2, O3).equal_detour_point()[0]
    c123 = g.Circle(center, radius)
    sphere = pretty_sphere(radius, center)
    plotter.add_mesh(sphere, smooth_shading=True, specular=10, **kwargs)
    return (
      ("ccc", c123, c1, c2),
      ("ccc", c123, c2, c3),
      ("ccc", c123, c1, c3)
    )

def side_circle_circle(plotter, A, B, cA, cB, **kwargs):
    if A[1] > B[1]:
        return side_circle_circle(plotter, B, A, cB, cA, **kwargs)
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
    sphere = pretty_sphere(soddyR, center)
    plotter.add_mesh(sphere, smooth_shading=True, specular=10, **kwargs)
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
    tangent = g.Line(P, P + np.array([-OP[1], OP[0]]))
    lineAC = g.Line(A, C)
    lineAB = g.Line(A, B)
    P1 = g.intersection_line_line(tangent, lineAC)
    P2 = g.intersection_line_line(tangent, lineAB)
    incircle = g.Triangle(A, P1, P2).incircle()
    sphere = pretty_sphere(
        incircle.radius, center = incircle.center
    )
    plotter.add_mesh(sphere, smooth_shading=True, specular=10, **kwargs)
    return (
        ("cll", A, B, C, incircle),
        ("ccl", circle, incircle, A, B),
        ("ccl", circle, incircle, A, C)        
    )

def startingCircles(c0, n, phi, shift): # c0 : exterior circle
	invphi = 1/phi
	I = (c0.radius * invphi, 0) + c0.center
	r = c0.radius * sqrt(invphi*invphi-1)
	circ0 = g.Circle(I, r)
	sine = sin(pi/n)
	Cradius = c0.radius / (1+sine)
	Cside = Cradius*sine
	circles0 = [None]*(n+1)
	K = [None]*n
	iota = g.Inversion.on_circle(circ0)
	for i in range(n):
		beta = (i+shift)*2*pi/n
		K[i] = iota.invert(c0.point_from_angle(beta, degrees=False))
		pti = Cradius*np.array([cos(beta), sin(beta)]) + c0.center
		circ1 = g.Circle(pti, Cside)
		circles0[i] = iota.invert_circle(circ1)
	circles0[n] = iota.invert_circle(g.Circle(c0.center, c0.radius-2*Cside))
	J = [None]*n
	invCircles = [None]*(n+1)
	for i in range(n):
		ip1 = (i+1) % n
		ip2 = (i+2) % n
		A = circles0[ip2].center
		rad = circles0[ip2].radius
		B = circles0[ip1].center
		J[i] = A + rad*(B-A)/distance(A, B)
		invCircles[i] = g.Triangle(J[i], K[ip1], K[ip2]).circumcircle()
	invCircles[n] = g.Triangle(J[0],J[1],J[2]).circumcircle(); 
	return (circles0, invCircles)

def new_holes(plotter, holes, color):
    newholes = ()
    for i in range(3):
        hole = holes[i]
        htype = hole[0]
        if htype == "ccc":
            newholes = newholes + (
                innerSoddyCircle(
                    plotter, hole[1], hole[2], hole[3], color=color
                ),
            )
        elif htype == "ccl":
            newholes = newholes + (
                side_circle_circle(
                    plotter, hole[3], hole[4], hole[1], hole[2], color=color
                ),
            )
        elif htype == "cll":
            newholes = newholes + (
                side_side_circle(
                    plotter, hole[1], hole[2], hole[3], hole[4], color=color
                ),
            )
    return newholes


# first generation of children ------------------------------------------
def Circles1(circles0, invCircles):
    circles1 = []
    invs = []
    n = len(circles0) - 1
    for i in range(n):
        ip1 = (i+1) % n
        ip2 = (i+2) % n;
        for j in range(n+1):
            if j != ip1 and j != ip2:
                iota = g.Inversion.on_circle(invCircles[i])
                circle = iota.invert_circle(circles0[j])
                circles1.append(circle)
                invs.append(i)
    return {"circles": circles1, "inversions": invs}

# calculate children ----------------------------------------------------
def children(invCircles, icircles1):
    m = len(invCircles)
    circles1 = icircles1["circles"]
    inversions = icircles1["inversions"]
    n = len(circles1)
    circles2 = []
    invs2 = []
    for i in range(n):
        for j in range(m):
            if inversions[i] != j:
                iota = g.Inversion.on_circle(invCircles[j])
                circle = iota.invert_circle(circles1[i])
                circles2.append(circle)
                invs2.append(j)
    return {"circles": circles2, "inversions": invs2}


def drawCircularGasket(plotter, circle0, n, depth, phi, shift):
    scircles = startingCircles(circle0, n, phi, n*shift)
    circles0 = scircles[0]
    invCircles = scircles[1]
    # construct the children ----------------------------------------------
    allCircles = [None]*depth
    allCircles[0] = {"circles": circles0};
    if depth > 1:
        allCircles[1] = Circles1(circles0, invCircles)
        for i in range(depth)[2:]:
            allCircles[i] = children(invCircles, allCircles[i - 1])
    colors = sb.color_palette(palette="dark", n_colors=depth)
    for i in range(depth):
        zcircles = allCircles[i]["circles"]
        for zcircle in zcircles:
            sphere = pretty_sphere(zcircle.radius, zcircle.center)
            plotter.add_mesh(
                sphere, smooth_shading=True, specular=10, color= colors[i]
            )


def c_list(list1, list2):
    return list1 + list2

def my_triangle(A, B, C):
    points = np.array([
        list(pair_to_triplet(A)), 
        list(pair_to_triplet(B)), 
        list(pair_to_triplet(C))    
    ])
    faces = [3, 0, 1, 2]
    return pv.PolyData(points, faces = faces)

def drawTriangularGasket(plotter, A, B, C, colors, depth):
    ABC = g.Triangle(A, B, C) 
    mcircles = ABC.malfatti_circles()["circles"]
    C1 = mcircles["cA"]
    C2 = mcircles["cB"]
    C3 = mcircles["cC"]
    triangle = my_triangle(A, B, C)
    plotter.add_mesh(
        triangle, show_edges=True, edge_color="yellow", 
        line_width=5, color="yellow", render_lines_as_tubes=True
    )
    holes = (
      side_circle_circle(plotter, A, B, C1, C2, color = colors[1]),
      side_circle_circle(plotter, B, C, C2, C3, color = colors[1]),
      side_circle_circle(plotter, C, A, C3, C1, color = colors[1]),
      innerSoddyCircle(plotter, C1, C2, C3, color = colors[1]),
      side_side_circle(plotter, A, B, C, C1, color = colors[1]),
      side_side_circle(plotter, B, A, C, C2, color = colors[1]),
      side_side_circle(plotter, C, A, B, C3, color = colors[1])
    )
  
    for d in range(depth):
        n_holes = len(holes)
        Holes = ()
        for i in range(n_holes):
            Holes = Holes + (new_holes(plotter, holes[i], colors[d+2]), )
        holes = reduce(c_list, Holes)

###############################################################################

A = (-5, -4)
B = (5, -2)
C = (0, 6)
depth = 3
colors = sb.color_palette(palette="Set1")#, n_colors=depth+2)

mcircles = g.Triangle(A, B, C).malfatti_circles()["circles"]
n1 = 3
n2 = 4
n3 = 5
depth2 = 4
phi = 0.25
shift = 0

pltr = pv.Plotter(window_size = [512, 512], off_screen=True)
pltr.set_background("#363940")
pltr.set_focus([0,0,0])
pltr.camera_position = "xy"
pltr.set_position([0,0,25])
drawTriangularGasket(pltr, A, B, C, colors, depth)
drawCircularGasket(pltr, mcircles["cA"], n1, depth2, phi, shift)
drawCircularGasket(pltr, mcircles["cB"], n2, depth2, phi, shift)
drawCircularGasket(pltr, mcircles["cC"], n3, depth2, phi, shift)
pltr.show(screenshot="MalfattiGaskets.png")
