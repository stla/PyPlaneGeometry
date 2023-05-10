from math import sqrt, pi, cos, sin
import matplotlib.pyplot as plt
import planegeometry.geometry as geom
import seaborn as sb
import numpy as np

# starting circles ####
c0 = geom.Circle((0,0), 1) # the exterior circle
n = 3
circles0 = geom.SteinerChain(c0, n, phi = 0.25, shift = 0)["circles"]

# construct the inversions ####
inversions = [None]*(n+1)
for i in range(n):
    inversions[i] = geom.Inversion.from_fixing_three_circles(
        c0, circles0[i], circles0[(i+1) % n]
    )
inversions[n] = geom.Inversion.from_swapping_two_circles(
    c0, circles0[n]
)

# first generation of children
circles1 = []
for i in range(n):
    ip1 = (i+1) % n
    for j in range(n+1): #(j in (1L:(n+1L))[-c(i,ip1)]){
        if j != i and j != ip1:
            circle = inversions[i].invert_circle(circles0[j])
            circles1.append((circle, i))


# function to construct the "children" ####
def children(inversions, circles1):
    m = len(inversions)
    n = len(circles1) 
    circles2 = [] 
    for i in range(n):
        k = circles1[i][1]
        for j in range(m):
            if j != k:
                circle = inversions[j].invert_circle(circles1[i][0])
                circles2.append((circle, j))
    return circles2

# construct children ####
depth = 5
allCircles = [None]*depth
allCircles[0] = circles0
allCircles[1] = circles1
for i in range(depth)[2:]:
    allCircles[i] = children(inversions, allCircles[i-1])
for i in range(depth)[1:]:
    allCircles[i] = [c[0] for c in allCircles[i]]

# Möbius transformation power t
def mobius(gamma, t):
    g = abs(gamma)
    h = sqrt(1-g*g)
    d2 = h**t * (cos(t*pi/2) + 1j*sin(t*pi/2))
    d1 = d2.conjugate()
    h11 = d1.real - 1j *d1.imag/h
    h12 = d2.imag * gamma/h
    h21 = h12.conjugate()
    h22 = h11.conjugate()
    return geom.Mobius([[h11,h12], [h21,h22]])

  
# plot ####
colors = sb.color_palette(palette="bright", n_colors=depth)
def draw_circle(C, color, fill=True):
    axes.add_artist(
        plt.Circle(
            C.center, C.radius, fill=fill, facecolor=color, 
            edgecolor="black", linewidth=2
        )
    )
gamma = 0.6 + 0.4j
t_ = np.linspace(0, 2, 11)[:10]
for j, t in enumerate(t_): 
    figure, axes = plt.subplots(figsize=(10, 10))
    axes.set_aspect(1)
    draw_circle(c0, "black", False)
    M = mobius(gamma,t)
    for i in range(depth):
        for circ in allCircles[i]:
            cc = M.transform_circle(circ)
            draw_circle(cc, colors[i])
    plt.xlim(-1.2, 1.2)
    plt.ylim(-1.2, 1.2)
    plt.axis("off")
    plt.show()
