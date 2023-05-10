# -*- coding: utf-8 -*-
from math import pi, cos, sin
import matplotlib.pyplot as plt
import planegeometry.geometry as geom
import seaborn as sb

def unique_with(L, f):
    size = len(L)
    for i in range(size-1):
        j = i + 1
        while j < size:
            if f(L[i], L[j]):
                del L[j]
                size -= 1
            else:
                j += 1
    return L[:size]


def draw_circle(axes, C, color):
    axes.add_artist(
        plt.Circle(
            C.center, C.radius, fill=False, edgecolor=color, linewidth=2
        )
    )

def draw_line(axes, L, color):
    axes.axline(L.A, L.B, linewidth=2, color=color)
    
def draw_gcircle(axes, G, color):
    if isinstance(G, geom.Circle):
        draw_circle(axes, G, color)
    else:
        draw_line(axes, G, color)


# generation 0
gen0 = [(geom.Circle((cos(beta), sin(beta)), 1), 0) for beta in [0, pi/2, pi, 3*pi/2]]
gen0.append((geom.Circle((0,0), 2), 0))

# generation 1
n0 = len(gen0)
n1 = n0*(n0-1)
gen1 = [None]*n1
k = 0
while k < n1:
    for j in range(n0):
        for i in range(n0):
            if i != j:
                circ = gen0[i][0]
                iota = geom.Inversion(circ.center, circ.radius**2)
                gen1[k] = (iota.invert_circle(gen0[j][0]), 1, i)
                k += 1

# generation 2
n2 = n0*n1-n1
gen2 = [None]*n2
k = 0
while k < n2:
    for j in range(n1):
        for i in range(n0):
            if gen1[j][2] != i:
                circ = gen0[i][0]
                iota = geom.Inversion(circ.center, circ.radius**2)
                G = gen1[j][0]
                if isinstance(G, geom.Circle):
                    gen2[k] = (iota.invert_circle(G), 2, i)
                else:
                    gen2[k] = (iota.invert_line(G), 2, i)
                k += 1

# generation 3
n3 = n0*n2-n2
gen3 = [None]*n3
k = 0
while k < n3:
    for j in range(n2):
        for i in range(n0):
            if gen2[j][2] != i:
                circ = gen0[i][0]
                iota = geom.Inversion(circ.center, circ.radius**2)
                G = gen2[j][0]
                if isinstance(G, geom.Circle):
                    gen3[k] = (iota.invert_circle(G), 3, i)
                else:
                    gen3[k] = (iota.invert_line(G), 3, i)
                k += 1
        
gcircles = gen0 + gen1 + gen2 + gen3
gcircles = unique_with(gcircles, lambda g1, g2: type(g1[0]) == type(g2[0]) and g1[0].is_equal(g2[0]))

colors = sb.color_palette(palette="viridis", n_colors=4)
        
figure, axes = plt.subplots(facecolor="black", figsize=(10, 10))
axes.set_aspect(1)
for gcircle in gcircles:
    draw_gcircle(axes, gcircle[0], colors[gcircle[1]])
plt.xlim(-2.3, 2.3)
plt.ylim(-2.3, 2.3)
plt.axis("off")
plt.show()
