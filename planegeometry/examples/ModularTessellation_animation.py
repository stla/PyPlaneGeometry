# -*- coding: utf-8 -*-
import os
import planegeometry.geometry as geom
import matplotlib.pyplot as plt
import numpy as np

# Möbius transformations
T = geom.Mobius(np.array([[0,-1], [1,0]]))
U = geom.Mobius(np.array([[1,1], [0,1]]))
R = U.compose(T)
# R**t, generalized power
def Rt(t):
    return R.gpower(t)

# starting circles
I = geom.Circle((0, 1.5), 0.5)
TI = T.transform_circle(I)

# modified Cayley transformation
Phi = geom.Mobius(np.array([[1j, 1], [1, 1j]]))

# function to draw a circle
def draw_circle(axes, C):
    axes.add_artist(
        plt.Circle(
            C.center, C.radius, fill=True, facecolor="magenta", 
            edgecolor="black"
        )
    )
    
def draw_pair(axes, M, u, compose=False):
    if compose:
        M = M.compose(T)
    A = M.compose(Rt(u)).compose(Phi)
    C = A.transform_circle(I)
    draw_circle(axes, C)
    C = A.transform_circle(TI)
    draw_circle(axes, C)
    if not compose:
        draw_pair(axes, M, u, compose=True)

n = 8
transfos = geom.unimodular_matrices(n)

u_ = np.linspace(0, 3, 181)[:180]

for i, u in enumerate(u_):
    figure, axes = plt.subplots(facecolor="black", figsize=(10, 10))
    axes.set_aspect(1)
    for transfo in transfos:
        M = geom.Mobius(transfo)
        draw_pair(axes, M, u)
        M = M.inverse()
        draw_pair(axes, M, u)
        np.fill_diagonal(transfo, -np.diag(transfo))    
        M = geom.Mobius(transfo)
        draw_pair(axes, M, u)
        M = M.inverse()
        draw_pair(axes, M, u)
        d = np.diag(transfo)
        if d[0] != d[1]:
            np.fill_diagonal(transfo, (d[1], d[0]))
            M = geom.Mobius(transfo)
            draw_pair(axes, M, u)
            M = M.inverse()
            draw_pair(axes, M, u)
    plt.xlim(-1.1, 1.1)
    plt.ylim(-1.1, 1.1)
    plt.axis("off")
    pngname = "zpic_%03d.png" % i
    plt.savefig(pngname, dpi=96)
        

os.system(
    "mogrify -resize 512x zpic_*.png"    
) 
os.system(
    "magick convert -dispose previous -loop 0 -delay 7 zpic_*.png ModularTessellation.gif"    
) 