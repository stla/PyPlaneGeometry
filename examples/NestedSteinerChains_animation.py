# -*- coding: utf-8 -*-
import os
import matplotlib.pyplot as plt
import planegeometry.geometry as g
import numpy as np

def draw_circle(axes, C, color):
    axes.add_artist(plt.Circle(
        C.center, C.radius, fill=False, color=color, linewidth=2
    ))

c0 = g.Circle((3,0), 3) # exterior circle
shift_ = np.linspace(0, 3, 101)[:100]
for i, shift in enumerate(shift_):
    circles = g.SteinerChain(c0, 3, -0.5, shift)["circles"]
    # Steiner chain for each circle, except the small one (it is too small)
    chains = [g.SteinerChain(c, 3, -0.4, -shift)["circles"] for c in circles[:3]]
    # draw the big Steiner chain
    figure, axes = plt.subplots(figsize=(10, 10))
    axes.set_aspect(1)
    [draw_circle(axes, circ, "blue") for circ in circles]
    draw_circle(axes, c0, "black")
    # draw the nested Steiner chain
    [[draw_circle(axes, c, "red") for c in chain] for chain in chains]
    plt.xlim(0, 6)
    plt.ylim(-4, 4)
    plt.axis("off")
    pngname = "zpic_%03d.png" % i
    plt.savefig(pngname, dpi=96)

os.system(
    "mogrify -resize 512x zpic_*.png"    
) 
os.system(
    "magick convert -dispose previous -loop 0 -delay 8 zpic_*.png NestedSteinerChains.gif"    
) 
