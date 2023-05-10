# -*- coding: utf-8 -*-
import os
import matplotlib.pyplot as plt
import planegeometry.geometry as geom
import numpy as np

c0 = geom.Circle((3,0), 3) # exterior circle
ell0 = geom.Ellipse((-4,0), 4, 2.5, 140)
f = geom.Affine.from_ellipse_to_ellipse(c0, ell0)

shift_ = np.linspace(0, 3, 101)[:100]

for i, shift in enumerate(shift_):
    chain = geom.SteinerChain(c0, 3, -0.2, shift=shift)
    circles = chain["circles"]
    ellipses = [f.transform_ellipse(c) for c in circles]
    figure, axes = plt.subplots(figsize=(10, 10))
    axes.set_aspect(1)
    def draw_ellipse(ell, color, fill):
        ellipse_path = ell.path()
        axes.add_artist(
            plt.Polygon(
                ellipse_path, closed=True, fill=fill, facecolor=color, 
                edgecolor="black", linewidth=2
            )
        )
    [draw_ellipse(ell, "blue", True) for ell in ellipses]
    draw_ellipse(ell0, "black", False)
    plt.xlim(-7.5, 0)
    plt.ylim(-3.5, 3.5)
    plt.axis("off")
    pngname = "zpic_%03d.png" % i
    plt.savefig(pngname, dpi=96)

os.system(
    "mogrify -resize 512x zpic_*.png"    
) 
os.system(
    "magick convert -dispose previous -loop 0 -delay 10 zpic_*.png EllipticalSteinerChain.gif"    
) 