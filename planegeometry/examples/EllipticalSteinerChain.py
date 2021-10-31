# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import planegeometry.geometry as geom


c0 = geom.Circle((3,0), 3) # exterior circle
ell0 = geom.Ellipse((-4,0), 4, 2.5, 140)
f = geom.Affine.from_ellipse_to_ellipse(c0, ell0)

chain = geom.SteinerChain(c0, 3, -0.2, shift=0)
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

plt.xlim(-8, 0)
plt.ylim(-4, 4)
plt.axis("off")
plt.show()
