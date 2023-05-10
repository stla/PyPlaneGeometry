# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import planegeometry.geometry as geom


c0 = geom.Circle((7,0), 3) # exterior circle
ell0 = geom.Ellipse((0,0), 4, 2.5, 140)
f = geom.Affine.from_ellipse_to_ellipse(c0, ell0)

circles = geom.SteinerChain(c0, 3, -0.15, 0.2)["circles"]
# Steiner chain for each circle, except the small one (it is too small)
chains = [geom.SteinerChain(c, 3, -0.4, -0.2)["circles"] for c in circles[:3]]
# draw the big Steiner chain
figure, axes = plt.subplots(figsize=(10, 10))
axes.set_aspect(1)

def draw_ellipse(ell, color, fill=False):
    ellipse_path = ell.path()
    axes.add_artist(
        plt.Polygon(
            ellipse_path, closed=True, fill=fill, color=color, 
            linewidth=2
        )
    )

[draw_ellipse(f.transform_ellipse(circ), "blue") for circ in circles]
draw_ellipse(ell0, "black")
# draw the nested Steiner chain
[[draw_ellipse(f.transform_ellipse(c), "red") for c in chain] for chain in chains]
plt.title("Elliptical nested Steiner chains", fontdict = {"fontsize": 30})
plt.xlim(-4, 4)
plt.ylim(-4, 4)
plt.axis("off")
plt.show()

