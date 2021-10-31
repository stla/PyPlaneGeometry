# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import planegeometry.geometry as g


c0 = g.Circle((3,0), 3) # exterior circle
circles = g.SteinerChain(c0, 3, -0.2, 0.5)["circles"]
# Steiner chain for each circle, except the small one (it is too small)
chains = [g.SteinerChain(c, 3, -0.2, 0.5)["circles"] for c in circles[:3]]
# draw the big Steiner chain
figure, axes = plt.subplots(figsize=(10, 10))
axes.set_aspect(1)
def draw_circle(C, color):
    axes.add_artist(plt.Circle(C.center, C.radius, fill=False, color=color))
[draw_circle(circ, "blue") for circ in circles]
draw_circle(c0, "black")
# draw the nested Steiner chain
[[draw_circle(c, "red") for c in chain] for chain in chains]
plt.title("Nested Steiner chains")
plt.xlim(0, 6)
plt.ylim(-4, 4)
plt.axis("off")
plt.show()

