# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import planegeometry.geometry as g


c0 = g.Circle((3,0), 3) # exterior circle
chain = g.SteinerChain(c0, 3, 0.6, 0.85)
circles = chain["circles"]
ellipse = chain["ellipse"]
ellipse_path = ellipse.path()


figure, axes = plt.subplots(figsize=(10, 10))
axes.set_aspect(1)

def draw_circle(C, color):
    axes.add_artist(plt.Circle(C.center, C.radius, fill=False, color=color))

[draw_circle(circ, "blue") for circ in circles]
draw_circle(c0, "black")

axes.add_artist(
    plt.Polygon(ellipse_path, closed=True, fill=False, color="red")
)

plt.title(
    "The red path is an ellipse which\npasses through the circle centers."
)
plt.xlim(0, 6)
plt.ylim(-4, 4)
plt.axis("off")
plt.show()

