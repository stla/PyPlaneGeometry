# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from planegeometry.geometry import (
    Circle, 
    Inversion, 
    Triangle, 
    Line, 
    intersection_circle_circle, 
    intersection_line_line,
    circleAB)


C1 = Circle((0,0), 2)
C2 = Circle((5,5), 3)
C3 = Circle((6,-2), 1)
# inversion swapping C1 and C3 with positive power
iota1 = Inversion.from_swapping_two_circles(C1, C3, positive = True)
# inversion swapping C2 and C3 with positive power
iota2 = Inversion.from_swapping_two_circles(C2, C3, positive = True)
# take an arbitrary point on C3
M = C3.point_from_angle(0)
# invert it with iota1 and iota2
M1 = iota1.invert(M)
M2 = iota2.invert(M)
# take the circle C passing through M, M1, M2
C = Triangle(M,M1,M2).circumcircle()
# take the line passing through the two inversion poles
cl = Line(iota1.pole, iota2.pole)
# take the radical axis of C and C3
L = C.radical_axis(C3)
# let H bet the intersection of these two lines
H = intersection_line_line(L, cl)
# take the circle Cp with diameter [HO3]
O3 = C3.center
Cp = circleAB(H, O3)
# get the two intersection points T0 and T1 of C3 with Cp
T0, T1 = intersection_circle_circle(C3, Cp)
# invert T0 with respect to the two inversions
T0p = iota1.invert(T0)
T0pp = iota2.invert(T0)
# the circle passing through T0 and its two images is a solution
Csolution0 = Triangle(T0, T0p, T0pp).circumcircle()
# invert T1 with respect to the two inversions
T1p = iota1.invert(T1)
T1pp = iota2.invert(T1)
# the circle passing through T1 and its two images is another solution
Csolution1 = Triangle(T1, T1p, T1pp).circumcircle()


def draw_circle(axes, C, bordercolor, fillcolor=None):
    if fillcolor is None:
        axes.add_artist(plt.Circle(
            C.center, C.radius, fill=False, edgecolor=bordercolor, linewidth=2
        ))
    else:
        axes.add_artist(plt.Circle(
            C.center, C.radius, fill=True, edgecolor=bordercolor, 
            facecolor=fillcolor, linewidth=2
        ))

figure, axes = plt.subplots(figsize=(10, 10))
axes.set_aspect(1)
draw_circle(axes, C1, "red", "yellow")
draw_circle(axes, C2, "red", "yellow")
draw_circle(axes, C3, "red", "yellow")
draw_circle(axes, Csolution0, "blue")
draw_circle(axes, Csolution1, "blue")
plt.xlim(-4, 9)
plt.ylim(-4, 9)
plt.axis("off")
plt.show()
