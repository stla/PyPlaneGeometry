# -*- coding: utf-8 -*-
from math import sqrt, cos, sin, tan
from sympy.combinatorics.free_groups import free_group
from itertools import product
import numpy as np
from functools import reduce
from planegeometry.geometry import Mobius, Circle, Line
from planegeometry.internal import dot_, distance_
import matplotlib.pyplot as plt

F, a, b = free_group("a, b")

A = a.inverse()
B = b.inverse()

elems = [a, A, b, B]

n = 6
points = product(elems, repeat=n) 
x = list(points)

def multiply(list_of_symbols):
    return reduce(lambda u, v: u*v, list_of_symbols)

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

transfos = unique_with(list(map(multiply, x)), lambda u, v: u == v)

def total(transfo):
    dec = transfo.array_form
    powers = list(map(lambda x: abs(x[1]), transfo))
    return sum(powers)

totals = list(map(total, transfos))

sizes = np.asarray(totals, dtype=int)
indices = np.where(np.equal(sizes, n))[0]
Gn = [transfos[i] for i in indices.tolist()]


# starting circles ####
Ca = Line((-1,0), (1,0))
Rc = sqrt(2)/4
yI = -3*sqrt(2)/4
CA = Circle((0,yI), Rc)
theta = -0.5
T = np.array([Rc*cos(theta), yI+Rc*sin(theta)])
P = np.array([T[0]+T[1]*tan(theta), 0])
PT = distance_(T, P)
xTprime = P[0] + PT
xPprime = -yI/tan(theta)
PprimeTprime = abs(xTprime-xPprime)
Rcprime = abs(yI*PprimeTprime/xPprime)
Cb = Circle((xTprime, -Rcprime), Rcprime)
CB = Circle((-xTprime, -Rcprime), Rcprime)

Moba = Mobius(np.array([[sqrt(2), 1j], [-1j, sqrt(2)]]))
Mobb = Mobius([
    [complex(*Cb.center), dot_(Cb.center)-Cb.radius**2],
    [1, -complex(*CB.center)]
])
MobA = Moba.inverse()
MobB = Mobb.inverse()

Mobs = {
    "a": Moba,
    "b": Mobb,
    "A": MobA,
    "B": MobB
}

Circles = {
    "a": Ca,
    "b": Cb,
    "A": CA,
    "B": CB
}

def transfo2seq(transfo):
    seq = []
    tup = transfo.array_form
    for t, i in tup:
        t = str(t)
        if i < 0:
            i = -i
            t = str.upper(t)
        seq = seq + [t for k in range(i)]
    return seq

def circle(transfo):
    seq = transfo2seq(transfo)
    mobs = list(map(lambda x: Mobs[x], seq))
    mobius = reduce(lambda M1, M2: M1.compose(M2), mobs[:(n-1)])
    if seq[-1] == "a":
        return mobius.transform_line(Circles["a"])
    return mobius.transform_circle(Circles[seq[-1]])


figure, axes = plt.subplots(facecolor="black", figsize=(10, 10))
axes.set_aspect(1)
axes.axline(Ca.A, Ca.B, linewidth=2, color="black")
for C in [CA, Cb, CB]:
    axes.add_artist(
        plt.Circle(
            C.center, C.radius, fill=False, edgecolor="black", linewidth=2
        )
    )
C1 = MobA.transform_circle(CA)
C2 = MobA.transform_circle(CB)
C3 = MobA.transform_circle(Cb)
for C in [C1, C2, C3]:
    axes.add_artist(
        plt.Circle(
            C.center, C.radius, fill=False, edgecolor="red", linewidth=2
        )
    )
C1 = Moba.transform_line(Ca)
C2 = Moba.transform_circle(CB)
C3 = Moba.transform_circle(Cb)
for C in [C1, C2, C3]:
    axes.add_artist(
        plt.Circle(
            C.center, C.radius, fill=False, edgecolor="green", linewidth=2
        )
    )
C1 = Mobb.transform_line(Ca)
C2 = Mobb.transform_circle(Cb)
C3 = Mobb.transform_circle(CA)
for C in [C1, C2, C3]:
    axes.add_artist(
        plt.Circle(
            C.center, C.radius, fill=False, edgecolor="blue", linewidth=2
        )
    )
C1 = MobB.transform_line(Ca)
C2 = MobB.transform_circle(CA)
C3 = MobB.transform_circle(CB)
for C in [C1, C2, C3]:
    axes.add_artist(
        plt.Circle(
            C.center, C.radius, fill=False, edgecolor="yellow", linewidth=2
        )
    )
for g in Gn:
    circ = circle(g)
    axes.add_artist(
        plt.Circle(
            circ.center, circ.radius, fill=False, edgecolor="orange", linewidth=2
        )
    )
plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.axis("off")
plt.show()
