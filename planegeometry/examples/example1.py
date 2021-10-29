import planegeometry.geometry as g

A = (0, 0)
B = (1, 0)
C = (2, 1)
t = g.Triangle(A, B, C)
#print(t.flatness())
t.show()
print(t.equal_detour_point())
print(t.incircle())
