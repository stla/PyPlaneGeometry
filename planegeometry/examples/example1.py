import planegeometry.geometry as g

A = (0, 0)
B = (1, 0)
C = (2, 1)
t = g.Triangle(A, B, C)
print("flatness", t.flatness)
t.show()
print(t.equal_detour_point())
print(t.incircle())

print("******************************")
ell1 = g.Ellipse((1,1), 5, 1, 30)
ell2 = g.Ellipse((4,-1), 3, 2, 50) 
f = g.Affine.from_ellipse_to_ellipse(ell1, ell2)
print(f.transform_ellipse(ell1)) # should be ell2
