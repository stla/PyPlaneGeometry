from planegeometry.geometry import Triangle, Circle, intersection_circle_circle
from numpy import allclose, array

def test_intersection_circle_circle():
      # two intersection points
      A = (0, 0) 
      B = (2, 0)
      C = (3, 2)
      D = (3,-3)
      ABC = Triangle(A, B, C)
      ABD = Triangle(A, B, D)
      c1 = ABC.circumcircle()
      c2 = ABD.circumcircle()
      Is = intersection_circle_circle(c1, c2)
      assert type(Is) == list
      assert (allclose(Is[0], A)) or (allclose(Is[0], B))
      assert (allclose(Is[1], A)) or (allclose(Is[1], B))
      # one intersection point
      c1 = Circle((0,0), 2)
      c2 = Circle((5,0), 3)
      I = intersection_circle_circle(c1, c2)
      assert allclose(array([2,0]), I)
      # no intersection
      ## c1 and c2 external
      c1 = Circle((0,0), 2)
      c2 = Circle((5,0), 1)
      assert intersection_circle_circle(c1, c2) is None
      ## c1 included in c2
      c1 = Circle((4,0), 2)
      c2 = Circle((5,0), 5)
      assert intersection_circle_circle(c1, c2) is None
