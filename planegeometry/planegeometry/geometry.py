from math import acos, sqrt, tan
import numpy as np
from .internal import distance_, vlength_, dot_


class Circle:
    def __init__(self, center, radius):
        self.center = np.asarray(center, dtype=float)
        self.radius = radius
        
    def __str__(self):
        return str(self.__dict__)

    def show(self):
        print("Circle:\n")
        print(" center: ", tuple(self.center), "\n")
        print(" radius: ", self.radius, "\n")


class Triangle:
    def __init__(self, A, B, C):
        self.A = np.asarray(A, dtype=float)
        self.B = np.asarray(B, dtype=float)
        self.C = np.asarray(C, dtype=float)
        
    def __str__(self):
        return str(self.__dict__)

    def show(self):
        print("Triangle:\n")
        print("       A: ", tuple(self.A), "\n")
        print("       B: ", tuple(self.B), "\n")
        print("       C: ", tuple(self.C), "\n")
        f = self.flatness()
        if f == 1:
            print("The triangle is flat.\n")
        elif f > 0.99:
            print("The triangle is almost flat (flatness: %s).\n" % f)
    
    def flatness(self):
        "Flatness, a number between 0 and 1; a triangle is flat when its flatness is 1."
        AB = self.B - self.A
        AC = self.C - self.A
        z = complex(*AB) * complex(*AC)
        re = z.real
        im = z.imag
        re2 = re * re
        return re2 / (re2 + im*im)

    def a(self):
        "Length of the side BC."
        return distance_(self.B, self.C)

    def b(self):
        "Length of the side AC."
        return distance_(self.A, self.C)

    def c(self):
        "Length of the side AB."
        return distance_(self.A, self.B)
    
    def edges(self):
        "Edge lengths of the triangle."
        return {"a": self.a(), "b": self.b(), "c": self.c()}
    
    def orientation(self):
        "Orientation of the triangle; 1 for counterclockwise, -1 for clockwise, 0 for collinear."
        A = self.A
        B = self.B
        C = self.C
        val = (B[1] - A[1])*(C[0] - B[0]) - (B[0] - A[0])*(C[1] - B[1])
        return 0 if val==0 else (1 if val>0 else -1)
    
    def contains(self, M):
        "Check whether a point lies inside the reference triangle."
        A0, A1 = self.A
        B0, B1 = self.B
        C0, C1 = self.C
        M0, M1 = M
        dsArea = -B1*C0 + A1*(-B0 + C0) + A0*(B1 - C1) + B0*C1
        s = (A1*C0 - A0*C1 + (C1 - A1)*M0 + (A0 - C0)*M1) / dsArea
        t = (A0*B1 - A1*B0 + (A1 - B1)*M0 + (B0 - A0)*M1) / dsArea
        return s > 0 and t > 0 and 1-s-t > 0
    
    def is_acute(self):
        "Check whether the triangle is acute."
        edges = [self.a(), self.b(), self.c()]
        edges.sort()
        edge0, edge1, edge2 = edges
        return edge0*edge0 + edge1*edge1 >= edge2*edge2
    
    def angleA(self):
        "The angle at the vertex A in radians."
        A = self.A
        B = self.B
        C = self.C
        AC = C - A
        AB = B - A
        b = vlength_(AC)
        c = vlength_(AB)
        return acos(dot_(AC, AB) / b / c)
        
    def angleB(self):
        "The angle at the vertex B in radians."
        A = self.A
        B = self.B
        C = self.C
        BC = C - B
        BA = A - B
        a = vlength_(BC)
        c = vlength_(BA)
        return acos(dot_(BC, BA) / a / c)
        
    def angleC(self):
        "The angle at the vertex C in radians."
        A = self.A
        B = self.B
        C = self.C
        CA = A - C
        CB = B - C
        b = vlength_(CA)
        a = vlength_(CB)
        return acos(dot_(CA, CB) / a / b)
    
    def incircle(self):
        "The incircle of the triangle."
        A = self.A
        B = self.B
        C = self.C
        a = distance_(B, C)
        b = distance_(A, C)
        c = distance_(B, A)
        p = a + b + c
        s = p / 2
        areaABC = sqrt(s*(s-a)*(s-b)*(s-c))
        center = (a*A + b*B + c*C) / p
        radius = areaABC / s
        return Circle(center, radius)
    
    def malfatti_circles(self):
        if self.flatness() == 1:
            return None
        A = self.A
        B = self.B
        C = self.C
        a = distance_(B, C)
        b = distance_(A, C)
        c = distance_(B, A)
        p = a + b + c
        s = p / 2
        smina = s - a
        sminb = s - b
        sminc = s - c
        areaABC = sqrt(s * smina * sminb * sminc)
        I = (a*A + b*B + c*C) / p # incenter
        r = areaABC / s # inradius
        # radii of Malfatti circles ####
        IA = vlength_(I - A)
        IB = vlength_(I - B)
        IC = vlength_(I - C)
        halfr = r / 2
        sminr = s - r
        r1 = halfr * (sminr - (IB+IC-IA)) / smina
        r2 = halfr * (sminr - (IC+IA-IB)) / sminb
        r3 = halfr * (sminr - (IA+IB-IC)) / sminc
        # centers of Malfatti circles ####
        d1 = r1 / tan(acos(dot_(C-A, B-A)/b/c)/2)
        d2 = r2 / tan(acos(dot_(C-B, A-B)/a/c)/2)
        d3 = r3 / tan(acos(dot_(A-C, B-C)/b/a)/2)
        w = d1 + d2 + 2*sqrt(r1*r2)
        u = d2 + d3 + 2*sqrt(r2*r3)
        v = d3 + d1 + 2*sqrt(r3*r1)
        d = sqrt((-u+v+w)*(u+v-w)*(u-v+w)*(u+v+w)) / 2
        x = d/r1 - (v+w)
        y = u
        z = u # trilinear coordinates
        O1 = (u*x*A + v*y*B + w*z*C) / (u*x + v*y + w*z)
        x = v
        y = d/r2 - (u+w)
        z = v # trilinear coordinates
        O2 = (u*x*A + v*y*B + w*z*C) / (u*x + v*y + w*z)
        x = w
        y = w
        z = d/r3 - (u+v) # trilinear coordinates
        O3 = (u*x*A + v*y*B + w*z*C) / (u*x + v*y + w*z)
        circles = {
          "cA": Circle(center = O1, radius = r1),
          "cB": Circle(center = O2, radius = r2),
          "cC": Circle(center = O3, radius = r3)
        }
        O2_O3 = O3 - O2
        O1_O3 = O3 - O1
        O1_O2 = O2 - O1
        tangency_points = {
          "TA": O2 + r2 / vlength_(O2_O3) * O2_O3,
          "TB": O1 + r1 / vlength_(O1_O3) * O1_O3,
          "TC": O1 + r1 / vlength_(O1_O2) * O1_O2
        }
        return {
            "circles": circles,
            "tangency_points": tangency_points
        }
    
    def equal_detour_point(self):
        "Equal detour point of the triangle, also known as the X(176) triangle center."
        A = self.A
        B = self.B
        C = self.C
        a = distance_(B, C)
        b = distance_(A, C)
        c = distance_(B, A)
        s = (a + b + c) / 2
        areaABC = sqrt(s*(s-a)*(s-b)*(s-c))
        abc = np.array([a, b, c])
        v = 2 * areaABC / np.array([b+c-a, c+a-b, a+b-c])
        tc = abc + v # triangular coordinates
        point = (tc[0]*A + tc[1]*B + tc[2]*C) / tc.sum()
        detour = distance_(A, point) + distance_(B, point) - c
        return (point, detour)



