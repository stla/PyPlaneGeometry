from math import acos, sqrt
import numpy as np
from .internal import distance_, vlength_, dot_


class Triangle:
    def __init__(self, A, B, C):
        self.A = np.asarray(A, dtype=float)
        self.B = np.asarray(B, dtype=float)
        self.C = np.asarray(C, dtype=float)

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
        return (center, radius)

    


        

