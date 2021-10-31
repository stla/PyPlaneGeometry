from math import acos, sqrt, tan, atan2, cos, sin, pi, inf, isinf
import numpy as np
from .internal import (
    distance_,
    vlength_,
    dot_,
    det2x2_,
    line_line_intersection_,
    unit_vector_,
    epsilon_,
    ellipse_points_,
    collinear_,
    circle_points_,
    circle_as_ellipse_
)


class Line:
    def __init__(self, A, B, extendA=True, extendB=True):
        self.A = np.asarray(A, dtype=float)
        self.B = np.asarray(B, dtype=float)
        self.extendA = extendA
        self.extendB = extendB
        
    def __str__(self):
        return str(self.__dict__)

    def show(self):
        print("Line:\n")
        print("       A: ", tuple(self.A), "\n")
        print("       B: ", tuple(self.B), "\n")
        print(" extendA: ", self.extendA, "\n")
        print(" extendB: ", self.extendB, "\n")
        if self.extendA and self.extendB:
            print("Infinite line passing through A and B.\n")
        elif self.extendA:
            print("Half-line with origin B and passing through A.\n")
        elif self.extendB:
            print("Half-line with origin A and passing through B.\n")
        else:
            print("Segment joining A and B.\n")
        
    def direction_offset(self):
        A = self.A
        B = self.B
        if A[0]==B[0]:
            if A[0]>0:
                direction = 0.0
                offset = A[0]
            else:
                direction = pi
                offset = -A[0]
        else:
            x = B[0] - A[0]
            y = B[1] - A[1]
            theta = -atan2(x, y)
            offset = A[0]*cos(theta) + A[1]*sin(theta)
            if offset < 0:
                theta = theta + pi
                offset = -offset
            direction = theta % (2*pi)
        return {
            "direction": direction,
            "offset": offset
        }
    
    def is_equal(self, line2):
        do1 = self.direction_offset()
        do2 = line2.direction_offset()
        do1 = (do1["direction"] % pi, do1["offset"])
        do2 = (do2["direction"] % pi, do2["offset"])
        return np.allclose(do1, do2)
    
    def is_parallel(self, line2):
        P1 = self.A
        P2 = self.B
        Q1 = line2.A
        Q2 = line2.B
        dx1 = P1[0] - P2[0]
        dx2 = Q1[0] - Q2[0]
        dy1 = P1[1] - P2[1]
        dy2 = Q1[1] - Q2[1]
        D = det2x2_((dx1, dy1), (dx2, dy2))
        return abs(D) < sqrt(epsilon_)
    
    def includes(self, M, strict = False, checkCollinear = True):
        """Check whether a point belongs to the line.
        
        :param M: the point for which we want to test whether it belongs to the line
        :param strict: Boolean, whether to take into account `extendA` and `extendB`
        :param checkCollinear: Boolean, whether to check the collinearity of `A`, `B`, and `M`; set to `False` only if you use `strict=True` and you are sure that `M` is on the line (AB)
        :returns: A Boolean value.        
        
        """
        A = self.A
        B = self.B
        if checkCollinear:
            if not collinear_(A, B, M):
                return False
        extendA = self.extendA
        extendB = self.extendB
        if (not strict) or (extendA and extendB):
            return collinear_(A, B, M)
        if (not extendA) and (not extendB):
            dotprod = dot_(A-M, B-M)
            if dotprod <= 0:
                return True
            print("The point is on the line (AB), but not on the segment [AB].")
            return False
        if extendA:
            if np.any((M-B)*(A-B)>0):
                return True
            print("The point is on the line (AB), but not on the half-line (AB].")
            return False
        # extendB
        if np.any((M-A)*(B-A)>0):
            return True
        print("The point is on the line (AB), but not on the half-line [AB).")
        return False
    
    def perpendicular(self, M, extendH = False, extendM = True):
        """Perpendicular line passing through a given point.
        
        :param M: the point through which the perpendicular passes
        :param extendH: Boolean, whether to extend the perpendicular line beyond the meeting point
        :param extendM: Boolean,  whether to extend the perpendicular line beyond the point `M`
        :returns: A `Line` object; its two points are the meeting point and the point `M`.
        
        """
        A = self.A
        B = self.B
        A_B = B - A
        v = np.array([-A_B[1], A_B[0]])
        if self.includes(M):
            print("M is on the line.")
            return Line(M, M + v, True, True)
        H = line_line_intersection_(A, B, M - v, M + v)
        return Line(H, M, extendH, extendM)



def intersection_line_line(line1, line2, strict=False):
    if line1.is_equal(line2):
        print("TO DO")
        return
    if line1.is_parallel(line2):
        print("Distinct parallel lines.")
        return None
    I = line_line_intersection_(line1.A, line1.B, line2.A, line2.B)
    if not strict:
        return I
    print("TO DO")
    return



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
        
    def orthogonalThroughTwoPointsWithinCircle(self, P1, P2, arc = False):
        """Orthogonal circle passing through two points within the reference circle.
        
        :param P1,P2: two distinct points in the interior of the reference circle
        :param arc: Boolean,  whether to return the arc joining the two points instead of the circle
        :returns: A `Circle` object or an `Arc` object, or a `Line` object if the two points are on a diameter.

        """
        if np.allclose(P1, P2):
            print("`P1` and `P2` must be distinct.")
            return
        I = self.center
        r = self.radius
        r2 = r * r
        if distance_(P1, I) >= r2:
            print("`P1` is not in the interior of the reference circle.")
            return
        if distance_(P2,I) >= r2:
            print("`P2` is not in the interior of the reference circle.")
            return
        if collinear_(I, P1, P2):
            return Line(P1, P2, not arc, not arc)
        iota = Inversion(I, r2)
        P1prime = iota.invert(P1)
        P2prime = iota.invert(P2)
        line1 = Line(P1, P1prime)
        line2 = Line(P2, P2prime)
        perp1 = line1.perpendicular((P1+P1prime)/2)
        perp2 = line2.perpendicular((P2+P2prime)/2)
        O = line_line_intersection_(perp1.A, perp1.B, perp2.A, perp2.B)
        if arc:
            theta1 = atan2(P1[1]-O[1], P1[0]-O[0]) % (2*pi)
            theta2 = atan2(P2[1]-O[1], P2[0]-O[0]) % (2*pi)
            return Arc(
                O, distance_(O, P1),
                min(theta1, theta2), max(theta1, theta2), False
            )
        return Circle(O, distance_(O, P1))
    
    def orthogonalThroughTwoPointsOnCircle(self, alpha1, alpha2, arc = False):
        """
        
        
        """
        I = self.center
        r = self.radius
        dalpha = alpha1 - alpha2
        if dalpha % pi == 0:
            eialpha1 = unit_vector_(alpha1)
            A = I + r*eialpha1
            B = I - r*eialpha1
            return Line(A, B, not arc, not arc)
        r0 = r * abs(tan(dalpha/2))
        IO = r / cos(dalpha/2)
        center = I + IO * unit_vector_((alpha1+alpha2)/2)
        if arc:
            dalpha = (alpha2-alpha1) % (2*pi)# - alpha1%%(2*pi)
            delta = pi if dalpha >= pi else 0
            beta1 = -pi/2 + delta
            beta2 = beta1 - pi + dalpha
            theta1 = beta1+alpha1 #%% (2*pi)
            theta2 = beta2+alpha1 #%% (2*pi)
            return Arc(
                center, r0, min(theta1,theta2), max(theta1,theta2), False
            )
        return Circle(center, r0)

class Arc:
    """Arc class.
    
    A circular arc is initialized by its center (array-like object of length 
    two), a radius, a starting angle and an ending angle. They are 
    respectively named `center`, `radius`, `alpha1` and `alpha2`.
    
    """
    def __init__(self, center, radius, alpha1, alpha2, degrees = True):
        self.center = np.asarray(center, dtype=float)
        self.radius = radius
        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.degrees = degrees

    def __str__(self):
        return str(self.__dict__)

    # def show(self):
    #     unit = "degree" if self.degrees else "radian"
    #     s = "" if alpha == 1 else "s"
    #     print("Ellipse:\n")
    #     print("       center: ", tuple(self.center), "\n")
    #     print(" major radius: ", self.rmajor, "\n")
    #     print(" minor radius: ", self.rminor, "\n")
    #     print("        angle: ", self.alpha, unit + s, "\n")
    
    def path(self, n_points=100):
        """Path that forms the arc.
        
        :param n_points: number of points of the path
        :returns: A matrix with two columns and `n_points` rows.
        
        """
        center = self.center
        radius = self.radius
        alpha1 = self.alpha1
        alpha2 = self.alpha2
        if self.degrees:
            alpha1 = alpha1 * pi/180
            alpha2 = alpha2 * pi/180
        alpha1 = alpha1 % (2*pi)
        alpha2 = alpha2 % (2*pi)
        dalpha = alpha2 - alpha1
        if dalpha < 0:
            dalpha += 2*pi
        theta = np.linspace(0, dalpha, n_points)
        return circle_points_(
            alpha1 + theta,
            center,
            radius
        )


class Ellipse:
    """Ellipse class.
    
    An ellipse is initialized by its center (array-like object of length two), 
    its major radius, its minor radius, and the angle `alpha` between the 
    x-axis and the major axis.
    
    """
    def __init__(self, center, rmajor, rminor, alpha, degrees = True):
        self.center = np.asarray(center, dtype=float)
        self.rmajor = rmajor
        self.rminor = rminor
        self.alpha = alpha
        self.degrees = degrees

    def __str__(self):
        return str(self.__dict__)

    def show(self):
        unit = "degree" if self.degrees else "radian"
        s = "" if self.alpha == 1 else "s"
        print("Ellipse:\n")
        print("       center: ", tuple(self.center), "\n")
        print(" major radius: ", self.rmajor, "\n")
        print(" minor radius: ", self.rminor, "\n")
        print("        angle: ", self.alpha, unit + s, "\n")
    
    def path(self, n_points=100):
        """Path that forms the ellipse.
        
        :param n_points: number of points of the path
        :returns: A matrix with two columns and `n_points` rows.
        
        """
        center = self.center
        alpha = self.alpha
        if self.degrees:
            alpha = alpha * pi/180
        return ellipse_points_(
            np.linspace(0, 2*pi, n_points+1)[:n_points],
            center,
            self.rmajor,
            self.rminor,
            alpha
        )
    
    def equation(self):
        """The coefficients of the implicit equation of the ellipse, 
        `Ax² + Bxy + Cy² + Dx + Ey + F = 0`.
        
        :returns: A dictionary giving the values of the coefficients.
        
        """
        a2 = self.rmajor**2
        b2 = self.rminor**2
        alpha = self.alpha
        if self.degrees:
            alpha *= pi/180
        x, y = self.center
        sine = sin(alpha)
        cosine = cos(alpha)
        sine2 = sine*sine
        cosine2 = 1 - sine2
        A = a2*sine2 + b2*cosine2
        B = 2*(b2-a2)*sine*cosine
        C = a2*cosine2 + b2*sine2
        return {
            "A": A,
            "B": B,
            "C": C,
            "D": -2*A*x - B*y,
            "E": -B*x - 2*C*y,
            "F": A*x*x + B*x*y + C*y*y - a2*b2
        }        



class Inversion:
    """Inversion class.
    
    An inversion is initialized by its pole and its power.
    
    """
    def __init__(self, pole, power):
        self.pole = np.asarray(pole, dtype=float)
        self.power = power
        
    def __str__(self):
        return str(self.__dict__)

    def show(self):
        print("Inversion:\n")
        print("      pole: ", tuple(self.pole), "\n")
        print("     power: ", self.power, "\n")
    
    def invert(self, M):
        pole = self.pole
        if np.isscalar(M) and isinf(M):
            return pole
        M = np.asarray(M)
        if np.allclose(pole, M):
            return inf
        k = self.power
        pole_M = M - pole
        return pole + k/dot_(pole_M) * pole_M
    
    def invert_circle(self, circ):
        c0 = self.pole
        k = self.power
        c1 = circ.center
        r1 = circ.radius
        D1 = (c1[0] - c0[0])**2 + (c1[1] - c0[1])**2 - r1*r1
        if abs(D1) > sqrt(epsilon_):
            s = k / D1
            return Circle(c0 + s*(c1-c0), abs(s)*r1)
        Ot = c0 - c1
        R180 = c1 - Ot
        R90 = np.array([-Ot[1], Ot[0]]) + c1
        return Line(self.invert(R180), self.invert(R90))
        


def SteinerChain_phi0_(c0, n, shift):
    R = c0.radius
    O = c0.center
    pi_over_n = pi / n
    sine = sin(pi_over_n)
    Cradius = R / (1 + sine)
    Cside = Cradius * sine
    circles0 = [None]*(n+1)
    for i in range(n):
        beta = (i + shift) * 2*pi_over_n
        pti = Cradius * unit_vector_(beta) + O
        circles0[i] = Circle(pti, Cside)
    circles0[n] = Circle(O, R - 2*Cside)
    return circles0

def SteinerChain(c0, n, phi, shift):
    circles = SteinerChain_phi0_(c0 = c0, n = n, shift = shift)
    R = c0.radius
    O = c0.center
    if phi != 0:
        invphi = 1 / phi
        I = np.array([R*invphi, 0]) + O
        r2 = R*R * (invphi*invphi - 1)
        iota = Inversion(I, r2)
        circles = [iota.invert_circle(circle) for circle in circles]
    # the ellipse of centers
    inner = circles[-1]
    O2 = inner.center
    r = inner.radius
    c = (O2[0] - O[0])/2
    a = (r + R)/2
    b = sqrt(a*a - c*c)
    ell = Ellipse((O+O2)/2, a, b, 0)
    return {
        "circles": circles,
        "ellipse": ell
    }



class Triangle:
    """Triangle class.
    
    A triangle is initialized by its three vertices, some array-like objects 
    of length two.
    
    """
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
        f = self.flatness
        if f == 1:
            print("The triangle is flat.\n")
        elif f > 0.99:
            print("The triangle is almost flat (flatness: %s).\n" % f)
    
    @property
    def flatness(self):
        "Flatness, a number between 0 and 1; a triangle is flat when its flatness is 1."
        AB = self.B - self.A
        AC = self.C - self.A
        z = complex(*AB) * complex(*AC)
        re = z.real
        im = z.imag
        re2 = re * re
        return re2 / (re2 + im*im)

    @property
    def a(self):
        "Length of the side BC."
        return distance_(self.B, self.C)

    @property
    def b(self):
        "Length of the side AC."
        return distance_(self.A, self.C)

    @property
    def c(self):
        "Length of the side AB."
        return distance_(self.A, self.B)

    @property    
    def edges(self):
        "Edge lengths of the triangle."
        return {"a": self.a(), "b": self.b(), "c": self.c()}

    @property    
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

    @property    
    def is_acute(self):
        "Check whether the triangle is acute."
        edges = [self.a(), self.b(), self.c()]
        edges.sort()
        edge0, edge1, edge2 = edges
        return edge0*edge0 + edge1*edge1 >= edge2*edge2

    @property    
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

    @property        
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

    @property        
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
        if self.flatness == 1:
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


class Affine:
    """A class for affine transformations.
    
    An affine transformation is initialized by a 2x2 matrix (a linear transformation), 
    and a length two vector (the 'intercept', an array-like object).
    
    """
    def __init__(self, A, b):
        self.A = np.asarray(A, dtype=float)
        self.b = np.asarray(b, dtype=float)
        
    def __str__(self):
        return str(self.__dict__)

    def show(self):
        print("Affine map Ax + b:\n")
        print("                 A: ", self.A, "\n")
        print("                 b: ", self.b, "\n")
    
    def get3x3matrix(self):
        """Get the 3x3 matrix corresponding to the affine transformation.
        
        """
        b = self.b.reshape(2,1) 
        return np.vstack((np.hstack((self.A, b)), np.array([0, 0, 1])))
    
    def inverse(self):
        """The inverse affine transformation if it exists.
        
        """
        if np.linalg.det(self.A) == 0:
            print("The affine map is singular.")
            return
        M = np.linalg.inv(self.get3x3matrix())
        return Affine(M[0:2, 0:2], M[0:2, -1])
    
    def compose(self, transfo, left=True):
        """Compose the reference affine map with another affine map.
        
        :param transfo: an `Affine` object
        :param left: Boolean, whether to compose at left or at right (i.e. returns `f1 o f0` or `f0 o f1`)
        :returns: An `Affine` object.
        
        """
        M0 = self.get3x3matrix()
        M1 = transfo.get3x3matrix()
        M = np.matmul(M1, M0) if left else np.matmul(M0, M1)
        return Affine(M[0:2, 0:2], M[0:2, -1])
    
    def transform(self, m):
        """Transform a point or several points by the affine map.
        
        :param m: a point or a two-column matrix of points, one point per row
        :returns: a matrix or a vector
        
        """
        b = self.b
        if m.ndim == 2:
            b = np.repeat(b.reshape(2,1), m.shape[1], axis=1)
        return np.transpose(np.matmul(self.A, np.transpose(m)) + b)
    
    def transform_line(self, line):
        """Transform a line by the affine map.
        
        
        """
        M = self.A
        if np.linalg.det(M) == 0:
            print("The affine map is singular.")
            return
        return Line(M.dot(line.A), M.dot(line.B), line.extendA, line.extendB)
    
    def transform_ellipse(self, ell):
        """Transform an ellipse by the reference affine transformation (only for an invertible affine map).
        
        :param ell: an `Ellipse` object or a `Circle` object
        :returns: An `Ellipse` object.
        
        """
        if np.linalg.det(self.A) == 0:
            print("The affine map is singular.")
            return
        if isinstance(ell, Circle):
            ell = circle_as_ellipse_(ell)
        A, B, C, D, E, F = ell.equation().values()
        X = np.array([
            [A, B/2, D/2],
            [B/2, C, E/2],
            [D/2, E/2, F]
        ])
        Mat = np.linalg.inv(self.get3x3matrix())
        Y = np.matmul(np.matmul(np.transpose(Mat), X), Mat) 
        A = Y[0,0]
        B = 2 * Y[0,1]
        C = Y[1,1]
        D = 2 * Y[0,2]
        E = 2 * Y[1,2]
        F = Y[2,2]
        Delta = B*B - 4*A*C
        s = sqrt((A-C)**2 + B*B)
        V = np.array([s, -s])
        a, b = - np.sqrt(
            2 * (A*E*E + C*D*D - B*D*E + Delta*F)*((A + C) + V)
        ) / Delta 
        x0 = (2*C*D - B*E) / Delta
        y0 = (2*A*E - B*D) / Delta
        theta = atan2(C - A - s, B)
        degrees = ell.degrees
        theta = (theta*180/pi) % 180 if degrees else (theta % pi)
        return Ellipse((x0, y0), a, b, theta, degrees = degrees)
    
    @classmethod
    def from_three_points(cls, P1, P2, P3, Q1, Q2, Q3):
        """Affine transformation mapping three given points to three given points.
        
        :param P1,P2,P3: three non-collinear points
        :param Q1,Q2,Q3: three non-collinear points
        :returns: An `Affine` object representing the transformation which maps Pi to Qi for each i=1,2,3.
        
        """
        P1 = np.asarray(P1, dtype=float)
        P2 = np.asarray(P2, dtype=float)
        P3 = np.asarray(P3, dtype=float)
        Q1 = np.asarray(Q1, dtype=float)
        Q2 = np.asarray(Q2, dtype=float)
        Q3 = np.asarray(Q3, dtype=float)
        if collinear_(P1, P2, P3):
            print("P1, P2 and P3 are collinear.")
            return
        if collinear_(Q1, Q2, Q3):
            print("Q1, Q2 and Q3 are collinear.")
            return
        f1 = Affine(np.transpose(np.column_stack((P2-P1, P3-P1))), P1)
        f2 = Affine(np.transpose(np.column_stack((Q2-Q1, Q3-Q1))), Q1)
        return f1.inverse().compose(f2)
    
    @classmethod
    def from_ellipse_to_ellipse(cls, ell1, ell2):
        """Affine transformation mapping a given ellipse to a given ellipse.
        
        :param ell1,ell2: `Ellipse` or `Circle` objects
        :returns: An `Affine` object representing the transformation which maps `ell1` to `ell2`.
        
        """
        if isinstance(ell1, Circle):
            a = b = ell1.radius
            costheta = 1
            sintheta = 0
        else:
            a = ell1.rmajor
            b = ell1.rminor
            theta = ell1.alpha
            if ell1.degrees:
                theta *= pi/180
            costheta = cos(theta)
            sintheta = sin(theta)
        col1 = a * np.array([costheta, sintheta])
        col2 = b * np.array([-sintheta, costheta])
        f1 = Affine(np.column_stack((col1, col2)), ell1.center)
        if isinstance(ell2, Circle):
            a = b = ell2.radius
            costheta = 1
            sintheta = 0
        else:
            a = ell2.rmajor
            b = ell2.rminor
            theta = ell2.alpha
            if ell2.degrees:
                theta *= pi/180
            costheta = cos(theta)
            sintheta = sin(theta)
        col1 = a * np.array([costheta, sintheta])
        col2 = b * np.array([-sintheta, costheta])
        f2 = Affine(np.column_stack((col1, col2)), ell2.center)
        return f1.inverse().compose(f2)
