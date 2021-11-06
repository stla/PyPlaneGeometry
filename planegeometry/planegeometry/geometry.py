from math import (
    acos, 
    sqrt, 
    tan, 
    atan2, 
    cos, 
    sin, 
    pi, 
    inf, 
    isinf, 
    isclose,
    atan
)
import numpy as np
from functools import reduce
from .internal import (
    distance_,
    vlength_,
    dot_,
    det2x2_,
    line_line_intersection_,
    circle_line_intersection_, 
    unit_vector_,
    sepsilon_,
    ellipse_points_,
    collinear_,
    circle_points_,
    det2x2_mat_,
    from_complex_,
    mod2_,
    farey_stack_,
    is_number_, 
    is_vector_, 
    is_real_vector_, 
    error_if_not_point_,
    error_if_not_number_,
    error_if_not_positive_,
    error_if_not_boolean_,
    runif_on_circle_,
    runif_in_circle_,
    runif_on_ellipse_,
    runif_in_ellipse_,
    ellipse_from_center_and_eigen_
)

def is_inf(x):
    return isinstance(x, float) and isinf(x)

class Line:
    """A class for lines. A line is initialized by two points it passes through, 
    and for each of these points a Boolean value to indicate whether the line 
    should be extended besides this point.
    
    """
    def __init__(self, A, B, extendA=True, extendB=True):
        _ = error_if_not_point_(A=A)
        _ = error_if_not_point_(B=B)
        self.A = np.asarray(A, dtype=float)
        self.B = np.asarray(B, dtype=float)
        if np.allclose(self.A, self.B):
            raise ValueError("`A` and `B` must be distinct.")
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
        """Direction and offset of the line.
        
        """
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
        """Check whether the reference line is equal to another line.
        
        """
        do1 = self.direction_offset()
        do2 = line2.direction_offset()
        do1 = (do1["direction"], do1["offset"])
        do2 = (do2["direction"], do2["offset"])
        return np.allclose(do1, do2)
    
    def is_parallel(self, line2):
        """Check whether the reference line is parallel to another line.
        
        """
        P1 = self.A
        P2 = self.B
        Q1 = line2.A
        Q2 = line2.B
        dx1 = P1[0] - P2[0]
        dx2 = Q1[0] - Q2[0]
        dy1 = P1[1] - P2[1]
        dy2 = Q1[1] - Q2[1]
        D = det2x2_((dx1, dy1), (dx2, dy2))
        return abs(D) < sepsilon_
    
    def includes(self, M, strict = False, checkCollinear = True):
        """Check whether a point belongs to the line.
        
        :param M: the point for which we want to test whether it belongs to the line
        :param strict: Boolean, whether to take into account `extendA` and `extendB`
        :param checkCollinear: Boolean, whether to check the collinearity of `A`, `B`, and `M`; set to `False` only if you use `strict=True` and you are sure that `M` is on the line (AB)
        :returns: A Boolean value.        
        
        """
        _ = error_if_not_point_(M=M)
        M = np.asarray(M, dtype=float)
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
        _ = error_if_not_point_(M=M)
        A = self.A
        B = self.B
        A_B = B - A
        v = np.array([-A_B[1], A_B[0]])
        if self.includes(M):
            print("M is on the line.")
            return Line(M, M + v, True, True)
        H = line_line_intersection_(A, B, M - v, M + v)
        return Line(H, M, extendH, extendM)

    def projection(self, M):
        """Orthogonal projection of a point to the reference line.
        
        :param M: a point
        :returns: A point on the reference line.
        
        """
        _ = error_if_not_point_(M=M)
        A = self.A
        B = self.B
        x, y = B - A
        v = np.array([-y, x])
        return line_line_intersection_(A, B, M, M+v)

    def distance(self, M):
        """Distance from a point to the reference line.
        
        :param M: a point
        :returns: A number.
        
        """
        P = self.projection(M)
        return distance_(M, P)
    
    def reflection(self, M):
        """Reflection of a point with respect to the reference line.
        
        :param M: a point
        :returns: A point.
        
        """
        R = Reflection(self)
        return R.reflect(M)
    
    def rotate(self, alpha, O, degrees=True):
        """Rotate the reference line.
        
        :param alpha: angle of rotation
        :param O: center of rotation
        :param degrees: Boolean, whether `alpha` is given in degrees
        :returns: A `Line` object.
        
        """
        _ = error_if_not_point(O=O)
        O = np.asarray(O, dtype=float)
        if degrees:
            alpha *= pi/180
        cosalpha, sinalpha = unit_vector_(alpha)
        x, y = self.A - O
        RAt <- np.array([
            cosalpha*x - sinalpha*y, sinalpha*x + cosalpha*y
        ])
        x, y = self.B - O
        RBt <- np.array([
            cosalpha*x - sinalpha*y, sinalpha*x + cosalpha*y
        ])
        return Line(RAt + O, RBt + O, self.extendA, self.extendB)

    def translate(self, v):
        """Translate the reference line.
        
        :param v: the vector of translation
        :returns: A `Line` object.
        
        """
        _ = error_if_not_point(v=v)
        v = np.asarray(v, dtype=float)
        return Line(self.A + v, self.B + v, self.extendA, self.extendB)
        
    def invert(self, iota):
        """Invert the reference line.
        
        :param iota: an `Inversion` object
        :returns: A `Line` object or a `Circle` object.
        
        """
        return iota.invert_line(self)

    def intersection_with_circle(self, circ):
        """Intersection(s) of the line with a circle.
        
        :param circ: a `Circle` object
        :returns: `None`, a point, or a list of two points.
        
        """
        return intersection_circle_line(circ, self)

    def intersection_with_ellipse(self, ell):
        """Intersection(s) of the line with an ellipse.
        
        :param ell: an `Ellipse` object
        :returns: `None`, a point, or a list of two points.
        
        """
        return intersection_ellipse_line(ell, self)

    def intersection_with_line(self, line2, strict=False):
        """Intersection(s) of the reference line with another line.
        
        :param line2: a `Line` object
        :returns: `None` (the lines are parallel), a point, or a `Line` object (the two lines are equal).
        
        """
        return intersection_line_line(self, line2, strict=strict)



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
    """A class for circles. A circle is given by its center and its radius.
    
    """
    def __init__(self, center, radius):
        _ = error_if_not_point_(center=center)
        _ = error_if_not_number_(radius=radius)
        _ = error_if_not_positive_(radius=radius)
        self.center = np.asarray(center, dtype=float)
        self.radius = radius
        
    def __str__(self):
        return str(self.__dict__)

    def show(self):
        print("Circle:\n")
        print(" center: ", tuple(self.center), "\n")
        print(" radius: ", self.radius, "\n")
        
    def is_equal(self, circ2):
        """Check whether the reference circle is equal to another circle.
        
        :param circ2: a `Circle` object
        :returns: A Boolean value.
        
        """
        if not isinstance(circ2, Circle):
            raise ValueError("`circ2` is not a `Circle` object.")
        return isclose(self.radius, circ2.radius) and np.allclose(self.center, circ2.center)
        
    def includes(self, P):
        """Check whether a point belongs to the reference circle.
        
        :param P: a point
        :returns: A Boolean value.
        
        """
        _ = error_if_not_point_(P=P)
        P = np.asarray(P, dtype=float)
        r = self.radius
        O = self.center
        d2 = dot_(P-O)
        return isclose(d2, r*r)

    def contains(self, P):
        """Check whether a point is contained in the reference circle.
        
        :param P: a point
        :returns: A Boolean value.
        
        """
        _ = error_if_not_point_(P=P)
        P = np.asarray(P, dtype=float)
        r = self.radius
        O = self.center
        d2 = dot_(P-O)
        return d2 <= r*r
    
    def random_points(self, n_points, where="in"):
        """Random points in the circle or on the circle.
        
        :param n_points: desired number of points
        :param where: either `"in"` or `"on"`
        :returns: A matrix with `n_points` rows and two columns; each row is 
        a random point inside the circle if `where="in"` or on the boundary of 
        the circle if `where="on"`.
        
        """
        if not isinstance(n_points, int):
            raise ValueError("`n_points` must be an integer.")
        _ = error_if_not_positive_(n_points=n_points)
        center = self.center
        radius = self.radius
        if where == "in":
            return runif_in_circle_(n_points, radius) + center
        if where == "on":
            return runif_on_circle_(n_points, radius) + center
        raise ValueError("The `where` argument must be `'in'` or `'on'`.")

    
    def point_from_angle(self, alpha, degrees=True):
        """Get a point on the reference circle from its polar angle.
        
        :param alpha: a number, the angle
        :param degrees: Boolean, whether the angle is given in degrees
        :returns: The point on the circle with polar angle `alpha`.
        
        """
        _ = error_if_not_number_(alpha=alpha)
        if degrees:
            alpha *= pi/180
        return self.center + self.radius * unit_vector_(alpha)
        
    def orthogonalThroughTwoPointsWithinCircle(self, P1, P2, arc = False):
        """Orthogonal circle passing through two points within the reference circle.
        
        :param P1,P2: two distinct points in the interior of the reference circle
        :param arc: Boolean, whether to return the arc joining the two points instead of the circle
        :returns: A `Circle` object or an `Arc` object, or a `Line` object if the two points are on a diameter.

        """
        _ = error_if_not_point_(P1=P1)
        _ = error_if_not_point_(P2=P2)
        P1 = np.asarray(P1, dtype=float)
        P2 = np.asarray(P2, dtype=float)
        if np.allclose(P1, P2):
            raise ValueError("`P1` and `P2` must be distinct.")
        I = self.center
        r = self.radius
        r2 = r * r
        if distance_(P1, I) >= r2:
            raise ValueError("`P1` is not in the interior of the reference circle.")
        if distance_(P2,I) >= r2:
            raise ValueError("`P2` is not in the interior of the reference circle.")
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
        """Orthogonal circle passing through two points on the reference circle.
        
        :param alpha1,alpha2: two angles defining two points on the reference circle
        :param arc: Boolean, whether to return only the arc at the interior of the reference circle
        :returns: A `Circle` object if `arc=False`, an `Arc` object if `arc=True, or a `Line` object: the diameter of the reference circle defined by the two points in case when the two angles differ by `pi`.

        """
        _ = error_if_not_number_(alpha1=alpha1)
        _ = error_if_not_number_(alpha2=alpha2)
        _ = error_if_not_boolean_(arc=arc)
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
    
    def power(self, M):
        """Power of a point with respect to the reference circle.
        
        :param M: a point (array-like of length two)
        :returns: A number, the power of `M` with respect to the circle.
        
        """
        _ = error_if_not_point_(M=M)
        M = np.asarray(M, dtype=float)
        radius = self.radius
        return dot_(M - self.center) - radius*radius
    
    def radical_center(self, circ2):
        """Radical center of two circles.
        
        :param circ2: a `Circle` object
        :returns: the radical center of the reference circle and `circ2`.
        
        """
        if not isinstance(circ2, Circle):
            raise ValueError("`circ2` is not a `Circle` object.")
        C1 = self.center
        C2 = circ2.center
        r1 = self.radius
        r2 = circ2.radius
        k = r1*r1 - r2*r2
        C1_C2 = C2 - C1
        C1C2sqr = dot_(C1_C2)
        if C1C2sqr == 0:
            return np.array([np.inf, np.inf])
        return (C1 + C2)/2 + k/2 * C1_C2/C1C2sqr
    
    def radical_axis(self, circ2):
        """Radical axis of two circles.
        
        :param circ2: a `Circle` object
        :returns: A `Line` object, the radical axis of the reference circle and `circ2`.
        
        """
        if not isinstance(circ2, Circle):
            raise ValueError("`circ2` is not a `Circle` object.")
        C1 = self.center
        C2 = circ2.center
        if np.allclose(C1, C2):
            print("The two circles must have distinct centers.")
            return
        x, y = C2 - C1
        v = np.array([-y, x])
        R = self.radical_center(circ2)
        return Line(R, R+v, True, True)
    
    def as_ellipse(self):
        """Converts the circle to an `Ellipse` object.
 
        
        """
        r = self.radius
        return Ellipse(self.center, r, r, 0)

    def intersection_with_line(self, line):
        """
        
        """
        return intersection_circle_line(self, line)


def radical_center(circ1, circ2, circ3):
    """Radical center of three circles.
    
    :param circ1,circ2,circ3: `Circle` objects
    :returns: A point, the radical center of the three circles.
    
    """
    if not isinstance(circ1, Circle):
        raise ValueError("`circ1` is not a `Circle` object.")
    if not isinstance(circ2, Circle):
        raise ValueError("`circ2` is not a `Circle` object.")
    if not isinstance(circ3, Circle):
        raise ValueError("`circ3` is not a `Circle` object.")
    l1 = circ1.radical_axis(circ2)
    l2 = circ1.radical_axis(circ3)
    return line_line_intersection_(l1.A, l1.B, l2.A, l2.B)
   

def intersection_circle_line(circ, line):
    """Intersection(s) of a circle and a line.

    """
    if not isinstance(circ, Circle):
        raise ValueError("`circ` is not a `Circle` object.")
    if not isinstance(line, Line):
        raise ValueError("`line` is not a `Line` object.")
    C = circ.center
    intersections = circle_line_intersection_(line.A - C, line.B - C, circ.radius)
    if intersections is None:
        return None
    if isinstance(intersections, list):
        I1, I2 = intersections
        return [I1 + C, I2 + C]
    return intersections + C
        

def intersection_ellipse_line(ell, line):
    """Intersection(s) of an ellipse and a line.

    """
    if not isinstance(ell, Ellipse) and not isinstance(ell, Circle):
        raise ValueError("`ell` is not an `Ellipse` object neither a `Circle` object.")
    if not isinstance(line, Line):
        raise ValueError("`line` is not a `Line` object.")
    if isinstance(ell, Circle):
        return intersection_circle_line(ell, line)
    a = ell.rmajor
    b = ell.rminor
    theta = ell.alpha
    if ell.degrees:
        theta *= pi/180
    v = unit_vector_(theta)
    w = np.array([-v[1], v[0]])
    # f maps the unit circle to ell:
    f = Affine(np.column_stack((a*v, b*w)), ell.center)
    invf = f.inverse() # maps ell to the unit circle
    line2 = invf.transform_line(line)
    Is = intersection_circle_line(Circle((0, 0), 1), line2)
    if Is is None:
        return None
    if isinstance(Is, list):
        I1, I2 = Is
        return [f.transform(I1), f.transform(I2)]
    return f.transform(Is)


class Arc:
    """Arc class.
    
    A circular arc is initialized by its center (array-like object of length 
    two), a radius, a starting angle and an ending angle. They are 
    respectively named `center`, `radius`, `alpha1` and `alpha2`.
    
    """
    def __init__(self, center, radius, alpha1, alpha2, degrees = True):
        _ = error_if_not_point_(center=center)
        _ = error_if_not_number_(radius=radius)
        _ = error_if_not_positive_(radius=radius)
        _ = error_if_not_number_(alpha1=alpha1)
        _ = error_if_not_number_(alpha2=alpha2)
        _ = error_if_not_boolean_(degrees=degrees)
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
        if not isinstance(n_points, int):
            raise ValueError("`n_points` must be an integer.")
        _ = error_if_not_positive_(n_points=n_points)
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
        _ = error_if_not_point_(center=center)
        _ = error_if_not_number_(rmajor=rmajor)
        _ = error_if_not_positive_(rmajor=rmajor)
        _ = error_if_not_number_(rminor=rminor)
        _ = error_if_not_positive_(rminor=rminor)
        _ = error_if_not_number_(alpha=alpha)
        _ = error_if_not_boolean_(degrees=degrees)
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
    
    def is_equal(self, ell2):
        """Check whether the reference ellipse equals another ellipse.
        
        
        """
        if isinstance(ell2, Circle):
            ell2 = ell2.as_ellipse()
        x1, y1 = self.center
        alpha1 = self.alpha
        if not self.degrees:
            alpha1 *= 180/pi
        x2, y2 = ell2.center
        alpha2 = ell2.alpha
        if not ell2.degrees:
            alpha2 *= 180/pi
        v1 = np.array([x1, y1, self.rmajor, self.rminor, alpha1 % 180])
        v2 = np.array([x2, y2, ell2.rmajor, ell2.rminor, alpha2 % 180])
        return np.allclose(v1, v2)

    def path(self, n_points=100):
        """Path that forms the ellipse.
        
        :param n_points: number of points of the path
        :returns: A matrix with two columns and `n_points` rows.
        
        """
        if not isinstance(n_points, int):
            raise ValueError("`n_points` must be an integer.")
        _ = error_if_not_positive_(n_points=n_points)
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

    def random_points(self, n_points, where="in"):
        """Random points in/on the ellipse.
        
        :param n_points: desired number of points
        :param where: either `"in"` or `"on"`
        :returns: A matrix with `n_points` rows and two columns; each row is 
        a random point inside the ellipse if `where="in"` or on the boundary
        of the ellipse if `where="on"`.
        
        """
        if not isinstance(n_points, int):
            raise ValueError("`n_points` must be an integer.")
        _ = error_if_not_positive_(n_points=n_points)
        center = self.center
        A = self.shape_matrix()
        if where == "in":
            return runif_in_ellipse_(n_points, A) + center
        if where == "on":
            return runif_on_ellipse_(n_points, A) + center
        raise ValueError("The `where` argument must be `'in'` or `'on'`.")
    
    def point_from_angle(self, theta, degrees=True):
        """
        
        """
        _ = error_if_not_boolean_(degrees=degrees)
        if not is_number_(theta):
            if not is_vector_(theta):
                raise ValueError("`theta` must be a number or a vector of numbers.")
            theta = np.asarray(theta, dtype=float)
            if not is_real_vector_(theta):
                raise ValueError("`theta` cannot be converted to a vector of real numbers.")
        else:
            theta = np.array([theta], dtype=float)
        O = self.center
        a = self.rmajor
        b = self.rminor
        alpha = self.alpha
        if self.degrees:
            alpha *= pi/180
        if degrees:
            theta *= pi/180
        theta_mod_2pi = theta % (2*pi)
        sgn = 2 * (theta_mod_2pi <= sepsilon_) - 1
        theta_eps = theta + sepsilon_ * sgn
        t = np.arctan2(a/b, 1/np.tan(theta_mod_2pi)) + theta_eps - theta_eps % pi
        out = ellipse_points_(t, O, a, b, alpha)
        if len(theta) == 1:
            out = out[0]
        return out

    def normal(self, t):
        """Normal unit vector to the ellipse.
    
        :param t: a number, the eccentric angle in radians of the point of the
        ellipse at which we want the normal unit vector
        :returns: The normal unit vector to the ellipse at the point given by
        eccentric angle `t`.

        """        
        _ = error_if_not_number_(t=t)
        O = self.center
        a = self.rmajor
        b = self.rminor
        alpha = self.alpha
        if self.degrees:
            alpha *= pi/180
        cosalpha = cos(alpha)
        sinalpha = sin(alpha)
        x = -a * sin(t)
        y = b * cos(t)
        v = np.array([
            sinalpha*x + cosalpha*y,
            -cosalpha*x + sinalpha*y                
        ])
        return v / vlength_(v)
        
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
    
    def includes(self, P):
        """Check whether a point belongs to the ellipse.
        
        :param P: a point
        :returns: A Boolean value.
        
        """
        _ = error_if_not_point_(P=P)
        A, B, C, D, E, F = self.equation().values()
        x, y = P
        zero = A*x*x + B*x*y + C*y*y + D*x + E*y + F
        return isclose(1.0, zero+1.0)

    def contains(self, P):
        """Check whether a point is contained in the ellipse.
        
        :param P: a point
        :returns: A Boolean value.
        
        """
        _ = error_if_not_point_(P=P)
        A, B, C, D, E, F = self.equation().values()
        x, y = P
        return A*x*x + B*x*y + C*y*y + D*x + E*y + F <= 0
    
    def shape_matrix(self):
        """The 2x2 symmetric matrix `S` associated to the reference ellipse. 
        The equation of the ellipse is `(M-O)' S (M-O) = 1`.
        
        """
        A, B, C, D, E, F = self.equation().values()
        X = np.array([
            [A, B/2, D/2],
            [B/2, C, E/2],
            [D/2, E/2, F]            
        ])
        K = - np.linalg.det(X) / (A*C - B*B/4)
        return X[0:2, 0:2] / K
    
    def theta2t(self, theta, degrees=True):
        """Convert angle to eccentric angle.
        
        :param theta: angle between the major axis and the half-line starting
        at the center of the ellipse and passing through the point of interest
        on the ellipse
        :param degrees: Boolean, whether `theta` is given in degrees
        :returns: The eccentric angle of the point of interest on the ellipse,
        in radians.
        
        """
        _ = error_if_not_number_(theta=theta)
        _ = error_if_not_boolean_(degrees=degrees)
        if degrees:
            theta *= pi/180
        theta_mod_2pi = theta % (2*pi)
        sgn = 1 if theta_mod_2pi < sepsilon_ else -1
        theta_eps = theta + sgn*sepsilon_
        return atan2(self.rmajor/self.rminor, 1/tan(theta_mod_2pi)) + theta_eps - theta_eps % pi

    def intersection_with_line(self, line):
        """
        
        """
        return intersection_ellipse_line(self, line)
    
    @classmethod
    def from_equation(cls, A, B, C, D, E, F):
        """Ellipse from its implicit equation.
        
        :param A,B,C,D,E,F: coefficients of the implicit equation of the ellipse
        :returns: An `Ellipse` object.
        
        """
        if A*C <= 0 or B*B-4*A*C >= 0 or D*D + E*E <= 4*(A+C)*F:
            raise ValueError("These parameters do not define an ellipse.")
        Q = np.array([[2*A, B, D], [B, 2*C, E], [D, E, 2*F]])
        if np.linalg.det(Q) == 0:
            raise ValueError("These parameters do not define an ellipse.")
        M0 = np.array([F, D/2, E/2, D/2, A, B/2, E/2, B/2, C]).reshape((3,3))
        M = np.array([[A, B/2], [B/2, C]])
        eigvals = np.linalg.eigvals(M)
        detM0 = np.linalg.det(M0)
        detM = np.linalg.det(M)
        a = sqrt(-detM0 / (detM * eigvals[0]))
        b = sqrt(-detM0 / (detM * eigvals[1]))
        xy = np.array([B*E - 2*C*D, B*D - 2*A*E])
        phi = 0.0 if A == C else (atan(B/(A-C))/2 if abs(C) > abs(A) else pi/2 - atan(-B/(A-C))/2)
        return Ellipse(xy/(4*A*C - B*B), max(a,b), min(a,b), (phi*180/pi) % 180)
    
    @classmethod
    def equation_from_five_points(cls, P1, P2, P3, P4, P5):
        """The implicit equation of the ellipse is
        `Ax² + Bxy + Cy² + Dx + Ey + F = 0`. This function returns
        A, B, C, D, E and F.
        
        :param P1,P2,P3,P4,P5: five points
        :returns: A dictionary giving A, B, C, D, E and F.
        
        """
        _ = error_if_not_point_(P1=P1)
        _ = error_if_not_point_(P2=P2)
        _ = error_if_not_point_(P3=P3)
        _ = error_if_not_point_(P4=P4)
        _ = error_if_not_point_(P5=P5)
        P1 = np.asarray(P1, dtype=float)
        P2 = np.asarray(P2, dtype=float)
        P3 = np.asarray(P3, dtype=float)
        P4 = np.asarray(P4, dtype=float)
        P5 = np.asarray(P5, dtype=float)
        P = np.array([P1, P2, P3, P4, P5])
        if np.unique(P, axis=0).shape != P.shape:
            raise ValueError("The five points are not distinct.")
        x = P[:, 0]
        y = P[:, 1]
        M = np.column_stack((x*x, x*y, y*y, x, y, np.ones(5)))
        A = np.linalg.det(M[:, [1, 2, 3, 4, 5]])
        B = -np.linalg.det(M[:, [0, 2, 3, 4, 5]])
        C = np.linalg.det(M[:, [0, 1, 3, 4, 5]])
        if B*B-4*A*C >= 0:
            raise ValueError("The five points do not lie on an ellipse.")
        D = -np.linalg.det(M[:, [0, 1, 2, 4, 5]])
        E = np.linalg.det(M[:, [0, 1, 2, 3, 5]])
        F = -np.linalg.det(M[:, [0, 1, 2, 3, 4]])
        return {
            "A": A,
            "B": B,
            "C": C,
            "D": D,
            "E": E,
            "F": F
        }
    
    @classmethod
    def from_five_points(cls, P1, P2, P3, P4, P5):
        """Ellipse from five points on this ellipse.
        
        :param P1,P2,P3,P4,P5: five points
        :returns: An `Ellipse` object.
        
        """
        coeffs = Ellipse.equation_from_five_points(P1, P2, P3, P4, P5)
        A, B, C, D, E, F = coeffs.values()
        return Ellipse.from_equation(A, B, C, D, E, F)
    
    @classmethod
    def from_center_and_matrix(cls, center, S):
        """
        
        """
        _ = error_if_not_point_(center=center)
        center = np.asarray(center, dtype=float)
        e = np.linalg.eigh(S)
        params = ellipse_from_center_and_eigen_(center, e)
        return Ellipse(
            params["center"], params["rmajor"], params["rminor"], params["alpha"]
        )
        
        



# def circle_as_ellipse_(C):
#     r = C.radius
#     return Ellipse(C.center, r, r, 0)


class Inversion:
    """Inversion class.
    
    An inversion is initialized by its pole and its power.
    
    """
    def __init__(self, pole, power):
        _ = error_if_not_point_(pole=pole)
        _ = error_if_not_number_(power=power)
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
        if is_inf(M):
            return pole
        _ = error_if_not_point_(M=M)
        M = np.asarray(M)
        if np.allclose(pole, M):
            return inf
        k = self.power
        pole_M = M - pole
        return pole + k/dot_(pole_M) * pole_M
    
    def invert_circle(self, circ):
        if not isinstance(circ, Circle):
            raise ValueError("`circ` must be a `Circle` object.")
        c0 = self.pole
        k = self.power
        c1 = circ.center
        r1 = circ.radius
        D1 = (c1[0] - c0[0])**2 + (c1[1] - c0[1])**2 - r1*r1
        if abs(D1) > sepsilon_:
            s = k / D1
            return Circle(c0 + s*(c1-c0), abs(s)*r1)
        Ot = c0 - c1
        R180 = c1 - Ot
        R90 = np.array([-Ot[1], Ot[0]]) + c1
        return Line(self.invert(R180), self.invert(R90))

    def invert_line(self, line):
        if not isinstance(line, Line):
            raise ValueError("`line` must be a `Line` object.")
        c0 = self.pole
        A = line.A
        B = line.B
        if collinear_(A, B, c0):
            return line
        Ap = self.invert(A)
        Bp = self.invert(B)
        return Triangle(c0, Ap, Bp).circumcircle()
    
    @classmethod
    def on_circle(cls, circ):
        """An inversion on a circle is the inversion whose pole is the center 
        of the circle and whose power is the squared radius of the circle.
        
        :param circ: `Circle` object
        :returns: An `Inversion` object.
        
        """
        return Inversion(circ.center, circ.radius**2)
    
    @classmethod
    def from_swapping_two_circles(cls, circ1, circ2, positive=True):
        """Inversion swapping two circles.
        
        :param circ1,circ2: `Circle` objects
        :param positive: Boolean, whether the sign of the desired inversion power must be positive or negative
        :returns: An `Inversion` object, which maps `circ1` to `circ2` and `circ2` to `circ1`, except in the case when `circ1` and `circ2` are congruent and tangent: in this case a `Reflection` object is returned (a reflection is an inversion on a line).
        
        """
        if not isinstance(circ1, Circle):
            raise ValueError("`circ1` must be a `Circle` object.")
        if not isinstance(circ2, Circle):
            raise ValueError("`circ2` must be a `Circle` object.")
        _ = error_if_not_boolean_(positive=positive)
        c1 = circ1.center
        r1 = circ1.radius
        c2 = circ2.center
        r2 = circ2.radius
        c1c2 = distance_(c1, c2)
        if r1 == r2:
            I = (c1 + c2) / 2
            if c1c2 < r1+r2: # they intersect at two points or are equal
                if not positive:
                    print("`positive = False` not possible; switching to `True`")
                k = abs(dot_(I - c2) - r2*r2)
                return Inversion(I, k)
            if c1c2 > r1+r2: # they do not intersect
                if positive:
                    print("`positive = True` not possible; switching to `False`")
                k = -abs(dot_(I - c2) - r2*r2)
                return Inversion(I, k)
            # they touch
            print("No possible inversion; returning a reflection")
            x, y = c2 - c1
            v = np.array([-y, x])
            print("TO DO")
            return#Reflection(Line(I+v, I-v))
        inside = inside_tangent = False
        if max(r1, r2) >= c1c2 + min(r1, r2):
            inside = True
            inside_tangent = max(r1, r2) == c1c2 + min(r1, r2)
        if (not positive) and ((not inside) or inside_tangent):
            they_intersect = c1c2 <= r1+r2
            if they_intersect:
                print("`positive = False` not possible; switching to `True`")
                positive = True
        a = r1 / r2
        if positive:
            if inside:
                O = c1 + a/(1+a) * (c2 - c1)
            else:
                O = c1 - a/(1-a) * (c2 - c1)
            return Inversion(O, a * abs(dot_(O - c2) - r2*r2))
        if inside:
            O = c1 - a/(1-a) * (c2 - c1)
        else:
            O = c1 + a/(1+a) * (c2 - c1)
        return Inversion(O, -a * abs(dot_(O - c2) - r2*r2))

    @classmethod
    def from_fixing_two_circles(cls, circ1, circ2):
        """Inversion fixing two circles.
        
        :param circ1,circ2: `Circle` objects
        :returns: An `Inversion` object representing an inversion which leaves each of the two circles invariant.
        
        """
        if not isinstance(circ1, Circle):
            raise ValueError("`circ1` must be a `Circle` object.")
        if not isinstance(circ2, Circle):
            raise ValueError("`circ2` must be a `Circle` object.")
        O = circ1.radical_center(circ2)
        return Inversion(O, circ1.power(O))

    @classmethod
    def from_fixing_three_circles(cls, circ1, circ2, circ3):
        """Inversion fixing three circles.
        
        :param circ1,circ2,circ3: `Circle` objects
        :returns: An `Inversion` object representing an inversion which leaves each of the three circles invariant.
        
        """
        if not isinstance(circ1, Circle):
            raise ValueError("`circ1` must be a `Circle` object.")
        if not isinstance(circ2, Circle):
            raise ValueError("`circ2` must be a `Circle` object.")
        if not isinstance(circ3, Circle):
            raise ValueError("`circ3` must be a `Circle` object.")
        Rc = radical_center(circ1, circ2, circ3)
        return Inversion(Rc, circ1.power(Rc))



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
    if not isinstance(c0, Circle):
        raise ValueError("`c0` must be a `Circle` object.")
    if not isinstance(n, int):
        raise ValueError("`n` must be an integer.")
    _ = error_if_not_positive_(n=n)
    _ = error_if_not_number(phi=phi)
    _ = error_if_not_number(shift=shift)
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
        _ = error_if_not_point_(A=A)
        _ = error_if_not_point_(B=B)
        _ = error_if_not_point_(C=C)
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
        """Flatness, a number between 0 and 1; a triangle is flat when its flatness is 1.
        
        """
        AB = self.B - self.A
        AC = self.C - self.A
        z = complex(*AB).conjugate() * complex(*AC)
        re = z.real
        im = z.imag
        re2 = re * re
        return re2 / (re2 + im*im)

    @property
    def a(self):
        """Length of the side BC.
        
        """
        return distance_(self.B, self.C)

    @property
    def b(self):
        """Length of the side AC.
        
        """
        return distance_(self.A, self.C)

    @property
    def c(self):
        """Length of the side AB.
        
        """
        return distance_(self.A, self.B)

    @property    
    def edges(self):
        """Edge lengths of the triangle.
        
        """
        return {"a": self.a, "b": self.b, "c": self.c}

    @property    
    def orientation(self):
        """Orientation of the triangle; 1 for counterclockwise, -1 for clockwise, 0 for collinear.
        
        """
        A = self.A
        B = self.B
        C = self.C
        val = (B[1] - A[1])*(C[0] - B[0]) - (B[0] - A[0])*(C[1] - B[1])
        return 0 if val==0 else (1 if val>0 else -1)
    
    def contains(self, M):
        """Check whether a point lies inside the reference triangle.
        
        """
        _ = error_if_not_point_(M=M)
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
        """Check whether the triangle is acute.
        
        """
        edges = [self.a, self.b, self.c]
        edges.sort()
        edge0, edge1, edge2 = edges
        return edge0*edge0 + edge1*edge1 >= edge2*edge2

    @property    
    def angleA(self):
        """The angle at the vertex A in radians.
        
        """
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
        """The angle at the vertex B in radians.
        
        """
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
        """The angle at the vertex C in radians.
        
        """
        A = self.A
        B = self.B
        C = self.C
        CA = A - C
        CB = B - C
        b = vlength_(CA)
        a = vlength_(CB)
        return acos(dot_(CA, CB) / a / b)
    
    def incircle(self):
        """The incircle of the triangle.
        
        """
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
    
    def circumcircle(self):
        if self.flatness == 1:
            print("The triangle is flat.")
            return None
        A = self.A
        B = self.B
        C = self.C
        q = np.array([dot_(A), dot_(B), dot_(C)])
        ABC = np.array([A, B, C])
        ones = np.array([[1.0], [1.0], [1.0]])
        Dx = np.linalg.det(
            np.hstack(
                (np.column_stack((q, ABC[:, 1])), ones)
            )
        )
        Dy = - np.linalg.det(
            np.hstack(
                (np.column_stack((q, ABC[:, 0])), ones)
            )
        )
        center = np.array([Dx, Dy]) / np.linalg.det(np.hstack((ABC, ones))) / 2
        return Circle(center, distance_(center, A))
    
    def excircles(self):
        """The excircles of the triangle.
        
        :returns: A dictionary of three `Circle` objects.
        
        """
        A = self.A
        B = self.B
        C = self.C
        a = distance_(B, C)
        b = distance_(A, C)
        c = distance_(B, A)
        s = (a + b + c) / 2
        JA = (-a*A + b*B + c*C) / (-a + b + c)
        JB = (a*A - b*B + c*C) / (a - b + c)
        JC = (a*A + b*B - c*C) / (a + b - c)
        rA = sqrt(s*(s-b)*(s-c)/(s-a))
        rB = sqrt(s*(s-a)*(s-c)/(s-b))
        rC = sqrt(s*(s-a)*(s-b)/(s-c))
        return {
          "A": Circle(center = JA, radius = rA),
          "B": Circle(center = JB, radius = rB),
          "C": Circle(center = JC, radius = rC)
        }
        
    def orthic_triangle(self):
        """Orthic triangle. Its vertices are the feet of the altitudes
        of the reference triangle.
        
        :returns: A `Triangle` object.

        """
        A = self.A
        B = self.B
        C = self.C
        a = distance_(B, C)
        b = distance_(A, C)
        c = distance_(B, A)
        AC = C-A
        AB = B-A
        BC = C-B
        BA = A-B
        CA = A-C
        CB = B-C
        x = 1 / (dot_(AC,AB) / b / c)
        y = 1 / (dot_(BC,BA) / a / c)
        z = 1 / (dot_(CA,CB) / a / b)
        HA = (b/z*B + c/y*C) / (b/z + c/y)
        HB = (a/z*A + c/x*C) / (a/z + c/x)
        HC = (a/y*A + b/x*B) / (a/y + b/x)
        return Triangle(HA, HB, HC)

    
    def malfatti_circles(self):
        """The Malfatti circles of the triangle.
        
        :returns: Three circles and three tangency points.
        
        """
        if self.flatness == 1:
            print("The triangle is flat.")
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
        """Equal detour point of the triangle, also known as the X(176) triangle center.
        
        """
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
    
    def steiner_ellipse(self):
        """The Steiner ellipse (or circumellipse) of the reference triangle. 
        This is the ellipse passing through the three vertices of the triangle 
        and centered at the centroid of the triangle.
        
        :returns: An `Ellipse` object.   
        
        """
        if self.flatness == 1:
            print("The triangle is flat.")
            return None
        P = (0, 0)
        Q = (10, 0)
        R = (5, 5*sqrt(3))
        circ = Triangle(P, Q, R).circumcircle()
        f = Affine.from_mapping_three_points(P, Q, R, self.A, self.B, self.C)
        return f.transform_ellipse(circ)

    def steiner_inellipse(self):
        """The Steiner inellipse (or midpoint ellipse) of the reference triangle. 
        This is the ellipse tangent to the sides of the triangle at their 
        midpoints, and centered at the centroid of the triangle.
        
        :returns: An `Ellipse` object.   
        
        """
        if self.flatness == 1:
            print("The triangle is flat.")
            return None
        P = (0, 0)
        Q = (10, 0)
        R = (5, 5*sqrt(3))
        circ = Triangle(P, Q, R).incircle()
        f = Affine.from_mapping_three_points(P, Q, R, self.A, self.B, self.C)
        return f.transform_ellipse(circ)


class Affine:
    """A class for affine transformations.
    
    An affine transformation is initialized by a 2x2 matrix (a linear transformation), 
    and a length two vector (the 'intercept', an array-like object).
    
    """
    def __init__(self, A, b):
        _ = error_if_not_point_(b=b)
        self.A = np.asarray(A, dtype=float)
        if self.A.ndim != 2 or self.A.shape != (2,2):
            raise ValueError("`A` cannot be converted to a 2x2 matrix.")
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
            return None
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
        m = np.asarray(m, dtype=float)
        if m.ndim == 2:
            b = np.repeat(b.reshape(2,1), m.shape[0], axis=1)
        return np.transpose(np.matmul(self.A, np.transpose(m)) + b)
    
    def transform_line(self, line):
        """Transform a line by the affine map.
        
        
        """
        M = self.A
        b = self.b
        if np.linalg.det(M) == 0:
            print("The affine map is singular.")
            return
        return Line(M.dot(line.A) + b, M.dot(line.B) + b, line.extendA, line.extendB)
    
    def transform_ellipse(self, ell):
        """Transform an ellipse by the reference affine transformation (only for an invertible affine map).
        
        :param ell: an `Ellipse` object or a `Circle` object
        :returns: An `Ellipse` object.
        
        """
        if np.linalg.det(self.A) == 0:
            print("The affine map is singular.")
            return
        if isinstance(ell, Circle):
            ell = Circle.as_ellipse(ell)
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
    def from_mapping_three_points(cls, P1, P2, P3, Q1, Q2, Q3):
        """Affine transformation mapping three given points to three given points.
        
        :param P1,P2,P3: three non-collinear points
        :param Q1,Q2,Q3: three non-collinear points
        :returns: An `Affine` object representing the transformation which maps Pi to Qi for each i=1,2,3.
        
        """
        _ = error_if_not_point_(P1=P1)
        _ = error_if_not_point_(P2=P2)
        _ = error_if_not_point_(P3=P3)
        _ = error_if_not_point_(Q1=Q1)
        _ = error_if_not_point_(Q2=Q2)
        _ = error_if_not_point_(Q3=Q3)
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
        f1 = Affine(np.column_stack((P2-P1, P3-P1)), P1)
        f2 = Affine(np.column_stack((Q2-Q1, Q3-Q1)), Q1)
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


class Mobius:
    """A class for Möbius transformations.
    
    A Möbius transformation is initialized by a complex 2x2 matrix with a  
    non-zero determinant.
    
    """
    def __init__(self, M):
        self.M = np.asarray(M, dtype=complex)
        
    def __str__(self):
        return str(self.__dict__)

    def show(self):
        print("Möbius transformation:\n")
        print("                     M: ", self.M, "\n")
        
    @property
    def a(self):
        return self.M[0, 0]
    
    @property
    def b(self):
        return self.M[0, 1]

    @property
    def c(self):
        return self.M[1, 0]
    
    @property
    def d(self):
        return self.M[1, 1]

    def compose(self, M1, left=True):
        """Compose the reference Möbius transformation with another Möbius transformation.
        
        :param M1: a `Mobius` object
        :param left: Boolean, whether to compose at left or at right (i.e. returns `M1 o M0` or `M0 o M1`)
        :returns: A `Mobius` object.
        
        """
        A = self.M
        B = M1.M
        return Mobius(np.matmul(B, A)) if left else Mobius(np.matmul(A, B)) 

    def inverse(self):
        """Inverse of the Möbius transformation.
        
        """
        return Mobius(np.array([
                [self.d, -self.b],
                [-self.c, self.a]
            ]))

    def power(self, k):
        """Power of the Möbius transformation.

        :param k: an integer, possibly negative
        :returns: A `Mobius` object corresponding to the Möbius transformation raised to the power `k`.
        
        """
        t = int(k)
        if t != k:
            raise ValueError("The power `k` must be an integer.")
        M = self.M
        if t < 0:
            M = M.inverse()
            t = -t
        return Mobius(np.linalg.matrix_power(M, t))
    
    def gpower(self, t):
        """Generalized power of the Möbius transformation.

        :param t: a float, possibly negative
        :returns: A `Mobius` object corresponding to the Möbius transformation raised to the power `t`.
        
        """
        M = self.M
        detM = det2x2_mat_(M)
        trM = M[0, 0] + M[1, 1]
        if abs(trM*trM - 4*detM) < sepsilon_:
            alpha = trM / 2
            D = np.diag((alpha, alpha))
            if np.allclose(M, D):
                alpha_t = alpha**t
                return Mobius(np.diag(alpha_t, alpha_t))
            N = M - D
            if np.allclose(N[:, 1], np.zeros((2,))):
                v2 = np.array([1.0, 0.0])
                v1 = N[:, 0]
            else:
                v2 = np.array([0.0, 1.0])
                v1 = N[:, 1]
            P = np.column_stack((v1, v2))
            alpha_t = alpha**t
            Jt = np.array([
                    [alpha_t, t*alpha**(t-1)],
                    [0.0, alpha_t]
                ])
            return Mobius(np.matmul(np.matmul(P, Jt), np.linalg.inv(P)))
        eigen_values, eigen_vectors = np.linalg.eig(M)
        Dt = np.diag(eigen_values**t)
        return Mobius(
            np.matmul(
                np.matmul(eigen_vectors, Dt), np.linalg.inv(eigen_vectors)
            )
        )

    def transform(self, P):
        """Transform a point by the Möbius transformation.
        
        :param P: a point (array-like of length two) or `inf`
        :returns: The image of `P` by the Möbius transformation (can be `inf`).
        
        """
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        if is_inf(P):
            return inf if c == 0 else from_complex_(a / c)
        _ = error_if_not_point_(P=P)
        P = np.asarray(P)
        z = complex(*P)
        condition = c != 0 and z == -d/c
        return inf if condition else from_complex_((a*z+b)/(c*z+d))
    
    def transform_circle(self, circ):
        """Transform a circle by the Möbius transformation.
        
        :param circ: a `Circle` object
        :returns: A `Circle` object or a `Line` object.
        
        """
        if not isinstance(circ, Circle):
            raise ValueError("`circ` must be a `Circle` object.")
        c = self.c
        d = self.d
        R = circ.radius
        z0 = complex(*circ.center)
        x1 = mod2_(d + c*z0)
        x2 = R*R*mod2_(c)
        if x1 != x2: # we are in this case if c=0
            if x1 > 0:
                z = z0 if c == 0 else z0 - R*R / (d/c + z0).conjugate()
                w0 = self.transform(from_complex_(z))
            else:
                w0 = self.transform(inf)
            v = w0 - self.transform(from_complex_(z0 + R))
            return Circle(w0, abs(complex(*v)))
        M = from_complex_(-d/c)
        if c != 0 and circ.includes(M):
            alpha = 0.0
            while True:
                P = circ.point_from_angle(alpha, degrees=False)
                if not np.allclose(P, M):
                    break
                alpha += 0.1
            while True:
                Q = circ.point_from_angle(alpha+3, degrees=False)
                if not np.allclose(Q, M):
                    break
                alpha += 0.1
            A = self.transform(P)
            B = self.transform(Q)
        else:
            A = self.transform(circ.point_from_angle(0))
            B = self.transform(circ.point_from_angle(180, degrees=True))
        return Line(A, B)
    
    def transform_line(self, line):
        """Transform a line by the Möbius transformation.
        
        :param line: a `Line` object
        :returns: A `Circle` object or a `Line` object.
        
        """
        if not isinstance(line, Line):
            raise ValueError("`line` must be a `Line` object.")
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        do = line.direction_offset()
        theta = do["direction"]
        gamma0 = complex(*unit_vector_(theta))
        D0 = 2 * do["offset"]
        A = -2 * (gamma0 * c * d.conjugate()).real - D0 * mod2_(c)
        gamma = gamma0.conjugate() * b * c.conjugate() + gamma0 * d.conjugate() * a + D0 * c.conjugate() *a
        D = D0 * mod2_(a) + 2 * (gamma0 * b.conjugate() * a).real
        if abs(A) > sepsilon_:
            return Circle(from_complex_(gamma/A), sqrt(mod2_(gamma)/A/A + D/A))
        return Line(self.transform(line.A), self.transform(line.B))
    
    @classmethod
    def from_mapping_three_points(cls, P1, P2, P3, Q1, Q2, Q3):
        """Möbius transformation mapping three given points to three given points.
        
        :param P1,P2,P3: three distinct points, `inf` allowed
        :param Q1,Q2,Q3: three distinct points, `inf` allowed
        :returns: A `Mobius` object, representing the Möbius transformation which sends `Pi` to `Qi` for each i=1,2,3.
            
        """
        if is_inf(P1):
            P1 = [inf]
        else:
            _ = error_if_not_point_(P1=P1)
        if is_inf(P2):
            P2 = [inf]
        else:
            _ = error_if_not_point_(P2=P2)
        if is_inf(P3):
            P3 = [inf]
        else:
            _ = error_if_not_point_(P3=P3)
        if is_inf(Q1):
            Q1 = [inf]
        else:
            _ = error_if_not_point_(Q1=Q1)
        if is_inf(Q2):
            Q2 = [inf]
        else:
            _ = error_if_not_point_(Q2=Q2)
        if is_inf(Q3):
            Q3 = [inf]
        else:
            _ = error_if_not_point_(Q3=Q3)
        z1 = complex(*np.asarray(P1))
        z2 = complex(*np.asarray(P2))
        z3 = complex(*np.asarray(P3))
        if z1 == z2 or z1 == z3 or z2 == z3:
            print("`P1`, `P2` and `P3` must be distinct.")
            return
        Mob1 = MobiusMappingThreePoints2ZeroOneInf_(z1, z2, z3)
        w1 = complex(*np.asarray(Q1))
        w2 = complex(*np.asarray(Q2))
        w3 = complex(*np.asarray(Q3))
        if w1 == w2 or w1 == w3 or w2 == w3:
            print("`Q1`, `Q2` and `Q3` must be distinct.")
            return
        Mob2 = MobiusMappingThreePoints2ZeroOneInf_(w1, w2, w3)
        return Mob1.compose(Mob2.inverse())


def MobiusMappingThreePoints2ZeroOneInf_(z1, z2, z3):
    if isinf(z1.real):
        M = np.array([
            [0, z2 - z3],
            [1, z3]
        ])
        return Mobius(M)
    if isinf(z2.real):
        M = np.array([
            [1, -z1],
            [1, -z3]
        ])
        return Mobius(M)
    if isinf(z3.real):
        K = 1 / (z2 - z1)
        M = np.array([
            [K,  K * z1],
            [0, 1]
        ])
        return Mobius(M)
    K = (z2 - z3) / (z2 - z1)
    M = np.array([
        [K,  -K * z1],
        [1, -z3]
    ])
    return Mobius(M)

def unimodular_matrices(n):
    """Generates unimodular matrices.
    
    :param n: integer, the maximum size of entries of matrices, at least 1
    :returns: List of unimodular matrices.
    
    """
    if not isinstance(n, int):
        raise ValueError("`n` must be an integer.")
    if n < 1:
        raise ValueError("`n` must be at least one.")
    I2 = np.eye(2, dtype=int)
    if n == 1:
        return [I2]
    out = [farey_stack_(i) for i in range(1, n+1)]
    out = np.unique(np.array(reduce(lambda x, y: x + y, out)), axis = 0)
    return [I2] + [arr.reshape((2,2), order="F") for arr in out] 


class Reflection:
    """A class for reflections.
    
    A reflection is initialized by a line.
    
    """
    def __init__(self, line):
        if not isinstance(line, Line):
            raise ValueError("You must supply a `Line` object to initialize a reflection.")
        self.line = Line(line.A, line.B, True, True)
        
    def __str__(self):
        return str(self.__dict__)

    def show(self):
        print("Reflection with respect to the line passing through A and B.\n")
        print("                 A: ", self.line.A, "\n")
        print("                 B: ", self.line.B, "\n")
    
    def get3x3matrix(self):
        """Augmented matrix of the reflection.
        
        """
        line = self.line
        Q = line.A
        w = line.B - line.A
        wt = np.array([-w[1], w[0]])
        v = np.array([0.0, 0.0, 1.0])
        v_reshaped = v.reshape((3,1))
        M1 = np.hstack(
            (
                np.array([w, wt, Q]), 
                v_reshaped
            )    
        )
        M2 = np.hstack(
            (
                np.array([w, -wt, Q]), 
                v_reshaped
            )    
        )
        M = np.matmul(np.linalg.inv(M1), M2)
        M[:, 2] = M[2, :]
        M[2, :] = v
        return M
    
    def transform(self, P):
        """Transform a point by the refection.
        
        :param P: a point, `inf` allowed
        :returns: The image of `P`.
        
        """
        if is_inf(P):
            return inf
        _ = error_if_not_point_(P=P)
        P = np.asarray(P, dtype=float)
        line = self.line
        if line.includes(P):
            return P
        perp = line.perpendicular(P, False, False)
        return P + 2 * (perp.A - perp.B)
    
    def reflect(self, P):
        """An alias of `transform`.
        
        """
        return self.transform(P)
        
    def transform_circle(self, circ):
        """Reflect a circle.
        
        :param circ: a `Circle` object
        :returns: A `Circle` object.
        
        """
        return Circle(self.transform(circ.center), circ.radius)

    def transform_line(self, line):
        """Reflect a line.
        
        :param line: a `Line` object
        :returns: A `Line` object.
        
        """
        return Line(
            self.transform(line.A), self.transform(line.B),
            line.extendA, line.extendB
        )
    
    def as_affine(self):
        """Convert the reflection to an `Affine` object.
        
        """
        M = self.get3x3matrix()
        return Affine(M[0:2, 0:2], M[0:2, 2])


class Projection:
    """A class for projections. A projection on a line is given by the line of 
    projection `D` and the directrix line `Delta`. 
    
    For an orthogonal projection, you can also use the `projection` method of 
    the `Line` class.
    
    """
    def __init__(self, D, Delta):
        if not isinstance(D, Line):
            raise ValueError("`D` must be a `Line` object to initialize the projection.")
        if not Delta is None and not isinstance(Delta, Line):
            raise ValueError("The directrix `Delta` must be a `Line` object or `None` (for an orthogonal projection).")
        self.D = Line(D.A, D.B, True, True)
        if Delta is None:
            v = D.B - D.A
            self.Delta = Line((0,0), (-v[1], v[0]), True, True)
        else:
            self.Delta = Line(Delta.A, Delta.B, True, True)
        self.orthogonal = Delta is None

    def __str__(self):
        return str(self.__dict__)

    def show(self):
        print("Projection onto the line `D` passing through `A` and `B` parallel to the line `Delta` passing through `P` and `Q`.\n")
        print("                 A: ", self.D.A, "\n")
        print("                 B: ", self.D.B, "\n")
        print("                 P: ", self.Delta.P, "\n")
        print("                 Q: ", self.Delta.Q, "\n")
        if self.orthogonal:
            print("This is an orthogonal projection.\n")

    def project(self, M):
        """Projection of a point.
        
        :param M: a point
        :returns: a point on `D`, the projection of `M`
        
        """
        _ = error_if_not_point_(M=M)
        M = np.asarray(M)
        D = self.D
        if D.includes(M):
            return M
        Delta = self.Delta
        u = Delta.B - Delta.A
        do = D.direction_offset()
        ab = unit_vector_(do["direction"])
        k = - (dot_(ab, M) - do["offset"]) / dot_(ab, u)
        return k*u + M

    def transform(self, M):
        """An alias of `project`
        
        """
        return self.project(M)

    def get3x3matrix(self):
        """Augmented matrix of the projection.
        
        :returns: A 3x3 matrix.
        
        """
        b = self.project((0, 0))
        col1 = self.project((1, 0)) - b
        col2 = self.project((0, 1)) - b
        return np.hstack(
          (
            np.vstack((np.column_stack((col1, col2)), np.array([[0.0, 0.0]]))),
            np.array([[b[0]], [b[1]], [1.0]])
          )
        )

    def as_affine(self):
        """Convert the projection to an `Affine` object.
        
        """
        M = self.get3x3matrix()
        return Affine(M[0:2, 0:2], M[0:2, 2])

