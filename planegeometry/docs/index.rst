
Welcome to PyPlaneGeometry’s documentation!
*******************************************


Indices and tables
******************

*  `Index <genindex.rst>`_

*  `Module Index <py-modindex.rst>`_

*  `Search Page <search.rst>`_

**class planegeometry.geometry.Affine(A, b)**

   A class for affine transformations.

   An affine transformation is initialized by a 2x2 matrix (a linear
   transformation),  and a length two vector (the ‘intercept’, an
   array-like object).

   **compose(transfo, left=True)**

      Compose the reference affine map with another affine map.

      :Parameters:
         *  **transfo** – an *Affine* object

         *  **left** – Boolean, whether to compose at left or at right
            (i.e. returns *f1 o f0* or *f0 o f1*)

      :Returns:
         An *Affine* object.

   **classmethod from_ellipse_to_ellipse(ell1, ell2)**

      Affine transformation mapping a given ellipse to a given
      ellipse.

      :Parameters:
         **ell1****,****ell2** – *Ellipse* or *Circle* objects

      :Returns:
         An *Affine* object representing the transformation which maps
         *ell1* to *ell2*.

   **classmethod from_mapping_three_points(P1, P2, P3, Q1, Q2, Q3)**

      Affine transformation mapping three given points to three given
      points.

      :Parameters:
         *  **P1****,****P2****,****P3** – three non-collinear points

         *  **Q1****,****Q2****,****Q3** – three non-collinear points

      :Returns:
         An *Affine* object representing the transformation which maps
         Pi to Qi for each i=1,2,3.

   **get3x3matrix()**

      Get the 3x3 matrix corresponding to the affine transformation.

   **inverse()**

      The inverse affine transformation if it exists.

   **transform(m)**

      Transform a point or several points by the affine map.

      :Parameters:
         **m** – a point or a two-column matrix of points, one point
         per row

      :Returns:
         A matrix or a vector.

   **transform_ellipse(ell)**

      Transform an ellipse by the reference affine transformation
      (only for an invertible affine map).

      :Parameters:
         **ell** – an *Ellipse* object or a *Circle* object

      :Returns:
         An *Ellipse* object.

   **transform_line(line)**

      Transform a line by the affine map.

      :Returns:
         A *Line* object.

**class planegeometry.geometry.Arc(center, radius, alpha1, alpha2,
degrees=True)**

   Arc class.

   A circular arc is initialized by its center (array-like object of
   length  two), a radius, a starting angle and an ending angle. They
   are  respectively named *center*, *radius*, *alpha1* and *alpha2*.

   **ending_point()**

      Ending point of the arc.

   **path(n_points=100)**

      Path that forms the arc.

      :Parameters:
         **n_points** – number of points of the path

      :Returns:
         A matrix with two columns and *n_points* rows.

   **starting_point()**

      Starting point of the arc.

**class planegeometry.geometry.Circle(center, radius)**

   A class for circles. A circle is given by its center and its
   radius.

   **angle(circ2)**

      Angle between the reference circle and a given circle, if they
      intersect.

      :Parameters:
         **circ2** – a *Circle* object

      :Returns:
         The angle in radians.

   **as_ellipse()**

      Converts the circle to an *Ellipse* object.

      :Returns:
         An *Ellipse* object.

   **contains(P)**

      Check whether a point is contained in the reference circle.

      :Parameters:
         **P** – a point

      :Returns:
         A Boolean value.

   **includes(P)**

      Check whether a point belongs to the reference circle.

      :Parameters:
         **P** – a point

      :Returns:
         A Boolean value.

   **intersection_with_circle(circ2)**

      Intersection(s) of the reference circle with another circle.

      :Parameters:
         **circ2** – a *Circle* object

      :Returns:
         *None* (no intersection), a point (the two circles are
         tangent), or two points.

   **intersection_with_line(line)**

      Intersection(s) of the reference circle with a line.

      :Parameters:
         **line** – a *Line* object

      :Returns:
         *None* (no intersection), a point (the line is tangent to the
         circle), or two points.

   **is_equal(circ2)**

      Check whether the reference circle is equal to another circle.

      :Parameters:
         **circ2** – a *Circle* object

      :Returns:
         A Boolean value.

   **is_orthogonal(circ2)**

      Check whether the reference circle is orthogonal to a given
      circle.

      :Parameters:
         **circ2** – a *Circle* object

      :Returns:
         A Boolean value.

   **orthogonal_through_two_points_in_circle(P1, P2, arc=False)**

      Orthogonal circle passing through two points within the
      reference circle.

      :Parameters:
         *  **P1****,****P2** – two distinct points in the interior of
            the reference circle

         *  **arc** – Boolean, whether to return the arc joining the
            two points instead of the circle

      :Returns:
         A *Circle* object or an *Arc* object, or a *Line* object if
         the two points are on a diameter.

   **orthogonal_through_two_points_on_circle(alpha1, alpha2,
   arc=False)**

      Orthogonal circle passing through two points on the reference
      circle.

      :Parameters:
         *  **alpha1****,****alpha2** – two angles defining two points
            on the reference circle

         *  **arc** – Boolean, whether to return only the arc at the
            interior of the reference circle

      :Returns:
         A *Circle* object if *arc=False*, an *Arc* object if
         *arc=True, or a `Line* object: the diameter of the reference
         circle defined by the two points in case when the two angles
         differ by *pi*.

   **point_from_angle(alpha, degrees=True)**

      Get a point on the reference circle from its polar angle.

      :Parameters:
         *  **alpha** – a number, the angle

         *  **degrees** – Boolean, whether the angle is given in
            degrees

      :Returns:
         The point on the circle with polar angle *alpha*.

   **power(M)**

      Power of a point with respect to the reference circle.

      :Parameters:
         **M** – a point (array-like of length two)

      :Returns:
         A number, the power of *M* with respect to the circle.

   **radical_axis(circ2)**

      Radical axis of two circles.

      :Parameters:
         **circ2** – a *Circle* object

      :Returns:
         A *Line* object, the radical axis of the reference circle and
         *circ2*.

   **radical_center(circ2)**

      Radical center of two circles.

      :Parameters:
         **circ2** – a *Circle* object

      :Returns:
         The radical center of the reference circle and *circ2*.

   **random_points(n_points, where='in')**

      Random points in the circle or on the circle.

      :Parameters:
         *  **n_points** – desired number of points

         *  **where** – either *“in”* or *“on”*

      :Returns:
         A matrix with *n_points* rows and two columns; each row is a
         random point inside the circle if *where=”in”* or on the
         boundary of the circle if *where=”on”*.

**class planegeometry.geometry.Ellipse(center, rmajor, rminor, alpha,
degrees=True)**

   Ellipse class.

   An ellipse is initialized by its center (array-like object of
   length two),  its major radius, its minor radius, and the angle
   *alpha* between the  x-axis and the major axis.

   **classmethod LownerJohnEllipse(ipoints, bpoints=None)**

      Minimum area ellipse containing a set of points (ellipse hull).

      :Parameters:
         *  **ipoints** – an array of shape (n,2) containing points as
            row vectors, which are inside the desired ellipse; it must
            have at  least three distinct rows

         *  **bpoints** – an array of shape (n,2) containing points as
            row vectors, which are on the boundary of the desired
            ellipse;  could be *None* (the default)

      :Returns:
         An *Ellipse* object, the Löwner-John ellipse.

      :Author:
         Dor Shaviv

   **contains(P)**

      Check whether a point is contained in the ellipse.

      :Parameters:
         **P** – a point

      :Returns:
         A Boolean value.

   **equation()**

      The coefficients of the implicit equation of the ellipse,  *Ax²
      + Bxy + Cy² + Dx + Ey + F = 0*.

      :Returns:
         A dictionary giving the values of the coefficients.

   **classmethod equation_from_five_points(P1, P2, P3, P4, P5)**

      The implicit equation of the ellipse is *Ax² + Bxy + Cy² + Dx +
      Ey + F = 0*. This function returns A, B, C, D, E and F.

      :Parameters:
         **P1****,****P2****,****P3****,****P4****,****P5** – five
         points

      :Returns:
         A dictionary giving A, B, C, D, E and F.

   **classmethod from_boundary3(points)**

      Compute the smallest ellipse that passes through 3 boundary
      points.

      :Parameters:
         **points** – an array of shape (3,2) containing points as row
         vectors, which are on the boundary of the desired ellipse.

      :Returns:
         An *Ellipse* object.

      :Author:
         Dor Shaviv

   **classmethod from_boundary4(points)**

      Compute the smallest ellipse that passes through 4 boundary
      points, based on the algorithm by: B. W. Silverman and D. M.
      Titterington, “Minimum covering ellipses,” SIAM Journal on
      Scientific and Statistical Computing 1, no. 4 (1980): 401-409.

      :Parameters:
         **points** – an array of shape (4,2) containing points as row
         vectors, which are on the boundary of the desired ellipse.

      :Returns:
         An *Ellipse* object.

      :Author:
         Dor Shaviv

   **classmethod from_center_and_matrix(center, S)**

   **classmethod from_equation(A, B, C, D, E, F)**

      Ellipse from its implicit equation.

      :Parameters:
         **A****,****B****,****C****,****D****,****E****,****F** –
         coefficients of the implicit equation of the ellipse

      :Returns:
         An *Ellipse* object.

   **classmethod from_five_points(P1, P2, P3, P4, P5)**

      Ellipse from five points on this ellipse.

      :Parameters:
         **P1****,****P2****,****P3****,****P4****,****P5** – five
         points

      :Returns:
         An *Ellipse* object.

   **includes(P)**

      Check whether a point belongs to the ellipse.

      :Parameters:
         **P** – a point

      :Returns:
         A Boolean value.

   **intersection_with_line(line)**

   **is_equal(ell2)**

      Check whether the reference ellipse equals another ellipse.

      :Parameters:
         **ell2** – an *Ellipse* object

      :Returns:
         A Boolean value.

   **normal(t)**

      Normal unit vector to the ellipse.

      :Parameters:
         **t** – a number, the eccentric angle in radians of the point
         of the ellipse at which we want the normal unit vector

      :Returns:
         The normal unit vector to the ellipse at the point given by
         eccentric angle *t*.

   **path(n_points=100)**

      Path that forms the ellipse.

      :Parameters:
         **n_points** – number of points of the path

      :Returns:
         A matrix with two columns and *n_points* rows.

   **point_from_angle(theta, degrees=True)**

   **random_points(n_points, where='in')**

      Random points in/on the ellipse.

      :Parameters:
         *  **n_points** – desired number of points

         *  **where** – either *“in”* or *“on”*

      :Returns:
         A matrix with *n_points* rows and two columns; each row is a
         random point inside the ellipse if *where=”in”* or on the
         boundary of the ellipse if *where=”on”*.

   **shape_matrix()**

      The 2x2 symmetric matrix *S* associated to the reference
      ellipse.  The equation of the ellipse is *(M-O)’ S (M-O) = 1*.

   **theta2t(theta, degrees=True)**

      Convert angle to eccentric angle.

      :Parameters:
         *  **theta** – angle between the major axis and the half-line
            starting at the center of the ellipse and passing through
            the point of interest on the ellipse

         *  **degrees** – Boolean, whether *theta* is given in degrees

      :Returns:
         The eccentric angle of the point of interest on the ellipse,
         in radians.

**class planegeometry.geometry.Homothety(center, scale)**

   A homothety is given by a center and a scale factor.

   **as_affine()**

      Converting to an *Affine* object.

   **get3x3matrix()**

      Get the augmented matrix of the homothety.

   **transform(M)**

      Transform one or more points.

      :Parameters:
         **M** – a point or a matrix of points

      :Returns:
         A point or a matrix of points.

   **transform_circle(circ)**

      Transform a circle by the homothety.

      :Parameters:
         **circ** – a *Circle* object

      :Returns:
         A *Circle* object.

**class planegeometry.geometry.Inversion(pole, power)**

   Inversion class.

   An inversion is initialized by its pole and its power (a number,
   possibly negative).

   **compose(iota2, left=True)**

      Compose the reference inversion with another inversion. The
      result  is a Möbius transformation.

      :Parameters:
         *  **iota2** – an *Inversion* object

         *  **left** – Boolean, whether to compose at left or at right
            (i.e. returns *iota2 o iota1* or *iota1 o iota2*)

      :Returns:
         A *Mobius* object.

   **classmethod from_fixing_three_circles(circ1, circ2, circ3)**

      Inversion fixing three circles.

      :Parameters:
         **circ1****,****circ2****,****circ3** – *Circle* objects

      :Returns:
         An *Inversion* object representing an inversion which leaves
         each of the three circles invariant.

   **classmethod from_fixing_two_circles(circ1, circ2)**

      Inversion fixing two circles.

      :Parameters:
         **circ1****,****circ2** – *Circle* objects

      :Returns:
         An *Inversion* object representing an inversion which leaves
         each of the two circles invariant.

   **classmethod from_swapping_two_circles(circ1, circ2,
   positive=True)**

      Inversion swapping two circles.

      :Parameters:
         *  **circ1****,****circ2** – *Circle* objects

         *  **positive** – Boolean, whether the sign of the desired
            inversion power must be positive or negative

      :Returns:
         An *Inversion* object, which maps *circ1* to *circ2* and
         *circ2* to *circ1*, except in the case when *circ1* and
         *circ2* are congruent and tangent: in this case a
         *Reflection* object is returned (a reflection is an inversion
         on a line).

   **invert_gcircle(gcircle)**

      Invert a generalized circle, that is, a circle or a line.

      :Params gcircle`:
         a *Circle* object or a *Line* object

      :Returns:
         A *Circle* object or a *Line* object.

   **classmethod on_circle(circ)**

      An inversion on a circle is the inversion whose pole is the
      center  of the circle and whose power is the squared radius of
      the circle.

      :Parameters:
         **circ** – *Circle* object

      :Returns:
         An *Inversion* object.

**class planegeometry.geometry.Line(A, B, extendA=True,
extendB=True)**

   A class for lines. A line is initialized by two points it passes
   through,  and for each of these points a Boolean value to indicate
   whether the line  should be extended besides this point.

   **direction_offset()**

      Direction and offset of the line. The equation of the line is
      *cos(direction)x + sin(direction)y = offset*.

      :Returns:
         The direction and the offset in a dictionary.

   **distance(M)**

      Distance from a point to the reference line.

      :Parameters:
         **M** – a point

      :Returns:
         A number.

   **includes(M, strict=False, checkCollinear=True)**

      Check whether a point belongs to the line.

      :Parameters:
         *  **M** – the point for which we want to test whether it
            belongs to the line

         *  **strict** – Boolean, whether to take into account
            *extendA* and *extendB*

         *  **checkCollinear** – Boolean, whether to check the
            collinearity of *A*, *B*, and *M*; set to *False* only if
            you use *strict=True* and you are sure that *M* is on the
            line (AB)

      :Returns:
         A Boolean value.

   **intersection_with_circle(circ)**

      Intersection(s) of the line with a circle.

      :Parameters:
         **circ** – a *Circle* object

      :Returns:
         *None*, a point, or a list of two points.

   **intersection_with_ellipse(ell)**

      Intersection(s) of the line with an ellipse.

      :Parameters:
         **ell** – an *Ellipse* object

      :Returns:
         *None*, a point, or a list of two points.

   **intersection_with_line(line2, strict=False)**

      Intersection(s) of the reference line with another line.

      :Parameters:
         **line2** – a *Line* object

      :Returns:
         *None* (the lines are parallel), a point, or a *Line* object
         (the two lines are equal).

   **invert(iota)**

      Invert the reference line.

      :Parameters:
         **iota** – an *Inversion.fro* object

      :Returns:
         A *Line* object or a *Circle* object.

   **is_equal(line2)**

      Check whether the reference line is equal to another line.

      :Parameters:
         **line2** – a *Line* object

      :Returns:
         A Boolean value.

   **is_parallel(line2)**

      Check whether the reference line is parallel to another line.

      :Parameters:
         **line2** – a *Line* object

      :Returns:
         A Boolean value.

   **perpendicular(M, extendH=False, extendM=True)**

      Perpendicular line passing through a given point.

      :Parameters:
         *  **M** – the point through which the perpendicular passes

         *  **extendH** – Boolean, whether to extend the perpendicular
            line beyond the meeting point

         *  **extendM** – Boolean, whether to extend the perpendicular
            line beyond the point *M*

      :Returns:
         A *Line* object; its two points are the meeting point and the
         point *M*.

   **projection(M)**

      Orthogonal projection of a point to the reference line.

      :Parameters:
         **M** – a point

      :Returns:
         A point on the reference line.

   **reflection(M)**

      Reflection of a point with respect to the reference line.

      :Parameters:
         **M** – a point

      :Returns:
         A point.

   **rotate(alpha, O, degrees=True)**

      Rotate the reference line.

      :Parameters:
         *  **alpha** – angle of rotation

         *  **O** – center of rotation

         *  **degrees** – Boolean, whether *alpha* is given in degrees

      :Returns:
         A *Line* object.

   **translate(v)**

      Translate the reference line.

      :Parameters:
         **v** – the vector of translation

      :Returns:
         A *Line* object.

**class planegeometry.geometry.Mobius(M)**

   A class for Möbius transformations.

   A Möbius transformation is initialized by a complex 2x2 matrix with
   a   non-zero determinant.

   **compose(M1, left=True)**

      Compose the reference Möbius transformation with another Möbius
      transformation.

      :Parameters:
         *  **M1** – a *Mobius* object

         *  **left** – Boolean, whether to compose at left or at right
            (i.e. returns *M1 o M0* or *M0 o M1*)

      :Returns:
         A *Mobius* object.

   **classmethod from_mapping_three_points(P1, P2, P3, Q1, Q2, Q3)**

      Möbius transformation mapping three given points to three given
      points.

      :Parameters:
         *  **P1****,****P2****,****P3** – three distinct points,
            *inf* allowed

         *  **Q1****,****Q2****,****Q3** – three distinct points,
            *inf* allowed

      :Returns:
         A *Mobius* object, representing the Möbius transformation
         which sends *Pi* to *Qi* for each i=1,2,3.

   **gpower(t)**

      Generalized power of the Möbius transformation.

      :Parameters:
         **t** – a float, possibly negative

      :Returns:
         A *Mobius* object corresponding to the Möbius transformation
         raised to the power *t*.

   **inverse()**

      Inverse of the Möbius transformation.

      :Returns:
         A Möbius transformation.

   **power(k)**

      Power of the Möbius transformation.

      :Parameters:
         **k** – an integer, possibly negative

      :Returns:
         A *Mobius* object corresponding to the Möbius transformation
         raised to the power *k*.

   **transform(P)**

      Transform a point by the Möbius transformation.

      :Parameters:
         **P** – a point (array-like of length two) or *inf*

      :Returns:
         The image of *P* by the Möbius transformation (can be *inf*).

   **transform_circle(circ)**

      Transform a circle by the Möbius transformation.

      :Parameters:
         **circ** – a *Circle* object

      :Returns:
         A *Circle* object or a *Line* object.

   **transform_line(line)**

      Transform a line by the Möbius transformation.

      :Parameters:
         **line** – a *Line* object

      :Returns:
         A *Circle* object or a *Line* object.

**class planegeometry.geometry.Projection(D, Delta)**

   A class for projections. A projection on a line is given by the
   line of  projection *D* and the directrix line *Delta*.

   For an orthogonal projection, you can also use the *projection*
   method of  the *Line* class.

   **as_affine()**

      Convert the projection to an *Affine* object.

   **get3x3matrix()**

      Augmented matrix of the projection.

      :Returns:
         A 3x3 matrix.

   **project(M)**

      Projection of a point.

      :Parameters:
         **M** – a point

      :Returns:
         A point on *D*, the projection of *M*.

   **transform(M)**

      An alias of *project*

**class planegeometry.geometry.Reflection(line)**

   A class for reflections.

   A reflection is initialized by a line.

   **as_affine()**

      Convert the reflection to an *Affine* object.

   **get3x3matrix()**

      Augmented matrix of the reflection.

   **reflect(P)**

      An alias of *transform*.

   **reflect_circle(circ)**

      Reflect a circle.

      :Parameters:
         **circ** – a *Circle* object

      :Returns:
         A *Circle* object.

   **reflect_line(line)**

      Reflect a line.

      :Parameters:
         **line** – a *Line* object

      :Returns:
         A *Line* object.

   **transform(P)**

      Transform a point by the refection.

      :Parameters:
         **P** – a point, *inf* allowed

      :Returns:
         The image of *P*.

   **transform_circle(circ)**

      An alias of *reflect_circle*.

   **transform_line(line)**

      An alias of *reflect_line*.

**class planegeometry.geometry.Rotation(center, theta, degrees=True)**

   A class for rotations. A rotation is given by its center and its
   angle.

   :Parameters:
      *  **center** – a point

      *  **theta** – a number, the angle of the rotation

      *  **degrees** – Boolean, whether the angle is given in degrees

   **as_affine()**

      Convert the rotation to an *Affine* object.

   **get3x3matrix()**

      Augmented matrix of the rotation.

   **rotate(M)**

      Rotate a point.

      :Parameters:
         **M** – a point or a matrix of points

      :Returns:
         A point or a matrix of points.

   **rotate_circle(circ)**

      Rotate a circle.

      :Parameters:
         **circ** – a *Circle* object

      :Returns:
         A *Circle* object.

   **rotate_ellipse(ell)**

      Rotate an ellipse.

      :Parameters:
         **ell** – an *Ellipse* object

      :Returns:
         An *Ellipse* object.

   **rotate_line(line)**

      Rotate a line.

      :Parameters:
         **line** – a *Line* object

      :Returns:
         A *Line* object.

   **transform(M)**

      An alias of *rotate*.

   **transform_circle(circ)**

      An alias of *rotate_circle*.

   **transform_ellipse(ell)**

      An alias of *rotate_ellipse*.

   **transform_line(line)**

      An alias of *rotate_line*.

**class planegeometry.geometry.ScalingXY(center, sx, sy)**

   A class for axis-scalings. An axis-scaling is given by a center,
   and  two scale factors *sx* and *sy*, one for the x-axis and one
   for the y-axis.

   **as_affine()**

      Converting to an *Affine* object.

   **get3x3matrix()**

      Get the augmented matrix of the axes-scaling.

   **transform(M)**

      Transform one or more points.

      :Parameters:
         **M** – a point or a matrix of points

      :Returns:
         A point or a matrix of points.

**class planegeometry.geometry.Shear(vertex, vector, ratio, angle,
degrees=True)**

   A class for shear transformations. A shear is given by a vertex,
   two perpendicular vectors, and an angle.

   :Example:
   
   >>> P = np.array([0,0]); w = np.array([1,0]); ratio = 1; angle = 45
   >>> shear = Shear(P, w, ratio, angle)
   >>> wt = ratio * np.array([-w[1], w[0]])
   >>> Q = P + w; R = Q + wt; S = P + wt
   >>> A = shear.transform(P); B = shear.transform(Q)
   >>> C = shear.transform(R); D = shear.transform(S)
   >>> import matplotlib.pyplot as plt
   >>> figure, axes = plt.subplots(figsize=(10, 10))
   >>> axes.set_aspect(1)
   >>> unit_square = plt.Polygon([P,Q,R,S], fill=False)
   >>> axes.add_artist(unit_square)
   >>> image = plt.Polygon([A,B,C,D], fill=False, color="red")
   >>> axes.add_artist(image)
   >>> plt.xlim(0, 1)
   >>> plt.ylim(0, 2)
   >>> plt.show()

   **as_affine()**

      Convert the shear to an *Affine* object.

   **get3x3matrix()**

      Get the augmented matrix of the shear.

   **transform(M)**

      Transform one or more points.

      :Parameters:
         **M** – a point or a matrix of points

      :Returns:
         A point or a matrix of points.

**class planegeometry.geometry.Triangle(A, B, C)**

   Triangle class.

   A triangle is initialized by its three vertices, some array-like
   objects  of length two.

   ``property a``

      Length of the side BC.

   ``property angleA``

      The angle at the vertex A in radians.

   ``property angleB``

      The angle at the vertex B in radians.

   ``property angleC``

      The angle at the vertex C in radians.

   ``property b``

      Length of the side AC.

   ``property c``

      Length of the side AB.

   **contains(M)**

      Check whether a point lies inside the reference triangle.

   ``property edges``

      Edge lengths of the triangle.

   **equal_detour_point()**

      Equal detour point of the triangle, also known as the X(176)
      triangle center.

      :Returns:
         A pair, the equal detour point and the detour.

   **excircles()**

      The excircles of the triangle.

      :Returns:
         A dictionary of three *Circle* objects.

   ``property flatness``

      Flatness, a number between 0 and 1; a triangle is flat when its
      flatness is 1.

   **incircle()**

      The incircle of the triangle.

   ``property is_acute``

      Check whether the triangle is acute.

   **malfatti_circles()**

      The Malfatti circles of the triangle.

      :Returns:
         Three circles and three tangency points.

   ``property orientation``

      Orientation of the triangle; 1 for counterclockwise, -1 for
      clockwise, 0 for collinear.

   **orthic_triangle()**

      Orthic triangle. Its vertices are the feet of the altitudes of
      the reference triangle.

      :Returns:
         A *Triangle* object.

   **random_points(n_points, where='in')**

      Random points inside the triangle or on the boundary of the
      triangle.

      :Parameters:
         *  **n_points** – desired number of points

         *  **where** – either *“in”* or *“on”*

      :Returns:
         A matrix with *n_points* rows and two columns; each row is a
         random point inside the triangle if *where=”in”* or on the
         boundary of the triangle if *where=”on”*.

   **steiner_ellipse()**

      The Steiner ellipse (or circumellipse) of the reference
      triangle.  This is the ellipse passing through the three
      vertices of the triangle  and centered at the centroid of the
      triangle.

      :Returns:
         An *Ellipse* object.

   **steiner_inellipse()**

      The Steiner inellipse (or midpoint ellipse) of the reference
      triangle.  This is the ellipse tangent to the sides of the
      triangle at their  midpoints, and centered at the centroid of
      the triangle.

      :Returns:
         An *Ellipse* object.

**planegeometry.geometry.circleAB(A, B)**

   Circle with diameter AB.

   :Parameters:
      **A****,****B** – two distinct points

   :Returns:
      A *Circle* object.

**planegeometry.geometry.intersection_circle_circle(circ1, circ2)**

   Intersection(s) of two circles.

   :Parameters:
      **circ1****,****circ2** – *Circle* objects

   :Returns:
      A *Circle* object if the two circles are equal, *None* if the
      two circles do not intersect, a point if the two circles are
      tangent, or a list of two points.

**planegeometry.geometry.intersection_circle_line(circ, line)**

   Intersection(s) of a circle and a line.

   :Parameters:
      *  **circ** – a *Circle* object

      *  **line** – a *Line* object

   :Returns:
      *None* if the intersection is empty, otherwise either one point
      (the line is tangent to the circle) or a list of two points.

**planegeometry.geometry.intersection_ellipse_line(ell, line)**

   Intersection(s) of an ellipse and a line.

   :Parameters:
      *  **ell** – an *Ellipse* object

      *  **line** – a *Line* object

   :Returns:
      *None* if the intersection is empty, otherwise either one point
      (the line is tangent to the ellipse) or a list of two points.

**planegeometry.geometry.mid_circles(circ1, circ2)**

   Return the mid-circle(s) of two circles. A mid-circle of two
   circles is  a generalized circle (i.e. a circle or a line) such
   that the Inversion.fro on  this circle swaps the two circles. The
   case of a line appears only when  the two circles have equal radii.

   :Parameters:
      **circ1****,****circ2** – *Circle* objects

   :Returns:
      A *Circle* object, or a *Line* object, or a list of two such
      objects.

**planegeometry.geometry.radical_center(circ1, circ2, circ3)**

   Radical center of three circles.

   :Parameters:
      **circ1****,****circ2****,****circ3** – *Circle* objects

   :Returns:
      A point, the radical center of the three circles.

**planegeometry.geometry.unimodular_matrices(n)**

   Generates unimodular matrices.

   :Parameters:
      **n** – integer, the maximum size of entries of matrices, at
      least 1

   :Returns:
      List of unimodular matrices.
