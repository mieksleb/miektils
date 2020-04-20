module geom
  implicit none

contains


subroutine ij_next ( i, j, n )

!*****************************************************************************80
!
!! IJ_NEXT returns the next matrix index.
!
!  Discussion:
!
!    For N = 3, the sequence of indices returned is:
!
!      (1,1), (1,2), (1,3), (2,1), (2,2), (2,3), (3,1), (3,2), (3,3), (0,0).
!
!    Note that once the value (N,N) is returned, the next value returned
!    will be (0,0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On input, the current pair of
!    indices.  On output, the next pair of indices.  If either index is illegal
!    on input, the output value of (I,J) will be (1,1).
!
!    Input, integer ( kind = 4 ) N, the maximum value for I and J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n

  if ( n < 1 ) then
    i = 0
    j = 0
    return
  end if

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j ) then
    i = 1
    j = 1
    return
  end if

  if ( j < n ) then
    j = j + 1
  else if ( i < n ) then
    i = i + 1
    j = 1
  else
    i = 0
    j = 0
  end if

  return
end
subroutine ij_next_gt ( i, j, n )

!*****************************************************************************80
!
!! IJ_NEXT_GT returns the next matrix index, with the constraint that I < J.
!
!  Discussion:
!
!    For N = 5, the sequence of indices returned is:
!
!      (1,2), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5), (3,4), (3,5), (4,5).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) I, J.  On input, the current pair of
!    indices.  On output, the next pair of indices.  If either index is illegal
!    on input, the output value of (I,J) will be (1,2).
!
!    Input, integer ( kind = 4 ) N, the maximum value for I and J.
!    A value of N less than 2 is nonsense.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n

  if ( n < 2 ) then
    i = 0
    j = 0
    return
  end if

  if ( i < 1 .or. n < i .or. j < 1 .or. n < j .or. j <= i ) then
    i = 1
    j = 2
    return
  end if

  if ( j < n ) then
    j = j + 1
  else if ( i < n - 1 ) then
    i = i + 1
    j = i + 1
  else
    i = 0
    j = 0
  end if

  return
end
subroutine line_exp2imp_2d ( x1, y1, x2, y2, a, b, c )

!*****************************************************************************80
!
!! LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2.  (X1,Y1) and (X2,Y2) are
!    two points on the line. (X1,Y1) must be different
!    from (X2,Y2).
!
!    Output, real ( kind = 8 ) A, B, C, three coefficients which describe
!    the line that passes through (X1,Y1) and (X2,Y2).
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
!
!  Take care of degenerate cases.
!
  if ( x1 == x2 .and. y1 == y2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LINE_EXP2IMP_2D - Fatal error!'
    write ( *, '(a)' ) '  (X1,Y1) = (X2,Y2)'
    write ( *, '(a,g14.6,2x,g14.6)' ) '  (X1,Y1) = ', x1, y1
    write ( *, '(a,g14.6,2x,g14.6)' ) '  (X2,Y2) = ', x2, y2
    stop
  end if

  a = y2 - y1
  b = x1 - x2
  c = x2 * y1 - x1 * y2

  return
end
subroutine line_exp_point_dist_2d ( x1, y1, x2, y2, x, y, dist )

!*****************************************************************************80
!
!! LINE_EXP_POINT_DIST_2D: distance ( explicit line, point ) in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2.  (X1,Y1) and (X2,Y2) are two
!    points on the line.
!
!    Input, real ( kind = 8 ) X, Y, the point whose distance from the line is
!    to be measured.
!
!    Output, real ( kind = 8 ) DIST, the distance from the point to the line.
!
  implicit none

  real ( kind = 8 ) bot
  real ( kind = 8 ) dist
  real ( kind = 8 ) dot
  real ( kind = 8 ) t
  real ( kind = 8 ) x
  real ( kind = 8 ) xn
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y
  real ( kind = 8 ) yn
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2

  bot = ( x1 - x2 )**2 + ( y1 - y2 )**2

  if ( bot == 0.0D+00 ) then

    xn = x1
    yn = y1
!
!  (P-P1) dot (P2-P1) = Norm(P-P1) * Norm(P2-P1) * Cos(Theta).
!
!  (P-P1) dot (P2-P1) / Norm(P-P1)**2 = normalized coordinate T
!  of the projection of (P-P1) onto (P2-P1).
!
  else

    dot = ( x - x1 ) * ( x2 - x1 ) + ( y - y1 ) * ( y2 - y1 )

    t = dot / bot

    xn = x1 + t * ( x2 - x1 )
    yn = y1 + t * ( y2 - y1 )

  end if

  dist = sqrt ( ( xn - x )**2 + ( yn - y )**2 )

  return
end
subroutine line_exp_point_dist_signed_2d ( x1, y1, x2, y2, x, y, dist_signed )

!*****************************************************************************80
!
!! LINE_EXP_POINT_DIST_SIGNED_2D: signed distance ( explicit line, point ).
!
!  Discussion:
!
!    The signed distance has two interesting properties:
!
!    *  The absolute value of the signed distance is the
!        usual (Euclidean) distance.
!
!    *  Points with signed distance 0 lie on the line,
!       points with a negative signed distance lie on one side
!         of the line,
!       points with a positive signed distance lie on the
!         other side of the line.
!
!    Assuming that C is nonnegative, then if a point is a positive
!    distance away from the line, it is on the same side of the
!    line as the point (0,0), and if it is a negative distance
!    from the line, it is on the opposite side from (0,0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, define the two points
!    (X1,Y1) and (X2,Y2) that determine the line.
!
!    Input, real ( kind = 8 ) X, Y, the point (X,Y) whose signed distance
!    is desired.
!
!    Output, real ( kind = 8 ) DIST_SIGNED, the signed distance from the
!    point to the line.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) dist_signed
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
!
!  Convert the line to implicit form.
!
  call line_exp2imp_2d ( x1, y1, x2, y2, a, b, c )
!
!  Compute the signed distance from the point to the line.
!
  dist_signed = ( a * x + b * y + c ) / sqrt ( a * a + b * b )

  return
end
subroutine line_seg_contains_point_2d ( x1, y1, x2, y2, x3, y3, u, v )

!*****************************************************************************80
!
!! LINE_SEG_CONTAINS_POINT_2D reports if a line segment contains a point in 2D.
!
!  Discussion:
!
!    In exact arithmetic, point P3 = (X3,Y3) is on the line segment between
!    P1=(X1,Y1) and P2=(X2,Y2) if and only if 0 <= U <= 1 and V = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, the endpoints P1 and P2 of
!    a line segment.
!
!    Input, real ( kind = 8 ) X3, Y3, a point P3 to be tested.
!
!    Output, real ( kind = 8 ) U, the coordinate of (X3,Y3) along the axis from
!    with origin at P1 and unit at P2.
!
!    Output, real ( kind = 8 ) V, the magnitude of the off-axis portion of the
!    vector P3-P1, measured in units of (P2-P1).
!
  implicit none

  real ( kind = 8 ) u
  real ( kind = 8 ) unit
  real ( kind = 8 ) v
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3

  unit = sqrt ( ( x2 - x1 )**2 + ( y2 - y1 )**2 )

  if ( unit == 0.0D+00 ) then

    if ( x3 == x1 .and. y3 == y1 ) then
      u = 0.5D+00
      v = 0.0D+00
    else
      u = 0.5D+00
      v = huge ( 1.0D+00 )
    end if

  else

    u = ( ( x3 - x1 ) * ( x2 - x1 ) + ( y3 - y1 ) * ( y2 - y1 ) ) / unit**2

    v = sqrt ( ( ( u - 1.0D+00 ) * x1 - u * x2 + x3 )**2 &
             + ( ( u - 1.0D+00 ) * y1 - u * y2 + y3 )**2 ) / unit

  end if

  return
end
subroutine line_seg_vec_int_2d ( n, x1, y1, x2, y2, i, j, flag, xint, yint )

!*****************************************************************************80
!
!! LINE_SEG_VEC_INT_2D computes intersections of a set of line segments.
!
!  Discussion:
!
!    This is an implementation of the relatively inefficient algorithm
!    for computing all intersections of a set of line segments.
!
!    This is an "incremental" code, which returns as soon as it finds
!    a single intersection.  To find the next intersection, simply call
!    again.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of line segments.
!
!    Input, real ( kind = 8 ) X1(N), Y1(N), X2(N), Y2(N), the coordinates of the
!    endpoints of the line segments.
!
!    Input/output, integer ( kind = 4 ) I, J, used to keep track of the
!    computation.  On first call with a given set of line segments,
!    set I = J = 0.  On return with FLAG = 1, I and J record the indices of the
!    line segments whose intersection has just been found.  To find the
!    next intersection, simply call again, but do not alter I and J.
!
!    Output, integer ( kind = 4 ) FLAG:
!    0, no more intersections, the computation is done.
!    1, an intersection was detected between segments I and J.
!
!    Output, real ( kind = 8 ) XINT, YINT, the location of an intersection
!    of line segments I and J, if FLAG is 1.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) flag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) x1(n)
  real ( kind = 8 ) x2(n)
  real ( kind = 8 ) xint
  real ( kind = 8 ) y1(n)
  real ( kind = 8 ) y2(n)
  real ( kind = 8 ) yint

  do

    call ij_next_gt ( i, j, n )

    if ( i == 0 ) then
      flag = 0
      exit
    end if

    call lines_seg_int_2d ( x1(i), y1(i), x2(i), y2(i), x1(j), y1(j), &
      x2(j), y2(j), flag, xint, yint )

    if ( flag == 1 ) then
      return
    end if

  end do

  return
end
subroutine lines_exp_int_2d ( ival, x, y, x1, y1, x2, y2, x3, y3, x4, y4 )

!*****************************************************************************80
!
!! LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
!
!  Discussion:
!
!    The explicit form of a line in 2D is:
!
!      (X1,Y1), (X2,Y2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IVAL, reports on the intersection.
!
!     0, no intersection, the lines may be parallel or degenerate.
!     1, one intersection point, returned in X, Y.
!     2, infinitely many intersections, the lines are identical.
!
!    Output, real ( kind = 8 ) X, Y, if IVAl = 1, then X, Y contains
!    the intersection point.  Otherwise, X = 0, Y = 0.
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, define the first line.
!
!    Input, real ( kind = 8 ) X3, Y3, X4, Y4, define the second line.
!
  implicit none

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  integer ( kind = 4 ) ival
  logical point_1
  logical point_2
  real ( kind = 8 ) x
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4

  ival = 0
  x = 0.0D+00
  y = 0.0D+00
!
!  Check whether either line is a point.
!
  if ( x1 == x2 .and. y1 == y2 ) then
    point_1 = .true.
  else
    point_1 = .false.
  end if

  if ( x3 == x4 .and. y3 == y4 ) then
    point_2 = .true.
  else
    point_2 = .false.
  end if
!
!  Convert the lines to ABC format.
!
  if ( .not. point_1 ) then
    call line_exp2imp_2d ( x1, y1, x2, y2, a1, b1, c1 )
  end if

  if ( .not. point_2 ) then
    call line_exp2imp_2d ( x3, y3, x4, y4, a2, b2, c2 )
  end if
!
!  Search for intersection of the lines.
!
  if ( point_1 .and. point_2 ) then
    if ( x1 == x3 .and. y1 == y3 ) then
      ival = 1
      x = x1
      y = y1
    end if
  else if ( point_1 ) then
    if ( a2 * x1 + b2 * y1 == c2 ) then
      ival = 1
      x = x1
      y = y1
    end if
  else if ( point_2 ) then
    if ( a1 * x3 + b1 * y3 == c1 ) then
      ival = 1
      x = x3
      y = y3
    end if
  else
    call lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, x, y )
  end if

  return

end

subroutine lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, x, y )

!*****************************************************************************80
!
!! LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
!
!  Discussion:
!
!    The implicit form of a line in 2D is:
!
!      A * X + B * Y + C = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A1, B1, C1, define the first line.
!    At least one of A1 and B1 must be nonzero.
!
!    Input, real ( kind = 8 ) A2, B2, C2, define the second line.
!    At least one of A2 and B2 must be nonzero.
!
!    Output, integer ( kind = 4 ) IVAL, reports on the intersection.
!
!    -1, both A1 and B1 were zero.
!    -2, both A2 and B2 were zero.
!     0, no intersection, the lines are parallel.
!     1, one intersection point, returned in X, Y.
!     2, infinitely many intersections, the lines are identical.
!
!    Output, real ( kind = 8 ) X, Y, if IVAL = 1, then X, Y contains
!    the intersection point.  Otherwise, X = 0, Y = 0.
!
  implicit none

  real ( kind = 8 ) a(2,2)
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b(2,2)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) det
  integer ( kind = 4 ) ival
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  x = 0.0D+00
  y = 0.0D+00
!
!  Refuse to handle degenerate lines.
!
  if ( a1 == 0.0D+00 .and. b1 == 0.0D+00 ) then
    ival = -1
    return
  else if ( a2 == 0.0D+00 .and. b2 == 0.0D+00 ) then
    ival = -2
    return
  end if
!
!  Set up a linear system, and compute its inverse.
!
  a(1,1) = a1
  a(1,2) = b1
  a(2,1) = a2
  a(2,2) = b2

  call r8mat2_inverse ( a, b, det )
!
!  If the inverse exists, then the lines intersect.
!  Multiply the inverse times -C to get the intersection point.
!
  if ( det /= 0.0D+00 ) then

    ival = 1
    x = - b(1,1) * c1 - b(1,2) * c2
    y = - b(2,1) * c1 - b(2,2) * c2
!
!  If the inverse does not exist, then the lines are parallel
!  or coincident.  Check for parallelism by seeing if the
!  C entries are in the same ratio as the A or B entries.
!
  else

    ival = 0

    if ( a1 == 0.0D+00 ) then
      if ( b2 * c1 == c2 * b1 ) then
        ival = 2
      end if
    else
      if ( a2 * c1 == c2 * a1 ) then
        ival = 2
      end if
    end if

  end if

  return
end


subroutine r8mat2_inverse ( a, b, det )

!*****************************************************************************80
!
!! R8MAT2_INVERSE inverts a 2 by 2 real matrix using Cramer's rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(2,2), the matrix to be inverted.
!
!    Output, real ( kind = 8 ) B(2,2), the inverse of the matrix A.
!
!    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
!
!    If DET is zero, then A is singular, and does not have an
!    inverse.  In that case, B is simply set to zero, and a
!    message is printed.
!
!    If DET is nonzero, then its value is roughly an estimate
!    of how nonsingular the matrix A is.
!
  implicit none

  real ( kind = 8 ) a(2,2)
  real ( kind = 8 ) b(2,2)
  real ( kind = 8 ) det
!
!  Compute the determinant.
!
  det = a(1,1) * a(2,2) - a(1,2) * a(2,1)
!
!  If the determinant is zero, bail out.
!
  if ( det == 0.0D+00 ) then

    b(1:2,1:2) = 0.0D+00

    return
  end if
!
!  Compute the entries of the inverse matrix using an explicit formula.
!
  b(1,1) = + a(2,2) / det
  b(1,2) = - a(1,2) / det
  b(2,1) = - a(2,1) / det
  b(2,2) = + a(1,1) / det

  return
end

subroutine lines_seg_dist_2d ( x1, y1, x2, y2, x3, y3, x4, y4, flag, x5, y5 )

!*****************************************************************************80
!
!! LINES_SEG_DIST_2D computes the distance of two line segments in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, the endpoints of the first
!    segment.
!
!    Input, real ( kind = 8 ) X3, Y3, X4, Y4, the endpoints of the second
!    segment.
!
!    Output, integer ( kind = 4 ) FLAG, records the results.
!    0, the line segments do not intersect.
!    1, the line segments intersect.
!
!    Output, real ( kind = 8 ) X5, Y5.
!    If FLAG = 0, X5 = Y5 = 0.
!    If FLAG = 1, then (X5,Y5) is a point of intersection.
!
  implicit none

  integer ( kind = 4 ) flag
  integer ( kind = 4 ) ival
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) y5

  x5 = 0.0D+00
  y5 = 0.0D+00
!
!  Find the intersection of the two lines.
!
  call lines_exp_int_2d ( ival, x5, y5, x1, y1, x2, y2, x3, y3, x4, y4 )

  if ( ival == 0 ) then
    flag = 0
    return
  end if
!
!  Is the point on the first segment?
!
  call line_seg_contains_point_2d ( x1, y1, x2, y2, x5, y5, u, v )

  if ( u < 0 .or. 1.0D+00 < u .or. 0.001D+00 < v ) then
    flag = 0
    return
  end if
!
!  Is the point on the second segment?
!
  call line_seg_contains_point_2d ( x3, y3, x4, y4, x5, y5, u, v )

  if ( u < 0 .or. 1.0D+00 < u .or. 0.001D+00 < v ) then
    flag = 0
    return
  end if

  flag = 1

  return
end
subroutine lines_seg_int_1d ( x1, x2, x3, x4, flag, x5, x6 )

!*****************************************************************************80
!
!! LINES_SEG_INT_1D computes the intersection of two line segments in 1D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, X2, the endpoints of the first segment.
!
!    Input, real ( kind = 8 ) X3, X4, the endpoints of the second segment.
!
!    Output, integer ( kind = 4 ) FLAG, records the results.
!    0, the line segments do not intersect.
!    1, the line segments intersect.
!
!    Output, real ( kind = 8 ) X5, X6, the endpoints of the intersection
!    segment.  If FLAG = 0, X5 = X6 = 0.
!
  implicit none

  integer ( kind = 4 ) flag
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) x6
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4

  y1 = min ( x1, x2 )
  y2 = max ( x1, x2 )
  y3 = min ( x3, x4 )
  y4 = max ( x3, x4 )

  flag = 0
  x5 = 0.0D+00
  x6 = 0.0D+00

  if ( y4 < y1 ) then
    return
  else if ( y2 < y3 ) then
    return
  end if

  flag = 1
  x5 = max ( y1, y3 )
  x6 = min ( y2, y4 )

  return
end


subroutine lines_seg_int_2d(x1, y1, x2, y2, x3, y3, x4, y4, flag, x5, y5 )

!*****************************************************************************80
!
!! LINES_SEG_INT_2D computes the intersection of two line segments in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 September 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X1, Y1, X2, Y2, the endpoints of the first
!    segment.
!
!    Input, real ( kind = 8 ) X3, Y3, X4, Y4, the endpoints of the second
!    segment.
!
!    Output, integer ( kind = 4 ) FLAG, records the results.
!    0, the line segments do not intersect.
!    1, the line segments intersect.
!
!    Output, real ( kind = 8 ) X5, Y5.
!    If FLAG = 0, X5 = Y5 = 0.
!    If FLAG = 1, then (X5,Y5) is a point of intersection.
!
  implicit none

  integer ( kind = 4 ) flag
  integer ( kind = 4 ) ival
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ) x4
  real ( kind = 8 ) x5
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) y3
  real ( kind = 8 ) y4
  real ( kind = 8 ) y5

  x5 = 0.0
  y5 = 0.0
!
!  Find the intersection of the two lines.
!
  call lines_exp_int_2d ( ival, x5, y5, x1, y1, x2, y2, x3, y3, x4, y4 )

  if ( ival == 0 ) then
    flag = 0
    return
  end if
!
!  Is the point on the first segment?
!
  call line_seg_contains_point_2d ( x1, y1, x2, y2, x5, y5, u, v )

  if ( u < 0 .or. 1.0D+00 < u .or. 0.001D+00 < v ) then
    flag = 0
    return
  end if
!
!  Is the point on the second segment?
!
  call line_seg_contains_point_2d ( x3, y3, x4, y4, x5, y5, u, v )

  if ( u < 0 .or. 1.0D+00 < u .or. 0.001D+00 < v ) then
    flag = 0
    return
  end if

  flag = 1

  return
end

end module geom
