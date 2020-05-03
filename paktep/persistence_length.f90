module persistence_length
  use :: spline
  implicit none

contains

subroutine tangent_correlation(bp,npoints,tx1,cx1,nx1,ty1,cy1,ny1,tz1,cz1,nz1,tx2,cx2,nx2,ty2,cy2,ny2,tz2&
                                              &,cz2,nz2,circular,corr,ss)

  real, dimension(:), allocatable   :: tx1,ty1,tz1,tx2,ty2,tz2,cx1,cy1,cz1,cx2,cy2,cz2,corr
  integer                           :: ii,ier,k=3,nx1,ny1,nz1,nx2,ny2,nz2,bp
  integer                           :: npoints,m,nsx,nsy,nsz
  logical                           :: circular
  real                              :: bpinc, norm, norm_0
  real                              :: delta_s
  real, dimension(:), allocatable   :: m1xx,m1yy,m1zz,dmxx,dmyy,dmzz,xx,yy,zz
  real, dimension(:), allocatable   :: csx,csy,csz
  real, dimension(:), allocatable   :: tsx,tsy,tsz,ss,contour,bpi
  real, dimension(:), allocatable   :: sx1,sy1,sz1,sx2,sy2,sz2,dum_x2

  m = npoints

  allocate(contour(npoints))
  allocate(ss(npoints))
  allocate(bpi(npoints))
  allocate(sx1(npoints))
  allocate(sy1(npoints))
  allocate(sz1(npoints))
  allocate(sx2(npoints))
  allocate(sy2(npoints))
  allocate(sz2(npoints))
 
 ! calculate the base pair index (bpi) common to both strands
  bpinc = (real(bp,8)-1)/(real(npoints,8)-1)
  do ii=1,npoints
    bpi(ii)=(ii-1)*bpinc
  end do

 ! evaluate the splines in bpi, output is sij for i=x,y,z j=1,2
  call evaluate_spline(bpi,sx1,tx1,cx1,k,nx1,npoints,circular,ier,0)
  call evaluate_spline(bpi,sy1,ty1,cy1,k,ny1,npoints,circular,ier,0)
  call evaluate_spline(bpi,sz1,tz1,cz1,k,nz1,npoints,circular,ier,0)
  call evaluate_spline(bpi,sx2,tx2,cx2,k,nx2,npoints,circular,ier,0)
  call evaluate_spline(bpi,sy2,ty2,cy2,k,ny2,npoints,circular,ier,0)
  call evaluate_spline(bpi,sz2,tz2,cz2,k,nz2,npoints,circular,ier,0)
 ! calculate the midpoint splines in bpi
 ! calculate the midpoint splines in bpi
  m1xx = (sx1+sx2)/2
  m1yy = (sy1+sy2)/2
  m1zz = (sz1+sz2)/2

  ! now contruct contour_len, a sequence which passes through all base pairs of lenght npoints
  contour(1)=0
  do ii=2,npoints
    contour(ii) = contour(ii-1) + ((m1xx(ii)-m1xx(ii-1))**2+(m1yy(ii)-m1yy(ii-1))**2+(m1zz(ii)-m1zz(ii-1))**2)**0.5
  end do
  delta_s = contour(npoints)/(real(npoints,8)-1)

  ss(1)=0
  do ii=2,npoints
    ss(ii) = delta_s*(ii-1)
  end do

  allocate(dum_x2(npoints))
  do ii=1,npoints
    dum_x2(ii)=ii-1
  end do


  ! now we compute the spline objects t,c,k,n of the midpoint spline 
  call get_spline(contour,m1xx,tsx,csx,k,nsx,npoints,circular,ier)
  call get_spline(contour,m1yy,tsy,csy,k,nsy,npoints,circular,ier)
  call get_spline(contour,m1zz,tsz,csz,k,nsz,npoints,circular,ier)

  allocate(xx(npoints))
  allocate(yy(npoints))
  allocate(zz(npoints))
  allocate(dmxx(npoints))
  allocate(dmyy(npoints))
  allocate(dmzz(npoints))
  ! we now evaluate this midpoint spline and its derivatives
  call evaluate_spline(ss,xx,tsx,csx,k,nsx,npoints,circular,ier,0)
  call evaluate_spline(ss,yy,tsy,csy,k,nsy,npoints,circular,ier,0)
  call evaluate_spline(ss,zz,tsz,csz,k,nsz,npoints,circular,ier,0)
  call evaluate_spline(ss,dmxx,tsx,csx,k,nsx,npoints,circular,ier,1)
  call evaluate_spline(ss,dmyy,tsy,csy,k,nsy,npoints,circular,ier,1)
  call evaluate_spline(ss,dmzz,tsz,csz,k,nsz,npoints,circular,ier,1)

  ! the autocorrelation of the tangent vectors will be stored in the array corr(i)
  norm_0 = (dmxx(1)**2+dmyy(1)**2+dmzz(1)**2)**0.5
  do ii=1,npoints
    norm = (dmxx(ii)**2+dmyy(ii)**2+dmzz(ii)**2)**0.5
    corr(ii)=(dmxx(1)*dmxx(ii)+dmyy(1)*dmyy(ii)+dmzz(1)*dmzz(ii))/(norm_0*norm)
  end do
 



end subroutine tangent_correlation

! for n data points (x,y), find the exponential fit of the form y(x)=exp(a*x)
real function exp_fit(x,y,n) result(a)
  real, dimension(:), allocatable, intent(in) :: x,y
  integer, intent(in) :: n
  integer :: ii
  real :: a1,a2,sum_y=0,sum_x2y=0,sum_xylny=0,sum_ylny=0,sum_xy=0
  real :: sum_x=0,sum_xlny=0,sum_x2=0,sum_lny=0,diff1=0,diff2=0
  
  do ii=1,n
    sum_y = sum_y + y(ii)
    sum_x2y = sum_x2y + ((x(ii))**2)*y(ii)
    sum_xylny = sum_xylny + x(ii)*y(ii)*log(y(ii))
    sum_ylny = sum_ylny + y(ii)*log(y(ii))
    sum_xy = sum_xy + x(ii)*y(ii)

    sum_x = sum_x + x(ii)
    sum_xlny = sum_xlny + x(ii)*log(y(ii))
    sum_x2 = sum_x2 + x(ii)*x(ii)
    sum_lny = sum_lny + log(y(ii))
  end do
  
  ! calculates a from two different least square fits
  a1 = (n*sum_xlny-sum_x*sum_lny)/(n*sum_x2-(sum_x)**2)
  a2 = (sum_y*sum_xylny-sum_xy*sum_ylny)/((sum_y*sum_x2y)-(sum_xy)**2)

  ! check which gives better fit
  do ii=1,n
    diff1 = diff1 + (exp(a1*x(ii))-y(ii))**2
    diff2 = diff2 + (exp(a2*x(ii))-y(ii))**2
  end do

  
  if (diff1.le.diff2) then
    a = a1
  else
    a = a2
  end if

  return

end function exp_fit




end module persistence_length
