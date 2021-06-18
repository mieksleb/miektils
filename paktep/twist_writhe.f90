module twist_writhe_fortran
  use spline
  implicit none

contains

subroutine twist_writhe_xyz(bp,x1,y1,z1,x2,y2,z2,twist_integral,writhe_integral,circular)
  ! This takes in the 2 sets of vectors each of length bp (or bp+1 if circular)
  ! x1,y1,z1 are lists of the components of the vectors which describe the backbone of strand1
  ! x2,y2,z2 are the components of the the list of vectors which descrive the backbone of strand2
  ! reverse is a logical which is True if one strand is to be reversed and False otherwise

  real :: twist_integral,writhe_integral
  real, dimension(:), allocatable :: dum_x1,x1,y1,z1,x2,y2,z2
  real, dimension(:), allocatable   :: tx1,ty1,tz1,tx2,ty2,tz2,cx1,cy1,cz1,cx2,cy2,cz2
  integer :: bp,npoints=1000,nx1,ny1,nz1,nx2,ny2,nz2,k=3,ier,i,step,length
  logical :: circular,reverse,energy_out,reverse1=.False.
  character :: conf_file_name*20,top_file_name*20


  if (circular.eqv..False.) then
    length = bp
  else
    length = bp+1
  end if

  ! generate a linear sequence to feed as independent variable in spline fitting procedure
  allocate(dum_x1(length))
  do i=1,length
    dum_x1(i)=i-1
  end do

  ! calculate the splines of the x,y,z postions for both strands
  ! by default we retain the direction of one spline, but allow for the second one to be reversed subject to the reverse logical
  call get_spline(dum_x1,x1,tx1,cx1,k,nx1,length,circular,ier)
  call get_spline(dum_x1,y1,ty1,cy1,k,ny1,length,circular,ier)
  call get_spline(dum_x1,z1,tz1,cz1,k,nz1,length,circular,ier)
  call get_spline(dum_x1,x2,tx2,cx2,k,nx2,length,circular,ier)
  call get_spline(dum_x1,y2,ty2,cy2,k,ny2,length,circular,ier)
  call get_spline(dum_x1,z2,tz2,cz2,k,nz2,length,circular,ier)

  ! call main subroutine 
  call twist_writhe(bp,npoints,tx1,cx1,nx1,ty1,cy1,ny1,tz1,cz1,nz1,tx2,cx2,nx2,ty2,cy2,ny2,tz2,&
                                            &cz2,nz2,circular,twist_integral,writhe_integral)


end subroutine twist_writhe_xyz



subroutine twist_writhe(bp,npoints,tx1,cx1,nx1,ty1,cy1,ny1,tz1,cz1,nz1,tx2,cx2,nx2,ty2,cy2,ny2,tz2&
                                              &,cz2,nz2,circular,twist_integral,writhe_integral)

  real :: pi = 3.1415926535897932
  real, dimension(:), allocatable   :: tx1,ty1,tz1,tx2,ty2,tz2,cx1,cy1,cz1,cx2,cy2,cz2
  integer                           :: ii,jj,srange,ier,k=3,nx1,ny1,nz1,nx2,ny2,nz2,bp,nuxx,nuyy,nuzz
  integer                           :: npoints,m,nnuxx,nnuyy,nnuzz,nsx,nsy,nsz
  logical                           :: circular,reverse
  real                              :: twist_integral,writhe_integral
  real                              :: bpinc
  real                              :: tsp,tsp2,delta_s,ds
  real, dimension(:), allocatable   :: uxx,uyy,uzz
  real, dimension(:,:), allocatable :: uu,tt
  real                              :: dot_tu, norm
  real, dimension(:), allocatable   :: uxx_bpi,uyy_bpi,uzz_bpi
  real, dimension(:), allocatable   :: m1xx,m1yy,m1zz,dmxx,dmyy,dmzz,duxx,duyy,duzz,xx,yy,zz
  real, dimension(:,:), allocatable :: duu
  real, dimension(:), allocatable   :: diff,tuxx,tuyy,tuzz,cnuxx,cnuyy,cnuzz,csx,csy,csz,cuxx,cuyy,cuzz
  real, dimension(:), allocatable   :: tnuxx,tnuyy,tnuzz,tsx,tsy,tsz,ss,contour,bpi
  real, dimension(:), allocatable   :: sx1,sy1,sz1,sx2,sy2,sz2,snuxx,snuyy,snuzz,dum_x2


  m = npoints

  allocate(contour(npoints-1))
  allocate(ss(npoints))
  allocate(bpi(npoints))
  allocate(tt(3,npoints))
  allocate(diff(3))
  allocate(uu(3,npoints))
  allocate(sx1(npoints))
  allocate(sy1(npoints))
  allocate(sz1(npoints))
  allocate(sx2(npoints))
  allocate(sy2(npoints))
  allocate(sz2(npoints))
  
 ! calculate the base pair index (bpi) common to both strands
  if (circular.eqv..True.) then  
    bpinc = real(bp,4)/(real(npoints,4)-1)
  else 
    bpinc = (real(bp,4)-1)/(real(npoints,4)-1)
  end if

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
  m1xx = (sx1+sx2)/2
  m1yy = (sy1+sy2)/2
  m1zz = (sz1+sz2)/2

  ! now contruct contour_len, a sequence which passes through all base pairs of length npoints
  contour(1)=0
  do ii=2,npoints
    contour(ii) = contour(ii-1) + ((m1xx(ii)-m1xx(ii-1))**2+(m1yy(ii)-m1yy(ii-1))**2+(m1zz(ii)-m1zz(ii-1))**2)**0.5
  end do
  delta_s = contour(npoints)/(real(npoints,8)-1)


  ! from the average arclength parameter delta_s, we construct ss, the average arclength array
  ! this will be used as integration domain for both twist and writhe
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
  tt(1,:)=dmxx(:)
  tt(2,:)=dmyy(:)
  tt(3,:)=dmzz(:)


  ! we also need the vector u(s) which points between the base pairs
  uxx_bpi = sx2-sx1
  uyy_bpi = sy2-sy1
  uzz_bpi = sz2-sz1
  ! we then the spline of this normal and evaluate it
  call get_spline(contour,uxx_bpi,tuxx,cuxx,k,nuxx,npoints,circular,ier)
  call get_spline(contour,uyy_bpi,tuyy,cuyy,k,nuyy,npoints,circular,ier)
  call get_spline(contour,uzz_bpi,tuzz,cuzz,k,nuzz,npoints,circular,ier)
  allocate(uxx(npoints))
  allocate(uyy(npoints))
  allocate(uzz(npoints))
  call evaluate_spline(ss,uxx,tuxx,cuxx,k,nuxx,npoints,circular,ier,0)
  call evaluate_spline(ss,uyy,tuyy,cuyy,k,nuyy,npoints,circular,ier,0)
  call evaluate_spline(ss,uzz,tuzz,cuzz,k,nuzz,npoints,circular,ier,0)
  do ii=1,npoints
    uu(1,ii) = uxx(ii)
    uu(2,ii) = uyy(ii)
    uu(3,ii) = uzz(ii)
    dot_tu = tt(1,ii)*uu(1,ii)+tt(2,ii)*uu(2,ii)+tt(3,ii)*uu(3,ii)
    uu(1,ii) = uu(1,ii)-dot_tu*tt(1,ii)
    uu(2,ii) = uu(2,ii)-dot_tu*tt(2,ii)
    uu(3,ii) = uu(3,ii)-dot_tu*tt(3,ii)
    norm = ((uu(1,ii))**2+(uu(2,ii))**2+(uu(3,ii))**2)**0.5
    uu(:,ii)=uu(:,ii)/norm
  end do
  
  allocate(snuxx(npoints))
  allocate(snuyy(npoints))
  allocate(snuzz(npoints))
  snuxx(:)=uu(1,:)
  snuyy(:)=uu(2,:)
  snuzz(:)=uu(3,:)
  

  ! we now do a spline fit to the components of uu
  call get_spline(ss,snuxx,tnuxx,cnuxx,k,nnuxx,npoints,circular,ier)
  call get_spline(ss,snuyy,tnuyy,cnuyy,k,nnuyy,npoints,circular,ier)
  call get_spline(ss,snuzz,tnuzz,cnuzz,k,nnuzz,npoints,circular,ier)
  ! evalauate thje derivative of this new spline
  allocate(duxx(npoints))
  allocate(duyy(npoints))
  allocate(duzz(npoints))
  call evaluate_spline(ss,duxx,tnuxx,cnuxx,k,nnuxx,npoints,circular,ier,1)
  call evaluate_spline(ss,duyy,tnuyy,cnuyy,k,nnuyy,npoints,circular,ier,1)
  call evaluate_spline(ss,duzz,tnuzz,cnuzz,k,nnuzz,npoints,circular,ier,1)

  allocate(duu(3,npoints))
  duu(1,:)=duxx(:)
  duu(2,:)=duyy(:)
  duu(3,:)=duzz(:)

  ds = (contour(npoints)-contour(1))/(npoints-1)
  twist_integral = 0.0
  writhe_integral = 0.0
  if(circular) then
    srange = size(ss) - 1
  else
    srange = size(ss)
  end if
  do ii=1,srange
    tsp=tt(1,ii)*(uu(2,ii)*duu(3,ii)-uu(3,ii)*duu(2,ii)) + tt(2,ii)*(uu(3,ii)*duu(1,ii)-uu(1,ii)*duu(3,ii)) &
          &+ tt(3,ii)*(uu(1,ii)*duu(2,ii)-uu(2,ii)*duu(1,ii))
    twist_integral = twist_integral + tsp
    do jj = 1, srange
      if (ii > jj) then
        diff = (/ xx(ii)-xx(jj), yy(ii)-yy(jj), zz(ii)-zz(jj)/)
        diff = diff / (sqrt(sum(diff**2)))**3
        tsp2 = tt(1,ii)*(tt(2,jj)*diff(3)-tt(3,jj)*diff(2)) + tt(2,ii)*(tt(3,jj)*diff(1)-tt(1,jj)*diff(3)) &
              &+ tt(3,ii)*(tt(1,jj)*diff(2)-tt(2,jj)*diff(1))
        writhe_integral = writhe_integral + tsp2
      end if
    end do
  end do

  twist_integral = twist_integral * ds / (2 * pi)
  writhe_integral = writhe_integral * ds**2 / (2 * pi) 
  

end subroutine twist_writhe 

  real function triple_scalar_product(a, b, c)
    implicit none
    real, dimension(3) :: a, b, c
    triple_scalar_product = a(1)*(b(2)*c(3)-b(3)*c(2)) + a(2)*(b(3)*c(1)-b(1)*c(3)) + a(3)*(b(1)*c(2)-b(2)*c(1))
    return
  end function triple_scalar_product


end module twist_writhe_fortran


