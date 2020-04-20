program twist
  implicit none

  real :: twist_integral = 0.0

  real :: pi = 3.1415926535897932
  integer                         :: ii, srange
  real, dimension(3,4)            :: tt = 2.0, uu =4.0, duu =7.0
  real, dimension(4)              :: ss, xx, yy, zz 
  logical                         :: circular = .True.
  real                            :: ds = 1.0
  real                                             :: tsp, tsp2
  real, dimension(3)                               :: diff

  call twist_writhe(tt, uu, duu, xx, yy, zz, ss, ds, circular, twist_integral, writhe_integral)

  write(*,*) twist_integral, writhe_integral


contains

subroutine twist_writhe(npoints, tt, uu, duu, xx, yy, zz, ss, ds, circular, twist_integral, writhe_integral)

  implicit none 

  real :: pi = 3.1415926535897932
  integer                          :: ii, jj, srange
  integer              :: npoints
  real, dimension(3, npoints) :: tt, uu, duu
  real, dimension(npoints)   :: ss, xx, yy, zz 
  logical             :: circular
  real              :: twist_integral, writhe_integral
  real                 :: ds
  real                             :: tsp, tsp2, triple_scalar_product
  real, dimension(3)               :: diff

  !allocate(tt(3, npoints))
  !allocate(uu(3, npoints))
  !allocate(duu(3, npoints))
  !allocate(xx(3, npoints))
  !allocate(yy(3, npoints))
  !allocate(zz(3, npoints))
  !allocate(ss(npoints))

  if (circular) then
    srange = size(ss) - 1
  else
    srange = size(ss)
  end if
 
  do ii = 1, srange
    tsp = tt(1,ii)*(uu(2,ii)*duu(3,ii)-uu(3,ii)*duu(2,ii)) + tt(2,ii)*(uu(3,ii)*duu(1,ii)-uu(1,ii)*duu(3,ii)) &
          + tt(3,ii)*(uu(1,ii)*duu(2,ii)-uu(2,ii)*duu(1,ii))
    twist_integral = twist_integral + tsp
    do jj = 1, srange
      if (ii > jj) then
        diff = (/ xx(ii)-xx(jj), yy(ii)-yy(jj), zz(ii)-zz(jj)/)
        diff = diff / (sqrt(sum(diff**2)))**3
        tsp2 = tt(1,ii)*(tt(2,jj)*diff(3)-tt(3,jj)*diff(2)) + tt(2,ii)*(tt(3,jj)*diff(1)-tt(1,jj)*diff(3)) &
              + tt(3,ii)*(tt(1,jj)*diff(2)-tt(2,jj)*diff(1))
        writhe_integral = writhe_integral + tsp2
      end if
    end do
  end do

  twist_integral = twist_integral * ds / (2 * pi)
  writhe_integral = writhe_integral * ds**2 / (2 * pi) 
  
  return 

end subroutine twist_writhe 

  real function triple_scalar_product(a, b, c)
    implicit none
    real, dimension(3) :: a, b, c
    triple_scalar_product = a(1)*(b(2)*c(3)-b(3)*c(2)) + a(2)*(b(3)*c(1)-b(1)*c(3)) + a(3)*(b(1)*c(2)-b(2)*c(1))
    return
  end function triple_scalar_product


end program twist_writhe_fortran


