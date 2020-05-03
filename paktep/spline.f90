module spline
  implicit none


contains


! input m data points (x,y) and returns spline object t,c,k,n. If data is periodic then set circular=.True.
subroutine get_spline(x,y,t,c,k,n,m,circular,ier)

  real :: xb,xe,s=0,fp
  integer :: iopt=0,m,k,nest,n,ier,i 
  integer :: lwrk
  logical :: circular
  real, dimension(:), allocatable :: x,y,w
  real, dimension(:), allocatable :: t,c
  real, dimension(:), allocatable :: wrk
  integer, dimension(:), allocatable :: iwrk
  nest=m+k+1
  n=nest
  lwrk=m*(k+1)+nest*(7+3*k)
  allocate(w(m))
  allocate(t(nest))
  allocate(c(nest))
  allocate(wrk(lwrk))
  allocate(iwrk(nest))
  
  do i=1,m
    w(i)=1
  end do
  
  xb=x(1)
  xe=x(m)

  if(circular.eqv..False.) then 
    call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
  else 
    call percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
    write(*,*) "ya boi a coircle"
  end if

  if (ier.eq.10) then
     print *,"something went wrong"
  else
  end if

end subroutine get_spline



! evaluate a spline object t,c,k,n at positions x and get values y. m is length of evalutation points.
subroutine evaluate_spline(x,y,t,c,k,n,m,circular,ier,der)

  integer :: k,n,m,ier,der
  real, dimension(:), allocatable :: x,y
  real, dimension(:), allocatable :: t,c
  real, dimension(:), allocatable :: wrk
  logical :: circular
!  allocate(y(m))
 

  if(der.eq.0) then
    call splev(t,n,c,k,x,y,m,ier)
  else
    allocate(wrk(n))
    call splder(t,n,c,k,der,x,y,m,wrk,ier)
  end if

end subroutine evaluate_spline


end module spline

