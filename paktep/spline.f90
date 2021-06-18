module spline
  implicit none
contains

! input m data points (x,y) and returns spline object t,c,k,n. If data is periodic then set circular=.True.
subroutine get_spline(x,y,t,c,k,n,m,circular,ier)

  real :: xb,xe,s=0,fp
  integer :: iopt=0,m,k,nest,n,ier,i
  integer :: lwrk,array_len
  logical :: circular,reverse
  real, dimension(:), allocatable :: x,y,w,y_rev
  real, dimension(:), allocatable :: t,c
  real, dimension(:), allocatable :: wrk
  integer, dimension(:), allocatable :: iwrk

  ! difference in size of arrays for periodic splines, importantly the arrays themselves are 
  ! one component longer such that x(1)=x(m)
  if(circular.eqv..False.) then
    nest=m+k+1
    lwrk=m*(k+1)+nest*(7+3*k)
  else 
    nest=m+2*k
    lwrk=m*(k+1)+nest*(8+5*k)
  end if
  n=nest
  allocate(w(m))
  allocate(t(nest))
  allocate(c(nest))
  allocate(wrk(lwrk))
  allocate(iwrk(nest))
  allocate(y_rev(m))

!  REVERSE LOGICAL REMOVED, for now I will comment out redundant lines until I am sure it can be deleted
 
  ! if a strand is reversed, we must reverse the strand first and then ensure x(1)=x(m)
  ! doing to the other way round will give a different answer
!  if (circular.eqv..True.) then
!    if (reverse.eqv..True.) then
!      do i=1,m-1 ! only loop to the penultimate element before reversal, as array is now one longer!
!        y_rev(i)=y(m-i)
!      end do
!    else
!    end if
!    y_rev(m)=y(1) ! now the array is reversed, we insert the first element at the end
!    y(m)=y_rev(1) ! for unreversed array, we do the same
!  else
!    if (reverse.eqv..True.) then
!      do i=1,m ! now strand is not circular, we can simply reverse it!
!        y_rev(i)=y(m+1-i)
!      end do
!    else
!    end if
!  end if
      
  ! use unit weights uniformly spread
  do i=1,m
    w(i)=1
  end do
  
  xb=x(1)
  xe=x(m)
  

  ! now we call the main routines from CURFIT found in dierckx library
!  if (circular.eqv..False.) then
!    if (reverse.eqv..True.) then
!      call curfit(iopt,m,x,y_rev,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
!    else
!      call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
!    end if
!  else
!    if (reverse.eqv..True.) then
!      call percur(iopt,m,x,y_rev,w,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
!    else 
!      call percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
!    end if

  if (circular.eqv..False.) then
      call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
  else
      call percur(iopt,m,x,y,w,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
  end if

  if (ier.eq.10) then
     print *,"something went wrong, spline conditions not met"
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

