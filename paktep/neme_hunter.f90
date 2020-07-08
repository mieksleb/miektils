module nemehunter
  use readers
  use spline
  use geom
  implicit none

contains 

subroutine neme_hunter_conf(conf_file_name,top_file_name,neme_pos,neme_len,step)
! callable subroutine for the plectoneme position of a single configuration
  integer :: neme_pos
  character :: conf_file_name*20, top_file_name*20


  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2,tx1,ty1,tz1,tx2,ty2,tz2,cx1,cy1,cz1,cx2,cy2,cz2
  real, dimension(:), allocatable :: dum_x1
  integer :: bp,npoints=1000,nx1,ny1,nz1,nx2,ny2,nz2,k=3,ier,i,step,array_len
  logical :: circular,reverse=.False.,energy_out
  real :: neme_len
  
  ! call reader to load in positions
  call reader(conf_file_name,top_file_name,step,bp,x1,y1,z1,x2,y2,z2,reverse.eqv..False.,circular,energy_out)

  ! generate a linear sequence to feed as independent variable in spline fitting procedure
 if (circular.eqv..True.) then
    array_len = bp+1
  else
    array_len = bp
  end if

  allocate(dum_x1(array_len))
   do i=1,array_len
    dum_x1(i) = i-1
  end do

  ! calculate the splines of the x,y,z postions for both strands
  call get_spline(dum_x1,x1,tx1,cx1,k,nx1,array_len,circular,reverse,ier)
  call get_spline(dum_x1,y1,ty1,cy1,k,ny1,array_len,circular,reverse,ier)
  call get_spline(dum_x1,z1,tz1,cz1,k,nz1,array_len,circular,reverse,ier)
  call get_spline(dum_x1,x2,tx2,cx2,k,nx2,array_len,circular,.True.,ier)
  call get_spline(dum_x1,y2,ty2,cy2,k,ny2,array_len,circular,.True.,ier)
  call get_spline(dum_x1,z2,tz2,cz2,k,nz2,array_len,circular,.True.,ier)
  
  call neme_hunter(bp,npoints,tx1,cx1,nx1,ty1,cy1,ny1,tz1,cz1,nz1,tx2,cx2,nx2,ty2,cy2,ny2,tz2&
                                              &,cz2,nz2,circular,neme_pos,neme_len)

end subroutine neme_hunter_conf

subroutine neme_hunter(bp,npoints,tx1,cx1,nx1,ty1,cy1,ny1,tz1,cz1,nz1,tx2,cx2,nx2,ty2,cy2,ny2,tz2&
                                              &,cz2,nz2,circular,neme_pos,neme_len)

  real, dimension(:), allocatable   :: tx1,ty1,tz1,tx2,ty2,tz2,cx1,cy1,cz1,cx2,cy2,cz2
  integer                           :: ii,ier,k=3,nx1,ny1,nz1,nx2,ny2,nz2,bp
  integer                           :: npoints,m,nsx,nsy,nsz,neme_pos,i_end_loop,j_end_loop,diff_min,diff_max,i,j
  integer                           :: i_plect_begin,j_plect_begin
  logical                           :: circular,reverse=.False.
  real                              :: bpinc,delta_s
  real, dimension(:), allocatable   :: m1xx,m1yy,m1zz,dmxx,dmyy,dmzz,xx,yy,zz
  real, dimension(:), allocatable   :: csx,csy,csz
  real, dimension(:), allocatable   :: tsx,tsy,tsz,ss,contour,bpi
  real, dimension(:), allocatable   :: sx1,sy1,sz1,sx2,sy2,sz2,dum_x2,curvature
  real :: x_spread,y_spread,z_spread, tol=0.01,neme_len
  logical :: x_proj,y_proj,z_proj
  integer :: flag, max_curv,coarse_factor=5,npoints_coarse,mid,neme_pos2,cutoff=10
  real(kind=8) :: x5,y5,alpha=0.1


  ! we begin with a much coarser spline fit as the self intersection algorith scales quadratically
  npoints_coarse = int(npoints/coarse_factor)
  m = npoints_coarse

  allocate(contour(npoints_coarse))
  allocate(ss(npoints_coarse))
  allocate(bpi(npoints_coarse))
  allocate(sx1(npoints_coarse))
  allocate(sy1(npoints_coarse))
  allocate(sz1(npoints_coarse))
  allocate(sx2(npoints_coarse))
  allocate(sy2(npoints_coarse))
  allocate(sz2(npoints_coarse))

! calculate the base pair index (bpi) common to both strands
  if (circular.eqv..True.) then
    bpinc = real(bp,4)/(real(npoints_coarse,4)-1)
  else
    bpinc = (real(bp,4)-1)/(real(npoints_coarse,4)-1)
  end if

  do ii=1,npoints_coarse
    bpi(ii)=(ii-1)*bpinc
  end do
  
 ! evaluate the splines in bpi, output is sij for i=x,y,z j=1,2
  call evaluate_spline(bpi,sx1,tx1,cx1,k,nx1,npoints_coarse,circular,ier,0)
  call evaluate_spline(bpi,sy1,ty1,cy1,k,ny1,npoints_coarse,circular,ier,0)
  call evaluate_spline(bpi,sz1,tz1,cz1,k,nz1,npoints_coarse,circular,ier,0)
  call evaluate_spline(bpi,sx2,tx2,cx2,k,nx2,npoints_coarse,circular,ier,0)
  call evaluate_spline(bpi,sy2,ty2,cy2,k,ny2,npoints_coarse,circular,ier,0)
  call evaluate_spline(bpi,sz2,tz2,cz2,k,nz2,npoints_coarse,circular,ier,0)
 ! calculate the midpoint splines in bpi
  m1xx = (sx1+sx2)/2
  m1yy = (sy1+sy2)/2
  m1zz = (sz1+sz2)/2

  ! now contruct contour_len, a sequence which passes through all base pairs of lenght npoints
  contour(1)=0
  do ii=2,npoints_coarse
    contour(ii) = contour(ii-1) + ((m1xx(ii)-m1xx(ii-1))**2+(m1yy(ii)-m1yy(ii-1))**2+(m1zz(ii)-m1zz(ii-1))**2)**0.5
  end do
  delta_s = contour(npoints_coarse)/(real(npoints_coarse,8)-1)

  ss(1)=0
  do ii=2,npoints_coarse
    ss(ii) = delta_s*(ii-1)
  end do

  allocate(dum_x2(npoints_coarse))
  do ii=1,npoints_coarse
    dum_x2(ii)=ii-1
  end do


  ! now we compute the spline objects t,c,k,n of the midpoint spline
  call get_spline(contour,m1xx,tsx,csx,k,nsx,npoints_coarse,circular,reverse,ier)
  call get_spline(contour,m1yy,tsy,csy,k,nsy,npoints_coarse,circular,reverse,ier)
  call get_spline(contour,m1zz,tsz,csz,k,nsz,npoints_coarse,circular,reverse,ier)

  allocate(xx(npoints_coarse))
  allocate(yy(npoints_coarse))
  allocate(zz(npoints_coarse))
  allocate(dmxx(npoints_coarse))
  allocate(dmyy(npoints_coarse))
  allocate(dmzz(npoints_coarse))

  ! we now evaluate this midpoint spline and its derivatives
  call evaluate_spline(ss,xx,tsx,csx,k,nsx,npoints_coarse,circular,ier,0)
  call evaluate_spline(ss,yy,tsy,csy,k,nsy,npoints_coarse,circular,ier,0)
  call evaluate_spline(ss,zz,tsz,csz,k,nsz,npoints_coarse,circular,ier,0)
  call evaluate_spline(ss,dmxx,tsx,csx,k,nsx,npoints_coarse,circular,ier,1)
  call evaluate_spline(ss,dmyy,tsy,csy,k,nsy,npoints_coarse,circular,ier,1)
  call evaluate_spline(ss,dmzz,tsz,csz,k,nsz,npoints_coarse,circular,ier,1)
  ! we now find which axis we have the most spread and choose to project out this dimension to obtain a 2D spline
  x_spread = maxval(xx)-minval(xx)
  y_spread = maxval(yy)-minval(yy)
  z_spread = maxval(zz)-minval(zz)

  if (x_spread.le.y_spread.and.x_spread.le.z_spread) then
    x_proj = .True.
    y_proj = .False.
    z_proj = .False.
  else if (y_spread.le.x_spread.and.y_spread.le.z_spread) then
    x_proj = .False. 
    y_proj = .True.
    z_proj = .False.
    
  else if(z_spread.le.x_spread.and.z_spread.le.y_spread) then
    x_proj = .False.
    y_proj = .False.
    z_proj = .True.
  end if

  ! now we project out the dimension with the smallest range
  ! we now call a routine from DUTCH which determines if our polygon has any self intersections (scales quadratically)
  diff_min=npoints_coarse  ! need a large upper bound for diff_min to detect end loop
  diff_max = 0            ! need lower bound for minimum diff_max to detect beginning and ending of plectoneme
  i_end_loop=0
  j_end_loop=0
  do i=1,npoints_coarse-1
    do j=1,npoints_coarse-1
      if (j>i+cutoff) then
        
        if (x_proj.eqv..True.) then
          call LINES_SEG_INT_2D(real(yy(i),8), real(zz(i),8), real(yy(i+1),8),real(zz(i+1),8),&
                      & real(yy(j),8),real(zz(j),8),real(yy(j+1),8),real(zz(j+1),8), flag,x5,y5)
          if (flag==1) then
!            write(*,*) "hoony"
!            write(*,*) i,j,bpi(npoints_coarse-i),bpi(npoints_coarse-j)
            if ((j-i).le.diff_min.and.(j-i).ge.cutoff) then ! if the differnece between i and j is smaller than the previous, then set i and j to be the new end loop indices
              diff_min=j-i
              i_end_loop = i
              j_end_loop = j
            end if
            if ((j-i).ge.diff_max.and.(j-i).ge.cutoff) then
              diff_max = j-i
              i_plect_begin = i
              j_plect_begin = j
            else
            end if
          else
          end if

        else if(y_proj.eqv..True.) then
          call LINES_SEG_INT_2D(real(xx(i),8), real(zz(i),8), real(xx(i+1),8),real(zz(i+1),8),&
                      & real(xx(j),8),real(zz(j),8),real(xx(j+1),8),real(zz(j+1),8), flag,x5,y5)
          if (flag==1) then
!            write(*,*) "hoony"
!            write(*,*) i,j,bpi(npoints_coarse-i),bpi(npoints_coarse-j)
            if ((j-i).le.diff_min.and.(j-i).ge.cutoff) then ! if the differnece between i and j is smaller than the previous, then set i and j to be the new end loop indices
              diff_min=j-i
              i_end_loop = i
              j_end_loop = j
            end if
            if ((j-i).ge.diff_max.and.(j-i).ge.cutoff) then
              diff_max = j-i
              i_plect_begin = i
              j_plect_begin = j
            end if
          else
          end if


        else if (z_proj.eqv..True.) then
          call LINES_SEG_INT_2D(real(xx(i),8), real(yy(i),8), real(xx(i+1),8),real(yy(i+1),8),&
                      & real(xx(j),8),real(yy(j),8),real(xx(j+1),8),real(yy(j+1),8), flag,x5,y5)
          if (flag==1) then
!            write(*,*) "hoony"
!            write(*,*) i,j,bpi(npoints_coarse-i),bpi(npoints_coarse-j)
            if ((j-i).le.diff_min.and.(j-i).ge.cutoff) then ! if the differnece between i and j is smaller than the previous, then set i and j to be the new end loop indices
              diff_min=j-i
              i_end_loop = i
              j_end_loop = j
            end if
            if ((j-i).ge.diff_max.and.(j-i).ge.cutoff) then
              diff_max = j-i
              i_plect_begin = i
              j_plect_begin = j
            end if
          else
          end if
        else
        end if
      else
      end if
    end do
  end do
   
!  we then do a curvature minimisation on the loop between the end loop indices to locate the plectoneme
!  write(*,*) "end loops"
!  write(*,*) i_end_loop,j_end_loop, bpi(i_end_loop),bpi(j_end_loop)
  if (i_end_loop.eq.0.and.j_end_loop.eq.0) then
    neme_pos=0
    neme_len=0
  else 
    allocate(curvature(npoints_coarse))
    mid = int((i_end_loop + j_end_loop)/2)
!    write(*,*) mid
    do i=1,npoints_coarse
      if (i.le.i_end_loop) then
        curvature(i) = 0
      else if (i.ge.j_end_loop) then
        curvature(i)=0
      else
        curvature(i) = (dmxx(i)**2+dmyy(i)**2+dmzz(i)**2)*exp(-alpha*abs(i-mid))
        curvature(1)=0
        curvature(npoints_coarse)=0
      end if
    end do
    i_end_loop= npoints_coarse-i_end_loop
    j_end_loop= npoints_coarse-j_end_loop
    neme_len = bpi(j_plect_begin)-bpi(i_plect_begin)
!    write(*,*) "end loops"
!    write(*,*) i_end_loop,j_end_loop, bpi(i_end_loop),bpi(j_end_loop)
    max_curv = maxloc(curvature,dim=1)
!    write(*,*) max_curv
    max_curv=npoints_coarse-max_curv
    neme_pos=int(bpi(max_curv))
    neme_pos2=int(bpi(mid))
!    write(*,*) neme_pos
  end if
end subroutine neme_hunter







end module nemehunter
