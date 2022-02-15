module readers
  implicit none 
  public :: oxdna_reader, amber_reader

contains

subroutine oxdna_reader(conf_file_name,top_file_name,step,bp,x1,y1,z1,x2,y2,z2,circular,energy_out,reverse)
  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2,x2_rev,y2_rev,z2_rev
  integer :: i,bp,nlines,io,step,bases,strands,strandid,circ_int,dum_int,length,start,fin
  logical :: energy_out,circular,reverse
  character :: conf_file_name*50,top_file_name*50, search_string*4, top_line*20, energy_out_string*4
  character :: base*1, second_line*20,bp_string*5
  real :: tol=0.01

  open(3,file=top_file_name)
  read(3,"(A20)") top_line
  read(top_line,*) bases,strands
  bp = int(bases/strands)
  
  write(*,"(i4,a21)") bp," base pairs detected"
  read(3,"(A20)") second_line
  read(second_line,"(i2,a1,i1,i3)") strandid,base,circ_int,dum_int
  if (circ_int.eq.-1) then
    circular=.False.
  else
    circular=.True.
    write(*,*) 'Circular DNA detected'
  end if
  close(3)

  search_string = 't = '
  open(1,file=conf_file_name)
  nlines=0

  ! searches for string pattern "t = " to set logical energy_out and obatin timestep
  read(1,"(A20)") top_line
  read(top_line,'(A4)') energy_out_string
  close(1)
  open(1,file=conf_file_name)
  if (energy_out_string==search_string) then
    read (1,"(A4,I8)") search_string,step
    energy_out=.True.
  else
    energy_out=.False.
  end if
  close(1)

  ! now we load the base pair positions of both strands
  open(1,file=conf_file_name) 
  if (circular.eqv..False.) then
    length = bp
  else
    length = bp+1
  end if


  allocate(x1(length))
  allocate(y1(length))
  allocate(z1(length))
  allocate(x2(length))
  allocate(y2(length))
  allocate(z2(length))
  if (energy_out.eqv..False.) then
    do i = 1,bp
      read(1,*) x1(i), y1(i), z1(i)
    end do
    do i = 1,bp
      read(1,*) x2(i), y2(i), z2(i)
    end do
  else
    do i = 1,3
      read(1,*,iostat=io)
    end do
    do i = 1,bp
      read(1,*,iostat=io) x1(i), y1(i), z1(i)
    end do
    do i = 1,bp
      read(1,*) x2(i), y2(i), z2(i)
    end do
  end if
  
  if (reverse.eqv..True.) then
    allocate(x2_rev(length))
    allocate(y2_rev(length))
    allocate(z2_rev(length))
    if (circular.eqv..True.) then
      do i=1,length-1 ! only loop to the penultimate element before reversal, as array is now one longer!
        x2_rev(i)=x2(length-i)
        y2_rev(i)=y2(length-i)
        z2_rev(i)=z2(length-i)
      end do
      x2_rev(length) = x2(length-1)
      y2_rev(length) = y2(length-1)
      z2_rev(length) = z2(length-1)
      x2(:)=x2_rev(:)
      y2(:)=y2_rev(:)
      z2(:)=z2_rev(:)
    else
      do i=1,length
        x2_rev(i)=x2(length-i)
        y2_rev(i)=y2(length-i)
        z2_rev(i)=z2(length-i)
      end do
      x2_rev(length) = x2_rev(1)
      y2_rev(length) = y2_rev(1)
      z2_rev(length) = z2_rev(1)
      x2(:)=x2_rev(:)
      y2(:)=y2_rev(:)
      z2(:)=z2_rev(:)
    end if
  else
    if (circular.eqv..True.) then
      x2(length)=x2(1)
      y2(length)=y2(1)
      z2(length)=z2(1)
    end if
  end if
    
        

  close(1)

end subroutine oxdna_reader


subroutine amber_reader(pdb_file_name,top_file_name,step,bp,x1,y1,z1,x2,y2,z2,reverse,circular)

  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2,x2_rev,y2_rev,z2_rev
  integer :: i,bp,nlines,io,step,bases,circ_int,dum_int,m,n,nres,natoms,k,s
  integer(kind = 4) :: bpn 
  logical :: circular,reverse
  character :: pdb_file_name*40,top_file_name*40,search_string*4,top_line*20,atom_name*4,a1*1,str*4
  character :: base*1, second_line*20,resname*10,atom*4,word*10,line*20,string*128,altloc,icode,chains,segid*4,charge*2,element*2
  real :: tol=0.01,x,y,z,occ,temp


  open(3,file=top_file_name)
  write(*,*) 'Reading amber topology'
  do i=1,6
    read(3,"(A20)") line
  end do
  read(3,"(I8)") natoms
  read(3,"(I8,I8)") k,nres
  
  bp = int(nres/2)
  write(*,"(i5,a21)") bp," base pairs detected"

  ! now we load the base pair positions of both strands
  open(1, file=pdb_file_name, iostat=io)
  write(*,*) 'Reading pdb'
  if (circular.eqv..False.) then
    allocate(x1(bp))
    allocate(y1(bp))
    allocate(z1(bp))
    allocate(x2(bp))
    allocate(y2(bp))
    allocate(z2(bp))
    s = 1
    do
      read (1,'(A)', iostat=io) string
      if (string=='END') then
        exit
      else if (string=="TER") then
        s = s + 1
      else
        read(string,"(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)",iostat=io) word,n,atom,altloc,resname,&
             &chains,bpn,icode,x,y,z,occ,temp,segid,element,charge
        if ( atom(1:1) == ' ' ) then
          atom = atom(2:)
        end if
        if (resname=="DT5") then
          str = "H05'"
        else
          str = "P"
        end if
        if (atom==str) then
          if(s.eq.1) then
            x1(bpn)=x
            y1(bpn)=y
            z1(bpn)=z
          else
            bpn = bpn - bp
            x2(bpn)=x
            y2(bpn)=y
            z2(bpn)=z
          end if
        end if
      end if
    end do

    if (reverse.eqv..True.) then
    allocate(x2_rev(bp))
    allocate(y2_rev(bp))
    allocate(z2_rev(bp))
      do i=1,bp
        x2_rev(i)=x2(bp+1-i)
        y2_rev(i)=y2(bp+1-i)
        z2_rev(i)=z2(bp+1-i)
      end do
    x2(:)=x2_rev(:)
    y2(:)=y2_rev(:)
    z2(:)=z2_rev(:)
    else
    end if

  else if (circular.eqv..True.) then
    allocate(x1(bp+1))
    allocate(y1(bp+1))
    allocate(z1(bp+1))
    allocate(x2(bp+1))
    allocate(y2(bp+1))
    allocate(z2(bp+1))
    
    s = 1
    do
      read (1,'(A)', iostat=io) string
      if (string=='END') then
        exit
      else if (string=="TER") then
        s = s + 1
      else
        read(string,"(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)",iostat=io) word,n,atom,altloc,resname,&
             chains,bpn,icode,x,y,z,occ,temp,segid,element,charge
      write(*,*) resname
        if ( atom(1:1) == ' ' ) then
          atom = atom(2:)
        end if
        if (atom=='P') then
          if(s.eq.1) then
            x1(bpn)=x
            y1(bpn)=y
            z1(bpn)=z
          else
            bpn = bpn - bp
            x2(bpn)=x
            y2(bpn)=y
            z2(bpn)=z
          end if
        end if
      end if
    end do

    x1(bp+1)=x1(1)
    y1(bp+1)=y1(1)
    z1(bp+1)=z1(1)
    x2(bp+1)=x2(1)
    y2(bp+1)=y2(1)
    z2(bp+1)=z2(1)
    
    if (reverse.eqv..True.) then
    allocate(x2_rev(bp+1))
    allocate(y2_rev(bp+1))
    allocate(z2_rev(bp+1))
      do i=1,bp+1
        x2_rev(i)=x2(bp+2-i)
        y2_rev(i)=y2(bp+2-i)
        z2_rev(i)=z2(bp+2-i)
      end do
    x2(:)=x2_rev(:)
    y2(:)=y2_rev(:)
    z2(:)=z2_rev(:)
    else
    end if
  end if
  close(1)
  
end subroutine amber_reader




subroutine oxdna_top_reader(top_file_name,bp,circular)

  ! retrieves the number of base pairs and the circular logical from an oxdna topology
  
  integer :: bp,bases,strands,circ_int,strandid,dum_int
  logical :: circular
  character :: top_line*30,second_line*30,top_file_name*30,base*1

  open(3,file=top_file_name)
  write(*,*) 'Reading oxdna topology'
  read(3,"(A20)") top_line
  read(top_line,*) bases,strands
  bp = int(bases/strands)

  read(3,"(A20)") second_line
  write(*,*) second_line
  read(second_line,"(i2,a1,i1,i3)") strandid,base,circ_int,dum_int
  if (circ_int.eq.-1) then
    circular=.False.
  else
    circular=.True.
  end if
  close(3)



end subroutine oxdna_top_reader

subroutine amber_top_reader(top_file_name,bp,circular)

  ! retrieves the number of base pairs and the circular logical from an amber topology

  integer :: k,nres,natoms,bp,circ_int,i
  logical :: circular
  character :: line*30,top_file_name*40

  open(3,file=top_file_name)
  write(*,*) 'Reading amber topology'
  do i=1,6
    read(3,"(A20)") line
  end do
  read(3,"(I8)") natoms
  read(3,"(I8,I8)") k,nres

  circ_int = 1
  if (circ_int.eq.1) then
    circular=.True.
  else
    circular=.False.
  end if
  bp = int(nres/2)
  write(*,"(i5,a21)") bp," base pairs detected"
  close(3)

end subroutine amber_top_reader


subroutine oxdna_block_reader(unit_num,step,bp,x1,y1,z1,x2,y2,z2,circular,energy_out,reverse)

  ! 
  ! This subroutine reads a block of the file (unit=unit_num) and loads the BDNA base-pair positions
  ! ESsentially a reworking of amber_reader but deisgned 

  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2,x2_rev,y2_rev,z2_rev
  integer :: i,bp,nlines,io,step,bases,strands,strandid,circ_int,dum_int,length,start,fin,unit_num
  logical :: energy_out,circular,reverse
  character :: search_string*4, top_line*20, energy_out_string*4
  character :: base*1, second_line*20,bp_string*5
  real :: tol=0.01

  search_string = 't = '
  nlines=0

  if (circular.eqv..False.) then
    length = bp
  else
    length = bp+1
  end if


  allocate(x1(length))
  allocate(y1(length))
  allocate(z1(length))
  allocate(x2(length))
  allocate(y2(length))
  allocate(z2(length))
  if (energy_out.eqv..False.) then
    do i = 1,bp
      read(unit_num,*) x1(i), y1(i), z1(i)
    end do
    do i = 1,bp
      read(unit_num,*) x2(i), y2(i), z2(i)
    end do
  else
    do i = 1,3
      read(unit_num,*,iostat=io)
    end do
    do i = 1,bp
      read(unit_num,*,iostat=io) x1(i), y1(i), z1(i)
    end do
    do i = 1,bp
      read(unit_num,*) x2(i), y2(i), z2(i)
    end do
  end if
  
  if (reverse.eqv..True.) then
    allocate(x2_rev(length))
    allocate(y2_rev(length))
    allocate(z2_rev(length))
    if (circular.eqv..True.) then
      do i=1,length-1 ! only loop to the penultimate element before reversal, as array is now one longer!
        x2_rev(i)=x2(length-i)
        y2_rev(i)=y2(length-i)
        z2_rev(i)=z2(length-i)
      end do
      x2_rev(length) = x2(length-1)
      y2_rev(length) = y2(length-1)
      z2_rev(length) = z2(length-1)
      x2(:)=x2_rev(:)
      y2(:)=y2_rev(:)
      z2(:)=z2_rev(:)
    else
      do i=1,length
        x2_rev(i)=x2(length-i)
        y2_rev(i)=y2(length-i)
        z2_rev(i)=z2(length-i)
      end do
      x2_rev(length) = x2_rev(1)
      y2_rev(length) = y2_rev(1)
      z2_rev(length) = z2_rev(1)
      x2(:)=x2_rev(:)
      y2(:)=y2_rev(:)
      z2(:)=z2_rev(:)
    end if
  else
    if (circular.eqv..True.) then
      x2(length)=x2(1)
      y2(length)=y2(1)
      z2(length)=z2(1)
    end if
  end if
        

end subroutine oxdna_block_reader


subroutine amber_block_reader(unit_num,step,bp,x1,y1,z1,x2,y2,z2,reverse,circular)
  
  ! 
  ! This subroutine reads a block of the file (unit=unit_num) and loads the BDNA base-pair positions
  ! ESsentially a reworking of amber_reader but deisgned to be used in a loop
  ! NOTE: x1,y1,z1,x2,y2,z2 are not-allocated here

  real, dimension(:), allocatable :: x1,y1,z1,x2,y2,z2,x_rev,y_rev,z_rev
  integer :: i,bp,nlines,io,step,bases,circ_int,dum_int,m,n,nres,natoms,k,s,unit_num,length
  integer(kind = 4) :: bpn 
  logical :: circular,reverse
  character :: search_string*4,top_line*20,atom_name*4,a1*1,str*4
  character :: base*1, second_line*20,resname*10,atom*4,word*10,line*20,string*128,altloc,icode,chains,segid*4,charge*2,element*2
  real :: tol=0.01,x,y,z,occ,temp

  if (circular.eqv..True.) then
    length = bp + 1
  else
    length = bp
  end if

  s = 1
  do
    read (unit_num,'(A)', iostat=io) string
    if (string(1:3)=='END') then
      exit
    else if (string(1:3)=="TER") then
      s = s + 1
    else if (string(1:4)=="ATOM") then
      read(string,"(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)",iostat=io) word,n,atom,altloc,resname,&
           &chains,bpn,icode,x,y,z,occ,temp,segid,element,charge
      if ( atom(1:1) == ' ' ) then
        atom = atom(2:)
      end if
      if (resname=="DT5") then
        str = "H05'"
      else
        str = "P"
      end if
      if (atom==str) then
        if(s.eq.1) then
          x1(bpn)=x
          y1(bpn)=y
          z1(bpn)=z
        else
          bpn = bpn - bp
          x2(bpn)=x
          y2(bpn)=y
          z2(bpn)=z
        end if
      end if
    else
    end if
  end do
  if (circular.eqv..True.) then
    x1(bp+1)=x1(1)
    y1(bp+1)=y1(1)
    z1(bp+1)=z1(1)
    x2(bp+1)=x2(1)
    y2(bp+1)=y2(1)
    z2(bp+1)=z2(1)
  end if
  
  if (reverse.eqv..True.) then
    allocate(x_rev(length))
    allocate(y_rev(length))
    allocate(z_rev(length))
    do i=1,length
      x_rev(i)=x1(length+1-i)
      y_rev(i)=y1(length+1-i)
      z_rev(i)=z1(length-i)
    end do
  x1(:)=x_rev(:)
  y1(:)=y_rev(:)
  z1(:)=z_rev(:)
  !deallocate(x2_rev(length))
  !deallocate(y2_rev(length))
  !deallocate(z2_rev(length))
  end if
  
end subroutine amber_block_reader

subroutine amber_init_reader(pdb_file_name,bp,circular)

  ! amber_init_reader obtains the logical circular via a search for terminal residues as well as the number of residues
  ! required for explicit solvent simulations where the number of base pairs is not obtainable from the topology

  character :: pdb_file_name*40,search_string*4,top_line*20,atom_name*4,a1*1,trunc_string*3
  integer :: i,bp,nlines,io,step,bases,circ_int,dum_int,m,n,nres,natoms,k,s
  integer(kind = 4) :: bpn
  logical :: circular
  character :: base*1, second_line*20,resname*10,atom*4,word*10,line*20,string*128,altloc,icode,chains,segid*4,charge*2,element*2
  
  ! strand assumed circular unless terminal residues can be found
  circular=.True.

  ! we must open the amber file to set the circular logical by looking for terminal residues
  open(1, file=pdb_file_name, iostat=io)
  s = 1
  do
    read (1,'(A)', iostat=io) string
    trunc_string = string(1:3)
    if (trunc_string=='END') then
      exit
    else if (trunc_string=="TER") then
      exit
    else
      read(string,"(a6,i5,1x,a4,a1,a3,1x,a1,i4)",iostat=io) word,n,atom,altloc,resname,chains,bpn
    end if

    if (resname=='DT5') then
      circular=.False.
    else if (resname=='DT3') then
      circular=.False.
    else if (resname=='DG5') then
      circular=.False.
    else if (resname=='DG3') then
      circular=.False.
    else if (resname=='DC5') then
      circular=.False.
    else if (resname=='DC3') then
      circular=.False.
    else if (resname=='DA5') then
      circular=.False.
    else if (resname=='DA3') then
      circular=.False.
    end if
  end do

  bp = bpn
  close(1)

end subroutine amber_init_reader

end module readers
