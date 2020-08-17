	program circularise
	implicit real*8 (a-h,o-z)
c
c reads in a DNA duplex made by nucgen, and
c makes it into a circle. The first and last
c residues in each strand are skipped, as the
c microcircle has no ends
c
	logical OK,getiarg,getrarg,findmove
	integer STDIN,STDOUT,STDERR
	real*4 twist

	parameter(MAXAT=100, pi=3.1415927, rise=3.6)
	parameter(STDIN=5,STDOUT=6,STDERR=0)
	dimension xyz(3,MAXAT),xyz2(3,MAXAT)
	dimension r(3,3),v(3)
	dimension xsug(3,8),xsugref(3,8)
	character*9 atomid(MAXAT)
	character*10 resid(MAXAT)

	OK=.true.
	findmove=.true.
	if(.not.getiarg('-n',nres)) then
	  write(STDERR,*) 'Specify total no. of residues with -n'
	  OK=.false.
	endif
	if(.not.getrarg('-t',twist)) then
	  write(STDERR,*) 'Specify target twist with -t'
	  OK=.false.
	endif
	if(.not.OK) call exit(1)

c
c adjust twist so that circle will be closed
c
	nbp=nres/2-2
	turns=real(nint(nbp*twist/360.0))
	twist2=turns*360.0/nbp

	write(STDERR,*) 'Actual twist will be ',twist2
	write(STDERR,*) 'Number of turns: ',turns
c
c calculate some parameters for the circle
c
	circum = nbp*rise
	rad = circum/(2*pi)
	theta = 360.0/nbp
	write(STDERR,*) 'Radius: ',rad
c
c initialise transformations
c
	tw=0.0
	th=0.0
c
c process first residue to obtain sugar coordinates
c
	call resread(STDIN,atomid,resid,xyz,nat)
	call sugget(atomid,xyz,xsugref,nat,nsugat)
c
c now process first strand
c
	do 10,ires=1,nbp
c
c read in the residue, fit to the reference position by
c overlaying the sugar atoms
c
	  call resread(STDIN,atomid,resid,xyz,nat)
	  call sugget(atomid,xyz,xsug,nat,nsugat)
	  call matfit(nsugat,xsugref,xsug,r,v,rmsd,findmove)
c	  if(rmsd.gt.0.5) stop 'suspicious sugar fit'
	  call transform(xyz,xyz2,r,v,nat)
c
c rotate about z to required twist, translate along x to required
c radius, then rotate about y to required theta
c
	  call makermat(r,3,tw)
	  v(1)=rad
	  v(2)=0.0
	  v(3)=0.0
	  call transform(xyz2,xyz,r,v,nat)
	  call makermat(r,2,th)
	  v(1)=0.0
	  call transform(xyz,xyz2,r,v,nat)
c
c  now write out residue with new coordinates
c
	  call reswrite(STDOUT,atomid,resid,xyz2,nat)
c
c next residue - update angles
c
	  tw=tw+twist2
	  th=th+theta
 10	continue
c
c end of strand - add TER record
c
	write(STDOUT,9001)
 9001	format('TER')
c
c now discard last residue of strand1 and first of strand2
c
	call resread(STDIN,atomid,resid,xyz,nat)
	call resread(STDIN,atomid,resid,xyz,nat)
	tw=tw-twist2
	th=th-theta
c
c adjust reference sugar coordinates for second strand
c
	do 33,i=1,8
	 do 34,j=2,3
	   xsugref(j,i)=-xsugref(j,i)
 34	 continue
 33	continue
c
c  now off we go again...
c
	do 20,ires=1,nbp
c
c read in the residue, fit to the reference position by
c overlaying the sugar atoms
c
	  call resread(STDIN,atomid,resid,xyz,nat)
	  call sugget(atomid,xyz,xsug,nat,nsugat)
	  call matfit(nsugat,xsugref,xsug,r,v,rmsd,findmove)
c	  if(rmsd.gt.0.5) stop 'suspicious sugar fit'
	  call transform(xyz,xyz2,r,v,nat)
c
c rotate about z to required twist, translate along x to required
c radius, then rotate about y to required theta
c
	  call makermat(r,3,tw)
	  v(1)=rad
	  v(2)=0.0
	  v(3)=0.0
	  call transform(xyz2,xyz,r,v,nat)
	  call makermat(r,2,th)
	  v(1)=0.0
	  call transform(xyz,xyz2,r,v,nat)
c
c  now write out residue with new coordinates
c
	  call reswrite(STDOUT,atomid,resid,xyz2,nat)
c
c next residue - update angles
c
	  tw=tw-twist2
	  th=th-theta
 20	continue
c
c end of strand - add TER record
c
	write(STDOUT,9001)
	end
