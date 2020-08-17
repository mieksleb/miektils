	logical function getiarg(flag,ivalue)
c
c  returns the integer argument specified by flag
c
	character*(*) flag
	character*80 arg

	l=lnblnk(flag)

	getiarg=.false.
	i=0
 	do 10,j=1,iargc()-1
	  call getarg(j,arg)
	  if(arg.eq.flag) i=j
 10	continue
	if(i.gt.0) then
	  call getarg(i+1,arg)
	  read(arg,*) ivalue
	  getiarg=.true.
	endif
	return
	end

	logical function getrarg(flag,rvalue)
c
c  returns the real argument specified by flag
c
	character*(*) flag
	character*80 arg

	l=lnblnk(flag)

	getrarg=.false.
	i=0
 	do 10,j=1,iargc()-1
	  call getarg(j,arg)
	  if(arg.eq.flag) i=j
 10	continue
	if(i.gt.0) then
	  call getarg(i+1,arg)
	  read(arg,*) rvalue
	  getrarg=.true.
	endif
	return
	end

	logical function getcarg(flag,cvalue)
c
c  returns the character argument specified by flag
c
	character*(*) flag,cvalue
	character*80 arg

	l=lnblnk(flag)

	getcarg=.false.
	i=0
 	do 10,j=1,iargc()-1
	  call getarg(j,arg)
	  if(arg.eq.flag) i=j
 10	continue
	if(i.gt.0) then
	  call getarg(i+1,cvalue)
	  getcarg=.true.
	endif
	return
	end

	logical function getnarg(flag)
c
c  returns true if flag is present
c
	character*(*) flag
	character*80 arg

	getnarg=.false.
	i=0
 	do 10,j=1,iargc()
	  call getarg(j,arg)
	  if(arg.eq.flag) i=j
 10	continue
	if(i.gt.0) then
	  getnarg=.true.
	endif
	return
	end
