	subroutine resread(iunit,atomid,resid,xyz,nat)
c
c  returns the next residue read from unit iunit 
c 
	character*80 line 
	character*9 atomid(*) 
	character*10 resid(*) 
	real*8 xyz(3,*)

	nat=0
 10	continue
	  read(iunit,1001,end=100) line
 1001	  format(a80)
	  if(line(1:4).ne.'ATOM') goto 100
	  nat=nat+1
	  read(line,1002) atomid(nat),resid(nat),(xyz(j,nat),j=1,3)
 1002	  format(7x,a9,a10,4x,3f8.3)
	  if(nat.gt.1) then
	    if(resid(nat).ne.resid(1)) then
	      nat=nat-1
	      backspace(iunit)
	      goto 100
	    endif
	  endif
	  goto 10
 100	continue
	return
	end

	subroutine reswrite(iunit,atomid,resid,xyz,nat)
c
c  write out the residue in PDB format
c
	character*9 atomid(*)
	character*10 resid(*)
	real*8 xyz(3,*)

	do 10,iat=1,nat
	  write(iunit,1001) atomid(iat),resid(iat),(xyz(j,iat),j=1,3)
 1001	  format('ATOM   ',a9,a10,4x,3f8.3)
 10	continue
	end

	subroutine sugget(atomid,xyz,xsug,nat,is)
c
c  identify sugar coordinates in xyz and transfer to xsug
c
	character*9 atomid(*)
	real*8 xyz(3,*),xsug(3,*)

	is=0
	do 10,iat=1,nat
	  if(atomid(iat)(9:9).eq.'''') then
	    is=is+1
	    do 20,j=1,3
	      xsug(j,is)=xyz(j,iat)
 20	    continue
	  endif
 10	continue
	if(is.eq.0) stop 'no sugar atoms!'
	return
	end

	subroutine makermat(r,i,t)
c
c  creates rotation matrix for rotation angle t about axis i
c  x=1, y=2, z=3
c
	real*8 r(3,3),st,ct,t

	st=sind(t)
	ct=cosd(t)

	if(i.eq.1) then
	  r(1,1)=1.0
	  r(1,2)=0.0
	  r(1,3)=0.0
	  r(2,1)=0.0
	  r(2,2)=ct
	  r(2,3)=st
	  r(3,1)=0.0
	  r(3,2)=-st
	  r(3,3)=ct
	else if(i.eq.2) then
	  r(1,1)=ct
	  r(1,2)=0.0
	  r(1,3)=-st
	  r(2,1)=0.0
	  r(2,2)=1.0
	  r(2,3)=0.0
	  r(3,1)=st
	  r(3,2)=0.0
	  r(3,3)=ct
	else
	  r(1,1)=ct
	  r(1,2)=st
	  r(1,3)=0.0
	  r(2,1)=-st
	  r(2,2)=ct
	  r(2,3)=0.0
	  r(3,1)=0.0
	  r(3,2)=0.0
	  r(3,3)=1.0
	endif
	return
	end

	subroutine transform(xyz1,xyz2,r,v,ntot)
c
c  applies the rotation matrix r and vector v to the coordinates in xyz1
c  and returns result in xyz2
c
	real*8 xyz1(3,ntot),xyz2(3,ntot),r(3,3),v(3)

	do 10,i=1,ntot
	  do 20,j=1,3
	    xyz2(j,i)=v(j)
	    do 30,k=1,3
	      xyz2(j,i)=xyz2(j,i)+r(j,k)*xyz1(k,i)
 30	    continue
 20	  continue
 10	continue
	return
	end

      SUBROUTINE MATFIT (N, XA, XB, R, V, RMSE, ENTRY )
	implicit real*8 (a-h,o-z)
      DIMENSION   XA(3,N), XB(3,N), R(3,3), V(3)
      LOGICAL ENTRY
C
C     SUBROUTINE TO FIT THE COORD SET XA(3,N) TO THE SET XB(3,N)
C     IN THE SENSE OF:
C            XA= R*XB +V
C     R IS A UNITARY 3.3 RIGHT HANDED ROTATION MATRIX
C     AND V IS THE OFFSET VECTOR. THIS IS AN EXACT SOLUTION
C
C     IF ENTRY IS LOGICALLY FALSE ONLY THE RMS COORDINATE ERROR
C     WILL BE RETURNED (BUT QUICKLY)
C
C    THIS SUBROUTINE IS A COMBINATION OF MCLACHLAN'S AND KABSCH'S
C    TECHNIQUES. SEE
C     KABSCH, W. ACTA CRYST A34, 827,1978
C     MCLACHAN, A.D., J. MOL. BIO. NNN, NNNN 1978
C     WRITTEN BY S.J. REMINGTON 11/78.
C
C     THIS SUBROUTINE USES THE IBM SSP EIGENVALUE ROUTINE 'EIGEN'
C
      DIMENSION CMA(3),CMB(3),UMAT(3,3)
      XN=N
      XASQ=0.0
      XBSQ=0.0
      XNI=1.0/XN
C
C     ACCUMULATE UNCORRECTED (FOR C.M.) SUMS AND SQUARES
C
      DO 40 I=1,3
      CMA(I) = 0.0
      CMB(I) = 0.0
      DO 10 J=1,3
   10 UMAT(I,J) = 0.0
C
      DO 30 J=1,N
      DO 20 K=1,3
   20 UMAT(I,K) = UMAT(I,K) + XA(I,J)*XB(K,J)
C
      T = XA(I,J)
      CMA(I) = CMA(I) + T
      XASQ = XASQ + T*T
      T = XB(I,J)
      CMB(I) = CMB(I) + T
      XBSQ = XBSQ + T*T
   30 CONTINUE
   40 CONTINUE
C
C     SUBTRACT CM OFFSETS
C
      DO 50 I=1,3
      XASQ = XASQ - CMA(I)*CMA(I)*XNI
      XBSQ = XBSQ - CMB(I)*CMB(I)*XNI
      DO 50 J=1,3
      UMAT(I,J) = (UMAT(I,J) - CMA(I)*CMB(J)*XNI)*XNI
   50 CONTINUE
C
C     FIT IT
C
      CALL QKFIT( UMAT, RTSUM, R, ENTRY )
      RMSE =(XASQ + XBSQ)*XNI - 2.0*RTSUM
      IF( RMSE .LT. 0.0 ) RMSE = 0.0
      RMSE = SQRT (RMSE)
C
C      CALCULATE OFFSET IF ENTRY=.TRUE.
C
      IF(.NOT. ENTRY) RETURN
C
      DO 70 I=1,3
      T = 0.0
      DO 60 J=1,3
   60 T = T + R(I,J)*CMB(J)
      V(I) = ( CMA(I) - T )*XNI
   70 CONTINUE
      RETURN
      END
      SUBROUTINE QKFIT (UMAT, RTSUM, R, ENTRY )
	IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION UMAT(3,3), R(3,3)
      LOGICAL ENTRY
C
C     THE 'EIGENVALUE ONLY' ENTRY WAS
C     ADAPTED FROM PROGRAM BY A.D. MCLACHAN 7/78
C
      DIMENSION USQMAT(3,3),ROOT(3),A(3,3),B(3,3),UTR(6)
C
      EQUIVALENCE (AAM,USQMAT(1,1)),(BAM,USQMAT(2,2)),
     * (CAM,USQMAT(3,3))
      EQUIVALENCE (FAM,USQMAT(2,3)),(GAM,USQMAT(1,3))
     *,(HAM,USQMAT(1,2))
      EQUIVALENCE( A(1), USQMAT(1) ), (UTR(1), B(1))
C
      DATA EPS/1.0D-12/
      DATA PI/3.14159265358979/
      DATA ONE/1.0/,ZERO/0.0/,TWO/2.0/,THREE/3.0/
      DATA HALF/0.5/,THIRD/0.333333333/,FORTHR/1.333333333/
      DATA USQMAT(2,1),USQMAT(3,1),USQMAT(3,2)/3*0.0/
      ISIG=1
C
C      IF ENTRY IS .TRUE. GET OUT THE ROTATION MATRIX
C
      IF (ENTRY) GO TO 200
C
C     CALC DET OF UMAT
C
      DU11=(UMAT(2,2)*UMAT(3,3))-(UMAT(2,3)*UMAT(3,2))
      DU21=(UMAT(2,3)*UMAT(3,1))-(UMAT(2,1)*UMAT(3,3))
      DU31=(UMAT(2,1)*UMAT(3,2))-(UMAT(2,2)*UMAT(3,1))
      DETU=(UMAT(1,1)*DU11)+(UMAT(1,2)*DU21)+(UMAT(1,3)*DU31)
C%    WRITE(6,999) DETU
C%999 FORMAT(/(3F12.5))
      IF(DETU.LT. 0.0) ISIG=-1
C
C     FORM USQMAT AS POSITIVE SEMI DEFINITE MATRIX
C
      DO 110 J=1,3
      DO 105 I=1,J
      USQMAT(I,J)=(UMAT(1,I)*UMAT(1,J))+(UMAT(2,I)*UMAT(2,J))+
     *   (UMAT(3,I)*UMAT(3,J))
105   CONTINUE
110   CONTINUE
C%    WRITE(6,999) USQMAT
C
C     REDUCE AVG OF DIAGONAL TERMS TO ZERO
C
      DIGAV=(AAM+BAM+CAM)*THIRD
C%    WRITE(6,999) DIGAV
      AAM=AAM-DIGAV
      BAM=BAM-DIGAV
      CAM=CAM-DIGAV
C
C     SETUP COEFFS OF SECULAR EQUATION OF MATRIX WITH TRACE ZERO
C
      CC=(FAM*FAM)+(GAM*GAM)+(HAM*HAM)-(AAM*BAM)-(BAM*CAM)-(CAM*AAM)
      DD=(AAM*BAM*CAM)+TWO*(FAM*GAM*HAM)-(AAM*FAM*FAM)
     *   -(BAM*GAM*GAM)-(CAM*HAM*HAM)
C
C     THE SECULAR EQN IS Y**3-CC*Y-DD=0  AND DD IS DET(USQMAT)
C     REDUCE THIS TO THE FORM COS**3-(3/4)COS-
C     (1/4)COS3THETA = 0
C     WITH SOLUTIONS COSTHETA.  SO Y=QQ*COSTHETA
C
      IF(CC.LE.EPS) GO TO 115
      QQ=SQRT(FORTHR*CC)
      COS3TH=(THREE*DD)/(CC*QQ)
      IF(ABS(COS3TH).GT.ONE) COS3TH=SIGN(ONE,COS3TH)
C
C     FUNCTION ARCOS
C
      IF(COS3TH.NE.0.0) GOTO 1200
1100   THETA= 1.570796327
      GOTO 1400
1200  ARGSQ=COS3TH*COS3TH
      THETA=ATAN(SQRT(1.0-ARGSQ)/COS3TH)
      IF(COS3TH.LT.0.0) THETA=PI-THETA
1400  CONTINUE
C
C     ROOTS IN ORDER OF SIZE GO 1,2,3 1 LARGEST
C
      THETA=THETA*THIRD
      ROOT(1)=QQ*COS(THETA)
      DIFF=HALF*SQRT(THREE*(QQ*QQ-ROOT(1)*ROOT(1)))
      ROOT(2)=-ROOT(1)*HALF+DIFF
      ROOT(3)=-ROOT(1)*HALF-DIFF
      GO TO 120
115   CONTINUE
C
C     SPECIAL FOR TRIPLY DEGENERATE
C
      ROOT(1)=0.0
      ROOT(2)=0.0
      ROOT(3)=0.0
120   CONTINUE
C     ADD ON DIGAV AND TAKE SQRT
      DO 125 J=1,3
      RT=ROOT(J)+DIGAV
      IF(RT.LT.EPS) RT=0.0
      ROOT(J)=SQRT(RT)
125   CONTINUE
C%    WRITE(6,999) ROOT
C     IF DETU IS <0 CHANGE SIGN OF ROOT(3)
      IF(ISIG.EQ.-1) ROOT(3)=-ROOT(3)
      RTSUM=ROOT(1)+ROOT(2)+ROOT(3)
C%    WRITE(6,999) RTSUM
      RETURN
C
C     THIS IS THE FANCY PART
C
200   CONTINUE
C
C     FORM USQ = (UT).U    (IN UPPER TRIANGULAR SYMMETRIC STORAGE MODE)
C
      DO 220 I=1,3
      DO 220 J=I,3
      T = 0.0
      DO 210 K=1,3
      T = T + UMAT(K,I)*UMAT(K,J)
  210 CONTINUE
      IA = I + (J*J-J)/2
      UTR(IA) = T
  220 CONTINUE
C%    WRITE(6,999) UTR
C
C     CALCULATE EIGENVALUES AND VECTORS
C
      CALL EIGEN (UTR, A, 3, 0)
	CALL ESORT(UTR,A,3,0)
C%    WRITE(6,999) UTR
C
      ROOT(1) = UTR(1)
      ROOT(2) = UTR(3)
      ROOT(3) = UTR(6)
C%    WRITE(6,999) ROOT
C%    WRITE(6,999) A
C
C     SET A3 = A1 CROSS A2
C     ROOTS ARE IN ORDER R(1) >= R(2) >= R(3) >= 0
C
      A(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
      A(2,3) = A(3,1)*A(1,2) - A(1,1)*A(3,2)
      A(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)
C%    WRITE(6,999) A
C
C     VECTOR SET B=U.A
C
      DO 240 I=1,3
      DO 240 J=1,3
      T = 0.0
      DO 230 K=1,3
  230 T = T + UMAT(J,K)*A(K,I)
      B(J,I) = T
  240 CONTINUE
C
C      NORMALIZE B1 AND B2 AND CALCULATE B3 = B1 CROSS B2
C
      B1 = SQRT( B(1,1)*B(1,1) + B(2,1)*B(2,1) + B(3,1)*B(3,1) )
      B2 = SQRT( B(1,2)*B(1,2) + B(2,2)*B(2,2) + B(3,2)*B(3,2) )
      DO 250 I=1,3
      B(I,1) = B(I,1)/B1
  250 B(I,2) = B(I,2)/B2
C
C      CHECK FOR LEFT HANDED ROTATION
C
      B13 = B(2,1)*B(3,2) - B(3,1)*B(2,2)
      B23 = B(3,1)*B(1,2) - B(1,1)*B(3,2)
      B33 = B(1,1)*B(2,2) - B(2,1)*B(1,2)
C
      S = B13*B(1,3) + B23*B(2,3) + B33*B(3,3)
      IF (S .LT. 0.0) ISIG = -1
      B(1,3) = B13
      B(2,3) = B23
      B(3,3) = B33
C%    WRITE(6,999) B
C
C     CALCULATE ROTATION MATRIX R
C
      DO 270 I=1,3
      DO 270 J=1,3
      T = 0.0
      DO 260 K=1,3
  260 T = T + B(I,K)*A(J,K)
      R(I,J) = T
  270 CONTINUE
C
C     RMS ERROR
C
      DO 280 I=1,3
      IF (ROOT(I) .LT. 0.0) ROOT(I) = 0.0
      ROOT(I) = SQRT( ROOT(I) )
  280 CONTINUE
C
C     CHANGE SIGN OF EVAL #3 IF LEFT HANDED
C
      IF (ISIG .LT. 0) ROOT(3)=-ROOT(3)
      RTSUM = ROOT(3) + ROOT(2) + ROOT(1)
      RETURN
      END
C---- SUBROUTINE TO COMPUTE EIGENVALUES & EIGENVECTORS OF A REAL SYMMETRIC
C---- MATRIX, STOLEN FROM IBM SSP MANUAL (SEE P165)
C---- DESCRIPTION OF PARAMETERS -
C---- A - ORIGINAL MATRIX STORED COLUMNWISE AS UPPER TRIANGLE ONLY,
C---- I.E. "STORAGE MODE" = 1.  EIGENVALUES ARE WRITTEN INTO DIAGONAL
C---- ELEMENTS OF A  I.E.  A(1)  A(3)  A(6)  FOR A 3*3 MATRIX.
C---- R - RESULTANT MATRIX OF EIGENVECTORS STORED COLUMNWISE IN SAME
C---- ORDER AS EIGENVALUES.
C---- N - ORDER OF MATRICES A & R.
C---- MV = 0 TO COMPUTE EIGENVALUES & EIGENVECTORS.
	SUBROUTINE EIGEN(A,R,N,MV)
        DOUBLE PRECISION A,R,ANORM,ANRMX,THR,X,Y,SINX,SINX2,COSX,
	1COSX2,SINCS,RANGE
C---- FOR DOUBLE PRECISION, SQRT IN STATEMENTS 40,68,75&78 MUST BE DSQRT,
C---- ABS IN 62 MUST BE DABS AND 1.E-6 IN 5 MUST BE 1.D-12 .
	DIMENSION A(1),R(1)
C5	RANGE=1.E-6
5	RANGE=1.D-12
	IF(MV-1)10,25,10
10	IQ=-N
	DO20J=1,N
	IQ=IQ+N
	DO20I=1,N
	IJ=IQ+I
	R(IJ)=0.
	IF(I-J)20,15,20
15	R(IJ)=1.
20	CONTINUE
C---- INITIAL AND FINAL NORMS (ANORM & ANRMX)
25	ANORM=0.
	DO35I=1,N
	DO35J=I,N
	IF(I-J)30,35,30
30	IA=I+(J*J-J)/2
	ANORM=ANORM+A(IA)**2
35	CONTINUE
	IF(ANORM)165,165,40
C40	ANORM=SQRT(2.*ANORM)
40	ANORM=DSQRT(2.*ANORM)
	ANRMX=ANORM*RANGE/N
C---- INITIALIZE INDICATORS AND COMPUTE THRESHOLD
	IND=0
	THR=ANORM
45	THR=THR/N
50	L=1
55	M=L+1
C---- COMPUTE SIN & COS
60	MQ=(M*M-M)/2
	LQ=(L*L-L)/2
	LM=L+MQ
62	IF(ABS(A(LM))-THR)130,65,65
65	IND=1
	LL=L+LQ
	MM=M+MQ
	X=.5*(A(LL)-A(MM))
C68	Y=-A(LM)/SQRT(A(LM)**2+X*X)
68	Y=-A(LM)/DSQRT(A(LM)**2+X*X)
	IF(X)70,75,75
70	Y=-Y
C75	SINX=Y/SQRT(2.*(1.+(SQRT(1.-Y*Y))))
75	SINX=Y/DSQRT(2.*(1.+(DSQRT(1.-Y*Y))))
	SINX2=SINX**2
C78	COSX=SQRT(1.-SINX2)
78	COSX=DSQRT(1.-SINX2)
	COSX2=COSX**2
	SINCS=SINX*COSX
C---- ROTATE L & M COLUMNS
	ILQ=N*(L-1)
	IMQ=N*(M-1)
	DO125I=1,N
	IQ=(I*I-I)/2
	IF(I-L)80,115,80
80	IF(I-M)85,115,90
85	IM=I+MQ
	GOTO95
90	IM=M+IQ
95	IF(I-L)100,105,105
100	IL=I+LQ
	GOTO110
105	IL=L+IQ
110	X=A(IL)*COSX-A(IM)*SINX
	A(IM)=A(IL)*SINX+A(IM)*COSX
	A(IL)=X
115	IF(MV-1)120,125,120
120	ILR=ILQ+I
	IMR=IMQ+I
	X=R(ILR)*COSX-R(IMR)*SINX
	R(IMR)=R(ILR)*SINX+R(IMR)*COSX
	R(ILR)=X
125	CONTINUE
	X=2.*A(LM)*SINCS
	Y=A(LL)*COSX2+A(MM)*SINX2-X
	X=A(LL)*SINX2+A(MM)*COSX2+X
	A(LM)=(A(LL)-A(MM))*SINCS+A(LM)*(COSX2-SINX2)
	A(LL)=Y
	A(MM)=X
C---- TESTS FOR COMPLETION
C---- TEST FOR M = LAST COLUMN
130	IF(M-N)135,140,135
135	M=M+1
	GOTO60
C---- TEST FOR L =PENULTIMATE COLUMN
140	IF(L-(N-1))145,150,145
145	L=L+1
	GOTO55
150	IF(IND-1)160,155,160
155	IND=0
	GOTO50
C---- COMPARE THRESHOLD WITH FINAL NORM
160	IF(THR-ANRMX)165,165,45
165	RETURN
C---- SORT EIGENVALUES AND EIGENVECTORS IN DESCENDING ORDER OF EIGENVALUES
	END
	SUBROUTINE ESORT(A,R,N,MV)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION A(1),R(1)
	IQ=-N
	DO185I=1,N
	IQ=IQ+N
	LL=I+(I*I-I)/2
	JQ=N*(I-2)
	DO185J=I,N
	JQ=JQ+N
	MM=J+(J*J-J)/2
	IF(A(LL)-A(MM))170,185,185
170	X=A(LL)
	A(LL)=A(MM)
	A(MM)=X
	IF(MV-1)175,185,175
175	DO180K=1,N
	ILR=IQ+K
	IMR=JQ+K
	X=R(ILR)
	R(ILR)=R(IMR)
180	R(IMR)=X
185	CONTINUE
	RETURN
	END
