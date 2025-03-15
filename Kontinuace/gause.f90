subroutine gause(n,n1,nmax,a,b,delta,tanv,kfix,ierr)

! solution of n linear equations for n1=n+1 unknowns
! based on gaussian elimination with pivoting.
!
! n      : number of equations
! a(,)   : n x n1 matrix of system
! b()    : right-hand sides
! delta(): on output vector of Newton corrections
! ierr   : if (ierr.ne.0) after return then a was numerically singular
! tanv() : tangent vector to the continuation curve,normalized on output
! pref(i): preference number for x(i) to be independent variable,
!          the lower is pref(i) the higher is preference of x(i)
! kfix    : resulting index of independent variable
!
! subprograms called: none
! fortran library functions: abs,sqrt
!
! ----------------------------------------------------------------------
  implicit double precision(a-h,o-z)
  dimension a(nmax,nmax+1),b(nmax+1),delta(nmax+1),tanv(nmax+1)
  dimension pref(nmax+1),bold(nmax+1)
  dimension y(nmax+1),x(nmax+1),irr(nmax+1),irk(nmax+1)
! ----------------------------------------------------------------------

  do i=1,n1
  bold(i)=b(i)
  enddo
  
  ierr=0

  do i=1,n1
    irk(i)=0
    irr(i)=0
    pref(i)=1.0
  enddo

  do id=1,n
    ir=1
    is=1
    amax=0.0

  do i=1,n
    if (irr(i).eq.0) then
      do j=1,n1
        p=pref(j)*abs(a(i,j))
        if (p.gt.amax) then
          ir=i
          is=j
          amax=p
        endif
      enddo
    endif
  enddo

  if (amax.eq.0.0) then
    ierr=1
    return
  endif

  irr(ir)=is
  
  do i=1,n
    if (i.ne.ir.and.a(i,is).ne.0.0) then
      p=a(i,is)/a(ir,is)
      do j=1,n1
        a(i,j)=a(i,j)-p*a(ir,j)
      enddo
      a(i,is)=0.0
      b(i)=b(i)-p*b(ir)
    endif
  enddo
  enddo
  
  do i=1,n
  ir=irr(i)
  x(ir)=b(i)/a(i,ir)
  irk(ir)=1
  enddo
  do i=1,n1
  if(irk(i).eq.0) kfix=i
  enddo
  do i=1,n
  ir=irr(i)
  y(ir)=-a(i,kfix)/a(i,ir)
  enddo
  do i=1,n1
  b(i)=x(i)
  delta(i)=y(i)
  enddo
  b(kfix)=0.0
  delta(kfix)=1.0
  s=0.0
  do i=1,n1
  s=s+delta(i)**2
  enddo
  do i=1,n1
  delta(i)=delta(i)/sqrt(s)
  enddo
  
  do i=1,n1
  tanv(i)=delta(i)
  delta(i)=b(i)
  b(i)=bold(i)
  enddo
  
  return
  end
