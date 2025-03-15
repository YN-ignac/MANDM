program kontinuace
    implicit none
    
    ! Declaration 
    integer :: nmax
    parameter (nmax=10) ! maximal dimension of y
    integer :: ndim, i
    real(8) :: f, J
    dimension f(1:nmax), J(1:nmax, 1:nmax+1)
    real(8) :: y
    dimension y(1:nmax+1) ! pre-made array of size 1:nmax+1

    real(8) :: gama, B, lambda, thetaC, beta
    COMMON /pars/ gama, B, lambda, thetaC, beta 

    ! Load data & assign data
    OPEN (10, file = 'input.txt')

    READ (10,*) ndim, (y(i), i=1,ndim+1)
    READ(10,*) gama, B, lambda, thetaC, beta



    call RHS(ndim, nmax, y, f, J)
    WRITE(*,*) 'y=', (y(i), i=1,ndim+1) 
    WRITE(*,*) 'f=', (f(i), i=1,ndim)
    WRITE(*,*) 'j=', (j(1,i), i=1,ndim+1)
    WRITE(*,*) '  ', (j(2,i), i=1,ndim+1)

    CLOSE (10)

end program kontinuace

subroutine RHS(ndim, nmax, y, f, J)

    implicit none

    integer :: ndim, nmax
    real(8) :: y, f, J
    dimension y(1:nmax+1), f(1:nmax), J(1:nmax, 1:nmax+1)
    real(8) :: x, theta, Da, exponent

    real(8) :: gama, B, lambda, thetaC, beta
    COMMON /pars/ gama, B, lambda, thetaC, beta 

    x = y(1)
    theta = y(2)
    Da = y(3)

    exponent = exp(theta/1+theta/gama)

    f(1) = -lambda*x + Da*(1-x)*exponent
    f(2) = -lambda*theta + Da*B*(1-x)*exponent - beta*(theta-thetaC)

    J(1,1) = -lambda-Da*exponent
    J(1,2) = Da*(1-X) * exponent * (1+theta/gama-theta/gama)/(1+theta/gama)**2
    J(1,3) = (1-x)*exponent

    J(2,1) = -Da*B*exponent
    J(2,2) = -lambda + B * Da * (1-X) * exponent * (1+theta/gama-theta/gama)/(1+theta/gama)**2 - beta
    J(2,3) = B*(1-x)*exponent
  

end subroutine RHS