program cstr1exo

    implicit none

    ! Deklarace
    real(8) :: theta_min, theta_max
    integer :: steps
    real(8) :: gama, B, lambda, thetaC, beta
    COMMON /pars/ gama, B, lambda, thetaC, beta
    integer :: i
    real(8) :: theta, X, Da
    real(8) :: J
    dimension J (1:2, 1:2)
    real(8) :: lambda1, lambda2
    dimension lambda1(1:2), lambda2(1:2)


    OPEN(10,file='input_cstr1exo.txt')
    OPEN(11,file='output_cstr1exo.txt')
    OPEN(12,file='output_stable.txt')
    OPEN(13,file='output_unstable.txt')


    READ(10,*) theta_min, theta_max, steps
    READ(10,*) gama, B, lambda, thetaC, beta

    theta = theta_min

    do i = 1, steps
        call map(theta, X, Da)
        call jacobi(theta, X, Da, J)
        call eigen(J, lambda1, lambda2)

        WRITE(11,*) X, theta, Da, lambda1(1), lambda1(2), lambda2(1), lambda2(2)

        if (lambda1(1) .lt. 0 .and. lambda2(1) .lt. 0) then
        WRITE(12,*) X, theta, Da
        else
        WRITE(13,*) X, theta, Da
        end if

        WRITE(*,*) "i   X   theta   Re1     Im1     Re2     Im2"
        WRITE(*,*) "----------------------------------------------------"
        WRITE(*,*) i, X, theta, Da, lambda1(1), lambda1(2), lambda2(1), lambda2(2)
        
        theta = theta + (theta_max - theta_min) / steps
    end do
        

    CLOSE(10)
    CLOSE(11)
    CLOSE(12)
    CLOSE(13)

end program cstr1exo

subroutine map(theta, X, Da)
    
    implicit none

    real(8) :: theta, X, Da
    real(8) :: gama, B, lambda, thetaC, beta
    COMMON /pars/ gama, B, lambda, thetaC, beta

    X = theta/B + beta*(theta - thetaC)/(B*lambda)
    Da = lambda*X/(1-X)/exp(theta/(1+theta/gama))

    
end subroutine map

subroutine jacobi(theta, X, Da, J)
    
    implicit none

    real(8) :: theta, X, Da
    real(8) :: gama, B, lambda, thetaC, beta
    COMMON /pars/ gama, B, lambda, thetaC, beta
    real(8) :: J, exponent
    dimension J (1:2, 1:2)

    exponent = exp(theta/1+theta/gama)

    J(1,1) = -lambda-Da*exponent
    J(1,2) = Da*(1-X) * exponent * (1+theta/gama-theta/gama)/(1+theta/gama)**2
    J(2,1) = -Da*B*exponent
    J(2,2) = -lambda + B * Da * (1-X) * exponent * (1+theta/gama-theta/gama)/(1+theta/gama)**2 - beta
  
end subroutine jacobi

subroutine eigen(J, lambda1, lambda2)
    
    implicit none

    real(8) :: J
    dimension J (1:2, 1:2)
    real(8) :: lambda1, lambda2
    dimension lambda1(1:2), lambda2(1:2)
    real(8) :: tr_J, det_J, disc

    tr_J = J(1,1) + J(2,2)
    det_J = J(1,1)*J(2,2) - J(1,2)*J(2,1)

    disc = tr_J**2 - 4 * det_J

    if (disc .ge. 0.0) then
        lambda1(1) = (tr_J + disc**0.5) / 2.0
        lambda1(2) = 0.0
        lambda2(1) = (tr_J - disc**0.5) / 2.0
        lambda2(2) = 0.0
    else
        lambda1(1) = tr_J / 2
        lambda1(2) = abs(disc)**0.5 / 2.0
        lambda2(1) = tr_J / 2
        lambda2(2) = -abs(disc)**0.5 / 2.
    end if    
  
end subroutine eigen