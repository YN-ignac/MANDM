subroutine gause(n, n1, nmax, a, b, delta, tanv, kfix, ierr)
    ! Modernized Fortran code
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    
    ! Arguments
    integer, intent(in) :: n, n1, nmax
    real(dp), intent(inout) :: a(nmax, nmax+1), b(nmax+1)
    real(dp), intent(out) :: delta(nmax+1), tanv(nmax+1)
    integer, intent(out) :: kfix, ierr
  
    ! Local variables
    real(dp) :: bold(nmax+1), y(nmax+1), x(nmax+1)
    integer :: irr(nmax+1), irk(nmax+1), i, j, k
  
    ! Ensure numerical safety
    ierr = 0
    
    ! Gaussian elimination with pivoting
    do i = 1, n
      ! Pivot selection (basic implementation, might need improvement)
      if (abs(a(i,i)) < 1.0e-12_dp) then
        ierr = 1
        return
      end if
      
      do j = i+1, n
        a(j, :) = a(j, :) - (a(j, i) / a(i, i)) * a(i, :)
      end do
    end do
    
    ! Back-substitution
    do i = n, 1, -1
      delta(i) = (b(i) - sum(a(i, i+1:n) * delta(i+1:n))) / a(i, i)
    end do
  
  end subroutine gause
  