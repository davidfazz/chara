! checkpot.f90 - elapid - david gallina - dg@physik.uni-kassel.de - 2022
!-------------------------------------------------------------------------------
! module to check the gradient and hessian of a given potential by taylor 
! expanding the expression of the energy.
!-------------------------------------------------------------------------------
! 1st order (gradient):
! taylor exp. > f(x + tv) = f(x) + t<grad[f(x)],v>
! error func. > e(t) = |f(x + tv)| - f(x) - t<grad[f(x)],v> ~ O(t^2)
! logarithm   > log(e(t)) ~ 2 log(t) + const.
!-------------------------------------------------------------------------------
! 2nd order (hessian):
! taylor exp. > f(x + tv) = f(x) + t<grad[f(x)],v> + t^2<hess[f(x)][v],v>/2
! error func. > e(t) = |f(x + tv)| - f(x) - t<grad[f(x)],v> 
!                    - t^2<hess[f(x)][v],v>/2 ~ O(t^3)
! logarithm   > log(e(t)) ~ 3 log(t) + const.
!-------------------------------------------------------------------------------

module checkpot
  use constants 
  use potential
  use random
  implicit none 
  
contains 

!...............................................................................

  subroutine checkpot_()
    implicit none 
    double precision, dimension(c_n,c_n)  :: hessfx
    double precision, dimension(c_n)      :: xopt
    double precision, dimension(c_n)      :: gradfx 
    double precision, dimension(c_n)      :: hessfxv
    double precision, dimension(c_n)      :: v 
    double precision, dimension(c_n)      :: newgrad 
    double precision                      :: hessfxvv
    double precision                      :: gradfxv
    double precision                      :: fx
    double precision                      :: t 
    double precision                      :: newfx 
    double precision                      :: dfrst 
    double precision                      :: dscnd
    integer                               :: i
    integer                               :: ival
    integer                               :: nvals


! we calculate the error function for 1st and 2nd order:

    nvals = 14
    open(unit=13,file='checkpotential.dat',action='write',status='replace')

! 1) compute f(x), <grad[f(x)],v>, <hess[f(x)][v],v>;

    call potential_random(xopt)
    call potential_random(v)
    call potential_(xopt,fx,gradfx)
    call hessian_(xopt,hessfx)

    gradfxv = dot_product(gradfx,v)
    hessfxv = matmul(hessfx,v)
    hessfxvv = dot_product(hessfxv,v)
    
! 2) compute e(t) for several logarithmically spaced values of t in[1.d-8,1.d0]:

    do ival=0,nvals 
      t = 1.d0
      do i=1,ival 
        t = t / sqrt(10.d0)
      enddo

      call potential_(xopt+t*v,newfx,newgrad)
      dfrst = abs(newfx - fx - t*gradfxv)
      dscnd = abs(newfx - fx - t*gradfxv - t*t*hessfxvv/2.d0)

      write(13,'(3(F12.5))') log10(t),log10(dfrst),log10(dscnd)
    enddo
    close(13)

! 3) check that e(t) on a log-log plot has a slope of 2 for first order and a 
!    slope of 3 for second order.

    return 
  end subroutine checkpot_

!...............................................................................

end module checkpot