! potential.f90 - elapid - david gallina - dg@physik.uni-kassel.de - 2022
!-------------------------------------------------------------------------------

module potential
  use constants 
  use functions
  use dpls1d
  use dpls2d
  implicit none
  save

contains

!-------------------------------------------------------------------------------
! initialization routine of potential:

  subroutine potential_init()
    implicit none

    if (dpls1d_log) then 
      call dpls1d_init()
    elseif (dpls2d_log) then
      call dpls2d_init()
    endif

    return
  end subroutine potential_init

!-------------------------------------------------------------------------------
! deinitialization:

  subroutine potential_deinit()
    implicit none

    if (dpls1d_log) then 
      call dpls1d_deinit()
    elseif (dpls2d_log) then
      call dpls2d_deinit()
    endif

    return
  end subroutine potential_deinit

!-------------------------------------------------------------------------------
! this routine calculates the energy and gradient for a given configuration.

  subroutine potential_(x,e,g)
    implicit none
    double precision, dimension(c_n), intent(in)  :: x 
    double precision, intent(out)                 :: e 
    double precision, dimension(c_n), intent(out) :: g

    if (dpls1d_log) then 
      call dpls1d_egrad(x,e,g)
    elseif (dpls2d_log) then
      call dpls2d_egrad(x,e,g)
    endif

    return
  end subroutine potential_

!-------------------------------------------------------------------------------
! this routine calculates the hessian matrix for a given magnetic configuration.

  subroutine hessian_(x,h)
    implicit none
    double precision, dimension(c_n)      :: x
    double precision, dimension(c_n,c_n)  :: h

    if (dpls1d_log) then 
      call dpls1d_hess(x,h)
    elseif (dpls2d_log) then
      call dpls2d_hess(x,h)
    endif

    return
  end subroutine hessian_

!-------------------------------------------------------------------------------
! due to the periodicity of spherical coordinates, it is necessary to project
! them back into their initial interval from time to time. we use the following
! interval: x in ([0,pi],[0,2pi]).

  subroutine potential_rescale(x)
    implicit none
    double precision, dimension(c_n) :: x
    
    if (dpls1d_log) then 
      call dpls1d_rescale(x)
    elseif (dpls2d_log) then
      call dpls2d_rescale(x)
    endif

    return
  end subroutine potential_rescale

!...............................................................................
! the distance between two vectors on the unit sphere is given by the angle 
! between the two vectors.

  double precision function potential_distance(x,y)
    implicit none
    double precision, dimension(c_n) :: x
    double precision, dimension(c_n) :: y

    if (dpls1d_log) then 
      potential_distance = dpls1d_distance(x,y)
    elseif (dpls2d_log) then
      potential_distance = dpls2d_distance(x,y)
    endif

    return
  end function potential_distance

!-------------------------------------------------------------------------------
! create a random configuration.

  subroutine potential_random(x)
    implicit none
    double precision, dimension(c_n), intent(inout)  :: x

    if (dpls1d_log) then 
      call dpls1d_random(x)
    elseif (dpls2d_log) then
      call dpls2d_random(x)
    endif

    return
  end subroutine potential_random

!-------------------------------------------------------------------------------
! calculate entropy:

  subroutine potential_entropy(x,s,nev)
    implicit none
    double precision, dimension(c_n), intent(in) :: x
    double precision, intent(out)                :: s
    integer, intent(out)                         :: nev 

    if (dpls1d_log) then 
      call dpls1d_entropy(x,s,nev)
    elseif (dpls2d_log) then
      call dpls2d_entropy(x,s,nev)
    endif

    return
  end subroutine potential_entropy

!-------------------------------------------------------------------------------
! calculate time-reversal vector:

  subroutine potential_reversal(x,y)
    implicit none 
    double precision, dimension(c_n), intent(in)  :: x
    double precision, dimension(c_n), intent(out) :: y

    if (dpls1d_log) then 
      call dpls1d_reversal(x,y)
    elseif (dpls2d_log) then
      call dpls2d_reversal(x,y)
    endif

    return 
  end subroutine potential_reversal

!-------------------------------------------------------------------------------
! calculate order parameters:

  subroutine potential_order(x,order,angle)
    implicit none 
    double precision, dimension(c_n), intent(in)  :: x
    double precision, intent(out)                 :: order
    double precision, intent(out)                 :: angle

    if (dpls1d_log) then 
      call dpls1d_order(x,order,angle)
    elseif (dpls2d_log) then
      call dpls2d_order(x,order,angle)
    endif

    return 
  end subroutine potential_order

!-------------------------------------------------------------------------------
! calculate time-reversal vector:

  subroutine potential_rates(xone,xts,xtwo,alpha,rone,rtwo)
    implicit none 
    double precision, dimension(c_n), intent(in)  :: xone
    double precision, dimension(c_n), intent(in)  :: xts
    double precision, dimension(c_n), intent(in)  :: xtwo
    double precision, intent(in)                  :: alpha
    double precision, intent(out)                 :: rone
    double precision, intent(out)                 :: rtwo

    rone = 0.d0
    rtwo = 0.d0

    return 
  end subroutine potential_rates

!-------------------------------------------------------------------------------

end module potential
