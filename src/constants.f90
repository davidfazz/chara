! constants.f90 - elapid - david gallina - university of kassel - 2022
!-------------------------------------------------------------------------------
! module containts all global constants and variables that are used in the 
! elapid code.
!-------------------------------------------------------------------------------

module constants
  implicit none 

!-------------------------------------------------------------------------------
! mathematical constants:

  double precision, parameter :: c_pi = 3.1415926535897932d0    ! pi 
  double precision, parameter :: c_e = 2.7182818284590452d0     ! exp 
  double precision, parameter :: c_squ2 = 1.4142135623730950d0  ! sqrt(2)

! physical constants:

  double precision, parameter :: c_hbar = 6.582119569d-16       ! [eV*s]
  double precision, parameter :: c_mu0 = 201335.4520789104d0    ! [meV*nm3/mub2]
  double precision, parameter :: c_mub = 5.7883818012d-5        ! [eV/T]
  double precision, parameter :: c_kb = 8.617333262d-5          ! [eV/K]
  double precision, parameter :: c_gamma = 1.760859644d2        ! [rad/mT*µs]

!-------------------------------------------------------------------------------
! global variables:

  double precision    :: c_temp       ! temperature [K]
  double precision    :: c_alpha      ! damping parameter
  integer             :: c_n          ! number of dimensions 
  integer             :: c_p          ! number of particles/atoms/molecules
  integer             :: c_seed       ! rng seed value
  character(len=15)   :: c_rngtype    ! type of rng
  character(len=15)   :: c_potential  ! potential

! potential variables:

  logical             :: dpls1d_log   ! dpls1d potential
  logical             :: dpls2d_log   ! dpls2d potential

! coordinate variables:

  logical             :: sys_pc_log   ! polar coordinates
  logical             :: sys_sc_log   ! spherical coordinates

!-------------------------------------------------------------------------------
! type definitions

  type t_atom
    double precision, dimension(3)  :: xyz    ! position [xyz]
    double precision                :: r      ! radius [nm]
    double precision                :: v      ! volume [nm3]
    double precision                :: m      ! magnetic moment [mub]
    integer                         :: knd    ! kind
  end type t_atom

!-------------------------------------------------------------------------------

end module constants                                                            