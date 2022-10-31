! random.f90 - elapid - david gallina - university of kassel - 2022
!-------------------------------------------------------------------------------
! module contains two implementations of random number generators: i) the 
! ziggurat generator by marsaglias adn the mersenne twister.
!-------------------------------------------------------------------------------

module random
  use constants
  implicit none

! general globals:

  integer                                 :: rng_seed
  logical                                 :: rng_set
  character(len=10)                       :: rng_type 

! mersenne twister parameters:

  integer, parameter                      :: mt_n = 624
  integer, parameter                      :: mt_n1 = 625
  integer, parameter                      :: mt_m = 397
  integer, parameter                      :: mt_mata =-1727483681
  integer, parameter                      :: mt_umask =-2147483647
  integer, parameter                      :: mt_lmask = 2147483647
  integer, parameter                      :: mt_tmaskb =-1658038656
  integer, parameter                      :: mt_tmaskc =-272236544 

! mersenne twister globals:

  integer                                 :: mti
  integer, dimension(mt_n)                :: mt 
  integer, dimension(2)                   :: mt_mag01 = (/ 0, mt_mata /)

! ziggurat parameters:

  double precision, parameter             :: zig_m1=2147483648.d0
  double precision, parameter             :: zig_m2=2147483648.d0

! ziggurat globals:

  double precision, dimension(128), save  :: zig_wn
  double precision, dimension(128), save  :: zig_fn
  double precision, dimension(256), save  :: zig_we
  double precision, dimension(256), save  :: zig_fe
  double precision                        :: zig_dn = 3.442619855899d0
  double precision                        :: zig_tn = 3.442619855899d0
  double precision                        :: zig_vn = 0.00991256303526217d0
  double precision                        :: zig_de = 7.697117470131487d0
  double precision                        :: zig_te =7.697117470131487d0
  double precision                        :: zig_ve = 0.003949659822581572d0
  double precision                        :: zig_q

  integer, dimension(128), save           :: zig_kn
  integer, dimension(256), save           :: zig_ke
  integer, save                           :: zig_iz
  integer, save                           :: zig_jz
  integer, save                           :: zig_hz
  integer, save                           :: zig_jsr = 123456789

  logical, save                           :: zig_init = .false.

contains

!-------------------------------------------------------------------------------
! initialization:

  subroutine rng_init(iseed)
    implicit none
    integer, intent(in) :: iseed 

    if (trim(rng_type) .eq. 'mersenne') then
    call mersenne_init(iseed)
    elseif (trim(rng_type) .eq. 'ziggurat') then 
    call ziggurat_init(iseed)
    endif 

    return 
  end subroutine rng_init

!-------------------------------------------------------------------------------
! create uniform random number in [0,1]:

  function rng_number()
    implicit none 
    double precision  :: rng_number

    if (.not. rng_set) then 
    stop '    rng is not set. program is aborted.'
    endif

    if (trim(rng_type) .eq. 'mersenne') then
    rng_number = mersenne_grnd()
    elseif (trim(rng_type) .eq. 'ziggurat') then 
    rng_number = ziggurat_uni()
    endif

    return 
  end function rng_number

!-------------------------------------------------------------------------------
! create gaussian distributed number:

  function rng_gauss(am,sd)
    implicit none 
    double precision              :: rng_gauss
    double precision, intent(in)  :: am 
    double precision, intent(in)  :: sd 
    double precision              :: uni
    integer                       :: i

    rng_gauss = 0.d0
    uni = 0.d0
    do i=1,12
    if (trim(rng_type) .eq. 'mersenne') then
        uni = uni + mersenne_grnd()
    elseif (trim(rng_type) .eq. 'ziggurat') then 
        uni = uni + ziggurat_uni()
    endif
    enddo
    rng_gauss = (uni - 6.d0)*sd + am

    return 
  end function rng_gauss

!-------------------------------------------------------------------------------
! ziggurat algorithm:

  subroutine ziggurat_init(iseed)
    implicit none 
    integer, intent(in) :: iseed 
    integer             :: i 

    zig_jsr = iseed

    zig_q = zig_vn * exp(-0.5d0*zig_dn*zig_dn)
    zig_kn(1) = int((zig_dn/zig_q)*zig_m1)
    zig_kn(2) = 0
    zig_wn(1) = zig_q/zig_m1
    zig_wn(128) = zig_dn/zig_m1
    zig_fn(1) = 1.d0
    zig_fn(128) = exp(-0.5d0*zig_dn*zig_dn)
    do i=127,1,-1
    zig_dn = sqrt(-2.d0*log(zig_vn/zig_dn + exp(-0.5d0*zig_dn*zig_dn)))
    zig_kn(i+1) = int((zig_dn/zig_tn)*zig_m1)
    zig_tn = zig_dn
    zig_fn(i) = exp(-0.5d0*zig_dn*zig_dn)
    zig_wn(i) = zig_dn/zig_m1
    enddo

    zig_q = zig_ve*exp(zig_de)
    zig_ke(1) = int((zig_de/zig_q)*zig_m2)
    zig_ke(2) = 0
    zig_we(1) = zig_q/zig_m2
    zig_we(256) = zig_de/zig_m2
    zig_fe(1) = 1.d0
    zig_fe(256) = exp(-zig_de)
    do i=255,1,-1
    zig_de = -log(zig_ve/zig_de + exp(-zig_de))
    zig_ke(i+1) = int(zig_m2 * (zig_de/zig_te))
    zig_te = zig_de
    zig_fe(i) = exp(-zig_de)
    zig_we(i) = zig_de/zig_m2
    enddo
    rng_set = .true.

    return
  end subroutine ziggurat_init

!-------------------------------------------------------------------------------

  function ziggurat_shr3()
    implicit none 
    integer             :: ziggurat_shr3

    zig_jz = zig_jsr 
    zig_jsr = ieor(zig_jsr,ishft(zig_jsr,13))
    zig_jsr = ieor(zig_jsr,ishft(zig_jsr,-17))
    zig_jsr = ieor(zig_jsr,ishft(zig_jsr,5))
    ziggurat_shr3 = zig_jz + zig_jsr
    
    return
  end function ziggurat_shr3

!-------------------------------------------------------------------------------

  function ziggurat_uni()
    implicit none 
    double precision    :: ziggurat_uni

    ziggurat_uni = 0.5d0 + 0.2328306d-9 * ziggurat_shr3()

    return
  end function ziggurat_uni

!-------------------------------------------------------------------------------
! mersenne twister:

  subroutine mersenne_init(iseed)
    implicit none 
    integer, intent(in) :: iseed
    integer :: seed 
    
    seed = iseed
    mt(1) = iand(seed,-1)
    do mti=2,mt_n
    mt(mti) = iand(69069 * mt(mti-1),-1)
    enddo
    mti = mt_n
    rng_set = .true.

    return
  end subroutine mersenne_init

!-------------------------------------------------------------------------------

  function mersenne_grnd()
    implicit none 
    double precision        :: mersenne_grnd
    integer    :: y
    integer    :: kk

    if (mti .ge. mt_n) then 
    if (mti .eq. mt_n + 1) then 
        stop '    rng is not set. program aborted.'
    endif

    do kk=1,mt_n-mt_m
        y = ior(iand(mt(kk),mt_umask), iand(mt(kk+1),mt_lmask))
        mt(kk) = ieor(ieor(mt(kk+mt_m), ishft(y,-1)),mt_mag01(iand(y,1)))
    enddo

    do kk=mt_n-mt_m+1,mt_n-1
        y = ior(iand(mt(kk),mt_umask), iand(mt(kk+1),mt_lmask))
        mt(kk) = ieor(ieor(mt(kk+(mt_m-mt_n)), ishft(y,-1)),mt_mag01(iand(y,1)))
    enddo

    y = ior(iand(mt(mt_n),mt_umask), iand(mt(1),mt_lmask))
    mt(mt_n) = ieor(ieor(mt(mt_m),ishft(y,-1)),mt_mag01(iand(y,1)))
    mti = 1
    endif

    y = mt(mti)
    mti = mti + 1
    y = ieor(y,ishft(y,-11))
    y = ieor(y,iand(ishft(y,7),mt_tmaskb))
    y = ieor(y,iand(ishft(y,15),mt_tmaskc))
    y = ieor(y,ishft(y,-18))

    if (y .lt. 0) then
    mersenne_grnd = (dble(y) + 2.d0**32) / (2.d0**32)
    else
    mersenne_grnd = dble(y) / (2.d0**32)
    endif

    return
  end function mersenne_grnd

!-------------------------------------------------------------------------------

end module random