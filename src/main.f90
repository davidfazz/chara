! main.f90 - elapid- david gallina - dg@physik.uni-kassel.de - 2022
!..............................................................................

program elapid
  use constants
  use functions
  use potential
  use explore
  use random
  use checkpot
  use omp_lib
  implicit none 
  double precision    :: stime
  double precision    :: etime
  double precision    :: dtime

  print *,'            _             _     _      __   ___    '
  print *,'           | |           (_)   | |    /_ | / _ \   '
  print *,'        ___| | __ _ _ __  _  __| |     | || | | |  '
  print *,'       / _ \ |/ _` | `_ \| |/ _` |     | || | | |  '
  print *,'      |  __/ | (_| | |_) | | (_| |     | || |_| |  '
  print *,'       \___|_|\__,_| .__/|_|\__,_|     |_(_)___/   '
  print *,'                   | |                             '            
  print *,'                   |_|                             '      
  
  stime = omp_get_wtime()

! initialize:

  call readparms_()

! intialise potential:  

  call potential_init()

! explore:

  !call checkpot_()
  call explore_()

  etime = omp_get_wtime()

  dtime = etime - stime

contains 

!..............................................................................
! read main parameters:

  subroutine readparms_()
    implicit none 
    integer   :: info

! random number generator:

    call f_read_('randomtype','params.dat',c_rngtype,info)
    if (trim(c_rngtype) .eq. 'ziggurat') then 
      rng_type = 'ziggurat'
    else 
      rng_type = 'mersenne'
    endif

    call f_read_('randomseed','params.dat',c_seed,info)
    call rng_init(c_seed)

! set potential:

    dpls1d_log = .false.
    dpls2d_log = .false.
    
    sys_pc_log = .false.
    sys_sc_log = .false.
    call f_read_('potential','params.dat',c_potential,info)
    if (trim(c_potential) .eq. 'dpls1d') then
      dpls1d_log = .true.
    elseif (trim(c_potential) .eq. 'dpls2d') then 
      dpls2d_log = .true.
    endif

    return 
  end subroutine readparms_

!...............................................................................

end program elapid

