! dpls1d.f90 - elapid - david gallina - university of kassel - 2022
!-------------------------------------------------------------------------------
! potential of n magnetic nanoparticles (mnp) located on a planar surface. each 
! mnp is characterized by a single superspin µ_i and they interact with each 
! other through the dipole-dipole interaction. other magnetic effects such as 
! magnetic anisotropies and higher order multipole corrections are not 
! considered.
! each superspin µ_i is characterized by a single polar angle phi_i, i.e., the 
! superspin is assumed to always lie in the plane of the surface.  
!-------------------------------------------------------------------------------
! settings:
! > p_lattice: current options are 'pargram', 'honeycomb', and 'kagome'.
! > p_omega:   lattice constant of the lattice;
! > p_kappa:   second lattice constant of the lattice;
! > p_sigma:   disorder of the square lattice. the deviations from the perfect 
!              lattice positions follow a gaussian distribution with the width 
!              p_omega*p_sigma.
! > p_rad:     radius of the magnetic nanoparticles associated with the dipoles.
! > p_mag:     magnetization of the magnetic nanoparticles in [µB / nm^3].
! > p_nx:      number of nanoparticles along the x-direction of the lattice.
! > p_ny:      number of nanoparticles along the y-direction of the lattice.
! > p_bnx:     number of identical cells along the x-direction to simulate 
!              periodic boundary conditions.
! > p_bny:     number of identical cells along the y-direction to simulate 
!              periodic boundary conditions.
! > p_hext:    size and direction of an externally applied magnetic field.
!-------------------------------------------------------------------------------
! references:
! > [Bae11]: baek et al., "Kosterlitz-Thouless transition of magnetic dipoles on 
!            the two-dimensional plane", phys. rev. b 83, 184409 (2011).
! > [Hol15]: holden et al., "monte carlo simulations of a kagome lattice with
!            magnetic dipolar interactions", phys. rev. b 91, 224425 (2015).
!-------------------------------------------------------------------------------
! implemented lattice geometries:
! > parallelogram
! > honeycomb
! > kagome
!-------------------------------------------------------------------------------
! important:
! the lattice dimensions in case of parallelogram lattices are not that 
! important. however, in the case of honeycomb and kagome lattices, one should
! pick the dimensions such that the ground state degeneracy is not lifted. in 
! honeycomb lattices the ground state is continuously degenerate and in kagome
! lattices it is 6-fold degenerate.
!------------------------------------------------------------------------------- 

module dpls1d
  use constants
  use functions
  use random
  implicit none 
    
!-------------------------------------------------------------------------------
! module arrays:

! interaction arrays:

  double precision, dimension(:,:), allocatable, private  :: xx3
  double precision, dimension(:,:), allocatable, private  :: xy3
  double precision, dimension(:,:), allocatable, private  :: yy3

! ensemble arrays:

  type(t_atom), dimension(:), allocatable, private        :: p_ens

!-------------------------------------------------------------------------------
! module variables:

! lattice type:

  character(len=15), private              :: p_lattice

! lattice constant and lattice disorder:

  double precision, private               :: p_omega
  double precision, private               :: p_kappa
  double precision, private               :: p_sigma 
  double precision, private               :: p_alpha

! particle radius and magnetization:

  double precision, private               :: p_rad
  double precision, private               :: p_mag

! number of particles in the unit cell in x and y direction:

  integer, private                        :: p_nx
  integer, private                        :: p_ny

! number of unit cell copies along x and y direction:

  integer, private                        :: p_bnx
  integer, private                        :: p_bny

! maximum distance from center of system:

  double precision, private               :: p_maxd 

! magnetic field:
  
  double precision, dimension(3), private :: p_hext

! optimization steps for disordered ensembles:

  integer, parameter, private             :: p_imaxopt = 50000

! debugging output;

  logical, private                        :: p_debug 
  integer, private                        :: p_copies                                 

contains

!-------------------------------------------------------------------------------
! initialization routine:

  subroutine dpls1d_init()
    implicit none 
    double precision, dimension(3)                        :: cvec
    double precision, dimension(3)                        :: cpos
    double precision, dimension(3)                        :: rvec
    double precision, dimension(3)                        :: xvec
    double precision, dimension(3)                        :: yvec
    double precision, dimension(:), allocatable           :: xgrd
    double precision                                      :: angle
    double precision                                      :: radius
    double precision                                      :: rabs
    double precision                                      :: rab5
    double precision                                      :: dx
    double precision                                      :: dy
    double precision                                      :: xpos
    double precision                                      :: ypos
    double precision                                      :: dgauss
    integer                                               :: isp
    integer                                               :: ix 
    integer                                               :: iy 
    integer                                               :: iz 
    integer                                               :: pi
    integer                                               :: ibx 
    integer                                               :: iby 
    integer                                               :: i 
    integer                                               :: j
    integer                                               :: k
    logical                                               :: update

!...............................................................................
! first, parameters are read from the parameter file. then, the number of 
! particles are calculated. one peculiarity of the kagome lattice is the factor 
! of 3 that enters the equation. in this case, the lattice basis contains three 
! particles. 
 
    call dpls1d_readparms()

    c_p = p_nx*p_ny
    if (trim(p_lattice) .eq. 'kagome') then 
      c_p = 3*c_p 
    endif
    c_n = c_p 

! we create some output of the parameters:

    print *,'    > dpls1d potential initializing..'
    print *,'      > lattice type:  ',trim(p_lattice)
    print *,'      > x dimensions:  ',f_trimmed_(p_nx)
    print *,'      > y dimensions:  ',f_trimmed_(p_ny)
    print *,'      > constant (nm): ',f_trimmed_(p_omega)
    print *,'      > sigma (%):     ',f_trimmed_(100.d0*p_sigma)
    print *,'      > radius (nm):   ',f_trimmed_(p_rad)
    print *,''
    
! we allocate the memory for the arrays that we use to store the configurations 
! and interaction parameters. moreover, some of the arrays are directly 
! initialized. 

    allocate(p_ens(c_p))
    allocate(xgrd(3*c_p))
    allocate(xx3(c_p,c_p))
    allocate(xy3(c_p,c_p))
    allocate(yy3(c_p,c_p))

    xgrd = 0.d0
    xx3 = 0.d0
    xy3 = 0.d0
    yy3 = 0.d0

!...............................................................................
! we create the magnetic nanoparticles and place them on the chosen lattice. 

    do i=1,c_p 
      p_ens(i)%r = p_rad
      p_ens(i)%v = 4.d0*c_pi/3.d0 * p_ens(i)%r**3
      p_ens(i)%m = p_ens(i)%v * p_mag * c_mub
    enddo

! we calculate the height and width of the unit cell:

    if (trim(p_lattice) .eq. 'pargram') then
      dx = dble(p_nx) * p_omega 
      dy = dble(p_ny) * p_kappa * sin(p_alpha)
    elseif (trim(p_lattice) .eq. 'honeycomb') then
      dx = dble(p_nx) * sqrt(3.d0) * p_omega
      dy = dble(p_ny) * 3.d0 / 4.d0 * p_omega
    elseif (trim(p_lattice) .eq. 'kagome') then
      dx = dble(p_nx) * p_omega * 2.d0
      dy = dble(p_ny) * p_omega * sqrt(3.d0)
    endif

! maximum distance from center of unit cell:

    if (p_bnx .gt. 100) then
      p_maxd = dble(p_bnx - 25) * dy
    else 
      p_maxd = 1.d15
    endif

! oblique lattice structures:

    isp = 0
    if (trim(p_lattice) .eq. 'pargram') then
      do ix=0,p_nx-1
        do iy=0,p_ny-1
          isp = isp + 1
          p_ens(isp)%xyz(3) = 0.d0
          p_ens(isp)%xyz(2) = (dble(iy) + 0.5d0) * p_kappa * sin(p_alpha)
          p_ens(isp)%xyz(1) = (dble(ix) + 0.5d0) * p_omega &
                            + p_ens(isp)%xyz(2) * cos(p_alpha)
        enddo
      enddo

! honeycomb lattice structures:

    elseif (trim(p_lattice) .eq. 'honeycomb') then
      do ix=1,p_nx 
        xpos = dble(ix) * sqrt(3.d0) * p_omega 
        ypos = p_omega / 2.d0
        do iy=1,p_ny 
          isp = isp + 1
  
          if (mod(iy,4) .eq. 1) then
            ypos = ypos + p_omega / 2.d0
            xpos = xpos - sqrt(3.d0) * p_omega / 2.d0
          elseif (mod(iy,4) .eq. 2) then 
            ypos = ypos + p_omega 
          elseif (mod(iy,4) .eq. 3) then 
            ypos = ypos + p_omega / 2.d0
            xpos = xpos + sqrt(3.d0) * p_omega /2.d0
          elseif (mod(iy,4) .eq. 0) then 
            ypos = ypos + p_omega
          endif 
          p_ens(isp)%xyz = [ xpos, ypos, 0.d0]
        enddo
      enddo

! kagome lattice structures:

    elseif (trim(p_lattice) .eq. 'kagome') then
      do ix=0,p_nx-1
        do iy=0,p_ny-1     
          xvec = dble(ix) * p_omega * (/ 2.d0, 0.d0, 0.d0 /) &
               + dble(iy) * p_omega * (/ 1.d0, sqrt(3.d0), 0.d0 /)  
          do iz=1,3
            if (iz .eq. 1) then 
              yvec = (/ 0.d0, 0.d0, 0.d0 /)
            elseif (iz .eq. 2) then
              yvec = p_omega * (/ 1.d0, 0.d0, 0.d0 /)
            elseif (iz .eq. 3) then
              yvec = p_omega * (/ 1.d0, sqrt(3.d0) , 0.d0 /) / 2.d0
            endif
            isp = isp + 1
            p_ens(isp)%xyz = xvec + yvec 
          enddo
        enddo
      enddo
    endif

!...............................................................................
! in order to create disorder, we displace all particles from their initial 
! lattice nodes according to a gaussian distribution.

    do i=1,c_p 
      angle = rng_number()*c_pi*2.d0
      dgauss = sqrt(cos(angle)*p_omega*cos(angle)*p_omega &
             + sin(p_alpha)*p_kappa*sin(p_alpha)*p_kappa*sin(angle)*sin(angle))
      radius = abs(rng_gauss(0.d0,dgauss*p_sigma))
      p_ens(i)%xyz = p_ens(i)%xyz + radius*(/ cos(angle), sin(angle), 0.d0 /)
    enddo

! we optimize the initial structure with a simple hard shell potential until 
! no particles touch each other. overall, a maximum of p_imaxopt steps are 
! taken. 

    do k=1,p_imaxopt  
      xgrd = 0.d0
      update = .false.
    
! calculate the force acting on each particle:

      do i=1,c_p 
        pi = 3*(i-1)
        do j=1,c_p 
          if (i .eq. j) then
            cycle 
          endif

! calculate distance between particles:
          
          rvec = p_ens(j)%xyz - p_ens(i)%xyz 
          rabs = sqrt(dot_product(rvec,rvec))

! if particles are touching, there is a repulsive force:

          if (rabs .lt. p_ens(i)%r + p_ens(j)%r) then
            xgrd(pi+1) = xgrd(pi+1) - rvec(1)
            xgrd(pi+2) = xgrd(pi+2) - rvec(2)
            xgrd(pi+3) = 0.d0
            update = .true.
          endif
        enddo

! we calculate the force acting on a particle from collisions with the left and 
! right unit cell boundaries:

        rvec = p_ens(i)%xyz 
        if (rvec(1) .lt. p_ens(i)%r + rvec(2)*cos(p_alpha)) then
          xgrd(pi+1) = xgrd(pi+1) + 1.d0 
          update = .true.
        elseif (rvec(1) .gt. dx + rvec(2)*cos(p_alpha) - p_ens(i)%r) then 
          xgrd(pi+1) = xgrd(pi+1) - 1.d0
          update = .true.
        endif

! colissions with the top and bottom unit cell boundaries:
        
        if (rvec(2) .lt. p_ens(i)%r) then
          xgrd(pi+2) = xgrd(pi+2) + 1.d0
          update = .true.
        elseif (rvec(2) .gt. dy - p_ens(i)%r) then
          xgrd(pi+2) = xgrd(pi+2) - 1.d0 
          update = .true.
        endif
      enddo

     ! print *,k,xgrd(1:10)

      if (update) then
        do i=1,c_p 
          pi = 3*(i-1)
          p_ens(i)%xyz(1) = p_ens(i)%xyz(1) + xgrd(pi+1)/500.d0
          p_ens(i)%xyz(2) = p_ens(i)%xyz(2) + xgrd(pi+2)/500.d0
        enddo
      else 
        exit
      endif 
    enddo

!...............................................................................
! we calculate the interaction matrices - checked!:

    print *,''
    print *,'    > periodic boundaries..'
    print *,''

! calculate the center of the unit cell:

    cpos(1) = dx + dy / tan(p_alpha)
    cpos(2) = dy 
    cpos(3) = 0.d0
    cpos = cpos / 2.d0
    
! calculate the interactions:
    
    if (p_debug) then
      open(unit=16,file='debug.dat',action='write',status='replace')
    endif

    do ibx=-p_bnx,p_bnx
      do iby=-p_bny,p_bny
        cvec(1) = dble(ibx)*dx + dble(iby)*dy/tan(p_alpha) 
        cvec(2) = dble(iby)*dy
        cvec(3) = 0.d0
        do i=1,c_p 
          if ((p_debug) .and. (abs(ibx) .le. 1) .and. (abs(iby) .le. 1)) then
            write(16,*) p_ens(i)%xyz - cvec
          endif
          do j=1,c_p 
            if (f_norm(p_ens(j)%xyz + cvec - cpos) .gt. p_maxd) then
              cycle 
            endif
            rvec = p_ens(i)%xyz - p_ens(j)%xyz - cvec
            rabs = f_norm(rvec)
            rab5 = rabs**5

            if ((ibx .eq. 0) .and. (iby .eq. 0)) then
              if (i .eq. j) then
                cycle 
              endif
              xx3(i,j) = xx3(i,j) + rvec(1)*rvec(1) / rab5
              xy3(i,j) = xy3(i,j) + rvec(1)*rvec(2) / rab5
              yy3(i,j) = yy3(i,j) + rvec(2)*rvec(2) / rab5 
            else 
              xx3(i,j) = xx3(i,j) + rvec(1)*rvec(1) / rab5
              xy3(i,j) = xy3(i,j) + rvec(1)*rvec(2) / rab5
              yy3(i,j) = yy3(i,j) + rvec(2)*rvec(2) / rab5 
            endif 
          enddo
        enddo
      enddo
    enddo 

    if (p_debug) then 
      close(16)
    endif
  
!...............................................................................
! deallocate memory - checked!:

    deallocate(xgrd)

! output:
    
    open(unit=11,file='ensemble.data',action='write',status='replace')
    do i=1,c_p 
      write(11,'(5(F15.5))') p_ens(i)%xyz,p_ens(i)%r,p_ens(i)%m
    enddo 
    close(11)
    
    return
  end subroutine dpls1d_init

!-------------------------------------------------------------------------------
! deinitialization:

  subroutine dpls1d_deinit()
    implicit none

    deallocate(p_ens)
    deallocate(xx3)
    deallocate(xy3)
    deallocate(yy3)

    return 
  end subroutine dpls1d_deinit

!-------------------------------------------------------------------------------
! energy and gradient of input configuration x;

  subroutine dpls1d_egrad(x,e,g)
    implicit none
    double precision, dimension(c_n), intent(in)  :: x 
    double precision, intent(out)                 :: e 
    double precision, dimension(c_n), intent(out) :: g 
    double precision, dimension(c_n)              :: cosx
    double precision, dimension(c_n)              :: sinx
    double precision, dimension(3)                :: si
    double precision, dimension(3)                :: sj
    double precision, dimension(3)                :: di
    double precision, dimension(3)                :: term
    integer                                       :: i
    integer                                       :: j

! initialize variables:

    sinx = sin(x)
    cosx = cos(x)

! calculate dipolar energy and gradient;

    e = 0.d0
    g = 0.d0
    do i=1,c_n
      si = p_ens(i)%m * [ cosx(i), sinx(i), 0.d0 ]
      di = p_ens(i)%m * [-sinx(i), cosx(i), 0.d0 ]
      
      term = 0.d0
      do j=1,c_n
        sj = p_ens(j)%m * [ cosx(j), sinx(j), 0.d0 ]

        term(1) = term(1) + yy3(i,j)*sj(1) - 2.d0*xx3(i,j)*sj(1) &
                - 3.d0*xy3(i,j)*sj(2)
        term(2) = term(2) + xx3(i,j)*sj(2) - 2.d0*yy3(i,j)*sj(2) &
                - 3.d0*xy3(i,j)*sj(1)
        term(3) = 0.d0
      enddo
      e = e + c_mu0/(4.d0*c_pi) * dot_product(si,term)
      g(i) = c_mu0/(4.d0*c_pi) * dot_product(di,term)
    enddo
    e = e / 2.d0

! calculate zeeman energy and gradient;

    do i=1,c_n 
      si = p_ens(i)%m * [ cosx(i), sinx(i), 0.d0 ]
      di = p_ens(i)%m * [-sinx(i), cosx(i), 0.d0 ]

      e = e - dot_product(si,p_hext)
      g(i) = g(i) - dot_product(di,p_hext)
    enddo

    return 
  end subroutine dpls1d_egrad

!-------------------------------------------------------------------------------
! calculate hessian matrix of input configuration x;
  
  subroutine dpls1d_hess(x,h)
    implicit none 
    double precision, dimension(c_n), intent(in)      :: x 
    double precision, dimension(c_n,c_n), intent(out) :: h
    double precision, dimension(c_n)                  :: sinx
    double precision, dimension(c_n)                  :: cosx
    double precision, dimension(3)                    :: sj
    double precision, dimension(3)                    :: dj
    double precision, dimension(3)                    :: di
    double precision, dimension(3)                    :: ti
    double precision, dimension(3)                    :: oterm
    double precision, dimension(3)                    :: dterm
    integer                                           :: i
    integer                                           :: j

! initialize variables:

    sinx = sin(x)
    cosx = cos(x)

! calculate hessian of dipole interaction;

    h = 0.d0
    do i=1,c_n
      di = p_ens(i)%m * [ -sinx(i),  cosx(i), 0.d0 ]
      ti = p_ens(i)%m * [ -cosx(i), -sinx(i), 0.d0 ]

      dterm = 0.d0
      do j=1,c_n
        sj = p_ens(j)%m * [ cosx(j), sinx(j), 0.d0 ]
        dj = p_ens(j)%m * [-sinx(j), cosx(j), 0.d0 ]

! off-diagonal terms:

        oterm(1) = yy3(i,j)*dj(1) - 2.d0*xx3(i,j)*dj(1) - 3.d0*xy3(i,j)*dj(2)
        oterm(2) = xx3(i,j)*dj(2) - 2.d0*yy3(i,j)*dj(2) - 3.d0*xy3(i,j)*dj(1)
        oterm(3) = xx3(i,j)*dj(3) + yy3(i,j)*dj(3)
        h(i,j) = h(i,j) + c_mu0/(4.d0*c_pi) * dot_product(di,oterm)

! diagonal terms:

        dterm(1) = dterm(1) + yy3(i,j)*sj(1) - 2.d0*xx3(i,j)*sj(1) - 3.d0*xy3(i,j)*sj(2)
        dterm(2) = dterm(2) + xx3(i,j)*sj(2) - 2.d0*yy3(i,j)*sj(2) - 3.d0*xy3(i,j)*sj(1)
        dterm(3) = 0.d0
      enddo

      h(i,i) = h(i,i) + c_mu0/(4.d0*c_pi) * dot_product(ti,dterm)
    enddo

! calculate hessian of zeeman interaction;

    do i=1,c_n
      ti = p_ens(i)%m * [ -cosx(i), -sinx(i), 0.d0 ] 
      
      h(i,i) = h(i,i) - dot_product(ti,p_hext)
    enddo

    return
  end subroutine dpls1d_hess

!-------------------------------------------------------------------------------
! rescaling of the vector x into the interval [0,2pi];

  subroutine dpls1d_rescale(x)
    implicit none
    double precision, dimension(c_n), intent(inout) :: x
    integer                                         :: i

    do i=1,c_n 
      do while (x(i) .lt. 0.d0)
        x(i) = x(i) + 2.d0*c_pi 
      enddo
      do while (x(i) .gt. 2.d0*c_pi)
        x(i) = x(i) - 2.d0*c_pi 
      enddo 
    enddo

    return 
  end subroutine dpls1d_rescale

!-------------------------------------------------------------------------------
! create random magnetic configuration:

  subroutine dpls1d_random(x)
    implicit none
    double precision, dimension(c_n)  :: x 
    integer                           :: i

    do i=1,c_n 
      x(i) = rng_number() * c_pi * 2.d0
    enddo

    return 
  end subroutine dpls1d_random

!-------------------------------------------------------------------------------
! calculate distance between vectors:

  function dpls1d_distance(x1,x2)
    implicit none 
    double precision                  :: dpls1d_distance
    double precision, dimension(c_n)  :: x1
    double precision, dimension(c_n)  :: x2
    double precision                  :: dx
    integer                           :: i

    dpls1d_distance = 0.d0
    do i=1,c_n 

! find angle between x1(i) und x2(i);

      dx = min(abs(x1(i) - x2(i)),2.d0*c_pi - abs(x1(i) - x2(i)))

! exceptions for numerical inaccuracies;

      dx = min(dx,1.d0)
      dx = max(dx,0.d0)

! add square of dx;

      dpls1d_distance = dpls1d_distance + dx*dx 
    enddo
    dpls1d_distance = sqrt(dpls1d_distance)

    return 
  end function dpls1d_distance

!-------------------------------------------------------------------------------
! find state having the opposite magnetization of x1;

  subroutine dpls1d_reversal(x1,x2)
    implicit none 
    double precision, dimension(c_n), intent(in)  :: x1
    double precision, dimension(c_n), intent(out) :: x2

    x2 = x1 + c_pi 
    call dpls1d_rescale(x2)

    return 
  end subroutine dpls1d_reversal

!-------------------------------------------------------------------------------
! vibrational entropy of magnetic configuration x:

  subroutine dpls1d_entropy(x,s,nev)
    implicit none 
    double precision, dimension(c_n), intent(in)  :: x 
    double precision, intent(out)                 :: s 
    integer, intent(out)                          :: nev
    double precision, dimension(c_n,c_n)          :: h 
    double precision, dimension(c_n,c_n)          :: evec
    double precision, dimension(c_n)              :: eval 
    integer                                       :: i
    integer                                       :: nzero

    call dpls1d_hess(x,h)
    call f_dsyev(h,evec,eval)

    nev = 0
    s = 0.d0
    do i=1,c_n
      if (abs(eval(i)) .lt. 1.d-5) then
        nzero = nzero + 1 
      elseif (eval(i) .lt. -1.d-5) then
        nev = nev + 1 
      else 
        s = s + log(eval(i))
      endif 
    enddo

    return 
  end subroutine dpls1d_entropy

!-------------------------------------------------------------------------------
! order parameters - see [Bae11],[Hol15];

  subroutine dpls1d_order(x,order,angle)
    implicit none
    double precision, dimension(c_n), intent(in)  :: x 
    double precision, intent(out)                 :: order
    double precision, intent(out)                 :: angle
    double precision, dimension(c_n)              :: y 
    double precision, dimension(c_n)              :: yt 
    double precision, dimension(6)                :: svec
    double precision, dimension(3)                :: avec
    double precision, dimension(2)                :: vec
    integer                                       :: ix
    integer                                       :: iy
    integer                                       :: wmod3
    integer                                       :: umod4
    integer                                       :: i

! square lattice;

    if (abs(p_alpha - c_pi/2.d0) .lt. 1.d-2) then      
      
! create gauge vector:

      y = x 
      i = 0
      do ix=1,p_nx 
        do iy=1,p_ny 
          i = i + 1
          if (mod(ix,2) .eq. 1) then 
            if (mod(iy,2) .eq. 1) then 
              y(i) = x(i)
            elseif (mod(iy,2) .eq. 0) then 
              y(i) = c_pi - x(i)
            endif
          elseif (mod(ix,2) .eq. 0) then
            if (mod(iy,2) .eq. 1) then 
              y(i) =-x(i)
            elseif (mod(iy,2) .eq. 0) then 
              y(i) = c_pi + x(i)
            endif
          endif
        enddo
      enddo
      call dpls1d_rescale(y)

! calculate order parameter [0,1]:

      vec = 0.d0
      do i=1,c_n 
        vec(1) = vec(1) + cos(y(i))
        vec(2) = vec(2) + sin(y(i))
      enddo
      order = sqrt(vec(1)*vec(1) + vec(2)*vec(2)) / dble(c_p)

! calculate angle [0,2pi]:

      vec = vec / sqrt(dot_product(vec,vec))
      angle = atan2(vec(2),vec(1))
      if (angle .lt. 0.d0) then
        angle = angle + 2.d0*c_pi 
      endif
    
! triangular lattice;

    elseif (abs(p_alpha - c_pi/3.d0) .lt. 1.d-2) then

! calculate magnetization;

      vec = 0.d0
      do i=1,c_n 
        vec(1) = vec(1) + cos(x(i))
        vec(2) = vec(2) + sin(x(i))
      enddo
      order = sqrt(vec(1)*vec(1) + vec(2)*vec(2)) / dble(c_n)

! calculate angle;

      vec(1:2) = vec(1:2) / sqrt(dot_product(vec(1:2),vec(1:2)))
      angle = atan2(vec(2),vec(1))

! honeycomb lattice;

    elseif (trim(p_lattice) .eq. 'honeycomb') then

      y = x 

! create gauge vector:
      
      i = 0
      do ix=1,p_nx 
        wmod3 = mod(ix-1,3)
        do iy=1,p_ny 
          umod4 = mod(iy-1,4)
  
! make correct increment based on coordinate system:
  
          i = i + 1
          y(i) = x(i) - c_pi/2.d0
  
! perform gauge transformation:
  
! wmod3 = 0:
  
          if ((umod4 .eq. 0) .and. (wmod3 .eq. 0)) then
            y(i) = y(i)
          elseif ((umod4 .eq. 1) .and. (wmod3 .eq. 0)) then
            y(i) =-y(i)
          elseif ((umod4 .eq. 2) .and. (wmod3 .eq. 0)) then   
            y(i) = y(i) + 2.d0*c_pi/3.d0 
          elseif ((umod4 .eq. 3) .and. (wmod3 .eq. 0)) then  
            y(i) = 2.d0*c_pi/3.d0 - y(i)        
  
! wmod3 = 1:
  
          elseif ((umod4 .eq. 0) .and. (wmod3 .eq. 1)) then    
            y(i) = y(i) - 2.d0*c_pi/3.d0      
          elseif ((umod4 .eq. 1) .and. (wmod3 .eq. 1)) then 
            y(i) =-y(i) - 2.d0*c_pi/3.d0         
          elseif ((umod4 .eq. 2) .and. (wmod3 .eq. 1)) then    
            y(i) = y(i)      
          elseif ((umod4 .eq. 3) .and. (wmod3 .eq. 1)) then 
            y(i) =-y(i)         
  
  ! wmod = 2:
  
          elseif ((umod4 .eq. 0) .and. (wmod3 .eq. 2)) then 
            y(i) = y(i) + 2.d0*c_pi/3.d0         
          elseif ((umod4 .eq. 1) .and. (wmod3 .eq. 2)) then  
            y(i) =-y(i) + 2.d0*c_pi/3.d0        
          elseif ((umod4 .eq. 2) .and. (wmod3 .eq. 2)) then  
            y(i) = y(i) - 2.d0*c_pi/3.d0        
          elseif ((umod4 .eq. 3) .and. (wmod3 .eq. 2)) then 
            y(i) =-y(i) - 2.d0*c_pi/3.d0         
          endif
        enddo
      enddo
  
      call dpls1d_rescale(y)
  
  ! create shifted vector:
  
      do i=1,c_n 
        if (y(i) .gt. c_pi) then 
          yt(i) = y(i) - 2.d0*c_pi 
        else
          yt(i) = y(i)
        endif 
      enddo
  
  ! calculate order parameter [0,1]:
  
      vec = 0.d0
      do i=1,c_n 
        vec(1) = vec(1) + cos(y(i))
        vec(2) = vec(2) + sin(y(i))
      enddo
  
      order = sqrt(vec(1)*vec(1) + vec(2)*vec(2)) / dble(c_p)
  
  ! calculate angle [0,2pi]:
  
      vec(1:2) = vec(1:2) / sqrt(dot_product(vec(1:2),vec(1:2)))
      angle = atan2(vec(2),vec(1))
      if (angle .lt. 0.d0) then
        angle = angle + 2.d0*c_pi 
      endif

! kagome lattice;

    elseif (trim(p_lattice) .eq. 'kagome') then

! create 3-fold magnetization vector;

      svec = 0.d0
      do i=1,c_n,3
        svec(1) = svec(1) + cos(x(i))
        svec(2) = svec(2) + sin(x(i)) 
        svec(3) = svec(3) + cos(x(i+1))
        svec(4) = svec(4) + sin(x(i+1))
        svec(5) = svec(5) + cos(x(i+2))
        svec(6) = svec(6) + sin(x(i+2)) 
      enddo

! calculate order for each sublattice;

      avec(1) = sqrt(dot_product(svec(1:2),svec(1:2)))
      avec(2) = sqrt(dot_product(svec(3:4),svec(3:4)))
      avec(3) = sqrt(dot_product(svec(5:6),svec(5:6)))

      order = (avec(1) + avec(2) + avec(3)) / dble(c_p)
      angle = 0.d0
    endif

    return
  end subroutine dpls1d_order

!-------------------------------------------------------------------------------
! read parameters:

  subroutine dpls1d_readparms()
    implicit none 
    integer   :: info

! lattice type:
! > parallelogram; honeycomb; kagome;

    call f_read_store_('lattice_type','params.dat',p_lattice,'pargram',info)
    
! lattice parameters:

    call f_read_store_('lattice_omega','params.dat',p_omega,10.d0,info)
    call f_read_store_('lattice_kappa','params.dat',p_kappa,10.d0,info)
    call f_read_store_('lattice_sigma','params.dat',p_sigma,0.d0,info)  
    call f_read_store_('lattice_alpha','params.dat',p_alpha,0.5d0,info)
    p_alpha = p_alpha / 180.d0 * c_pi

    if (trim(p_lattice) .eq. 'kagome') then 
      p_alpha = c_pi / 3.d0 
    elseif (trim(p_lattice) .eq. 'honeycomb') then 
      p_alpha = c_pi / 2.d0 
    endif 
    
    call f_read_store_('lattice_xrows','params.dat',p_nx,8,info)
    call f_read_store_('lattice_ycols','params.dat',p_ny,8,info)
    call f_read_store_('lattice_xperiod','params.dat',p_bnx,0,info)
    call f_read_store_('lattice_yperiod','params.dat',p_bny,0,info)

! particle parameters:

    call f_read_store_('particle_radius','params.dat',p_rad,4.d0,info)
    call f_read_store_('particle_moment','params.dat',p_mag,180.d0,info)
    
! field parameters:

    p_hext = 0.d0
    call f_read_('field_vector','params.dat',p_hext,info)

! debugging:

    call f_read_store_('lattice_debug','params.dat',p_debug,.false.,info)
    call f_read_store_('lattice_copies','params.dat',p_copies,1,info)

    return 
  end subroutine dpls1d_readparms 

!...............................................................................

end module dpls1d
