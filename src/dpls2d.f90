! dpls2d.f90 - elapid - david gallina - university of kassel - 2022
!-------------------------------------------------------------------------------
! potential of n magnetic nanoparticles (mnp) located on a planar surface. each 
! mnp is characterized by a single superspin µ_i and they interact with each 
! other through the dipole-dipole interaction. other magnetic effects such as 
! magnetic anisotropies and higher order multipole corrections are not 
! considered.
! each superspin µ_i is characterized by the polar angle phi_i and azimuthal 
! angle theta_i.
!-------------------------------------------------------------------------------
! settings:
! > p_lattice: current options are 'square', 'triangular', 'honeycomb', and 
!              'kagome'.
! > p_lambda:  lattice constant of the underlying perfectly periodic square 
!              lattice.
! > p_sigma:   disorder of the square lattice. the deviations from the perfect 
!              lattice positions follow a gaussian distribution with the width 
!              p_lambda*p_sigma.
! > p_urec:    unit cell is rectangular. standard for square/honeycomb. 
!              available for square/triangular/honeycomb.
! > p_urho:    unit cell is rhomboidal. standard for triangular/kagome. 
!              available for triangular/kagome.
!
! > p_radius:  radius of the magnetic nanoparticles associated with the dipoles.
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
! > [Gra11]: graef et al., "On the computation of spherical designs by a new 
!            optimization approach based on fast spherical Fourier transforms", 
!            num. math. 119, 699 (2011).
! > [Hol15]: holden et al., "monte carlo simulations of a kagome lattice with
!            magnetic dipolar interactions", phys. rev. b 91, 224425 (2015).
!-------------------------------------------------------------------------------
! important:
! the lattice dimensions in case of square and triangular lattices are not that 
! important. however, in the case of honeycomb and kagome lattices, one should
! pick the dimensions such that the ground state degeneracy is not lifted. in 
! honeycomb lattices the ground state is continuously degenerate and in kagome
! lattices it is 6-fold degenerate.
!------------------------------------------------------------------------------- 

module dpls2d
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
  type(t_atom), dimension(:), allocatable, private        :: p_dbl

!-------------------------------------------------------------------------------
! module variables:

! lattice type:

  character(len=15), private                              :: p_lattice

! lattice constant and lattice disorder:

  double precision, private                               :: p_lambda
  double precision, private                               :: p_sigma 

! unit cell shape:

  logical, private                                        :: p_urec
  logical, private                                        :: p_urho

! particle radius and magnetization:

  double precision, private                               :: p_radius
  double precision, private                               :: p_mag

! number of particles in the unit cell in x and y direction:

  integer, private                                        :: p_nx
  integer, private                                        :: p_ny

! number of unit cell copies along x and y direction:

  integer, private                                        :: p_bnx
  integer, private                                        :: p_bny

! disorder with or without relaxing the particle positions:

  logical, private                                        :: p_rel_clc

! maximum distance from center of system:

  double precision, private                               :: p_maxd 

! magnetic field:
  
  double precision, dimension(3), private                 :: p_hext

contains

!-------------------------------------------------------------------------------
! initialization routine:

  subroutine dpls2d_init()
    implicit none 
    double precision, dimension(3)                        :: cvec
    double precision, dimension(3)                        :: cpos
    double precision, dimension(3)                        :: rvec
    double precision, dimension(3)                        :: xj
    double precision, dimension(3)                        :: xi
    double precision, dimension(3)                        :: xvec
    double precision, dimension(3)                        :: yvec
    double precision, dimension(:,:), allocatable         :: xgrd
    double precision                                      :: angle
    double precision                                      :: radius
    double precision                                      :: rabs
    double precision                                      :: rab5
    double precision                                      :: dx
    double precision                                      :: dy
    double precision                                      :: xpos
    double precision                                      :: ypos
    integer                                               :: isp
    integer                                               :: ix 
    integer                                               :: iy 
    integer                                               :: iz
    integer                                               :: ibx 
    integer                                               :: iby 
    integer                                               :: iout
    integer                                               :: i 
    integer                                               :: j
    integer                                               :: k
    logical                                               :: overlap
    logical                                               :: update

! read parameters:

    call dpls2d_readparms()

! set dimension and particle parameters. the kagome lattice has a base of 3 that
! has to be accounted for:

    if (trim(p_lattice) .eq. 'kagome') then 
      c_p = 3*p_nx*p_ny 
    else 
      c_p = p_nx*p_ny 
    endif
    c_n = c_p 

! allocate the required memory:

    allocate(p_ens(c_p))
    allocate(p_dbl(c_p))
    allocate(xgrd(c_p,2))

! output of settings:

    print *,'    > dpls2d potential initializing..'
    print *,'      > lattice type: ',trim(p_lattice)
    print *,'      > x dimensions: ',f_trimmed_(p_nx)
    print *,'      > y dimensions: ',f_trimmed_(p_ny)
    print *,'      > lambda (nm):  ',f_trimmed_(p_lambda)
    print *,'      > sigma (%):    ',f_trimmed_(100.d0*p_sigma)
    print *,'      > radius (nm):  ',f_trimmed_(p_radius)
    print *,''

! create magnetic nanoparticles:

    do i=1,c_p 
      p_ens(i)%r = p_radius 
      p_ens(i)%v = 4.d0*c_pi/3.d0 * p_ens(i)%r**3
      p_ens(i)%m = p_ens(i)%v * p_mag * c_mub
    enddo

! height and width of unit cell:

    if (trim(p_lattice) .eq. 'square') then
      dx = dble(p_nx) * p_lambda
      dy = dble(p_ny) * p_lambda
    elseif (trim(p_lattice) .eq. 'triangular') then
      dx = dble(p_nx) * p_lambda
      dy = dble(p_ny) * p_lambda * sqrt(3.d0) / 2.d0
    elseif (trim(p_lattice) .eq. 'honeycomb') then
      dx = dble(p_nx) * sqrt(3.d0) * p_lambda
      dy = dble(p_ny) * 3.d0 / 4.d0 * p_lambda
    elseif (trim(p_lattice) .eq. 'kagome') then
      dx = dble(p_nx) * p_lambda * 2.d0
      dy = dble(p_ny) * p_lambda * sqrt(3.d0)
    endif

! distance from center of unit cell:

    if (p_bnx .gt. 100) then
      p_maxd = dble(p_bnx - 25) * dy
    else 
      p_maxd = 1.d15
    endif

! create particle positions:

    isp = 0
    if (trim(p_lattice) .eq. 'square') then

! square lattice - rectangular unit cell;

      if (p_urec) then
        do ix=0,p_nx-1
          do iy=0,p_ny-1
            isp = isp + 1
            p_ens(isp)%xyz(1) = (dble(ix) + 0.5d0)*p_lambda
            p_ens(isp)%xyz(2) = (dble(iy) + 0.5d0)*p_lambda
            p_ens(isp)%xyz(3) = 0.d0
          enddo
        enddo
      else 

! raise critical error;

        call f_critical_('dpls2d',234)
      endif
    elseif (trim(p_lattice) .eq. 'triangular') then

! triangular lattice - rectangular unit cell;

      if (p_urec) then
        do ix=0,p_nx-1
          do iy=0,p_ny-1
            isp = isp + 1
            if (mod(iy,2) .eq. 1) then
              p_ens(isp)%xyz(1) = (dble(ix) + 1.d0)*p_lambda
              p_ens(isp)%xyz(2) = (dble(iy) + 0.5d0)*sqrt(3.d0)*p_lambda / 2.d0
              p_ens(isp)%xyz(3) = 0.d0
            else 
              p_ens(isp)%xyz(1) = (dble(ix) + 0.5d0)*p_lambda
              p_ens(isp)%xyz(2) = (dble(iy) + 0.5d0)*sqrt(3.d0)*p_lambda / 2.d0
              p_ens(isp)%xyz(3) = 0.d0    
            endif 
          enddo
        enddo

! triangular lattice - rhomboidal unit cell;

      elseif (p_urho) then 
        do ix=0,p_nx-1
          do iy=0,p_ny-1
            isp = isp + 1
            if (mod(iy,2) .eq. 1) then
              p_ens(isp)%xyz(1) = (dble(ix) + 1.d0)*p_lambda &
                                + dble(iy-1)*p_lambda / 2.d0 
              p_ens(isp)%xyz(2) = (dble(iy) + 0.5d0)*sqrt(3.d0)*p_lambda / 2.d0
              p_ens(isp)%xyz(3) = 0.d0
            else 
              p_ens(isp)%xyz(1) = (dble(ix) + 0.5d0)*p_lambda &
                                + dble(iy)*p_lambda / 2.d0 
              p_ens(isp)%xyz(2) = (dble(iy) + 0.5d0)*sqrt(3.d0)*p_lambda / 2.d0
              p_ens(isp)%xyz(3) = 0.d0    
            endif 
          enddo
        enddo
      else 

! raise critical error;

        call f_critical_('dpls2d',279)
      endif
    elseif (trim(p_lattice) .eq. 'honeycomb') then

! honeycomb lattice - rectangular unit cell;

      if (p_urec) then
        do ix=1,p_nx 
          xpos = dble(ix) * sqrt(3.d0) * p_lambda 
          ypos = p_lambda / 2.d0
          do iy=1,p_ny 
            isp = isp + 1
    
            if (mod(iy,4) .eq. 1) then
              ypos = ypos + p_lambda / 2.d0
              xpos = xpos - sqrt(3.d0) * p_lambda / 2.d0
            elseif (mod(iy,4) .eq. 2) then 
              ypos = ypos + p_lambda 
            elseif (mod(iy,4) .eq. 3) then 
              ypos = ypos + p_lambda / 2.d0
              xpos = xpos + sqrt(3.d0) * p_lambda /2.d0
            elseif (mod(iy,4) .eq. 0) then 
              ypos = ypos + p_lambda
            endif 
            p_ens(isp)%xyz = [ xpos, ypos, 0.d0]
          enddo
        enddo
      else 

! raise critical error;

        call f_critical_('dpls2d',310)
      endif
    elseif (trim(p_lattice) .eq. 'kagome') then

! kagome lattice - rhomboidal unit cell;

      if (p_urho) then 
        do ix=0,p_nx-1
          do iy=0,p_ny-1
    
    ! outside triangular lattice:
    
            xvec = dble(ix) * p_lambda * (/ 2.d0, 0.d0, 0.d0 /) &
                 + dble(iy) * p_lambda * (/ 1.d0, sqrt(3.d0), 0.d0 /) 
    
    ! inside triangular lattice:
    
            do iz=1,3
              if (iz .eq. 1) then 
                yvec = (/ 0.d0, 0.d0, 0.d0 /)
              elseif (iz .eq. 2) then
                yvec = p_lambda * (/ 1.d0, 0.d0, 0.d0 /)
              elseif (iz .eq. 3) then
                yvec = p_lambda * (/ 1.d0, sqrt(3.d0) , 0.d0 /) / 2.d0
              endif
              isp = isp + 1
              p_ens(isp)%xyz = xvec + yvec 
            enddo
          enddo
        enddo
      else 

! raise critical error;

        call f_critical_('dpls2d',344)
      endif
    endif

! create disorder - 1) classic way:

    do i=1,c_p 
      p_dbl(i)%xyz = p_ens(i)%xyz 
      p_dbl(i)%r = p_ens(i)%r
    enddo

    if (.not. p_rel_clc) then
      iout = 1
      overlap = .true.
      do while ((overlap) .and. (iout .lt. huge(ix)))
        iout = iout + 1
        do i=1,c_p
          p_dbl(i)%xyz = p_ens(i)%xyz
        enddo

  ! create offset:

        do i=1,c_p 
          angle = rng_number()*c_pi*2.d0
          radius = abs(rng_gauss(0.d0,p_lambda*p_sigma))
          p_dbl(i)%xyz = p_dbl(i)%xyz + radius*(/ cos(angle), sin(angle), 0.d0 /) 
        enddo

  ! check for overlap:

        overlap = .false.
        do ibx=-1,1
          do iby=-1,1
            do i=1,c_p 
              do j=1,c_p 
                if ((ibx .eq. 0) .and. (iby .eq. 0) .and. (i .eq. j)) then
                  cycle 
                endif
                xi = p_dbl(i)%xyz 
                xj = p_dbl(j)%xyz + [ dble(ibx)*dx , dble(iby)*dy,  0.d0 ]
                if (f_distance(xi,xj) .lt. p_ens(i)%r + p_ens(j)%r) then
                  overlap = .true.
                endif
              enddo
            enddo
          enddo
        enddo

        if (.not. overlap) then
          do i=1,c_p
            p_ens(i)%xyz = p_dbl(i)%xyz 
          enddo
          exit
        endif
      enddo

      if (overlap) then
        print *,''
        print *,'    dpls2d > no valid configuration could be created. aborted!'
        print *,''
        stop
      endif

! 2) we create initial configuration and relaxe with ball potential until there
! is no more overlap.

    elseif (p_rel_clc) then

! create initial configuration:

      do i=1,c_p 
        angle = rng_number()*c_pi*2.d0
        radius = abs(rng_gauss(0.d0,p_lambda*p_sigma))
        p_dbl(i)%xyz = p_dbl(i)%xyz + radius*(/ cos(angle), sin(angle), 0.d0 /)
      enddo

! update initial configuration by displacing particles until there is no 
! overlap. for this purpose, we use a simple ball potential:

      do k=1,50000  
        xgrd = 0.d0
        update = .false.
     
! calculate force:

        do i=1,c_p 
          do j=1,c_p 
            if (i .eq. j) then
              cycle 
            endif

! calculate distance between particles:
            
            rvec = p_dbl(j)%xyz - p_dbl(i)%xyz 
            rabs = sqrt(dot_product(rvec,rvec))

! if particles are not touching, there is no force:

            if (rabs .lt. p_dbl(i)%r + p_dbl(j)%r) then
              xgrd(i,1) = xgrd(i,1) - rvec(1)
              xgrd(i,2) = xgrd(i,2) - rvec(2)
              update = .true.
            endif
          enddo

! calculate colissions with unit cell boundaries:

          rvec = p_dbl(i)%xyz 

! x-bounds;

          if (p_urec) then
            if (rvec(1) .lt. p_dbl(i)%r) then
              xgrd(i,1) = xgrd(i,1) + 1.d0
              update = .true.
            elseif (rvec(1) .gt. dx - p_dbl(i)%r) then
              xgrd(i,1) = xgrd(i,1) - 1.d0
              update = .true.
            endif
          elseif (p_urho) then 
            if (rvec(1) .lt. p_dbl(i)%r + rvec(2)/sqrt(3.d0)) then 
              xgrd(i,1) = xgrd(i,1) + 1.d0
              update = .true.
            elseif (rvec(1) .gt. dx + rvec(2)/sqrt(3.d0) - p_dbl(i)%r) then 
              xgrd(i,1) = xgrd(i,1) - 1.d0
              update = .true.
            endif
          endif

! y-bounds:
          
          if (rvec(2) .lt. p_dbl(i)%r) then
            xgrd(i,2) = xgrd(i,2) + 1.d0
            update = .true.
          elseif (rvec(2) .gt. dy - p_dbl(i)%r) then
            xgrd(i,2) = xgrd(i,2) - 1.d0 
            update = .true.
          endif
        enddo

        if (update) then
          do i=1,c_p 
            p_dbl(i)%xyz(1) = p_dbl(i)%xyz(1) + xgrd(i,1)/500.d0
            p_dbl(i)%xyz(2) = p_dbl(i)%xyz(2) + xgrd(i,2)/500.d0
          enddo
        else 
          exit
        endif 
      enddo
    endif

    do i=1,c_p 
      p_ens(i)%xyz = p_dbl(i)%xyz 
    enddo

! create interaction arrays to increase efficiency of code:

    print *,''
    print *,'    > periodic boundaries..'
    print *,''

! allocate memory and initialize:

    allocate(xx3(c_p,c_p))
    allocate(xy3(c_p,c_p))
    allocate(yy3(c_p,c_p))
    xx3 = 0.d0
    xy3 = 0.d0
    yy3 = 0.d0

! store values in arrays:

    cpos(1) = dx / 2.d0
    cpos(2) = dy / 2.d0
    cpos(3) = 0.d0
    do ibx=-p_bnx,p_bnx
      do iby=-p_bny,p_bny
        if (p_urec) then
          cvec(1) = dble(ibx) * dx
          cvec(2) = dble(iby) * dy
          cvec(3) = 0.d0
        elseif (p_urho) then 
          cvec(1) = dble(ibx)*dx + dble(iby)*dx / 2.d0
          cvec(2) = dble(iby)*dy
          cvec(3) = 0.d0
        endif
        do i=1,c_p 
          do j=1,c_p 
            if (f_norm(p_ens(j)%xyz + cvec - cpos) .gt. p_maxd) then
              cycle 
            endif
            rvec =  p_ens(i)%xyz - p_ens(j)%xyz - cvec
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
    
! deallocate memory:

    deallocate(p_dbl)
    deallocate(xgrd)

! output:
    
    open(unit=11,file='ensemble.data',action='write',status='replace')
    do i=1,c_p 
      write(11,'(5(F15.5))') p_ens(i)%xyz,p_ens(i)%r,p_ens(i)%m
    enddo 
    close(11)
    
    return
  end subroutine dpls2d_init

!-------------------------------------------------------------------------------
! deinitialization:

  subroutine dpls2d_deinit()
    implicit none

    deallocate(p_ens)
    deallocate(xx3)
    deallocate(xy3)
    deallocate(yy3)

    return 
  end subroutine dpls2d_deinit

!-------------------------------------------------------------------------------
! energy and gradient of configuration x;

  subroutine dpls2d_egrad(x,e,g)
    implicit none
    double precision, dimension(c_n), intent(in)  :: x 
    double precision, intent(out)                 :: e 
    double precision, dimension(c_n), intent(out) :: g 
    double precision, dimension(c_n)              :: cosx
    double precision, dimension(c_n)              :: sinx
    double precision, dimension(3)                :: si
    double precision, dimension(3)                :: sj
    double precision, dimension(3)                :: dti
    double precision, dimension(3)                :: dpi
    double precision, dimension(3)                :: term
    integer                                       :: i
    integer                                       :: j
    integer                                       :: pi
    integer                                       :: pj

! initialize variables:

    sinx = sin(x)
    cosx = cos(x)

! calculate dipolar energy and gradient - the derivatives with respect to phi 
! are divided by sin(theta);

    e = 0.d0
    g = 0.d0
    do i=1,c_n,2
      pi = i/2 + 1
      si = p_ens(pi)%m * [ sinx(i)*cosx(i+1), sinx(i)*sinx(i+1), cosx(i) ]
      dti = p_ens(pi)%m * [ cosx(i)*cosx(i+1), cosx(i)*sinx(i+1),-sinx(i) ]
      dpi = p_ens(pi)%m * [-sinx(i+1), cosx(i+1), 0.d0 ]
      
      term = 0.d0
      do j=1,c_n,2
        pj = j/2 + 1
        sj = p_ens(pj)%m * [ sinx(j)*cosx(j+1), sinx(j)*sinx(j+1), cosx(j) ]

        term(1) = term(1) + yy3(pi,pj)*sj(1) - 2.d0*xx3(pi,pj)*sj(1) &
                - 3.d0*xy3(pi,pj)*sj(2)
        term(2) = term(2) + xx3(pi,pj)*sj(2) - 2.d0*yy3(pi,pj)*sj(2) &
                - 3.d0*xy3(pi,pj)*sj(1)
        term(3) = term(3) + xx3(pi,pj)*sj(3) + yy3(pi,pj)*sj(3)
      enddo
      e = e + c_mu0/(4.d0*c_pi) * dot_product(si,term)
      g(i+0) = c_mu0/(4.d0*c_pi) * dot_product(dti,term)
      g(i+1) = c_mu0/(4.d0*c_pi) * dot_product(dpi,term)
    enddo
    e = e / 2.d0

! calculate zeeman energy and gradient;

    do i=1,c_n,2 
      pi = i/2 + 1
      si = p_ens(pi)%m * [ sinx(i)*cosx(i+1), sinx(i)*sinx(i+1), cosx(i) ]
      dti = p_ens(pi)%m * [ cosx(i)*cosx(i+1), cosx(i)*sinx(i+1),-sinx(i) ]
      dpi = p_ens(pi)%m * [-sinx(i+1), cosx(i+1), 0.d0 ]     

      e = e - dot_product(p_hext,si) 
      g(i+0) = g(i+0) - dot_product(p_hext,dti)
      g(i+1) = g(i+1) - dot_product(p_hext,dpi)
    enddo

    return 
  end subroutine dpls2d_egrad

!-------------------------------------------------------------------------------
! hessian of configuration x;

  subroutine dpls2d_hess(x,h)
    implicit none 
    double precision, dimension(c_n), intent(in)      :: x 
    double precision, dimension(c_n,c_n), intent(out) :: h
    double precision, dimension(c_n)                  :: sinx
    double precision, dimension(c_n)                  :: cosx
    double precision, dimension(c_n)                  :: g
    double precision, dimension(3)                    :: sj
    double precision, dimension(3)                    :: dti
    double precision, dimension(3)                    :: dpi
    double precision, dimension(3)                    :: dtj
    double precision, dimension(3)                    :: dpj
    double precision, dimension(3)                    :: dt2
    double precision, dimension(3)                    :: dtp
    double precision, dimension(3)                    :: dp2
    double precision, dimension(6)                    :: oterm
    double precision, dimension(3)                    :: dterm
    double precision                                  :: e
    integer                                           :: i
    integer                                           :: j
    integer                                           :: pi
    integer                                           :: pj

! initialize variables:

    sinx = sin(x)
    cosx = cos(x)

! calculate hessian:

    h = 0.d0
    do i=1,c_n,2
      pi = i/2 + 1 
      dti = p_ens(pi)%m * [ cosx(i)*cosx(i+1), cosx(i)*sinx(i+1),-sinx(i) ]
      dt2 = p_ens(pi)%m * [-sinx(i)*cosx(i+1),-sinx(i)*sinx(i+1),-cosx(i) ]
      dpi = p_ens(pi)%m * [-sinx(i+1), cosx(i+1), 0.d0 ]
      dp2 = p_ens(pi)%m * [-sinx(i)*cosx(i+1),-sinx(i)*sinx(i+1), 0.d0 ]
      dtp = p_ens(pi)%m * [-cosx(i)*sinx(i+1), cosx(i)*cosx(i+1), 0.d0 ]

      dterm = 0.d0
      do j=1,c_n,2 
        pj = j/2 + 1
        sj = p_ens(pj)%m * [ sinx(j)*cosx(j+1), sinx(j)*sinx(j+1), cosx(j) ]
        dtj = p_ens(pj)%m * [-cosx(j)*sinx(j+1), cosx(j)*cosx(j+1),-sinx(j) ]
        dpj = p_ens(pj)%m * [-sinx(j+1), cosx(j+1), 0.d0 ]

! off-diagonal terms:

        oterm(1) = yy3(pi,pj)*dtj(1) - 2.d0*xx3(pi,pj)*dtj(1) &
                 - 3.d0*xy3(pi,pj)*dtj(2)
        oterm(2) = xx3(pi,pj)*dtj(2) - 2.d0*yy3(pi,pj)*dtj(2) &
                 - 3.d0*xy3(pi,pj)*dtj(1)
        oterm(3) = xx3(pi,pj)*dtj(3) + yy3(pi,pj)*dtj(3)
        oterm(4) = yy3(pi,pj)*dpj(1) - 2.d0*xx3(pi,pj)*dpj(1) &
                 - 3.d0*xy3(pi,pj)*dpj(2)
        oterm(5) = xx3(pi,pj)*dpj(2) - 2.d0*yy3(pi,pj)*dpj(2) &
                 - 3.d0*xy3(pi,pj)*dpj(1)
        oterm(6) = xx3(pi,pj)*dpj(3) + yy3(pi,pj)*dpj(3)

        h(i+0,j+0) = h(i+0,j+0) &
                   + c_mu0/(4.d0*c_pi) * dot_product(dti,oterm(1:3))
        h(i+0,j+1) = h(i+0,j+1) &
                   + c_mu0/(4.d0*c_pi) * dot_product(dti,oterm(4:6))
        h(i+1,j+0) = h(i+1,j+0) &
                   + c_mu0/(4.d0*c_pi) * dot_product(dpi,oterm(1:3))
        h(i+1,j+1) = h(i+1,j+1) &
                   + c_mu0/(4.d0*c_pi) * dot_product(dpi,oterm(4:6))

! diagonal terms:

        dterm(1) = dterm(1) + yy3(pi,pj)*sj(1) - 2.d0*xx3(pi,pj)*sj(1) &
                 - 3.d0*xy3(pi,pj)*sj(2)
        dterm(2) = dterm(2) + xx3(pi,pj)*sj(2) - 2.d0*yy3(pi,pj)*sj(2) &
                 - 3.d0*xy3(pi,pj)*sj(1)
        dterm(3) = dterm(3) + xx3(pi,pj)*sj(3) + yy3(pi,pj)*sj(3)
      enddo

      h(i+0,i+0) = h(i+0,i+0) + c_mu0/(4.d0*c_pi) * dot_product(dt2,dterm)
      h(i+0,i+1) = h(i+0,i+1) + c_mu0/(4.d0*c_pi) * dot_product(dtp,dterm)
      h(i+1,i+0) = h(i+1,i+0) + c_mu0/(4.d0*c_pi) * dot_product(dtp,dterm)
      h(i+1,i+1) = h(i+1,i+1) + c_mu0/(4.d0*c_pi) * dot_product(dp2,dterm)
    enddo

! zeeman terms:

    do i=1,c_n,2       
      dt2 = p_ens(pi)%m * [-sinx(i)*cosx(i+1),-sinx(i)*sinx(i+1),-cosx(i) ]
      dp2 = p_ens(pi)%m * [-sinx(i)*cosx(i+1),-sinx(i)*sinx(i+1), 0.d0 ]
      dtp = p_ens(pi)%m * [-cosx(i)*sinx(i+1), cosx(i)*cosx(i+1), 0.d0 ]

      h(i+0,i+0) = h(i+0,i+0) - dot_product(dt2,p_hext)
      h(i+0,i+1) = h(i+0,i+1) - dot_product(dtp,p_hext)
      h(i+1,i+0) = h(i+1,i+0) - dot_product(dtp,p_hext)
      h(i+1,i+1) = h(i+1,i+1) - dot_product(dp2,p_hext)
    enddo

! correct hessian according to [Gra11]]:

    call dpls2d_egrad(x,e,g)
    do i=1,c_n,2
      h(i+0,i+0) = h(i+0,i+0)
      h(i+0,i+1) = h(i+0,i+1)/sinx(i) - cosx(i)/sinx(i)*g(i+1)
      h(i+1,i+0) = h(i+1,i+0)/sinx(i) - cosx(i)/sinx(i)*g(i+1)
      h(i+1,i+1) = h(i+1,i+1)/sinx(i)/sinx(i) + cosx(i)/sinx(i)*g(i+0)
    enddo

    return
  end subroutine dpls2d_hess

!-------------------------------------------------------------------------------
! rescaling of each element of vector x into ([0,pi],[0,2pi]):

  subroutine dpls2d_rescale(x)
    implicit none
    double precision, dimension(c_n), intent(inout) :: x
    integer                                         :: i

    do i=1,c_n,2 
      do while (x(i) .lt. 0.d0)
        x(i) = x(i) + 2.d0*c_pi
      enddo
      
      do while (x(i) .gt. 2.d0*c_pi)
        x(i) = x(i) - 2.d0*c_pi 
      enddo

      if (x(i) .gt. c_pi) then
        x(i) = 2.d0*c_pi - x(i) 
        x(i+1) = x(i+1) + c_pi 
      endif

      do while(x(i+1) .lt. 0.d0)
        x(i+1) = x(i+1) + 2.d0*c_pi 
      enddo

      do while (x(i+1) .gt. 2.d0*c_pi)
        x(i+1) = x(i+1) - 2.d0*c_pi 
      enddo
    enddo

    return 
  end subroutine dpls2d_rescale

!-------------------------------------------------------------------------------
! create random configuration:

  subroutine dpls2d_random(x)
    implicit none
    double precision, dimension(c_n)  :: x 
    integer                           :: i

    do i=1,c_n,2
      x(i) = rng_number() * c_pi
      x(i+1) = rng_number() * c_pi * 2.d0
    enddo
    
    return 
  end subroutine dpls2d_random

!-------------------------------------------------------------------------------
! calculate distance between configurations x1 and x2;

  function dpls2d_distance(x1,x2)
    implicit none 
    double precision                  :: dpls2d_distance
    double precision, dimension(c_n)  :: x1
    double precision, dimension(c_n)  :: x2
    double precision                  :: dx
    integer                           :: i

    dpls2d_distance = 0.d0
    do i=1,c_n,2
      
! calculate dot product between x1 and x2;

      dx = dot_product(x1(i:i+1),x2(i:i+1))

! exceptions due to possible numerical inaccuracies such that arccosine is 
! always defined;

      dx = min(dx,1.d0)
      dx = max(dx,0.d0)

! calculate arccosine;

      dx = acos(dx)
      dpls2d_distance = dpls2d_distance + dx*dx
    enddo
    dpls2d_distance = sqrt(dpls2d_distance)

    return 
  end function dpls2d_distance

!-------------------------------------------------------------------------------
! find state x2 having the opposite magnetization of state x1;

  subroutine dpls2d_reversal(x1,x2)
    implicit none 
    double precision, dimension(c_n), intent(in)  :: x1
    double precision, dimension(c_n), intent(out) :: x2
    integer                                       :: i

    x2 = x1
    do i=1,c_n,2
      x2(i) = x2(i) + c_pi 
    enddo
    call dpls2d_rescale(x2)

    return 
  end subroutine dpls2d_reversal

!-------------------------------------------------------------------------------
! vibrational entropy:

  subroutine dpls2d_entropy(x,s,nev)
    implicit none 
    double precision, dimension(c_n), intent(in)  :: x 
    double precision, intent(out)                 :: s 
    integer, intent(out)                          :: nev
    double precision, dimension(c_n,c_n)          :: h 
    double precision, dimension(c_n,c_n)          :: evec
    double precision, dimension(c_n)              :: eval 
    integer                                       :: i
    integer                                       :: nzero

    call dpls2d_hess(x,h)
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
  end subroutine dpls2d_entropy

!-------------------------------------------------------------------------------
! order parameters - see [Bae11],[Hol15];

  subroutine dpls2d_order(x,order,angle)
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

    if (trim(p_lattice) .eq. 'square') then      
      
! create gauge vector:

      y = x 
      i = 0
      do ix=1,p_nx 
        do iy=1,p_ny 
          i = i + 2
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
  
      call dpls2d_rescale(y)
  
  ! calculate order parameter [0,1]:
  
      vec = 0.d0
      do i=2,c_n,2
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

    elseif (trim(p_lattice) .eq. 'triangular') then
      
! calculate magnetization;

      vec = 0.d0
      do i=1,c_n,2 
        vec(1) = vec(1) + sin(x(i))*cos(x(i+1))
        vec(2) = vec(2) + sin(x(i))*sin(x(i+1))
      enddo
      order = sqrt(vec(1)*vec(1) + vec(2)*vec(2)) / dble(c_n)
  
  ! calculate angle;
  
      vec = vec / sqrt(dot_product(vec,vec))
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
  
          i = i + 2
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
  
      call dpls2d_rescale(y)
        
! create shifted vector:

      do i=1,c_n,2
        y(i) = c_pi/2.d0
        if (y(i+1) .gt. c_pi) then 
          yt(i+1) = y(i+1) - 2.d0*c_pi 
        else 
          yt(i+1) = y(i+1)
        endif
      enddo

! calculate order parameter [0,1]:

      vec = 0.d0
      do i=2,c_n,2
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
      
! kagome lattice;

    elseif (trim(p_lattice) .eq. 'kagome') then

! create 3-fold magnetization vector;

      svec = 0.d0
      do i=1,c_n,6
        svec(1) = svec(1) + cos(x(i+1))
        svec(2) = svec(2) + sin(x(i+1)) 
        svec(3) = svec(3) + cos(x(i+3))
        svec(4) = svec(4) + sin(x(i+3))
        svec(5) = svec(5) + cos(x(i+5))
        svec(6) = svec(6) + sin(x(i+6)) 
      enddo

! calculate order for each sublattice;

      avec(1) = sqrt(dot_product(svec(1:2),svec(1:2)))
      avec(2) = sqrt(dot_product(svec(3:4),svec(3:4)))
      avec(3) = sqrt(dot_product(svec(5:6),svec(5:6)))

      order = (avec(1) + avec(2) + avec(3)) / dble(c_p)
      angle = 0.d0
    endif

    return
  end subroutine dpls2d_order

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! read parameters:

  subroutine dpls2d_readparms()
    implicit none 
    integer   :: info

! lattice parameters:

    call f_read_store_('lattice_type','params.dat',p_lattice,'square',info)
    call f_read_store_('lattice_urec','params.dat',p_urec,.false.,info)
    call f_read_store_('lattice_urho','params.dat',p_urho,.false.,info)
    if (trim(p_lattice) .eq. 'square') then
      if ((.not. p_urec) .and. (p_urho)) then 
        p_urec = .true.
        p_urho = .false.
        print *,''
        print *,'    > warning: urho option not available for square lattices.'
        print *,'               switched to urec option.' 
        print *,''
      elseif ((.not. p_urec) .and. (.not. p_urho)) then
        p_urec = .true.
        print *,''
        print *,'    > warning: neither urho option nor urec option chosen.'
        print *,'               switched to urec option.' 
        print *,''
      elseif ((p_urec) .and. (p_urho)) then
        p_urho = .false.
        print *,''
        print *,'    > warning: both urho and urec option chosen. switched to'
        print *,'               urec option.' 
        print *,''
      endif 
    elseif (trim(p_lattice) .eq. 'triangular') then
      if ((.not. p_urec) .and. (.not. p_urho)) then
        p_urho = .true.
        print *,''
        print *,'    > warning: neither urho option nor urec option chosen.'
        print *,'               switched to urho option.' 
        print *,''
      elseif ((p_urec) .and. (p_urho)) then
        p_urec = .false.
        print *,''
        print *,'    > warning: both urho and urec option chosen. switched to'
        print *,'               urho option.' 
        print *,''
      endif 
    elseif (trim(p_lattice) .eq. 'honeycomb') then
      if ((.not. p_urec) .and. (p_urho)) then 
        p_urec = .true.
        p_urho = .false.
        print *,''
        print *,'    > warning: urho option not available for honeycomb lattices.'
        print *,'               switched to urec option.' 
        print *,''
      elseif ((.not. p_urec) .and. (.not. p_urho)) then
        p_urec = .true.
        print *,''
        print *,'    > warning: neither urho option nor urec option chosen.'
        print *,'               switched to urec option.' 
        print *,''
      elseif ((p_urec) .and. (p_urho)) then
        p_urho = .false.
        print *,''
        print *,'    > warning: both urho and urec option chosen. switched to'
        print *,'               urec option.' 
        print *,''
      endif 
    elseif (trim(p_lattice) .eq. 'kagome') then
      if ((p_urec) .and. (.not. p_urho)) then 
        p_urec = .false.
        p_urho = .true.
        print *,''
        print *,'    > warning: urec option not available for kagome lattices.'
        print *,'               switched to urec option.' 
        print *,''
      elseif ((.not. p_urec) .and. (.not. p_urho)) then
        p_urho = .true.
        print *,''
        print *,'    > warning: neither urho option nor urec option chosen.'
        print *,'               switched to urec option.' 
        print *,''
      elseif ((p_urec) .and. (p_urho)) then
        p_urec = .false.
        print *,''
        print *,'    > warning: both urho and urec option chosen. switched to'
        print *,'               urec option.' 
        print *,''
      endif 
    endif

    call f_read_store_('lattice_constant','params.dat',p_lambda,10.d0,info)
    call f_read_store_('lattice_disorder','params.dat',p_sigma,0.05d0,info)
    call f_read_store_('lattice_xrows','params.dat',p_nx,8,info)
    call f_read_store_('lattice_ycols','params.dat',p_ny,8,info)
    call f_read_store_('lattice_xperiod','params.dat',p_bnx,0,info)
    call f_read_store_('lattice_yperiod','params.dat',p_bny,0,info)

! particle parameters:

    call f_read_store_('particle_radius','params.dat',p_radius,4.d0,info)
    call f_read_store_('particle_moment','params.dat',p_mag,180.d0,info)

! field parameters:

    p_hext = 0.d0
    call f_read_('field_vector','params.dat',p_hext,info)

! system parameters:

    call f_read_store_('particle_relaxe','params.dat',p_rel_clc,.false.,info)

    return 
  end subroutine dpls2d_readparms 

!-------------------------------------------------------------------------------

end module dpls2d
