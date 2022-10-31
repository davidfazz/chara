! explore.f90 - elapid - david gallina - dg@physik.uni-kassel.de - 2022
!-------------------------------------------------------------------------------
! algorithm to create a connected network of local minima and adjacent first-
! order saddle points for a given potential energy surface. more information on
! the routine can be found in the faq.
!-------------------------------------------------------------------------------
! double parameters:
! p_perturb - perturbation in random ef searches.
! p_offset - offset at the saddle point along the single negative eigenvalue.
! p_deps - accuracy of ef and lbfgs searches - rmsgrad < p_deps.
! p_dvec - criterion to check if two states are different - based on configs.
! p_dnrg - criterion to check if two states are different - based on energy.
!
! integer parameters:
! p_nmodes - number of eigenmodes to follow in ef searches.
! p_nrands - number of random directions to follow in ef searches.
! p_npool - pool of local minima to start the algorithm.
! p_iter - number of iterations.
! p_nm - maximum number of storable local minima.
! p_ns - maximum number of storable saddle points.
!
! boolean parameters:
! p_timerev_add - add states found through time-reversal symmetry.
! p_print_rates - print transition rates obtained through TST in 'rates.data'.
! p_print_order - print order parameters of all local minima in 'order.data'
! p_print_cfgs - print configurations of all states in '*.pts' files.
!-------------------------------------------------------------------------------

module explore
  use constants 
  use functions
  use potential
  use approx
  use lbfgs
  use efol
  use omp_lib
  implicit none

!-------------------------------------------------------------------------------
! explore parameters:

  double precision, private   :: p_perturb
  double precision, private   :: p_deps   
  double precision, private   :: p_dvec   
  double precision, private   :: p_dnrg   
  double precision, private   :: p_offset  
  double precision, private   :: p_alpha  
  integer, private            :: p_nmodes  
  integer, private            :: p_nrands  
  integer, private            :: p_npool  
  integer, private            :: p_iter    
  integer, private            :: p_nm    
  integer, private            :: p_ns     
  logical, private            :: p_timerev_add
  logical, private            :: p_print_rates
  logical, private            :: p_print_order
  logical, private            :: p_print_cfgs
  logical, private            :: p_contdeg = .true.

!-------------------------------------------------------------------------------

  type t_lm
    double precision, dimension(:), allocatable  :: x      ! magnetic configuration
    double precision                :: order
    double precision                :: angle
    double precision                :: m      ! magnetic order
    double precision                :: v      ! vortex order
    double precision                :: e      ! energy
    double precision                :: s      ! vibrational entropy
    integer                         :: pos    ! position
  end type t_lm

  type t_ts
    double precision, dimension(:), allocatable  :: x      ! magnetic configuration
    double precision                :: e                ! energy
    double precision                :: s      ! vibrational entropy
    double precision, dimension(2)  :: rates  ! rates
    integer, dimension(2)           :: trip   ! lm-ts-lm triplet
  end type t_ts

!-------------------------------------------------------------------------------
! storage arrays:

  type(t_lm), dimension(:), allocatable  :: mvec
  type(t_ts), dimension(:), allocatable  :: svec

! control variables:

  integer                                       :: nm
  integer                                       :: ns

contains

!-------------------------------------------------------------------------------

  subroutine explore_()
    implicit none
    double precision, dimension(c_n,c_n)  :: hess
    double precision, dimension(c_n)  :: xopt
    double precision, dimension(c_n)  :: xone
    double precision, dimension(c_n)  :: xts
    double precision, dimension(c_n)  :: xtwo
    double precision                  :: energy
    double precision                  :: rone
    double precision                  :: rtwo
    integer                           :: i
    integer                           :: j
    integer                           :: pos
    integer                           :: iter
    logical                           :: conv
    
! read parameters:

    call explore_readparms()

! output:

    print *,''
    print *,'    explore > max.lm:   ',trim(f_trimmed_(p_nm))
    print *,'            > max.ts:   ',trim(f_trimmed_(p_ns))

! allocate minima and saddle point arrays:

    allocate(mvec(p_nm))
    do i=1,p_nm
      allocate(mvec(i)%x(c_n))
    enddo
    nm = 0

    allocate(svec(p_ns))
    do i=1,p_ns
      allocate(svec(i)%x(c_n))
    enddo
    ns = 0

! create local minima pool:

    do i=1,p_npool
      call potential_random(xopt)
      call lbfgs_(c_n,5,xopt,p_deps,conv,energy,5000,j,.false.)
      if (.not. conv) then
        cycle
      endif
      call explore_add_lm(xopt,pos) 
    enddo

! start of the main algorithm:

    iter = 1
    do while (iter .le. nm)
      if (mod(iter-1,p_iter) .eq. 0) then
        print *,''
        print *,'    explore > explored lm:   ',trim(f_trimmed_(iter))
        print *,'            > total lm:      ',trim(f_trimmed_(nm))
      endif

      call explore_search(mvec(iter)%x)
      iter = iter + 1
    enddo

! calculate rates:

    if (p_print_rates) then 

!$omp parallel do default(shared) private(xone,xtwo,xts,rone,rtwo)

      do i=1,ns 
        xone = mvec(svec(i)%trip(1))%x
        xtwo = mvec(svec(i)%trip(2))%x
        xts = svec(i)%x 

        call potential_rates(xone,xts,xtwo,p_alpha,rone,rtwo)
        svec(i)%rates(1) = rone 
        svec(i)%rates(2) = rtwo
      enddo

!$omp end parallel do

    endif

! output:

    open(unit=15,file='min.data',action='write',status='replace')
    do i=1,nm
      write(15,*) mvec(i)%e,mvec(i)%s,i,mvec(i)%order,mvec(i)%angle,0
    enddo
    close(15)

    open(unit=15,file='ts.data',action='write',status='replace')
    do i=1,ns
      write(15,*) svec(i)%e,svec(i)%s,0,svec(i)%trip,0,0,0
    enddo
    close(15)

    if (p_print_rates) then 
      open(unit=15,file='rates.data',action='write',status='replace')
      do i=1,ns
        write(15,*) svec(i)%e,svec(i)%s,svec(i)%trip,svec(i)%rates
      enddo
      close(15)
    endif

    if (p_print_order) then
      open(unit=15,file='order.data',action='write',status='replace')
      do i=1,nm
        write(15,*) mvec(i)%e,mvec(i)%s,mvec(i)%m,mvec(i)%v
      enddo
      close(15)
    endif

    if (p_print_cfgs) then 
      open(unit=15,file='min.pts',action='write',status='replace')
      do i=1,nm
        write(15,*) mvec(i)%x
      enddo
      close(15)

      open(unit=15,file='ts.pts',action='write',status='replace')
      do i=1,ns
        write(15,*) svec(i)%x
      enddo
      close(15)
    endif

! deallocate:

    deallocate(mvec)
    deallocate(svec)

    return
  end subroutine explore_

!...............................................................................
! routine to add a local minimum to the current database of local minima:

  subroutine explore_add_lm(x,pos)
    implicit none
    double precision, dimension(c_n), intent(in)  :: x
    integer, intent(out)                          :: pos
    double precision, dimension(c_n)              :: gx
    double precision, dimension(c_n)              :: gy    
    double precision, dimension(c_n)              :: y
    double precision                              :: sx
    double precision                              :: sy
    double precision                              :: ex
    double precision                              :: ey
    double precision                              :: order
    double precision                              :: angle
    integer                                       :: i
    integer                                       :: nev

! at first, we calculate the energy and vibrational entropy:

    call potential_(x,ex,gx)
    call potential_entropy(x,sx,nev)
    call potential_order(x,order,angle)
    if (nev .gt. 0) then
      pos = -1
      return
    endif

! now, we search through the minima array. if continuous degeneracy is set, than 
! continiuously degenerate states are grouped as one state.
    
    pos = 0
    if (p_contdeg) then
      if (nm .gt. 0) then
        do i=1,nm 
          if (abs(mvec(i)%e - ex) .lt. p_dnrg) then
            if (abs(order) .gt. 0.99d0) then
              pos = i
              exit  
            else
              if (potential_distance(mvec(i)%x,x) .lt. p_dvec) then   
                pos = i 
                exit
              endif 
            endif 
          endif
        enddo
      endif

! no continuous degeneracy:

    elseif (.not. p_contdeg) then
      if (nm .gt. 0) then
        do i=1,nm 
          if (abs(mvec(i)%e - ex) .lt. p_dnrg) then
            if (potential_distance(mvec(i)%x,x) .lt. p_dvec) then   
              pos = i 
              exit 
            endif
          endif
        enddo
      endif
    endif

    if (pos .eq. 0) then
      if (nm + 1 .le. p_nm) then
        nm = nm + 1
        mvec(nm)%x = x 
        mvec(nm)%e = ex 
        mvec(nm)%s = sx 
        mvec(nm)%pos = nm 
        mvec(nm)%order = order
        mvec(nm)%angle = angle
        if (p_timerev_add) then 
          if (nm + 1 .le. p_nm) then
            nm = nm + 1
            call potential_reversal(x,y)
            call potential_(y,ey,gy)
            call potential_entropy(y,sy,nev)
            call potential_order(y,order,angle)
            mvec(nm)%x = y 
            mvec(nm)%e = ey 
            mvec(nm)%s = sy
            mvec(nm)%pos = nm 
            mvec(nm)%order = order
            mvec(nm)%angle = angle
          else 
            pos =-1
          endif
        endif
      else 
        pos =-1
      endif
    endif

    return
  end subroutine explore_add_lm

!...............................................................................
! routine to add a first-order saddle point to the current database:

  subroutine explore_add_ts(x,evec)
    implicit none
    double precision, dimension(c_n), intent(in) :: x
    double precision, dimension(c_n), intent(in) :: evec
    double precision, dimension(c_n)             :: y
    double precision, dimension(c_n)             :: gx
    double precision, dimension(c_n)             :: xp
    double precision, dimension(c_n)             :: xm
    double precision                             :: sx
    double precision                             :: ex
    double precision                             :: ep
    double precision                             :: em
    integer                                      :: i
    integer                                      :: j
    integer                                      :: mpos
    integer                                      :: ppos
    integer                                      :: nev
    logical                                      :: conv

! at first, we calculate the energy and moment:

    call potential_(x,ex,gx)
    if (ns .gt. 0) then
      do i=1,ns 
        if (abs(svec(i)%e - ex) .lt. p_dnrg) then
          if (potential_distance(svec(i)%x,x) .lt. p_dvec) then
            return 
          endif 
        endif 
      enddo
    endif

! check if array is full:

    if (ns+1 .gt. p_ns) then
      return
    endif

! if the state is new, we calculate its vibrational entropy:

    call potential_entropy(x,sx,nev)

! we calculate the adjacent minima:

    xp = x + p_offset * evec
    call lbfgs_(c_n,5,xp,p_deps,conv,ep,5000,j,.false.)
    if (.not. conv) then
      return
    endif
    call explore_add_lm(xp,ppos)

    xm = x - p_offset * evec
    call lbfgs_(c_n,5,xm,p_deps,conv,em,5000,j,.false.)
    if (.not. conv) then
      return
    endif
    call explore_add_lm(xm,mpos)

    if ((mpos .le. 0) .or. (ppos .le. 0)) then
      return
    endif

! we update the transition state array:

    ns = ns + 1
    svec(ns)%x = x
    svec(ns)%e = ex
    svec(ns)%s = sx
    svec(ns)%trip = (/ mpos, ppos /)

! time-reversal symmetry:

    if (p_timerev_add) then
      if (ns + 1 .le. p_ns) then
        call potential_reversal(x,y)

! find position of degenerate minima:

        if (mod(ppos,2) .eq. 0) then
          ppos = ppos - 1
        else 
          ppos = ppos + 1
        endif
        if (mod(mpos,2) .eq. 0) then
          mpos = mpos - 1
        else 
          mpos = mpos + 1
        endif

! add degenerate saddle point:
    
        ns = ns + 1
        svec(ns)%x = y
        svec(ns)%e = ex
        svec(ns)%s = sx 
        svec(ns)%trip = (/ mpos, ppos /)
      endif
    endif

    return
  end subroutine explore_add_ts

!...............................................................................
! routine to find first-order saddle points from a given local minimum 'x'. the 
! loop is optimized for openmp use.

  subroutine explore_search(x)
    implicit none
    double precision, dimension(c_n), intent(in) :: x
    double precision, dimension(c_n)             :: xopt
    double precision, dimension(c_n)             :: evec
    double precision                             :: prod
    integer                                      :: i
    integer                                      :: mode
    logical                                      :: conv

! searching along the modes of the hessian of the local minimum:

!$omp parallel default(shared) private(mode,xopt,evec,prod,conv)
!$omp do
    
    do i=1,2*p_nmodes
      if (i .le. p_nmodes) then
        mode =-i
      else 
        mode = i - p_nmodes 
      endif

      xopt = x
      call efol_(xopt,c_n,mode,evec,prod,conv)
      if (.not. conv) then
        cycle 
      endif

!$omp critical 
!$omp flush(ns)

      call explore_add_ts(xopt,evec)

!$omp end critical
    
    enddo

!$omp end do
!$omp end parallel

    return
  end subroutine explore_search

!...............................................................................

  
!...............................................................................

  subroutine explore_readparms()
    implicit none
    integer           :: info

! double parameters:

    call f_read_store_('explore_deps','params.dat',p_deps,1.d-7,info) 
    call f_read_store_('explore_dvec','params.dat',p_dvec,1.d-1,info)
    call f_read_store_('explore_dnrg','params.dat',p_dnrg,1.d-4,info)
    call f_read_store_('explore_ef_perturb','params.dat',p_perturb,1.d0,info)
    call f_read_store_('explore_alpha','params.dat',p_alpha,1.d-1,info)
    call f_read_store_('explore_offset','params.dat',p_offset,1.d-2,info)

! integer parameters:

    call f_read_store_('explore_max_pool','params.dat',p_npool,1000,info)
    call f_read_store_('explore_ef_modes','params.dat',p_nmodes,5,info)
    call f_read_store_('explore_ef_rands','params.dat',p_nrands,10,info)
    call f_read_store_('explore_max_lm','params.dat',p_nm,5000,info)
    call f_read_store_('explore_max_ts','params.dat',p_ns,20000,info)
    call f_read_store_('explore_output','params.dat',p_iter,500,info)

! boolean parameters:

    call f_read_store_('explore_timerev_add','params.dat',p_timerev_add,.false.,info)
    call f_read_store_('explore_print_cfgs','params.dat',p_print_cfgs,.false.,info)
    call f_read_store_('explore_print_rates','params.dat',p_print_rates,.false.,info)
    call f_read_store_('explore_print_order','params.dat',p_print_order,.false.,info)

    return
  end subroutine explore_readparms

!...............................................................................

end module explore





