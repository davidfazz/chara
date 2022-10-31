! functions.f90 - elapid - david gallina - university of kassel - 2022
!-------------------------------------------------------------------------------
! module containts a number of helpful functions and subroutines.
!-------------------------------------------------------------------------------

module functions
  use constants 
  implicit none 

!-------------------------------------------------------------------------------
! interface for f_read_ function:

  interface f_read_
    procedure f_read_int
    procedure f_read_dbl
    procedure f_read_chr
    procedure f_read_log
    procedure f_read_vec
  end interface f_read_

! interface for f_read_store_ function:

  interface f_read_store_
    procedure f_read_store_int
    procedure f_read_store_dbl
    procedure f_read_store_chr
    procedure f_read_store_log
    procedure f_read_store_vec
  end interface f_read_store_  

! interface for f_trimmed_ function:

  interface f_trimmed_
    procedure f_trimmed_int
    procedure f_trimmed_dbl
    procedure f_trimmed_chr
    procedure f_trimmed_log
    procedure f_trimmed_vec
  end interface f_trimmed_ 

contains 

!-------------------------------------------------------------------------------
! MATHEMATICAL routines
!-------------------------------------------------------------------------------
! outer vector product between two 3d vectors x and y:

  function f_outer(x,y)
    double precision, dimension(3)  :: f_outer 
    double precision, dimension(3)  :: x 
    double precision, dimension(3)  :: y 

    f_outer(1) = x(2)*y(3) - x(3)*y(2)
    f_outer(2) = x(3)*y(1) - x(1)*y(3)
    f_outer(3) = x(1)*y(2) - x(2)*y(1)

    return 
  end function f_outer 

!-------------------------------------------------------------------------------
! inner vector product between two vectors x and y:

  function f_inner(x,y)
    double precision                :: f_inner 
    double precision, dimension(:)  :: x 
    double precision, dimension(:)  :: y 
    integer                         :: i 

    f_inner = 0.d0
    do i=1,size(x)
      f_inner = f_inner + x(i)*y(i)
    enddo
    
    return 
  end function f_inner 

!------------------------------------------------------------------------------- 
! root mean square of an input vector x:

  function f_rms(x)
    implicit none
    double precision                :: f_rms
    double precision, dimension(:)  :: x

    f_rms = f_inner(x,x)
    f_rms = sqrt(f_rms/dble(size(x)))

    return
  end function f_rms

!------------------------------------------------------------------------------- 
! norm of an input vector x:

  function f_norm(x)
    implicit none 
    double precision                :: f_norm
    double precision, dimension(:)  :: x 

    f_norm = sqrt(f_inner(x,x))

    return 
  end function f_norm

!...............................................................................  
! euclidean distance between two vectors x and y:

  function f_distance(x,y)
    implicit none 
    double precision                :: f_distance
    double precision, dimension(:)  :: x
    double precision, dimension(:)  :: y

    f_distance = sqrt(f_inner(x-y,x-y))
    
    return 
  end function f_distance 

!-------------------------------------------------------------------------------
! transforms unit vector from spherical to euclidean coordinates:

  function f_sph2car(x)
    implicit none
    double precision, dimension(3)  :: f_sph2car 
    double precision, dimension(2)  :: x 

    f_sph2car(1) = sin(x(1))*cos(x(2))
    f_sph2car(2) = sin(x(1))*sin(x(2))
    f_sph2car(3) = cos(x(1))

    return 
  end function f_sph2car 

!-------------------------------------------------------------------------------
! convert unit vector from cartesian coordinates to spherical coordinates:

  function f_car2sph(x)
    implicit none
    double precision, dimension(2)  :: f_car2sph
    double precision, dimension(3)  :: x 

    f_car2sph(1) = acos(x(3))
    f_car2sph(2) = atan2(x(2),x(1))

    return 
  end function f_car2sph

!-------------------------------------------------------------------------------  
! convert unit vector from polar coordinates to cartesian coordinates:

  function f_pol2car(x)
    implicit none
    double precision, dimension(2)  :: f_pol2car 
    double precision                :: x 

    f_pol2car(1) = cos(x)
    f_pol2car(2) = sin(x)

    return 
  end function f_pol2car 

!------------------------------------------------------------------------------- 
! convert cartesian coordinates to polar coordinates:

  function f_car2pol(x)
    implicit none
    double precision                :: f_car2pol 
    double precision, dimension(2)  :: x 

    f_car2pol = atan2(x(2),x(1))

    return 
  end function f_car2pol

!-------------------------------------------------------------------------------
! angle between two cartesian vectors:

  function f_angle(x,y)
    implicit none
    double precision                :: f_angle
    double precision, dimension(3)  :: x 
    double precision, dimension(3)  :: y 
    
    f_angle = acos(dot_product(x/f_norm(x),y/f_norm(y)))

    return 
  end function f_angle

!------------------------------------------------------------------------------- 
! angle between cartesian vectors based on vincentys formula:

  function f_angle_vin(x,y)
    implicit none
    double precision                :: f_angle_vin
    double precision, dimension(3)  :: x 
    double precision, dimension(3)  :: y 

    f_angle_vin = atan2(f_norm(f_outer(x,y)),dot_product(x,y))

    return 
  end function f_angle_vin 

!-------------------------------------------------------------------------------
! eigenvalues and eigenvectors of a symmetric matrix hsym:

  subroutine f_dsyev(hsym,evec,eval)
    implicit none
    double precision, dimension(:,:), intent(in)    :: hsym
    double precision, dimension(:,:), intent(out)   :: evec
    double precision, dimension(:), intent(out)     :: eval
    double precision, dimension(:,:), allocatable   :: dhsym
    double precision, dimension(:), allocatable     :: deval
    double precision, dimension(:), allocatable     :: work
    integer                                         :: lwork
    integer                                         :: info
    integer                                         :: n 
    integer                                         :: i

  ! dimensionality of input matrix and size of work array:

    n = size(hsym,1)
    lwork = 33*n

  ! allocate arrays:

    allocate(dhsym(n,n))
    allocate(deval(n))
    allocate(work(lwork))
    dhsym = hsym
    deval = 0.d0
    work = 0.d0

  ! call dsyev routine:

    call dsyev('V','U',n,dhsym,n,deval,work,lwork,info)
    if (info .ne. 0) then
      print *,''
      print *,'    dsyev > diagonalisation failed!'
      print *,''
      return
    endif

  ! store eigenvalues and eigenvectors:  

    do i=1,n 
      eval(i) = deval(i)
      evec(i,:) = dhsym(:,i)
    enddo

    deallocate(dhsym)
    deallocate(deval)
    deallocate(work)

    return
  end subroutine f_dsyev  

!-------------------------------------------------------------------------------
! eigenvectors and eigenvalues of an antisymmetric matrix:

  subroutine f_dgeev(asym,evec,eval)
    implicit none 
    double precision, dimension(:,:), intent(in)  :: asym 
    double precision, dimension(:,:), intent(out) :: evec
    double precision, dimension(:), intent(out)   :: eval 
    double precision, dimension(:,:), allocatable :: dasym  
    double precision, dimension(:,:), allocatable :: devec
    double precision, dimension(:), allocatable   :: deval_r 
    double precision, dimension(:), allocatable   :: deval_i
    double precision, dimension(:), allocatable   :: dpwork
    double precision                              :: dp
    integer                                       :: n
    integer                                       :: i
    integer                                       :: info
    integer                                       :: lwork

! dimensionality of input matrix and size of work array:

    n = size(asym,1)
    lwork = 33*n

  ! allocate arrays:

    allocate(dasym(n,n))
    allocate(devec(n,n))
    allocate(deval_r(n))
    allocate(deval_i(n))
    allocate(dpwork(lwork))
    dasym = asym 
    devec = 0.d0
    deval_r = 0.d0
    deval_i = 0.d0
    dpwork = 0.d0

! call dgeev routine:

    call dgeev('N','V',n,dasym,n,deval_r,deval_i,dp,1,devec,n,dpwork,lwork,info)
    if (info .ne. 0) then
      print *,''
      print *,'    dgeev > diagonalization failed!'
      print *,''
      return
    endif

    do i=1,n 
      eval(i) = deval_r(i)
      evec(i,:) = devec(:,i)
    enddo

    deallocate(dasym)
    deallocate(devec)
    deallocate(deval_r)
    deallocate(deval_i)
    deallocate(dpwork)

  return
  end subroutine f_dgeev

!-------------------------------------------------------------------------------

  pure function f_mean(x)
    implicit none 
    double precision, dimension(:), intent(in)  :: x
    double precision                            :: f_mean 
    integer                                     :: i 

    f_mean = 0.d0 
    do i=1,size(x)
      f_mean = f_mean + x(i)
    enddo
    f_mean = f_mean / dble(size(x))

    return 
  end function f_mean

!-------------------------------------------------------------------------------

  function f_stdev(x,imean)
    implicit none 
    double precision, dimension(:), intent(in)  :: x
    double precision, intent(in), optional      :: imean
    double precision                            :: f_stdev 
    double precision                            :: mean
    integer                                     :: i 

    if (.not. present(imean)) then 
      mean = f_mean(x)
    else 
      mean = imean
    endif

    f_stdev = 0.d0
    do i=1,size(x)
      f_stdev = f_stdev + (x(i) - mean)*(x(i) - mean)
    enddo
    f_stdev = sqrt(f_stdev / dble(size(x)))

    return 
  end function f_stdev

!-------------------------------------------------------------------------------
! STRING routines
!-------------------------------------------------------------------------------
! removes all white spaces from a string:

  subroutine f_stripspaces(chr)
    implicit none
    character(len=*), intent(inout) :: chr
    integer                         :: stringlen
    integer                         :: last
    integer                         :: actual

    stringlen = len(chr)
    last = 1
    actual = 1
    do while (actual .lt. stringlen)
      if (chr(last:last) .eq. ' ') then
        actual = actual + 1
        chr(last:last) = chr(actual:actual)
        chr(actual:actual) = ' '
      else
        last = last + 1
        if (actual .lt. last) then
          actual = last
        endif
      endif
    enddo
    
    return
  end subroutine f_stripspaces

!-------------------------------------------------------------------------------
! read integers:

  subroutine f_read_int(chartag,charfile,value,info)
    implicit none
    character(len=*), intent(in)  :: chartag
    character(len=*), intent(in)  :: charfile
    integer, intent(out)          :: value
    integer, intent(out)          :: info
    character(len=25)             :: tag
    integer                       :: ioerr

    value = 0
    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_int > ',trim(charfile),' could not be openend. no value is'
      print *,'            returned.'
      return
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif
      
      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,value
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)

    return
  end subroutine f_read_int

!-------------------------------------------------------------------------------
! read doubles:

  subroutine f_read_dbl(chartag,charfile,value,info)
    implicit none
      character(len=*), intent(in)                :: chartag
      character(len=*), intent(in)                :: charfile
      double precision                            :: value
      integer, intent(out)                        :: info
      character(len=25)                           :: tag
      integer                                     :: ioerr

    value = 0.d0

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_dp > ',trim(charfile),' could not be openend. no value is '
      print *,'           returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,value
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)

    return
  end subroutine f_read_dbl

!-------------------------------------------------------------------------------
! read characters:

  subroutine f_read_chr(chartag,charfile,value,info)
    implicit none
      character(len=*), intent(in)  :: chartag
      character(len=*), intent(in)  :: charfile
      character(len=*), intent(out) :: value
      integer, intent(out)          :: info
      character(len=25)             :: tag
      integer                       :: ioerr

    value = ''

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_char > ',trim(charfile),' could not be openend. no value'
      print *,'             is returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,value
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    value = trim(value)
    close(12)

    return
  end subroutine f_read_chr

!-------------------------------------------------------------------------------
! read vectors:

  subroutine f_read_vec(chartag,charfile,value,info)
    implicit none
    character(len=*), intent(in)                :: chartag
    character(len=*), intent(in)                :: charfile
    double precision, dimension(:), intent(out) :: value
    integer, intent(out)                        :: info
    character(len=25)                           :: tag
    integer                                     :: ioerr

    value = 0.d0

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_vec > ',trim(charfile),' could not be openend. no value is '
      print *,'             returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,value
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)

    return
  end subroutine f_read_vec

!-------------------------------------------------------------------------------
! read logicals:

  subroutine f_read_log(chartag,charfile,value,info)
    implicit none
      character(len=*), intent(in)    :: chartag
      character(len=*), intent(in)    :: charfile
      logical, intent(out)            :: value
      integer, intent(out)            :: info
      character(len=30)               :: tag
      integer                         :: ioerr

    value = .false.

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_char > ',trim(charfile),' could not be openend. no value'
      print *,'             is returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        info = 0
        value = .true.
        exit
      endif
    enddo
    close(12)

    return
  end subroutine f_read_log

!-------------------------------------------------------------------------------
! read integers:

  subroutine f_read_store_int(chartag,charfile,values,default,info)
    implicit none
    character(len=*), intent(in)  :: chartag
    character(len=*), intent(in)  :: charfile
    integer, intent(out)          :: values
    integer, intent(in)           :: default
    integer, intent(out)          :: info
    character(len=25)             :: tag
    integer                       :: ioerr

    values = 0
    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_int > ',trim(charfile),' could not be openend. no value is'
      print *,'            returned.'
      return
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif
      
      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,values
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)

    if (info .ne. 0) then
      values = default 
    endif

    return
  end subroutine f_read_store_int

!-------------------------------------------------------------------------------
! read doubles:

  subroutine f_read_store_dbl(chartag,charfile,values,default,info)
    implicit none
      character(len=*), intent(in)                :: chartag
      character(len=*), intent(in)                :: charfile
      double precision, intent(out)               :: values
      double precision, intent(in)                :: default
      integer, intent(out)                        :: info
      character(len=25)                           :: tag
      integer                                     :: ioerr

    values = 0.d0

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_dp > ',trim(charfile),' could not be openend. no value is '
      print *,'           returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,values
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)

    if (info .ne. 0) then 
      values = default 
    endif

    return
  end subroutine f_read_store_dbl

!-------------------------------------------------------------------------------
! read characters:

  subroutine f_read_store_chr(chartag,charfile,values,default,info)
    implicit none
      character(len=*), intent(in)  :: chartag
      character(len=*), intent(in)  :: charfile
      character(len=*), intent(out) :: values
      character(len=*), intent(in)  :: default
      integer, intent(out)          :: info
      character(len=25)             :: tag
      integer                       :: ioerr

    values = ''

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_char > ',trim(charfile),' could not be openend. no value'
      print *,'             is returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,values
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)
    
    if (info .eq. 0) then
      values = trim(values)
    else 
      values = trim(default)
    endif

    return
  end subroutine f_read_store_chr

!-------------------------------------------------------------------------------
! read vectors:

  subroutine f_read_store_vec(chartag,charfile,values,default,info)
    implicit none
    character(len=*), intent(in)                :: chartag
    character(len=*), intent(in)                :: charfile
    double precision, dimension(:), intent(out) :: values
    double precision, dimension(:), intent(in)  :: default
    integer, intent(out)                        :: info
    character(len=25)                           :: tag
    integer                                     :: ioerr

    values = 0.d0

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_vec > ',trim(charfile),' could not be openend. no value is '
      print *,'             returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        backspace(12)
        read(12,*,iostat=ioerr) tag,values
        if (ioerr .ne. 0) then
          info = -1
          exit
        endif
        info = 0
        exit
      endif
    enddo
    close(12)

    if (info .ne. 0) then
      values = default 
    endif

    return
  end subroutine f_read_store_vec

!-------------------------------------------------------------------------------
! read logicals:

  subroutine f_read_store_log(chartag,charfile,values,default,info)
    implicit none
      character(len=*), intent(in)    :: chartag
      character(len=*), intent(in)    :: charfile
      logical, intent(out)            :: values
      logical, intent(in)             :: default
      integer, intent(out)            :: info
      character(len=30)               :: tag
      integer                         :: ioerr

    values = .false.

    open(12,file=trim(charfile),action='read',status='old',iostat=ioerr)
    if (ioerr .ne. 0) then
      print *,' read_char > ',trim(charfile),' could not be openend. no value'
      print *,'             is returned.'
    endif

    info = 1
    do
      read(12,*,iostat=ioerr) tag
      if (ioerr .ne. 0) then
        exit
      endif

      if (trim(chartag) .eq. trim(tag)) then
        info = 0
        values = .true.
        exit
      endif
    enddo
    close(12)
    if (info .ne. 0) then
      values = default 
    endif

    return
  end subroutine f_read_store_log

!-------------------------------------------------------------------------------
! trim integer string:
  
  function f_trimmed_int(int)
    implicit none
    character(len=15) :: f_trimmed_int
    integer           :: int
    character(len=30) :: dchr

    write(dchr,'(I15)') int
    call f_stripspaces(dchr)
    f_trimmed_int = trim(dchr)

    return
  end function f_trimmed_int

!-------------------------------------------------------------------------------
! trim character string:

  function f_trimmed_chr(chr)
    implicit none
    character(len=15) :: f_trimmed_chr
    character(len=30) :: chr
    character(len=30) :: dchr

    dchr = trim(chr)
    call f_stripspaces(dchr)
    f_trimmed_chr = trim(dchr)

    return
  end function f_trimmed_chr

!-------------------------------------------------------------------------------
! trim logical string:

  function f_trimmed_log(log)
    implicit none
    character(len=15) :: f_trimmed_log
    logical           :: log

    if (log) then
      f_trimmed_log = 'true'
    else
      f_trimmed_log = 'false'
    endif

    return
  end function f_trimmed_log

!-------------------------------------------------------------------------------
! trim double string:

  function f_trimmed_dbl(dbl)
    implicit none
    character(len=15) :: f_trimmed_dbl
    double precision  :: dbl
    character(len=30) :: dchr

    write(dchr,'(F15.5)') dbl
    call f_stripspaces(dchr)
    f_trimmed_dbl = trim(dchr)

    return
  end function f_trimmed_dbl
 
!-------------------------------------------------------------------------------
! trim vector:

  function f_trimmed_vec(dbl)
    implicit none
    character(len=15)               :: f_trimmed_vec
    double precision, dimension(3)  :: dbl
    character(len=15)               :: dchr

    write(dchr,100) dbl
100 format(F4.1,' ',F4.1,' ',F4.1)
    f_trimmed_vec = dchr 

    return
  end function f_trimmed_vec

!-------------------------------------------------------------------------------
! errors

  subroutine f_critical_(chr,line)
    implicit none 
    character(len=*), intent(in) :: chr 
    integer, intent(in)           :: line 

    print *,'    > critical error in line ',f_trimmed_(line)
    print *,'      of file ',trim(chr)
    stop 

    return 
  end subroutine f_critical_

!-------------------------------------------------------------------------------

end module functions