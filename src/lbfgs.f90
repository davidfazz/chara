! lbfgs.f90 - elapid - david gallina - university of kassel - 2022
!-------------------------------------------------------------------------------

module lbfgs
  use constants 
  use functions
  use random
  use potential
  implicit none

!-------------------------------------------------------------------------------
! parameters that are necessary for the algorithm.
! > maxbfgsstp: maximum step length of lbfgs step
! > emaxlimit: maximum energy
! > eminlimit: minimum energy
! > maxerise: maximum allowed energy increase
! > maxefall: maximum allowed energy decrease
! > rescale: rescaling at each step
! > rotate: rotating vectors that are close to the pole

  double precision, parameter, private :: maxbfgsstep = 0.2d0
  double precision, parameter, private :: emaxlimit   = 1.d6
  double precision, parameter, private :: eminlimit   =-1.d6
  double precision, parameter, private :: maxerise    = 1.d-10
  double precision, parameter, private :: maxefall    =-huge(1.d0)
  logical, parameter, private          :: lbfgs_rescale = .true.

contains

!-------------------------------------------------------------------------------

  subroutine lbfgs_(n,m,x,deps,converged,energy,itmax,itdone,output)
    implicit none
    integer, intent(in)                         :: n
    integer, intent(in)                         :: m
    double precision, dimension(n)              :: x
    double precision, intent(in)                :: deps
    logical, intent(out)                        :: converged
    double precision, intent(out)               :: energy
    integer, intent(in)                         :: itmax
    integer, intent(out)                        :: itdone
    logical, intent(in)                         :: output
    double precision, dimension(n*(2*m+1)+2*m)  :: w
    double precision, dimension(n)              :: grad
    double precision, dimension(n)              :: diag
    double precision, dimension(n)              :: gnew
    double precision, dimension(n)              :: wtemp
    double precision, dimension(n)              :: xsave
    double precision, dimension(n)              :: projstep
    double precision                            :: rmsgrad
    double precision                            :: ys
    double precision                            :: yy
    double precision                            :: yr
    double precision                            :: dummy
    double precision                            :: gnorm
    double precision                            :: slength
    double precision                            :: stp
    double precision                            :: enew
    double precision                            :: overlap
    double precision                            :: dot1
    double precision                            :: dot2
    double precision                            :: beta
    double precision                            :: sq
    double precision                            :: ddot
    integer                                     :: iter
    integer                                     :: nfail
    integer                                     :: ndecrease
    integer                                     :: i
    integer                                     :: cp
    integer                                     :: bound
    integer                                     :: point
    integer                                     :: ispt
    integer                                     :: npt
    integer                                     :: iypt
    integer                                     :: inmc
    integer                                     :: iscn
    integer                                     :: iycn
    logical                                     :: failed
    
! startup of lbfgs routine - setting initial values.

    iter   = 0
    nfail  = 0
    itdone = 0
    
! calculate gradient and energy aswell as the rms.


    call potential_(x,energy,grad)
    rmsgrad = f_rms(grad)

    if (rmsgrad .lt. deps) then
      converged = .true.
      return
    else
      converged = .false.
    endif
        
    do while ((.not. converged) .and. (iter .le. itmax))
      if (output) then
        print *,'    lbfgs > debugging - iteration:',iter
        print *,'                      - energy   :',energy
        print *,'                      - rms      :',rmsgrad
      endif
      ndecrease = 0
        
! on first entry, the parameters have to be initialised. also we check, if there 
! are improper input parameters.

      if (iter .eq. 0) then
        if ((n .le. 0) .or. (m .le. 0)) then
          print *,'    lbfgs > improper input parameters! n and/or m are not positive!'
          print *,'      '
          print *,'    chara > aborted!'
          stop
        endif
  
        point = 0
        converged = .false.

! setting up diagonal elements of inverse hessian.

        do i=1,n
          diag(i) = 1.d0
        enddo

! indices for storage and first step calculation:
!	 -> ispt: storage of search steps
!  -> iypt: storage of gradient differences
            
        ispt = n + 2*m
        iypt = ispt + n*m
        do i=1,n
          dummy = -grad(i)*diag(i)
          w(ispt+i) = dummy
          w(i) = dummy
        enddo
        gnorm = sqrt(ddot(n,grad,1,grad,1))

! first guess for step length is pretty cautious.

        stp = min(1.d0/gnorm,gnorm)

! iter greater than 0.

      else
        bound = iter
        if (iter .gt. m) then
          bound = m
        endif

! setting up diagonal elements of inverse hessian.

        ys = ddot(n,w(iypt+npt+1),1,w(ispt+npt+1),1)
        yy = ddot(n,w(iypt+npt+1),1,w(iypt+npt+1),1)
        if (yy .eq. 0.d0) then
          yy = 1.d0
          if (output) then
            print *,'    lbfgs > debugging - yy is set to 1.d0!'
          endif
        endif
        if (ys .eq. 0.d0) then
          ys = 1.d0
          if (output) then
            print *,'    lbfgs > debugging - ys is set to 1.d0!'
          endif
        endif
        do i=1,n
          diag(i) = ys/yy
        enddo
                
! computing -hess*grad using the formula given in nocedal(1980).

        cp = point
        if (point .eq. 0) then
          cp = m
        endif
        w(n+cp) = 1.d0/ys
        do i=1,n
          w(i) = -grad(i)
        enddo
        cp = point
        do i=1,bound
          cp = cp-1
          if (cp .eq. -1) then
            cp = m-1
          endif
          sq = ddot(n,w(ispt+cp*n+1),1,w,1)
          inmc = n + m + cp + 1
          iycn = iypt + cp*n
          w(inmc) = w(n+cp+1)*sq
          call daxpy(n,-w(inmc),w(iycn+1),1,w,1)
        enddo
        do i=1,n
          w(i) = diag(i)*w(i)
        enddo

        do i=1,bound
          yr = ddot(n,w(iypt+cp*n+1),1,w,1)
          beta = w(n+cp+1)*yr
          inmc = n + m + cp + 1
          beta = w(inmc) - beta
          iscn = ispt + cp*n
          call daxpy(n,beta,w(iscn+1),1,w,1)
          cp = cp + 1
          if (cp .eq. m) then
            cp = 0
          endif
          stp = 1.0
        enddo
      endif
    
! store the new search direction

      if (iter .gt. 0) then
        do i=1,n
          w(ispt+point*n+i) = w(i)
        enddo
      endif
      dot1 = sqrt(ddot(n,grad,1,grad,1))

      dummy = 1.d0
      do i=1,n
        if (abs(w(i)) .gt. dummy) then
          dummy = abs(w(i))
        endif
      enddo
      do i=1,n
        wtemp(i) = w(i)/dummy
      enddo
      dot2 = sqrt(ddot(n,wtemp,1,wtemp,1))
      overlap = 0.d0
      if (dot1*dot2 .ne. 0.d0) then
        overlap = ddot(n,grad,1,wtemp,1)/(dot1*dot2)
      endif

! if search direction and gradient have positive overlap, the last step is 
! reversed!

      if (overlap .gt. 0) then
        do i=1,n
          w(ispt+point*n+i) = -w(i)
        enddo
        if (output) then
          print *,'    lbfgs > debugging - positive overlap detected!'
          print *,'                      - step is reversed!'
        endif
      endif

      do i=1,n
        w(i) = grad(i)
      enddo
    
! finding appropriate step length.
    
      slength = 0.d0
      do i=1,n
        slength = slength + w(ispt+point*n+i)**2
      enddo
      slength = sqrt(slength)
      if (stp*slength .gt. maxbfgsstep) then
        stp = maxbfgsstep/slength
      endif

! we got our step length and search direction saving coordinates so that we can 
! undo the step reliably
                
      xsave(1:n) = x(1:n)
      do i=1,n
        projstep(i) = stp*w(ispt+point*n+i)
      enddo
      x = x + projstep
      failed = .false.

      do
        if (lbfgs_rescale) then
          call potential_rescale(x)
        endif
        call potential_(x,enew,gnew)

        rmsgrad = f_rms(gnew)
        if (rmsgrad .lt. deps) then
          if (output) then
            print *,'    lbfgs > debugging - converged after:',itdone,'iterations'
            print *,'                      - energy         :',energy
            print *,'                      - rms            :',rmsgrad
          endif
          converged = .true.
          return
        endif

  ! if the new calculated energy is outside of reasonable boundaries, we undo the 
  ! last step.

        if ((enew .gt. emaxlimit) .or. (enew .lt. eminlimit)) then
          energy = 1.d3
          enew = 1.d3
          rmsgrad = 1.d3
          converged = .false.
          if (output) then
            print *,'    lbfgs > debugging - energy is out of bounds!'
          endif

          return
        endif

  ! if the rise or fall of the energy is in reasonable boundaries, the step is 
  ! executed.

        if ((enew - energy .le. maxerise) .and. (enew - energy .ge. maxefall)) then
          iter = iter + 1
          itdone = itdone + 1
          energy = enew
          do i=1,n
            grad(i) = gnew(i)
          enddo
          exit
      
  ! energy decreased too much - trying again with smaller step size.

        elseif (enew - energy .le. maxefall) then
          if (output) then
            print *,'    lbfgs > debugging - decrease is not sufficient!'
            print *,'                      - step length is changed!'
          endif
          if (ndecrease .gt. 5) then
            if (output) then
              print *,'    lbfgs > debugging - run failed. resetting with new gradient!'
            endif
            failed = .true.
            nfail = nfail + 1
            x(1:n) = xsave(1:n)
            grad(1:n) = gnew(1:n)
            iter = 0

            if (nfail .gt. 20) then
              if (output) then
                print *,'    lbfgs > debugging - 20 runs failed! giving up!'
              endif
              converged = .false.
              return
            endif
            exit
          endif
          x(1:n) = xsave(1:n)
          do i=1,n
            x(i) = x(i) + 0.5d0*stp*w(ispt+point*n+i)
          enddo
          stp = stp/2.d0
          ndecrease = ndecrease + 1
          cycle

  ! energy increased too much - trying again with smaller step size.
      
        else
          if (output) then
            print *,'    lbfgs > debugging - no decrease detected!'
            print *,'                      - step length is changed!'
          endif
          if (ndecrease .gt. 10) then
            if (output) then
              print *,'    lbfgs > debugging - run failed. resetting with new gradient!'
            endif
            nfail = nfail + 1
            x(1:n) = xsave(1:n)
            grad(1:n) = gnew(1:n)
            iter = 0

            if (nfail .gt. 5) then
              if (output) then
                print *,'    lbfgs > debugging - 5 runs failed! giving up!'
              endif
              converged = .false.
              return
            endif
            exit
          endif

          x(1:n) = xsave(1:n)
          do i=1,n
            x(i) = x(i) + 0.1d0*stp*w(ispt+point*n+i)
          enddo
          stp = stp/1.d1
          ndecrease = ndecrease + 1
          cycle
        endif
      enddo
        
! check, if lfbgs routine failed.

      if (failed) then
        cycle
      endif
        
! compute the new step and gradient changes.
    
      npt = point*n
      do i=1,n
        w(ispt+npt+i) = stp*w(ispt+npt+i)
        w(iypt+npt+i) = grad(i) - w(i)
      enddo
      point = point + 1
      if (point .eq. m) then
        point = 0
      endif
    enddo

    return
  end subroutine lbfgs_

!-------------------------------------------------------------------------------

end module lbfgs
