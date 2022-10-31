! finitediffs.f90 - elapid - david gallina - university of kassel - 2022
!...............................................................................
! module contains a number of routines to calculate the gradient and hessian of
! a given configuration by means of finite energy differences. 
!...............................................................................

module approx
  use potential 
  implicit none 
  
contains

!-------------------------------------------------------------------------------
! routine to calculate forward the differences approximation of the hessian
! along a vector v;

  subroutine ap_forward_hess(x,v,hess,deps)
    implicit none 
    double precision, dimension(c_n), intent(in)  :: x 
    double precision, dimension(c_n), intent(in)  :: v
    double precision, dimension(c_n), intent(out) :: hess
    double precision, intent(in)                  :: deps 
    double precision, dimension(c_n)              :: gfx 
    double precision, dimension(c_n)              :: ngfx 
    double precision                              :: fx 
    double precision                              :: nfx

    call potential_(x,fx,gfx)
    call potential_(x+deps*v,nfx,ngfx)

    hess = (ngfx - gfx) / deps

    return 
  end subroutine ap_forward_hess

! !...............................................................................
! ! routine to calculate the central difference approximation of the gradient:

!   subroutine fd_grad(x,grad,deps,ctype)
!     implicit none
!     double precision, dimension(c_n), intent(in)  :: x 
!     double precision, dimension(c_n), intent(out) :: grad
!     double precision, intent(in)                  :: deps
!     character(len=*), intent(in)                  :: ctype
!     double precision, dimension(c_n)              :: gold
!     double precision, dimension(c_n)              :: xforw
!     double precision, dimension(c_n)              :: xback
!     double precision                              :: eback 
!     double precision                              :: eforw
!     integer                                       :: i 
  
! ! calculate the different entries of the gradient:
  
!     do i=1,c_n
!       xforw = x 
!       xback = x
!       xforw(i) = xforw(i) + deps 
!       xback(i) = xback(i) - deps 
!       call potential_(xforw,eforw,gold)
!       call potential_(xback,eback,gold)
!       grad(i) = (eforw - eback) / (2.d0 * deps)
!       if (trim(ctype) .eq. 'spherical') then 
!         if (mod(i,2) .eq. 0) then
!           grad(i) = grad(i) / sin(x(i-1))
!         endif
!       endif 
!     enddo
  
!     return
!   end subroutine fd_grad
  
! !...............................................................................
! ! routine to calculate the forward difference approximation of the hessian:
  
!   subroutine fd_hess(x,hess,deps,ctype)
!     implicit none
!     double precision, dimension(c_n), intent(in)      :: x 
!     double precision, dimension(c_n,c_n), intent(out) :: hess
!     double precision, intent(in)                      :: deps 
!     character(len=*), intent(in)                      :: ctype
!     double precision, dimension(c_n)                  :: gold 
!     double precision, dimension(c_n)                  :: xij 
!     double precision, dimension(c_n)                  :: xik 
!     double precision, dimension(c_n)                  :: xki 
!     double precision, dimension(c_n)                  :: xji 
!     double precision, dimension(c_n)                  :: xjj 
!     double precision, dimension(c_n)                  :: xii 
!     double precision, dimension(c_n)                  :: xi 
!     double precision, dimension(c_n)                  :: xj
!     double precision, dimension(c_n)                  :: xg
!     double precision                                  :: eold
!     double precision                                  :: eij
!     double precision                                  :: eik
!     double precision                                  :: eki
!     double precision                                  :: eji
!     double precision                                  :: eii
!     double precision                                  :: ejj
!     double precision                                  :: ei
!     double precision                                  :: ej 
!     integer                                           :: i 
!     integer                                           :: j

  
! ! calculate the entries of the hessian:
  
!     call potential_(x,eold,gold)
!     do i=1,c_n 
!       do j=1,c_n 
!         if (i .ne. j) then
!           xij = x 
!           xij(i) = xij(i) + deps
!           xij(j) = xij(j) + deps
!           xik = x 
!           xik(i) = xik(i) + deps
!           xik(j) = xik(j) - deps
!           xki = x 
!           xki(i) = xki(i) - deps
!           xki(j) = xki(j) + deps 
!           xji = x 
!           xji(i) = xji(i) - deps
!           xji(j) = xji(j) - deps
!           call potential_(xij,eij,gold)
!           call potential_(xik,eik,gold)
!           call potential_(xki,eki,gold)
!           call potential_(xji,eji,gold)
!           hess(i,j) = (eij - eik - eki + eji) / (4.d0*deps*deps)
!         elseif (i .eq. j) then
!           xii = x 
!           xii(i) = xii(i) + 2.d0*deps
!           xi = x 
!           xi(i) = xi(i) + deps 
!           xjj = x 
!           xjj(i) = xjj(i) - 2.d0*deps
!           xj = x 
!           xj(i) = xj(i) - deps
!           call potential_(xii,eii,gold)
!           call potential_(xi,ei,gold)
!           call potential_(xjj,ejj,gold)
!           call potential_(xj,ej,gold) 
!           hess(i,j) = (-eii + 16.d0*ei - 30.d0*eold + 16.d0*ej - ejj) &
!                     / (12.d0*deps*deps)
!         endif
!       enddo
!     enddo     
    
!     if (trim(ctype) .eq. 'spherical') then
!       call fd_grad(x,xg,deps,'cartesian')
!       do i=1,c_n,2 
!         do j=1,c_n,2
!           if (i .ne. j) then
!             hess(i+0,j+1) = hess(i+0,j+1) / sin(x(j))
!             hess(i+1,j+0) = hess(i+1,j+0) / sin(x(i))
!             hess(i+1,j+1) = hess(i+1,j+1) / (sin(x(i))*sin(x(j)))
!           elseif (i .eq. j) then
!             hess(i+0,i+1) = hess(i+0,i+1) / sin(x(i)) - cos(x(i)) &
!                           * xg(i+1) / sin(x(i))
!             hess(i+1,i+0) = hess(i+1,i+0) / sin(x(i)) - cos(x(i)) &
!                           * xg(i+1) / sin(x(i))
!             hess(i+1,i+1) = hess(i+1,i+1) / sin(x(i)) / sin(x(i)) &
!                           + cos(x(i)) / sin(x(i)) * xg(i)
!           endif
!         enddo
!       enddo
!     endif
    
!     return
!   end subroutine fd_hess
  
!...............................................................................
  
end module approx