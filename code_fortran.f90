!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      MERCURY.F90    (ErikSoft   3 May 2002)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Author: John E. Chambers (and Christophe Cossou)

! Mercury is a general-purpose N-body integration package for problems in
! celestial mechanics.


!------------------------------------------------------------------------------


program mercury
  
    use user_module
    implicit none
  
    integer :: j,nbod,nbig
    integer :: opflag,ngflag,ndump,nfun
    integer :: error

    real(double_precision) :: cefac,time,h0,tol,en(3),am(3),rcen,jcen(3)

    integer, dimension(:), allocatable :: stat ! (Number of bodies)
    real(double_precision), dimension(:), allocatable :: m,rho,rceh,epoch ! (Number of bodies)
    real(double_precision), dimension(:,:), allocatable :: xh,vh,s ! (3,Number of bodies)
    real(double_precision), dimension(:,:), allocatable :: ngf ! (4,Number of bodies)
    character(len=8), dimension(:), allocatable :: id ! (Number of bodies)
  
    !------------------------------------------------------------------------------

    call mfo_user()

end program mercury
