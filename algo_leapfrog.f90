!******************************************************************************
! MODULE: algo_radau
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Modules that gather various functions about the RADAU algorithm.\n\n
!! RADAU is described in E.Everhart (1985) in ``The Dynamics of Comets:
!! Their Origin and Evolution'' p185-202, eds. A.Carusi & G.B.Valsecchi,
!! pub. Reidel.
!
!******************************************************************************

module algo_radau

  use types_numeriques

  implicit none

  private
  ! Gauss-Radau spacings for substeps within a sequence, for the 15th order 
  ! integrator. The sum of the H values should be 3.733333333333333
  
  real(double_precision), dimension(8), parameter :: h = (/ 0.d0,.0562625605369221d0,.1802406917368924d0,&
    .3526247171131696d0,.5471536263305554d0,.7342101772154105d0,.8853209468390958d0,.9775206135612875d0/)
    
    ! Constant coefficients used in series expansions for X and V
    !  XC: 1/2,  1/6,  1/12, 1/20, 1/30, 1/42, 1/56, 1/72
    !  VC: 1/2,  1/3,  1/4,  1/5,  1/6,  1/7,  1/8
  real(double_precision), dimension(8), parameter :: xc = (/.5d0,.1666666666666667d0,.08333333333333333d0,.05d0,&
      .03333333333333333d0,.02380952380952381d0,.01785714285714286d0,.01388888888888889d0/)
      
  real(double_precision), dimension(7), parameter :: vc = (/.5d0,.3333333333333333d0,.25d0,.2d0,.1666666666666667d0,&
        .1428571428571429d0,.125d0/)

  ! values that need to be saved in mdt_ra15
  real(double_precision) :: c(21),d(21),r(28)
  real(double_precision), dimension(:,:), allocatable :: b, e !  (7,3*Number of bodies)

  public :: mdt_ra15, allocate_radau
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Christophe Cossou
!
!> @date 2011
!
! DESCRIPTION: 
!> @brief allocate various private variable of the module to avoid testing at each timestep
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine allocate_radau(nb_bodies)

  implicit none
  
  integer, intent(in) :: nb_bodies !< [in] number of bodies, that is, the size of the arrays
  
  if (.not. allocated(b)) then
    allocate(b(7,3*nb_bodies))
    b(1:7,1:3*nb_bodies) = 0.d0
  end if
  
  if (.not. allocated(e)) then
    allocate(e(7,3*nb_bodies))
    e(1:7,1:3*nb_bodies) = 0.d0
  end if
  
  
end subroutine allocate_radau

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> John E. Chambers
!
!> @date 2 March 2001
!
! DESCRIPTION: 
!> @brief Integrates NBOD bodies (of which NBIG are Big) for one timestep H0 using
!! Everhart's RA15 integrator algorithm. The accelerations are calculated
!! using the subroutine FORCE. The accuracy of the step is approximately 
!! determined by the tolerance parameter TOL.
!!\n\n
!! Based on RADAU by E. Everhart, Physics Department, University of Denver.
!! Comments giving equation numbers refer to Everhart (1985) ``An Efficient
!! Integrator that Uses Gauss-Radau Spacings'', in The Dynamics of Comets:
!! Their Origin and Evolution, p185-202, eds. A. Carusi & G. B. Valsecchi,
!! pub Reidel. (A listing of the original subroutine is also given in this 
!! paper.)
!!\n\n
!! DTFLAG = 0 implies first ever call to this subroutine, 
!!        = 1 implies first call since number/masses of objects changed.
!!        = 2 normal call
!
!> @note Input/output must be in coordinates with respect to the central body.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine mdt_ra15 (time,t,tdid,tol,jcen,nbod,nbig,mass,x1,v1,spin,rphys,rcrit,ngf,stat,dtflag,ngflag,nce,ice,jce,force)
  
  use physical_constant
  use mercury_constant
  use mercury_globals

  implicit none

  
  ! Input/Output
  integer, intent(in) :: nbod !< [in] current number of bodies (1: star; 2-nbig: big bodies; nbig+1-nbod: small bodies)
  integer, intent(in) :: nbig !< [in] current number of big bodies (ones that perturb everything else)
  integer, intent(inout) :: dtflag
  integer, intent(in) :: ngflag !< [in] do any bodies experience non-grav. forces?
!!\n                            ( 0 = no non-grav forces)
!!\n                              1 = cometary jets only
!!\n                              2 = radiation pressure/P-R drag only
!!\n                              3 = both
  integer, intent(in) :: stat(nbod) !< [in] status (0 => alive, <>0 => to be removed)
  integer, intent(in) :: nce
  integer, intent(in) :: ice(nce)
  integer, intent(in) :: jce(nce)
  real(double_precision), intent(in) :: time !< [in] current epoch (days)
  real(double_precision), intent(inout) :: t
  real(double_precision), intent(out) :: tdid
  real(double_precision), intent(in) :: tol !< [in] Integrator tolerance parameter (approx. error per timestep)
  real(double_precision), intent(in) :: jcen(3) !< [in] J2,J4,J6 for central body (units of RCEN^i for Ji)
  real(double_precision), intent(in) :: mass(nbod) !< [in] mass (in solar masses * K2)
  
  real(double_precision), intent(inout) :: x1(3*nbod)
  real(double_precision), intent(inout) :: v1(3*nbod)
  real(double_precision), intent(in) :: spin(3*nbod)
  real(double_precision), intent(in) :: ngf(4,nbod) !< [in] non gravitational forces parameters
  !! \n(1-3) cometary non-gravitational (jet) force parameters
  !! \n(4)  beta parameter for radiation pressure and P-R drag
  real(double_precision), intent(in) :: rphys(nbod)
  real(double_precision), intent(in) :: rcrit(nbod)
  
  external force
  
  ! Local
  integer :: nv,niter,j,k,n
  real(double_precision) :: x(3*nb_bodies_initial),v(3*nb_bodies_initial),a(3*nb_bodies_initial),a1(3*nb_bodies_initial)
  real(double_precision) :: g(7,3*nb_bodies_initial)
  real(double_precision) :: s(9)
  real(double_precision) :: q,q2,q3,q4,q5,q6,q7,temp,gk
  
  !------------------------------------------------------------------------------
  
  
  ! If this is first call to the subroutine, set values of the constant arrays
  ! (R = R21, R31, R32, R41, R42, R43 in Everhart's paper.)
  if (dtflag.eq.0) then
     n = 0
     do j = 2, 8
        do k = 1, j - 1
           n = n + 1
           r(n) = 1.d0 / (h(j) - h(k))
        end do
     end do
     
     ! Constants to convert between B and G arrays (C = C21, C31, C32, C41, C42...)
     c(1) = - h(2)
     d(1) =   h(2)
     n = 1
     do j = 3, 7
        n = n + 1
        c(n) = -h(j) * c(n-j+2)
        d(n) =  h(2) * d(n-j+2)
        do k = 3, j - 1
           n = n + 1
           c(n) = c(n-j+1)  -  h(j) * c(n-j+2)
           d(n) = d(n-j+1)  +  h(k) * d(n-j+2)
        end do
        n = n + 1
        c(n) = c(n-j+1) - h(j)
        d(n) = d(n-j+1) + h(j)
     end do
     
     dtflag = 1
  end if
  
  nv = 3 * nbod
do
  
  ! If this is first call to subroutine since number/masses of objects changed
  ! do 6 iterations and initialize B, E arrays, otherwise do 2 iterations.
  if (dtflag.eq.1) then
     niter = 6
     do j = 4, nv
        do k = 1, 7
           b (k,j) = 0.d0
           e (k,j) = 0.d0
        end do
     end do
  else
     niter = 2
  end if
  
  ! Calculate forces at the start of the sequence
  call force (time,jcen,nbod,nbig,mass,x1,v1,spin,rcrit,a1,stat,ngf,ngflag,nce,ice,jce)
  
  ! Find G values from B values predicted at the last call (Eqs. 7 of Everhart)
  do k = 4, nv
     g(1,k) = b(7,k)*d(16) + b(6,k)*d(11) + b(5,k)*d(7) + b(4,k)*d(4) + b(3,k)*d(2) + b(2,k)*d(1) + b(1,k)
     g(2,k) = b(7,k)*d(17) + b(6,k)*d(12) + b(5,k)*d(8) + b(4,k)*d(5) + b(3,k)*d(3) + b(2,k)
     g(3,k) = b(7,k)*d(18) + b(6,k)*d(13) + b(5,k)*d(9) + b(4,k)*d(6) + b(3,k)
     g(4,k) = b(7,k)*d(19) + b(6,k)*d(14) + b(5,k)*d(10) + b(4,k)
     g(5,k) = b(7,k)*d(20) + b(6,k)*d(15) + b(5,k)
     g(6,k) = b(7,k)*d(21) + b(6,k)
     g(7,k) = b(7,k)
  end do
  
  !------------------------------------------------------------------------------
  
  !  MAIN  LOOP  STARTS  HERE
  
  ! For each iteration (six for first call to subroutine, two otherwise)...
  do n = 1, niter
     
     ! For each substep within a sequence...
     do j = 2, 8
        
        ! Calculate position predictors using Eqn. 9 of Everhart
        s(1) = t * h(j)
        s(2) = s(1) * s(1) * .5d0
        s(3) = s(2) * h(j) * .3333333333333333d0
        s(4) = s(3) * h(j) * .5d0
        s(5) = s(4) * h(j) * .6d0
        s(6) = s(5) * h(j) * .6666666666666667d0
        s(7) = s(6) * h(j) * .7142857142857143d0
        s(8) = s(7) * h(j) * .75d0
        s(9) = s(8) * h(j) * .7777777777777778d0
        
        do k = 4, nv
           x(k) = s(9)*b(7,k) + s(8)*b(6,k) + s(7)*b(5,k) + s(6)*b(4,k) + s(5)*b(3,k) + s(4)*b(2,k) &
                + s(3)*b(1,k) + s(2)*a1(k)  + s(1)*v1(k) + x1(k)
        end do
        
        ! If necessary, calculate velocity predictors too, from Eqn. 10 of Everhart
        if (ngflag.ne.0) then
           s(1) = t * h(j)
           s(2) = s(1) * h(j) * .5d0
           s(3) = s(2) * h(j) * .6666666666666667d0
           s(4) = s(3) * h(j) * .75d0
           s(5) = s(4) * h(j) * .8d0
           s(6) = s(5) * h(j) * .8333333333333333d0
           s(7) = s(6) * h(j) * .8571428571428571d0
           s(8) = s(7) * h(j) * .875d0
           
           do k = 4, nv
              v(k) = s(8)*b(7,k) + s(7)*b(6,k) + s(6)*b(5,k) + s(5)*b(4,k) + s(4)*b(3,k)&
                   + s(3)*b(2,k) + s(2)*b(1,k) + s(1)*a1(k)  + v1(k)
           end do
        end if
        
        ! Calculate forces at the current substep
        call force (time,jcen,nbod,nbig,mass,x,v,spin,rcrit,a,stat,ngf,ngflag,nce,ice,jce)
        
        ! Update G values using Eqs. 4 of Everhart, and update B values using Eqs. 5
        select case (j)
        case (2)
           do k = 4, nv
              temp = g(1,k)
              g(1,k) = (a(k) - a1(k)) * r(1)
              b(1,k) = b(1,k) + g(1,k) - temp
           end do
        case (3)
           do k = 4, nv
              temp = g(2,k)
              gk = a(k) - a1(k)
              g(2,k) = (gk*r(2) - g(1,k))*r(3)
              temp = g(2,k) - temp
              b(1,k) = b(1,k)  +  temp * c(1)
              b(2,k) = b(2,k)  +  temp
           end do
        case (4)
           do k = 4, nv
              temp = g(3,k)
              gk = a(k) - a1(k)
              g(3,k) = ((gk*r(4) - g(1,k))*r(5) - g(2,k))*r(6)
              temp = g(3,k) - temp
              b(1,k) = b(1,k)  +  temp * c(2)
              b(2,k) = b(2,k)  +  temp * c(3)
              b(3,k) = b(3,k)  +  temp
           end do
        case (5)
           do k = 4, nv
              temp = g(4,k)
              gk = a(k) - a1(k)
              g(4,k) = (((gk*r(7) - g(1,k))*r(8) - g(2,k))*r(9)- g(3,k))*r(10)
              temp = g(4,k) - temp
              b(1,k) = b(1,k)  +  temp * c(4)
              b(2,k) = b(2,k)  +  temp * c(5)
              b(3,k) = b(3,k)  +  temp * c(6)
              b(4,k) = b(4,k)  +  temp
           end do
        case (6)
           do k = 4, nv
              temp = g(5,k)
              gk = a(k) - a1(k)
              g(5,k) =  ((((gk*r(11) - g(1,k))*r(12) - g(2,k))*r(13)- g(3,k))*r(14) - g(4,k))*r(15)
              temp = g(5,k) - temp
              b(1,k) = b(1,k)  +  temp * c(7)
              b(2,k) = b(2,k)  +  temp * c(8)
              b(3,k) = b(3,k)  +  temp * c(9)
              b(4,k) = b(4,k)  +  temp * c(10)
              b(5,k) = b(5,k)  +  temp
           end do
        case (7)
           do k = 4, nv
              temp = g(6,k)
              gk = a(k) - a1(k)
              g(6,k) = (((((gk*r(16) - g(1,k))*r(17) - g(2,k))*r(18)- g(3,k))*r(19) - g(4,k))*r(20) - g(5,k))*r(21)
              temp = g(6,k) - temp
              b(1,k) = b(1,k)  +  temp * c(11)
              b(2,k) = b(2,k)  +  temp * c(12)
              b(3,k) = b(3,k)  +  temp * c(13)
              b(4,k) = b(4,k)  +  temp * c(14)
              b(5,k) = b(5,k)  +  temp * c(15)
              b(6,k) = b(6,k)  +  temp
           end do
        case (8)
           do k = 4, nv
              temp = g(7,k)
              gk = a(k) - a1(k)
              g(7,k) = ((((((gk*r(22) - g(1,k))*r(23) - g(2,k))*r(24) - g(3,k))*r(25) &
                   - g(4,k))*r(26) - g(5,k))*r(27) - g(6,k))*r(28)
              temp = g(7,k) - temp
              b(1,k) = b(1,k)  +  temp * c(16)
              b(2,k) = b(2,k)  +  temp * c(17)
              b(3,k) = b(3,k)  +  temp * c(18)
              b(4,k) = b(4,k)  +  temp * c(19)
              b(5,k) = b(5,k)  +  temp * c(20)
              b(6,k) = b(6,k)  +  temp * c(21)
              b(7,k) = b(7,k)  +  temp
           end do
        end select
     end do
  end do
  
  !------------------------------------------------------------------------------
  
  !  END  OF  MAIN  LOOP
  
  ! Estimate suitable sequence size for the next call to subroutine (Eqs. 15, 16)
  temp = 0.d0
  do k = 4, nv
     temp = max( temp, abs( b(7,k) ) )
  end do
  temp = temp / (72.d0 * abs(t)**7)
  tdid = t
  if (temp.eq.0) then
     t = tdid * 1.4d0
  else
     t = sign( (tol/temp)**(1.d0/9.d0), tdid )
  end if
  
  ! If sequence size for the first subroutine call is too big, go back and redo
  ! the sequence using a smaller size.
  if ((dtflag.eq.1).and.(abs(t/tdid).lt.1)) then
     t = t * .8d0
  else
    exit
  end if
  end do
  
  ! If new sequence size is much bigger than the current one, reduce it
  if (abs(t/tdid).gt.1.4d0) t = tdid * 1.4d0
  
  ! Find new position and velocity values at end of the sequence (Eqs. 11, 12)
  temp = tdid * tdid
  do k = 4 , nv
     x1(k) = (xc(8)*b(7,k) + xc(7)*b(6,k) + xc(6)*b(5,k) + xc(5)*b(4,k) + xc(4)*b(3,k) + xc(3)*b(2,k)&
          + xc(2)*b(1,k) + xc(1)*a1(k))*temp + v1(k)*tdid + x1(k)
     
     v1(k) = (vc(7)*b(7,k) + vc(6)*b(6,k) + vc(5)*b(5,k) + vc(4)*b(4,k) + vc(3)*b(3,k) + vc(2)*b(2,k)&
          + vc(1)*b(1,k) + a1(k))*tdid + v1(k)
  end do
  
  ! Predict new B values to use at the start of the next sequence. The predicted
  ! values from the last call are saved as E. The correction, BD, between the
  ! actual and predicted values of B is applied in advance as a correction.
  q = t / tdid
  q2 = q  * q
  q3 = q  * q2
  q4 = q2 * q2
  q5 = q2 * q3
  q6 = q3 * q3
  q7 = q3 * q4
  
  do k = 4, nv
     s(1) = b(1,k) - e(1,k)
     s(2) = b(2,k) - e(2,k)
     s(3) = b(3,k) - e(3,k)
     s(4) = b(4,k) - e(4,k)
     s(5) = b(5,k) - e(5,k)
     s(6) = b(6,k) - e(6,k)
     s(7) = b(7,k) - e(7,k)
     
     ! Estimate B values for the next sequence (Eqs. 13 of Everhart).
     e(1,k) = q* (b(7,k)* 7.d0 + b(6,k)* 6.d0 + b(5,k)* 5.d0 + b(4,k)* 4.d0 + b(3,k)* 3.d0 + b(2,k)*2.d0 + b(1,k))
     e(2,k) = q2*(b(7,k)*21.d0 + b(6,k)*15.d0 + b(5,k)*10.d0 + b(4,k)* 6.d0 + b(3,k)* 3.d0 + b(2,k))
     e(3,k) = q3*(b(7,k)*35.d0 + b(6,k)*20.d0 + b(5,k)*10.d0 + b(4,k)*4.d0 + b(3,k))
     e(4,k) = q4*(b(7,k)*35.d0 + b(6,k)*15.d0 + b(5,k)*5.d0 + b(4,k))
     e(5,k) = q5*(b(7,k)*21.d0 + b(6,k)*6.d0 + b(5,k))
     e(6,k) = q6*(b(7,k)*7.d0 + b(6,k))
     e(7,k) = q7* b(7,k)
     
     b(1,k) = e(1,k) + s(1)
     b(2,k) = e(2,k) + s(2)
     b(3,k) = e(3,k) + s(3)
     b(4,k) = e(4,k) + s(4)
     b(5,k) = e(5,k) + s(5)
     b(6,k) = e(6,k) + s(6)
     b(7,k) = e(7,k) + s(7)
  end do
  dtflag = 2
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mdt_ra15
  
end module algo_radau

