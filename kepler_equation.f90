!******************************************************************************
! MODULE: kepler_equation
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Modules that solve the kepler equation
!
!******************************************************************************

module kepler_equation

  use types_numeriques

  implicit none
  
  contains

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   !> @author 
   !> John E. Chambers
   !
   !> @date 7 July 1999
   !
   ! DESCRIPTION: 
   !> @brief Solves Kepler's equation for eccentricities less than one.
   !! Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.
   !!\n\n
   !!  e = eccentricity\n
   !!  l = mean anomaly      (radians)\n
   !!  u = eccentric anomaly (   "   )\n
   !
   !> @return mco_kep blablabla
   !
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function mco_kep (e,oldl)

  implicit none

  
  ! Input/Outout
  real(double_precision), intent(in) :: oldl
  real(double_precision), intent(in) :: e
  real(double_precision) :: mco_kep
  
  ! Local
  real(double_precision) :: l,pi,twopi,piby2,u1,u2,ome,sign
  real(double_precision) :: x,x2,sn,dsn,z1,z2,z3,f0,f1,f2,f3
  real(double_precision) :: p,q,p2,ss,cc
  logical flag,big,bigg
  
  !------------------------------------------------------------------------------
  
  pi = 3.141592653589793d0
  twopi = 2.d0 * pi
  piby2 = .5d0 * pi
  
  ! Reduce mean anomaly to lie in the range 0 < l < pi
  if (oldl.ge.0) then
     l = mod(oldl, twopi)
  else
     l = mod(oldl, twopi) + twopi
  end if
  sign = 1.d0
  if (l.gt.pi) then
     l = twopi - l
     sign = -1.d0
  end if
  
  ome = 1.d0 - e
  
  if (l.ge..45d0.or.e.lt..55d0) then
     
     ! Regions A,B or C in Nijenhuis
     ! -----------------------------
     
     ! Rough starting value for eccentric anomaly
     if (l.lt.ome) then
        u1 = ome
     else
        if (l.gt.(pi-1.d0-e)) then
           u1 = (l+e*pi)/(1.d0+e)
        else
           u1 = l + e
        end if
     end if
     
     ! Improved value using Halley's method
     flag = u1.gt.piby2
     if (flag) then
        x = pi - u1
     else
        x = u1
     end if
     x2 = x*x
     sn = x*(1.d0 + x2*(-.16605 + x2*.00761) )
     dsn = 1.d0 + x2*(-.49815 + x2*.03805)
     if (flag) dsn = -dsn
     f2 = e*sn
     f0 = u1 - f2 - l
     f1 = 1.d0 - e*dsn
     u2 = u1 - f0/(f1 - .5d0*f0*f2/f1)
  else
     
     ! Region D in Nijenhuis
     ! ---------------------
     
     ! Rough starting value for eccentric anomaly
     z1 = 4.d0*e + .5d0
     p = ome / z1
     q = .5d0 * l / z1
     p2 = p*p
     z2 = exp( log( dsqrt( p2*p + q*q ) + q )/1.5 )
     u1 = 2.d0*q / ( z2 + p + p2/z2 )
     
     ! Improved value using Newton's method
     z2 = u1*u1
     z3 = z2*z2
     u2 = u1 - .075d0*u1*z3 / (ome + z1*z2 + .375d0*z3)
     u2 = l + e*u2*( 3.d0 - 4.d0*u2*u2 )
  end if
  
  ! Accurate value using 3rd-order version of Newton's method
  ! N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy
  
  ! First get accurate values for u2 - sin(u2) and 1 - cos(u2)
  bigg = (u2.gt.piby2)
  if (bigg) then
     z3 = pi - u2
  else
     z3 = u2
  end if
  
  big = (z3.gt.(.5d0*piby2))
  if (big) then
     x = piby2 - z3
  else
     x = z3
  end if
  
  x2 = x*x
  ss = 1.d0
  cc = 1.d0
  
  ss = x*x2/6.*(1. - x2/20.*(1. - x2/42.*(1. - x2/72.*(1. - x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))))
  cc =   x2/2.*(1. - x2/12.*(1. - x2/30.*(1. - x2/56.*(1. - x2/ 90.*(1. - x2/132.*(1. - x2/182.*(1. - x2/240.*(1. - x2/306.))))))))
  
  if (big) then
     z1 = cc + z3 - 1.d0
     z2 = ss + z3 + 1.d0 - piby2
  else
     z1 = ss
     z2 = cc
  end if
  
  if (bigg) then
     z1 = 2.d0*u2 + z1 - pi
     z2 = 2.d0 - z2
  end if
  
  f0 = l - u2*ome - e*z1
  f1 = ome + e*z2
  f2 = .5d0*e*(u2-z1)
  f3 = e/6.d0*(1.d0-z2)
  z1 = f0/f1
  z2 = f0/(f2*z1+f1)
  mco_kep = sign*( u2 + f0/((f3*z1+f2)*z2+f1) )
  
  !------------------------------------------------------------------------------
  
  return
end function mco_kep

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> M. Duncan 
!
!> @date May 26, 1992
!
! DESCRIPTION: 
!> @brief Solves Kepler's eqn. for hyperbola using hybrid approach. \n\n
!!
!!ALGORITHM: Uses power series for N in terms of F and Newton,s method
!
!> @warning REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function orbel_flon(e,capn0)

  use mercury_constant

  implicit none


  !...  Inputs Only: 
  real(double_precision), intent(in) :: e !< [in] 
  real(double_precision), intent(in) :: capn0 !< [in] 
  
  ! Output
  real(double_precision) :: orbel_flon !< [out] 

  !...  Internals:
  real(double_precision) :: capn !< copy of the input capn0 that is not modified
  integer :: iflag,i
  real(double_precision) :: a,b,sq,biga,bigb
  real(double_precision) :: x,x2
  real(double_precision) :: f,fp,dx
  real(double_precision) :: diff
  real(double_precision) :: a0,a1
  real(double_precision) :: b1
  
  integer, parameter :: IMAX = 10
  real(double_precision), parameter :: a3 = 1037836800.d0
  real(double_precision), parameter :: a5 = 51891840.d0
  real(double_precision), parameter :: a7 = 1235520.d0
  real(double_precision), parameter :: a9 = 17160.d0
  real(double_precision), parameter :: a11 = 156.d0
  
  real(double_precision), parameter :: b3 = 3.d0*a3
  real(double_precision), parameter :: b5 = 5.d0*a5
  real(double_precision), parameter :: b7 = 7.d0*a7
  real(double_precision), parameter :: b9 = 9.d0*a9
  real(double_precision), parameter :: b11 = 11.d0*a11

  !----
  !...  Executable code 

  ! We copy the input parameter to be able to modify it.
  capn = capn0
  
  ! Function to solve "Kepler's eqn" for F (here called
  ! x) for given e and CAPN. Only good for smallish CAPN 

  iflag = 0
  if( capn .lt. 0.d0) then
     iflag = 1
     capn = -capn
  endif

  a1 = 6227020800.d0 * (1.d0 - 1.d0/e)
  a0 = -6227020800.d0*capn/e
  b1 = a1

  !  Set iflag nonzero if capn < 0., in which case solve for -capn
  !  and change the sign of the final answer for F.
  !  Begin with a reasonable guess based on solving the cubic for small F  


  a = 6.d0*(e-1.d0)/e
  b = -6.d0*capn/e
  sq = sqrt(0.25*b*b +a*a*a/27.d0)
  biga = (-0.5*b + sq)**0.3333333333333333d0
  bigb = -(+0.5*b + sq)**0.3333333333333333d0
  x = biga + bigb
  !  write(6,*) 'cubic = ',x**3 +a*x +b
  orbel_flon = x
  ! If capn is tiny (or zero) no need to go further than cubic even for
  ! e =1.
  if( capn .lt. TINY) go to 100

  do i = 1,IMAX
     x2 = x*x
     f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))))
     fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.d0*x2)))))   
     dx = -f/fp
     !    write(6,*) 'i,dx,x,f : '
     !    write(6,432) i,dx,x,f
     orbel_flon = x + dx
     !   If we have converged here there's no point in going on
     if(abs(dx) .le. TINY) go to 100
     x = orbel_flon
  enddo

  ! Abnormal return here - we've gone thru the loop 
  ! IMAX times without convergence
  if(iflag .eq. 1) then
     orbel_flon = -orbel_flon
     capn = -capn
  endif
  write(6,*) 'FLON : RETURNING WITHOUT COMPLETE CONVERGENCE' 
  diff = e*sinh(orbel_flon) - orbel_flon - capn
  write(6,*) 'N, F, ecc*sinh(F) - F - N : '
  write(6,*) capn,orbel_flon,diff
  return

  !  Normal return here, but check if capn was originally negative
100 if(iflag .eq. 1) then
     orbel_flon = -orbel_flon
     capn = -capn
  endif

  return
end function orbel_flon

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> M. Duncan 
!
!> @date May 11, 1992.
!
! DESCRIPTION: 
!> @brief Solves Kepler's eqn. for hyperbola using hybrid approach. 
!!\n\n
!!ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
!!\n           Cel. Mech. ".  Quartic convergence from Danby's book.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function orbel_fget(e,capn)

  use mercury_constant

  implicit none


  !...  Inputs Only: 
  real(double_precision), intent(in) :: e !< [in] eccentricity
  real(double_precision), intent(in) :: capn !< [in] hyperbola mean anomaly
  
  ! Output
  real(double_precision) :: orbel_fget !< [out] eccentric anomaly

  !...  Internals:
  integer :: i,IMAX
  real(double_precision) :: tmp,x,shx,chx
  real(double_precision) :: esh,ech,f,fp,fpp,fppp,dx
  PARAMETER (IMAX = 10)

  !----
  !...  Executable code 

  ! Function to solve "Kepler's eqn" for F (here called
  ! x) for given e and CAPN. 

  !  begin with a guess proposed by Danby  
  if( capn .lt. 0.d0) then
     tmp = -2.d0*capn/e + 1.8d0
     x = -log(tmp)
  else
     tmp = +2.d0*capn/e + 1.8d0
     x = log( tmp)
  endif

  orbel_fget = x

  do i = 1,IMAX
    shx = sinh(x)
    chx = cosh(x)
!~      call mco_sinh (x,shx,chx)
     esh = e*shx
     ech = e*chx
     f = esh - x - capn
     !    write(6,*) 'i,x,f : ',i,x,f
     fp = ech - 1.d0  
     fpp = esh 
     fppp = ech 
     dx = -f/fp
     dx = -f/(fp + dx*fpp/2.d0)
     dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
     orbel_fget = x + dx
     !   If we have converged here there's no point in going on
     if(abs(dx) .le. TINY) RETURN
     x = orbel_fget
  enddo

  write(6,*) 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE' 
  return
end function orbel_fget

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> M. Duncan 
!
!> @date May 26,1992.
!
! DESCRIPTION: 
!> @brief Solves Kepler's eqn. for hyperbola using hybrid approach. \n\n
!!
!! ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON 
!!\n           For larger N, uses FGET
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function orbel_fhybrid(e,n)

  use mercury_constant

  implicit none


  !...  Inputs Only: 
  real(double_precision), intent(in) :: e !< [in] eccentricity
  real(double_precision), intent(in) :: n !< [in] hyperbola mean anomaly
  
  ! Output
  real(double_precision) :: orbel_fhybrid !< [out] eccentric anomaly
  
  !...  Internals:
  real(double_precision) :: abn

  !----
  !...  Executable code 

  abn = n
  if(n.lt.0.d0) abn = -abn

  if(abn .lt. 0.636d0*e -0.6d0) then
     orbel_fhybrid = orbel_flon(e,n)
  else 
     orbel_fhybrid = orbel_fget(e,n)
  endif

  return
end function orbel_fhybrid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> M. Duncan 
!
!> @date May 11, 1992
!
! DESCRIPTION: 
!> @brief Solves the equivalent of Kepler's eqn. for a parabola 
!! given Q (Fitz. notation.)\n\n
!!
!!ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
!
!> @remarks For a parabola we can solve analytically.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function orbel_zget(q0)

  use mercury_constant

  implicit none


  !...  Inputs Only: 
  real(double_precision), intent(in) :: q0 !< [in] parabola mean anomaly
  
  ! Output
  real(double_precision) :: orbel_zget !< [out] eccentric anomaly

  !...  Internals:
  real(double_precision) :: q !< copy of q0 to modify it without propagating the modified value
  integer :: iflag
  real(double_precision) :: x,tmp


  !----
  !...  Executable code 

  ! We copy the input parameter to be able to modify it.
  q = q0

  iflag = 0
  if(q.lt.0.d0) then
     iflag = 1
     q = -q
  endif

  if (q.lt.1.d-3) then
     orbel_zget = q*(1.d0 - (q*q/3.d0)*(1.d0 -q*q))
  else
     x = 0.5d0*(3.d0*q + sqrt(9.d0*(q**2) +4.d0))
     tmp = x**(1.d0/3.d0)
     orbel_zget = tmp - 1.d0/tmp
  endif

  if(iflag .eq.1) then
     orbel_zget = -orbel_zget
     q = -q
  endif

  return
end function orbel_zget    ! orbel_zget



end module kepler_equation
