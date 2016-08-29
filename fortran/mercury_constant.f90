!******************************************************************************
! MODULE: mercury_constant
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Global constant of the Mercury Code, numerical constant that fix
!! some numerical problems.
!
!******************************************************************************

module mercury_constant
  use types_numeriques

implicit none

integer, parameter :: CMAX = 50 !< CMAX  = maximum number of close-encounter minima monitored simultaneously
integer, parameter :: NMESS = 200 !< NMESS = maximum number of messages in message.in
integer, parameter :: NFILES = 50 !< NFILES = maximum number of files that can be open at the same time
real(double_precision), parameter :: HUGE = 9.9d29 !< HUGE  = an implausibly large number
real(double_precision), parameter :: TINY = 4.D-15 !< A small number

!...   convergence criteria for danby
real(double_precision), parameter :: DANBYAC= 1.0d-14
real(double_precision), parameter :: DANBYB = 1.0d-13

!...    loop limits in the Laguerre attempts
integer, parameter :: NLAG1 = 50 
integer, parameter :: NLAG2 = 400

end module mercury_constant
