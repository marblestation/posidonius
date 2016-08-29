!******************************************************************************
! MODULE: types_numeriques
!******************************************************************************
!
! DESCRIPTION: 
!> @brief To define several parameter types, in particular to define
!! real simple and double precision
!
!******************************************************************************
module types_numeriques

  implicit none
  
  integer, parameter :: simple_precision = selected_real_kind(6,37) !< to define real simple precision
  integer, parameter :: double_precision = selected_real_kind(15,307) !< to define real double precision
  integer, parameter :: sp = simple_precision !< alternate way to define simple precision
  integer, parameter :: dp = double_precision !< alternate way to define double precision

end module types_numeriques
