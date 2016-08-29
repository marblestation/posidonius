!******************************************************************************
! MODULE: mercury_globals
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Modules that contains all the globals variables of mercury
!
!******************************************************************************

module mercury_globals

  use types_numeriques
  use mercury_constant

  implicit none
  
  integer :: nb_bodies_initial !< number of bodies when we start the simulation. Used for local variables in several modules. 
  integer, dimension(8) :: opt = (/0,1,1,2,0,1,0,0/) !< Default options (can be overwritten later in the code) for mercury.\n
!!\n  OPT(1) = close-encounter option (0=stop after an encounter, 1=continue)
!!\n  OPT(2) = collision option (0=no collisions, 1=merge, 2=merge+fragment)
!!\n  OPT(3) = time style (0=days 1=Greg.date 2/3=days/years w/respect to start)
!!\n  OPT(4) = o/p precision (1,2,3 = 4,9,15 significant figures)
!!\n  OPT(5) = < Not used at present >
!!\n  OPT(6) = < Not used at present >
!!\n  OPT(7) = apply post-Newtonian correction? (0=no, 1=yes)
!!\n  OPT(8) = apply user-defined force routine mfo_user? (0=no, 1=yes)
  
  character(len=80), dimension(NMESS) :: mem !< Various messages and strings used by mercury
  integer, dimension(NMESS) :: lmem !< the length of each string of the 'mem' elements
  
  character(len=80), dimension(3) :: outfile !< filenames for output files (*.out)
!!\n  OUTFILE  (1) = osculating coordinates/velocities and masses
!!\n  OUTFILE  (2) = close encounter details
!!\n  OUTFILE  (3) = information file

  character(len=80), dimension(4) :: dumpfile !< filenames for dump files (*.dmp)
!!\n  DUMPFILE (1) = Big-body data
!!\n  DUMPFILE (2) = Small-body data
!!\n  DUMPFILE (3) = integration parameters
!!\n  DUMPFILE (4) = restart file  
  
  integer :: algor !< An index that represent the algorithm used. \n
!!\n  ALGOR = 1  ->  Mixed-variable symplectic
!!\n          2  ->  Bulirsch-Stoer integrator
!!\n          3  ->         "           "      (conservative systems only)
!!\n          4  ->  RA15 `radau' integrator
!!\n          10 ->  Hybrid MVS/BS (democratic-heliocentric coords)
!!\n          11 ->  Close-binary hybrid (close-binary coords)
!!\n          12 ->  Wide-binary hybrid (wide-binary coords)
  
  real(double_precision) :: tstart !< epoch of first required output (days)
  real(double_precision) :: tstop !< epoch final required output (days)
  real(double_precision) :: tdump !< time of next data dump (days)
  
  real(double_precision) :: dtout !< data output interval           (days)
  real(double_precision) :: dtdump !< data-dump interval             (days)
  real(double_precision) :: dtfun !< interval for other periodic effects (e.g. check for ejections)

  real(double_precision) :: rmax !< heliocentric distance at which objects are considered ejected (AU)
  
  
end module mercury_globals
