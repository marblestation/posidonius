!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      CODE_FORTRAN.F90    (ErikSoft   3 May 2002)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Mercury is a general-purpose N-body integration package for problems in
! celestial mechanics.


!------------------------------------------------------------------------------


program mercury
  
    use user_module
    use physical_constant
    use tides_constant_GR
    use orbital_elements

    implicit none
  
    integer :: j,nbod,nbig
    integer :: error

    real(double_precision) :: time,gm_tmp
    real(double_precision) :: q_tmp,e_tmp,i_tmp,p_tmp,n_tmp,l_tmp
    real(double_precision) :: x_tmp,y_tmp,z_tmp,u_tmp,v_tmp,w_tmp

    real(double_precision), dimension(:), allocatable :: m,sma,ecc,inc ! (Number of bodies)
    real(double_precision), dimension(:,:), allocatable :: xh,vh ! (3,Number of bodies)
    real(double_precision), dimension(:,:), allocatable :: x,v,a ! (3,Number of bodies)
  
    real(double_precision) :: time_step,time_limit
    real(double_precision) :: star_mass
    real(double_precision) :: planet_mass

    !------------------------------------------------------------------------------

    ! star + planets
    nbod = 2
    ! planet
    nbig = 1

    ! Initialization
    time        = 0.0d0
    allocate(m(nbod), stat=error)
    allocate(sma(nbod), stat=error)
    allocate(ecc(nbod), stat=error)
    allocate(inc(nbod), stat=error)
    m(1:nbod)   = 0.0d0
    sma(1:nbod) = 0.0d0
    ecc(1:nbod) = 0.0d0
    inc(1:nbod) = 0.0d0

    allocate(xh(3,nbod), stat=error)
    allocate(x(3,nbod), stat=error)
    xh(1,1:nbod) = 0.0d0
    xh(2,1:nbod) = 0.0d0
    xh(3,1:nbod) = 0.0d0
    allocate(vh(3,nbod), stat=error)
    allocate(v(3,nbod), stat=error)
    vh(1,1:nbod) = 0.0d0
    vh(2,1:nbod) = 0.0d0
    vh(3,1:nbod) = 0.0d0
    allocate(a(3,nbod), stat=error)

    !------------------------------------------------------------------------------
    !---- Simulation
    time_step = 0.08d0 ! in days
    time_limit = time_step * 365.25d0 * 10.0d2
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !  Initial conditions from CASE 3 in Bolmont et al. 2015
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !---- Star (central body)
    star_mass = 0.08d0 ! Solar masses

    !---------------------------------------------------------------------------
    !---- Planet
    planet_mass = 3.0d-6!/ Solar masses (3.0e-6 solar masses = 1 earth mass)

   !--------- Specify initial position and velocity for a stable orbit
   !----- Keplerian orbital elements, in the `asteroidal' format of Mercury code'
    sma(2) = 0.018d0        !/ semi-major axis (in AU)
    ecc(2) = 0.1d0          !/ eccentricity
    e_tmp = ecc(2)
    inc(2) = 5.d0*PI/180.d0           !/ inclination (degrees)
    i_tmp = inc(2)
    p_tmp = 0.d0              !/ argument of pericentre (degrees)
    n_tmp = 0.d0              !/ longitude of the ascending node (degrees)
    l_tmp = 0.d0              !/ mean anomaly (degrees)
    p_tmp = (p_tmp + n_tmp)             !/ Convert to longitude of perihelion !!
    q_tmp = sma(2) * (1.0d0 - e_tmp)       !/ perihelion distance
    gm_tmp = K2*(planet_mass+star_mass)
    !(x, y, z, vx, vy, vz) = posidonius::general::calculate_cartesian_coordinates(gm, q, e, i, p, n, l)

    ! reconciliation with mercury notations
    ! m is here actually gm 
    m(1) = K2*star_mass
    m(2) = K2*planet_mass

    ! Loop on planets
    call mco_el2x(gm_tmp,q_tmp,e_tmp,i_tmp,p_tmp,n_tmp,l_tmp,x_tmp,y_tmp,z_tmp,u_tmp,v_tmp,w_tmp)

    xh(1,2) = x_tmp
    xh(2,2) = y_tmp
    xh(3,2) = z_tmp
    vh(1,2) = u_tmp
    vh(2,2) = v_tmp
    vh(3,2) = w_tmp

    x(1,2) = xh(1,2)
    x(2,2) = xh(2,2)
    x(3,2) = xh(3,2)
    v(1,2) = vh(1,2)
    v(2,2) = vh(2,2)
    v(3,2) = vh(3,2)

    do j = 0, 5 
        call mfo_user(time,nbod,nbig,time_step,m,x,v,a)
        time = time + time_step
    enddo

    !time = time + dt
end program mercury
