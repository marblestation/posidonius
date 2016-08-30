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

    implicit none
  
    integer :: j,nbod,nbig

    real(double_precision) :: time,gm_tmp,dt
    real(double_precision) :: q_tmp,e_tmp,i_tmp,p_tmp,n_tmp,l_tmp
    real(double_precision) :: x_tmp,y_tmp,z_tmp,u_tmp,v_tmp,w_tmp

    real(double_precision), dimension(:), allocatable :: m,sma,ecc,inc ! (Number of bodies)
    real(double_precision), dimension(:,:), allocatable :: xh,vh ! (3,Number of bodies)
    real(double_precision), dimension(:,:), allocatable :: x,v,a ! (3,Number of bodies)
  
    real(double_precision) :: time_step,time_limit
    real(double_precision) :: star_mass, radius_factor, star_radius, star_love_number
    real(double_precision) :: star_dissipation_factor_scale,star_dissipation_factor
    real(double_precision) :: star_radius_of_gyration_2,planet_dissipation_factor_scale
    real(double_precision) :: star_position,star_velocity,star_acceleration
    real(double_precision) :: star_rotation_period, star_spin0, star_spin

    real(double_precision) :: planet_mass,planet_radius_factor,planet_radius
    real(double_precision) :: planet_love_number,k2pdelta,planet_dissipation_factor
    real(double_precision) :: planet_radius_of_gyration_2

    !------------------------------------------------------------------------------

    ! star + planets
    nbod = 2
    ! planet
    nbig = 1

    ! Initialization
    time        = 0.0d0
    m(1:nbod)   = 0.0d0
    sma(1:nbod) = 0.0d0
    ecc(1:nbod) = 0.0d0
    inc(1:nbod) = 0.0d0

    xh(1,1:nbod) = 0.0d0
    xh(2,1:nbod) = 0.0d0
    xh(3,1:nbod) = 0.0d0
    xh(1,1:nbod) = 0.0d0
    xh(2,1:nbod) = 0.0d0
    xh(3,1:nbod) = 0.0d0

    !------------------------------------------------------------------------------
    !---- Simulation
    time_step = 0.08 ! in days
    time_limit = time_step * 365.25 * 10.0e2
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !  Initial conditions from CASE 3 in Bolmont et al. 2015
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    !---- Star (central body)
    star_mass = 0.08 ! Solar masses
    radius_factor = 0.845649342247916
    star_radius = radius_factor * rsun
    star_love_number = 0.307!/ M Dwarf
    !---- Disipation factor (sigma)
    star_dissipation_factor_scale = 1.
    star_dissipation_factor = star_dissipation_factor_scale * 2.006*3.845764e4!/ -60+64
    !----- Radius of gyration
    star_radius_of_gyration_2 = 1.94e-1!/ Brown dwarf
    !- To calculate tidal forces, it is needed to have the central body/star at [0,0,0] and without velocity or acceleration (heliocentric)
    star_position     = 0.0d0 
    star_velocity     = 0.0d0 
    star_acceleration = 0.0d0 
    !----- Initialization of stellar spin
    star_rotation_period = 70.0!/ days
    star_spin0 = TWOPI/star_rotation_period!/ days^-1
    star_spin  = 0. 

    !---------------------------------------------------------------------------
    !---- Planet
    planet_mass = 3.0e-6!/ Solar masses (3.0e-6 solar masses = 1 earth mass)
    !----- Planetary radius in AU (rearth in AU) Rocky planet
    !  For other planet types, specify a factor:
    planet_radius_factor = 1.
    planet_radius = planet_radius_factor * rearth
    planet_love_number = 0.305!/ Earth

    !----- Disipation factor (sigma)
    planet_dissipation_factor_scale = 1.
    ! --- Hot Gas Giant:
    ! planet_dissipation_factor = planet_dissipation_factor_scale * 2.006*3.845764d4
    ! Terrestrial:
    k2pdelta = 2.465278e-3!/ Terrestrial planets (no gas)
    planet_dissipation_factor = planet_dissipation_factor_scale * 2.d0 * K2 * k2pdelta/(3.d0 * planet_radius**5)
    ! Radius of gyration
    planet_radius_of_gyration_2 = 3.308e-1!/ Earth type planet
    ! planet_radius_of_gyration_2 = 2.54e-1!/ Gas planet

   !--------- Specify initial position and velocity for a stable orbit
   !----- Keplerian orbital elements, in the `asteroidal' format of Mercury code'
    sma(2) = 0.018d0        !/ semi-major axis (in AU)
    ecc(2) = 0.1d0          !/ eccentricity
    e_tmp = ecc(2)
    inc(2) = 5.d0           !/ inclination (degrees)
    i_tmp = inc(j)
    p_tmp = 0.              !/ argument of pericentre (degrees)
    n_tmp = 0.              !/ longitude of the ascending node (degrees)
    l_tmp = 0.              !/ mean anomaly (degrees)
    p_tmp = (p_tmp + n_tmp)             !/ Convert to longitude of perihelion !!
    q_tmp = sma(2) * (1.0 - e_tmp)       !/ perihelion distance
    gm_tmp = K2*(planet_mass+star_mass)
    !(x, y, z, vx, vy, vz) = posidonius::general::calculate_cartesian_coordinates(gm, q, e, i, p, n, l)

    ! reconciliation with mercury notations
    ! m is here actually gm 
    m(1) = K2*star_mass
    m(2) = K2*planet_mass

    ! Loop on planets
    call mco_el2x(gm_tmp,q_tmp,e_tmp,i_tmp,p_tmp,n_tmp,l_tmp,x_tmp,y_tmp,z_tmp,u_tmp,v_tmp,w_tmp)

    xh(1,j) = x_tmp
    xh(2,j) = y_tmp
    xh(3,j) = z_tmp
    vh(1,j) = u_tmp
    vh(2,j) = v_tmp
    vh(3,j) = w_tmp

    x(1,j) = xh(1,j)
    x(2,j) = xh(2,j)
    x(3,j) = xh(3,j)
    v(1,j) = vh(1,j)
    v(2,j) = vh(2,j)
    v(3,j) = vh(3,j)

    call mfo_user(time,nbod,nbig,m,x,v,a)
    write(*,*) a(1,2), a(2,2), a(3,2)
    write(*,*) a(1,3), a(2,3), a(3,3)

    !time = time + dt
end program mercury
