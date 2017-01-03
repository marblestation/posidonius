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
    real(double_precision), dimension(:,:), allocatable :: x,v,spin,a,Torque ! (3,Number of bodies)
    real(double_precision), dimension(:,:), allocatable :: a_grav,a_tides

    real(double_precision) :: time_step,time_limit,half_time_step,time_write
    real(double_precision) :: output_tmp

    real(double_precision) :: star_mass,spin0,Pst
    real(double_precision) :: planet_mass,spinp0
    real(double_precision) :: rg2s,k2s,k2fs,Rsth,Rsth5,Rsth10,sigmast
    real(double_precision), dimension(:),allocatable :: Rp,Rp5,Rp10,sigmap
    real(double_precision), dimension(:),allocatable :: tintin,k2p,k2fp,k2pdeltap,rg2p
    real(double_precision) :: tmp, tmp1

    !------------------------------------------------------------------------------

    ! star + planets
    nbod = 2
    ! planet
    nbig = 1

    ! Initialization
    time        = 0.0d0
    time_write  = 0.0d0
    !------------------------------------------------------------------------------
    !---- Simulation
    time_step = 0.08d0 ! in days
    !time_limit = 365.25d0 * 1.d8
    time_limit = time_step*4.
    output_tmp = output
    !------------------------------------------------------------------------------


    allocate(Rp(nbod), stat=error)
    allocate(Rp5(nbod), stat=error)
    allocate(Rp10(nbod), stat=error)
    Rp(1:nbod)    = 0.0d0
    Rp5(1:nbod)   = 0.0d0
    Rp10(1:nbod)  = 0.0d0

    allocate(sigmap(nbod), stat=error)
    sigmap(1:nbod)= 0.0d0

    allocate(rg2p(nbod-1), stat=error)
    allocate(k2p(nbod-1), stat=error)
    allocate(k2fp(nbod-1), stat=error)
    allocate(k2pdeltap(nbod-1), stat=error)
    rg2p(1:nbod-1)  = 0.0d0
    k2p(1:nbod-1)   = 0.0d0
    k2fp(1:nbod-1)  = 0.0d0
    k2pdeltap(1:nbod-1) = 0.0d0

    allocate(m(nbod), stat=error)
    allocate(sma(nbod), stat=error)
    allocate(ecc(nbod), stat=error)
    allocate(inc(nbod), stat=error)
    m(1:nbod)   = 0.0d0
    sma(1:nbod) = 0.0d0
    ecc(1:nbod) = 0.0d0
    inc(1:nbod) = 0.0d0

    allocate(x(3,nbod), stat=error)
    x(1,1:nbod) = 0.0d0
    x(2,1:nbod) = 0.0d0
    x(3,1:nbod) = 0.0d0

    allocate(xh(3,nbod), stat=error)
    xh(1,1:nbod) = 0.0d0
    xh(2,1:nbod) = 0.0d0
    xh(3,1:nbod) = 0.0d0

    allocate(v(3,nbod), stat=error)
    v(1,1:nbod) = 0.0d0
    v(2,1:nbod) = 0.0d0
    v(3,1:nbod) = 0.0d0

    allocate(vh(3,nbod), stat=error)
    vh(1,1:nbod) = 0.0d0
    vh(2,1:nbod) = 0.0d0
    vh(3,1:nbod) = 0.0d0

    allocate(a_tides(3,nbod), stat=error)
    a_tides(1,1:nbod) = 0.0d0
    a_tides(2,1:nbod) = 0.0d0
    a_tides(3,1:nbod) = 0.0d0

    allocate(a_grav(3,nbod), stat=error)
    a_grav(1,1:nbod) = 0.0d0
    a_grav(2,1:nbod) = 0.0d0
    a_grav(3,1:nbod) = 0.0d0

    allocate(a(3,nbod), stat=error)
    a(1,1:nbod) = 0.0d0
    a(2,1:nbod) = 0.0d0
    a(3,1:nbod) = 0.0d0

    allocate(Torque(3,nbod), stat=error)
    Torque(1,1:nbod) = 0.0d0
    Torque(2,1:nbod) = 0.0d0
    Torque(3,1:nbod) = 0.0d0

    allocate(spin(3,nbod), stat=error)
    spin(1,1:nbod) = 0.0d0
    spin(2,1:nbod) = 0.0d0
    spin(3,1:nbod) = 0.0d0

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

    !---------------------------------------------------------------
    !-------------------------  HOST BODY  -------------------------
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    !-----------------------  NON EVOLVING -------------------------
    if (Rscst.eq.1) then 
        ! radius of gyration    
        rg2s  = rg2_what
        ! potential love number
        k2s   = k2st_what
        ! fluid love number
        k2fs  = k2fst_what

        ! Value of initial rotation period
        Pst = Period_st

        ! Value of the radius 
        Rsth   = radius_star*rsun
        Rsth5  = Rsth*Rsth*Rsth*Rsth*Rsth
        Rsth10 = Rsth5*Rsth5

        ! Dissipation
        if (sigma_what.gt.0) sigmast = dissstar*sigma_what
        if (k2sdeltats.gt.0) sigmast = dissstar*2.d0*K2*k2sdeltats/(3.d0*Rsth5)
    endif


    !---------------------------------------------------------------
    !--------------------------  PLANETS  --------------------------
    !---------------------------------------------------------------
    do j=2,ntid+1

        !-----------------------------------------------------------
        !------------------  Planetary radius in AU  ---------------

        ! If planet_type eq 0, rocky planet with prescription mass-radius
        if (planet_type(j-1).eq.0) Rp(j) = rearth*((0.0592d0*0.7d0+0.0975d0) &
             *(dlog10(m(j))+dlog10(m2earth)-dlog10(K2))**2+(0.2337d0*0.7d0+0.4938d0) &
             *(dlog10(m(j))+dlog10(m2earth)-dlog10(K2))+0.3102d0*0.7d0+0.7932d0)
        ! If planet_type ne 0, value given by radius_p of tides_constants
        if (planet_type(j-1).ne.0) Rp(j) = radius_p(j-1)*rearth

        Rp5(j)  = Rp(j)*Rp(j)*Rp(j)*Rp(j)*Rp(j)
        Rp10(j) = Rp5(j)*Rp5(j)

        !-----------------------------------------------------------
        !------------------  Planetary radius of gyration  ---------
        !---------------  love number and fluid love number  -------
        !------------------  planetary dissipation  ----------------

        ! planet_type = 0 or 1: terrestrial planet
        if ((planet_type(j-1).eq.0).or.(planet_type(j-1).eq.1)) then
            rg2p(j-1) = rg2p_terr
            k2p(j-1) = k2p_terr
            k2fp(j-1) = k2fp_terr
            if (tides.eq.1) then
                k2pdeltap(j-1) = k2pdeltap_terr
                sigmap(j) = dissplan(j-1)*2.d0*K2*k2pdeltap(j-1)/(3.d0*Rp5(j))
            endif
        endif

        ! planet_type = 2: gas giant (Jupiter)
        if (planet_type(j-1).eq.2) then
            rg2p(j-1) = rg2p_gg
            k2p(j-1) = k2p_gg
            k2fp(j-1) = k2p_gg
            if (tides.eq.1) then
                k2pdeltap(j-1) = k2pdeltap_gg
                sigmap(j) = dissplan(j-1)*sigma_gg
            endif
        endif

        ! planet_type = 3: you give the values you want in tides_constant_GR
        if (planet_type(j-1).eq.3) then
            rg2p(j-1) = rg2p_what(j-1)
            k2p(j-1) = k2tp_what(j-1)
            k2fp(j-1) = k2fp_what(j-1)
            if (tides.eq.1) then
                k2pdeltap(j-1) = k2pdeltap_what(j-1)
                sigmap(j) = dissplan(j-1)*2.d0*K2*k2pdeltap(j-1)/(3.d0*Rp5(j))
            endif
        endif
    enddo

    !---------------------------------------------------------------
    !-------------  Initialization of stellar spin (day-1)  --------
    if (crash.eq.0) then
        if (Rscst.eq.1) spin0 = TWOPI/Pst 
        spin(1,1) = 0.d0
        spin(2,1) = 0.d0
        spin(3,1) = spin0
    else
        spin(1,1) = rot_crash(1) 
        spin(2,1) = rot_crash(2) 
        spin(3,1) = rot_crash(3) 
    endif

    !--------------------------------------------------------------------
    !------------  Initialization of spin of planets (day-1)  -----------
    if (crash.eq.0) then

        ! rotation period given in
        ! tides_constants
        spinp0 = 24.d0*TWOPI/Pp0(1)

        ! If the planet has no inclination
        if (inc(2).eq.0.0) then
            spin(1,2) = 0.0d0
            spin(2,2) = -spinp0*sin(oblp(1))
            spin(3,2) = spinp0*cos(oblp(1)) 
        endif

        ! If the planet has an inclination
        if (inc(2).ne.0.0) then
            spin(1,2) = 0.0d0 
            spin(2,2) = -spinp0*sin(oblp(1)+inc(2))
            spin(3,2) = spinp0*cos(oblp(1)+inc(2)) 
        endif
    else

        spin(1,2) = rot_crashp1(1) !day-1
        spin(2,2) = rot_crashp1(2) !day-1
        spin(3,2) = rot_crashp1(3) !day-1 

    endif

    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    ! Integration !!!

    do while (time.le.time_limit)
        half_time_step = 0.5d0*time_step

        !write(*,*) 'Time',time
        call mfo_user(time,nbod,nbig,time_step,Rsth,Rsth5,Rsth10,rg2s,k2s,k2fs &
            ,sigmast,Rp,Rp5,Rp10,k2p,k2fp,rg2p,sigmap,m,x,v,spin,a_tides,Torque)
        !Torque(1,1:nbod) = 0.0d0
        !Torque(2,1:nbod) = 0.0d0
        !Torque(3,1:nbod) = 0.0d0
        !write(*,*) "a_tides",a_tides
        !write(*,*) "Torque ",Torque

        ! Leap frog, kick + spin
        ! Position of star
        !write(*,*) "star pos1",x(1,1),x(2,1),x(3,1)
        x(1,1) = x(1,1) + half_time_step * v(1,1)
        x(2,1) = x(2,1) + half_time_step * v(2,1)
        x(3,1) = x(3,1) + half_time_step * v(3,1)
        !write(*,*) "star pos2",x(1,1),x(2,1),x(3,1)

        ! Position of planet
        !write(*,*) "plan vel1",v(1,2),v(2,2),v(3,2)
        !write(*,*) "plan pos1",x(1,2),x(2,2),x(3,2)
        !write(*,*) "plan dis1",sqrt(x(1,2)*x(1,2)+x(2,2)*x(2,2)+x(3,2)*x(3,2))
        x(1,2) = x(1,2) + half_time_step * v(1,2)
        x(2,2) = x(2,2) + half_time_step * v(2,2)
        x(3,2) = x(3,2) + half_time_step * v(3,2)
        !write(*,*) "plan pos2",x(1,2),x(2,2),x(3,2)
        !write(*,*) "plan dis2",sqrt(x(1,2)*x(1,2)+x(2,2)*x(2,2)+x(3,2)*x(3,2))

        ! Spin of star
        !write(*,*) "star spin1",spin(1,1),spin(2,1),spin(3,1)
        tmp = -1.d0*K2/(m(1)*rg2s*Rsth*Rsth)
        spin(1,1) = spin(1,1) + half_time_step * tmp * Torque(1,1)
        spin(2,1) = spin(2,1) + half_time_step * tmp * Torque(2,1)
        spin(3,1) = spin(3,1) + half_time_step * tmp * Torque(3,1)
        !write(*,*) "star spin2",spin(1,1),spin(2,1),spin(3,1)
        !write(*,*) "dSpindt_s",tmp * Torque(1,1), tmp * Torque(2,1), tmp * Torque(3,1)

        ! Spin of planet
        !write(*,*) "plan spin1",spin(1,2),spin(2,2),spin(3,2)
        tmp = -1.d0*K2/(m(2)*rg2p(1)*Rp(2)*Rp(2))
        spin(1,2) = spin(1,2) + half_time_step * tmp * Torque(1,2)
        spin(2,2) = spin(2,2) + half_time_step * tmp * Torque(2,2)
        spin(3,2) = spin(3,2) + half_time_step * tmp * Torque(3,2)
        !write(*,*) "plan spin2",spin(1,2),spin(2,2),spin(3,2)
        !write(*,*) "dSpindt_p",tmp * Torque(1,2), tmp * Torque(2,2), tmp * Torque(3,2)

        time = time + half_time_step

        ! Calculate accelerations, gravity and tides
        call gravity_calculate_acceleration(nbod,m,x,v,a_grav)
        !write(*,*) "a_grav",a_grav

        call mfo_user(time,nbod,nbig,time_step,Rsth,Rsth5,Rsth10,rg2s,k2s,k2fs &
            ,sigmast,Rp,Rp5,Rp10,k2p,k2fp,rg2p,sigmap,m,x,v,spin,a_tides,Torque)
        !Torque(1,1:nbod) = 0.0d0
        !Torque(2,1:nbod) = 0.0d0
        !Torque(3,1:nbod) = 0.0d0
        !write(*,*) "a_tides",a_tides
        !write(*,*) "Torque ",Torque
        !write(*,*) 'end first part: K+D'

        ! Sum of accelerations
        a(1,1) = a_grav(1,1) + a_tides(1,1)
        a(2,1) = a_grav(2,1) + a_tides(2,1)
        a(3,1) = a_grav(3,1) + a_tides(3,1)

        a(1,2) = a_grav(1,2) + a_tides(1,2)
        a(2,2) = a_grav(2,2) + a_tides(2,2)
        a(3,2) = a_grav(3,2) + a_tides(3,2)
        !write(*,*) "acc tot",a(1,1), a(2,1),a(3,1),a(1,2),a(2,2),a(3,2)

        ! Velocity of star
        !write(*,*) "star vel1",v(1,1),v(2,1),v(3,1)
        v(1,1) = v(1,1) + time_step * a(1,1)
        v(2,1) = v(2,1) + time_step * a(2,1)
        v(3,1) = v(3,1) + time_step * a(3,1)
        !write(*,*) "star vel2",v(1,1),v(2,1),v(3,1)

        ! Velocity of planet
        !write(*,*) "plan vel1",v(1,2),v(2,2),v(3,2)
        v(1,2) = v(1,2) + time_step * a(1,2)
        v(2,2) = v(2,2) + time_step * a(2,2)
        v(3,2) = v(3,2) + time_step * a(3,2)
        !write(*,*) "plan vel2",v(1,2),v(2,2),v(3,2)

        ! Position of star
        !write(*,*) "star pos3",x(1,1),x(2,1),x(3,1)
        x(1,1) = x(1,1) + half_time_step * v(1,1)
        x(2,1) = x(2,1) + half_time_step * v(2,1)
        x(3,1) = x(3,1) + half_time_step * v(3,1)
        !write(*,*) "star pos4",x(1,1),x(2,1),x(3,1)

        ! Position of planet
        !write(*,*) "plan pos3",x(1,2),x(2,2),x(3,2)
        x(1,2) = x(1,2) + half_time_step * v(1,2)
        x(2,2) = x(2,2) + half_time_step * v(2,2)
        x(3,2) = x(3,2) + half_time_step * v(3,2)
        !write(*,*) "plan pos4",x(1,2),x(2,2),x(3,2)
        !write(*,*) "plan dis4",sqrt(x(1,2)*x(1,2)+x(2,2)*x(2,2)+x(3,2)*x(3,2))

        ! Spin of star
        tmp = -1.d0*K2/(m(1)*rg2s*Rsth*Rsth)
        spin(1,1) = spin(1,1) + half_time_step * tmp * Torque(1,1)
        spin(2,1) = spin(2,1) + half_time_step * tmp * Torque(2,1)
        spin(3,1) = spin(3,1) + half_time_step * tmp * Torque(3,1)

        ! Spin of planet
        tmp = -1.d0*K2/(m(2)*rg2p(1)*Rp(2)*Rp(2))
        spin(1,2) = spin(1,2) + half_time_step * tmp * Torque(1,2)
        spin(2,2) = spin(2,2) + half_time_step * tmp * Torque(2,2)
        spin(3,2) = spin(3,2) + half_time_step * tmp * Torque(3,2)
        !write(*,*) "spin star",spin(1,1),spin(2,1),spin(3,1)
        !write(*,*) "spin plan",spin(1,2),spin(2,2),spin(3,2)

        if (time.ge.time_write) then
            gm_tmp = m(1)+m(2)
            call mco_x2el (gm_tmp,x(1,2),x(2,2),x(3,2),v(1,2),v(2,2),v(3,2),q_tmp,e_tmp,i_tmp,p_tmp,n_tmp,l_tmp)
            open(13, file="aeipnl.out", access="append")
            write(13,'(7("  ", es20.10e3))') time/365.25d0,q_tmp/(1.d0-e_tmp),e_tmp,i_tmp,p_tmp,n_tmp,l_tmp
            close(13)
            !if (time/365.25d0.ge.10.d0*output_tmp) then
                !output_tmp = output_tmp * 10.d0
                !time_write = time_write + output_tmp*365.25d0
            !else 
                time_write = time_write + output_tmp*365.25d0
            !endif

        endif
        time = time + half_time_step
        !write(*,*) ""
    enddo

    !time = time + dt
end program mercury

