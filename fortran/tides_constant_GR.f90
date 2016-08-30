module tides_constant_GR
  use types_numeriques
  use physical_constant

  implicit none
  !
  ! Author: Emeline Bolmont
  ! Date: 04/04/13
  !
  !----------------------------------------------------------------------------- 
  ! Output of spin every 'output' years
  real(double_precision), parameter :: output = 1

  !---------------------------  effects  --------------------------------------- 
  ! If you want effect of rotational induced flattening or not
  integer, parameter :: rot_flat = 0
  ! If you want General Relativity or not
  integer, parameter :: GenRel = 0
  ! If you want Tides or not
  integer, parameter :: tides = 1

  ! Number of planets for which you want these effects
  integer, parameter :: ntid = 1

  !-----------------------  Nature of host body  ------------------------------- 
  integer, parameter :: brown_dwarf = 0
  integer, parameter :: M_dwarf = 0
  integer, parameter :: Sun_like_star = 0  
  integer, parameter :: Jupiter_host = 0  
  ! For an utilization of the code with non-evolving host body
  ! if Rscst = 1, Rs = cst, rg2s = cst
  integer, parameter :: Rscst = 1

  !---------------------  Integration parameters  ------------------------------ 

  ! If crash
  integer, parameter :: crash = 0

  ! if crash = 0, then t_init = initial integration time (day)
  !                           = age of the evolving body at beginning of simu
  ! if crash = 1, then t_init = time of last line of PLANETi.aei
  ! if crash = 1 and evolving body, then t_init = time of last line of PLANETi.aei + t_init(crash=0)
  real(double_precision), parameter :: t_init = 0.0d6 * 365.25d0

  ! if crash = 1, then t_crash = time of last line of spini.out (day)
  real(double_precision), parameter :: t_crash = 0.0d0 * 365.25d0
  ! if crash = 1, then give the values of last line of spins.out to rot_crash
  ! and do the same for the planets, spinpi.out corresponds to rot_crashpi
  real(double_precision), parameter, dimension(3) :: rot_crash   = (/0.0d0,0.0d0,0.0d0/)
  real(double_precision), parameter, dimension(3) :: rot_crashp1 = (/0.0d0,0.0d0,0.0d0/)
  real(double_precision), parameter, dimension(3) :: rot_crashp2 = (/0.0d0,0.0d0,0.0d0/)
  real(double_precision), parameter, dimension(3) :: rot_crashp3 = (/0.0d0,0.0d0,0.0d0/)
  real(double_precision), parameter, dimension(3) :: rot_crashp4 = (/0.0d0,0.0d0,0.0d0/)
  real(double_precision), parameter, dimension(3) :: rot_crashp5 = (/0.0d0,0.0d0,0.0d0/)
  real(double_precision), parameter, dimension(3) :: rot_crashp6 = (/0.0d0,0.0d0,0.0d0/)

  !----------------------------------------------------------------------------- 
  !------------------------------  PLANETS  ------------------------------------
  !----------------------------------------------------------------------------- 

  !--------------------------  planets radius  ---------------------------------
  ! Value of the radius of Jupiter in Earth radius
  real(double_precision), parameter :: Rjup = 10.9d0 !rearth	

  ! planet_type gives the type of planet
  ! 0: Earth-like (mass-radius relationship, Fortney 2007) 
  ! 1: Terrestrial (no mass-radius relationship)
  ! 2: Gas giant
  ! 3: user chooses everything
  integer, parameter, dimension(ntid) :: planet_type = (/3/)
  ! If planet_type ne 0, then indicate radius in Rearth
  ! for ex: 1 or 0.954d0*Rjup
  real(double_precision), parameter :: planet_radius_factor = 1.d0
  real(double_precision), parameter, dimension(ntid) :: radius_p = (/planet_radius_factor/)

  !---------------------------  initial spin  ----------------------------------
  ! If pseudo_rot eq 0 : initial period as given by Pp0 (in hr)
  ! If pseudo_rot eq toto : initial period = toto*pseudo_synchronization period 
  real(double_precision), parameter, dimension(ntid) :: pseudo_rot = (/0.d0/)
  real(double_precision), parameter :: planet_rotation_period = 24.d0
  real(double_precision), parameter, dimension(ntid) :: Pp0 = (/planet_rotation_period/)

  ! Planets obliquities in rad
  real(double_precision), parameter :: planet_obliquity = 11.459156d0*PI/180.d0
  real(double_precision), parameter, dimension(ntid) :: oblp = (/planet_obliquity/)

  !---------------------  for planets with planet_type=3  ----------------------
  !--------------  love number, radius of gyration, dissipation  ---------------
  ! Radius of gyration
  real(double_precision), parameter :: planet_radius_of_gyration_2 = 3.308e-1!/ Earth type planet
  real(double_precision), parameter, dimension(ntid) :: rg2p_what = (/planet_radius_of_gyration_2/)
  ! Tidal love number
  real(double_precision), parameter :: planet_love_number = 0.305d0 !/ Earth
  real(double_precision), parameter, dimension(ntid) :: k2tp_what = (/planet_love_number/)
  ! Fluid love number
  real(double_precision), parameter, dimension(ntid) :: k2fp_what = (/planet_love_number/)
  ! k2delta_t (day)
  real(double_precision), parameter :: k2pdelta = 2.465278e-3!/ Terrestrial planets (no gas)
  real(double_precision), parameter, dimension(ntid) :: k2pdeltap_what = (/k2pdelta/)

  !---------------------------  dissipation  -----------------------------------
  ! factor to multiply the value of the dissipation of the planets
  real(double_precision), parameter :: planet_dissipation_factor_scale = 1.
  real(double_precision), parameter, dimension(ntid) :: dissplan = (/planet_dissipation_factor_scale/)


  !----------------------------------------------------------------------------- 
  !-----------------------------  HOST BODY  -----------------------------------
  !----------------------------------------------------------------------------- 

  ! Star dissipation, and caracteristics in CGS
  ! factor to multiply the value of the dissipation of the host body
  real(double_precision), parameter :: star_dissipation_factor_scale = 1.0d0
  real(double_precision), parameter :: dissstar = star_dissipation_factor_scale

  ! For R=cst
  ! Dissipation, either you give sigma or k2deltat
  ! when you choose one, put the other at a negative value
  ! dissipation factor sigma
  real(double_precision), parameter :: sigma_what = 2.006*3.845764d4 !-60+64
  ! k2delta_t (day), here value of Earth for satellite systems for example 
  real(double_precision), parameter :: k2sdeltats = -2.465277778d-3
  ! Radius of gyration
  real(double_precision), parameter :: star_radius_of_gyration_2 = 1.94e-1!/ Brown dwarf
  real(double_precision), parameter :: rg2_what = star_radius_of_gyration_2
  ! Potential Love number
  real(double_precision), parameter :: star_love_number = 0.307!/ M Dwarf
  real(double_precision), parameter :: k2st_what = star_love_number
  ! Radius of star in Rsun
  real(double_precision), parameter :: radius_factor = 0.845649342247916
  real(double_precision), parameter :: radius_star = radius_factor

  ! For R=cst, or dM or Suns
  ! Initial period of rotation in day
  real(double_precision), parameter :: star_rotation_period = 70.d0
  real(double_precision), parameter :: Period_st   = star_rotation_period
  ! Fluid Love number
  real(double_precision), parameter :: k2fst_what = star_love_number 
  


  !----------------------------------------------------------------------------- 
  !----------------------------------------------------------------------------- 
  !--------------------  No Need to change stuff from here  -------------------- 
  !----------------------------------------------------------------------------- 
  ! Radius of gyration and love number for dM 
  real(double_precision), parameter :: rg2_dM = 2.0d-1
  real(double_precision), parameter :: k2st_dM = 0.307d0 
  ! Radius of gyration and love number for Suns
  real(double_precision), parameter :: rg2_Sun = 5.9d-2
  real(double_precision), parameter :: k2st_Sun = 0.03d0 
  ! Radius of gyration, love number and k2delta_t for terrestrial planets
  real(double_precision), parameter :: rg2p_terr = 3.308d-1
  real(double_precision), parameter :: k2fp_terr = 0.9532d0
  real(double_precision), parameter :: k2p_terr = 0.299d0
  real(double_precision), parameter :: k2pdeltap_terr = 2.465278d-3
  ! Radius of gyration, love number and k2delta_t for gas giants
  real(double_precision), parameter :: rg2p_gg = 2.54d-1
  real(double_precision), parameter :: k2p_gg = 0.38d0
  real(double_precision), parameter :: k2pdeltap_gg = 8.101852d-9

  ! Sigma for BD, dM, Suns, Jupiter
  ! BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
  real(double_precision), parameter :: sigma_BD = 2.006*3.845764d4 !-60+64
  real(double_precision), parameter :: sigma_dM = 2.006*3.845764d4 !-60+64
  ! Sun-like-star: sigmast = 4.992d-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
  real(double_precision), parameter :: sigma_Sun = 4.992*3.845764d-2 !-66+64
  ! If planet not terrestrial, dissipation factor Hot Gas Giant
  real(double_precision), parameter :: sigma_gg = 2.006*3.845764d4
  ! k2delta_t for Jupiter: 2-3d-2 s, here in day (Leconte)
  real(double_precision), parameter :: k2delta_jup = 2.893519d-7

  ! Radius of the Sun in AU 
  real(double_precision), parameter :: rsun = 4.649130365292d-3
  ! Radius of Earth in AU
  real(double_precision), parameter :: rearth = 4.25875047552248d-5
  ! Ratio between mass of the Sun and mass of Earth
  real(double_precision), parameter :: m2earth = 3.2890045607235d5
  ! meter in AU
  real(double_precision), parameter :: minau = 6.684587153547d-12
  ! Speed of light
  real(double_precision), parameter :: C2 = 1.731444830225d2

  ! Runge Kutta parameters
  real(double_precision), parameter, dimension(5) :: aa = (/0.2d0,0.3d0,0.6d0,1.d0,0.875d0/)
  real(double_precision), parameter :: bb2 = 0.2d0
  real(double_precision), parameter, dimension(2) :: bb3 = (/7.5d-2,2.25d-1/)
  real(double_precision), parameter, dimension(3) :: bb4 = (/0.3d0,-0.9d0,1.2d0/)
  real(double_precision), parameter, dimension(4) :: bb5 = (/-1.1d0/5.4d0,2.5d0,-7.d0/2.7d0 &
                                                             ,3.5d0/2.7d0/)
  real(double_precision), parameter, dimension(5) :: bb6 = (/1.631d0/5.5296d1,1.75d0/5.12d0 &
                                                              ,5.75d0/1.3824d2 &
                                                              ,4.4275d0/1.10592d1,2.53d0/4.096d1/)
  real(double_precision), parameter, dimension(6) :: cc = (/3.7d0/3.78d1,0.d0,2.5d0/6.21d0 &
                                                              ,1.25d0/5.94d0 &
                                                              ,0.d0,5.12d0/1.771d1/)                                                 
contains 

    !---------------------------------------------------------------------------
    ! function that return the initial timestep of the simulation
    function get_initial_timestep()

        use utilities, only : mio_spl
        implicit none
        !Input
        !None actually
        ! Outpout
        real(double_precision), dimension(3) :: get_initial_timestep 
        !Locals
        integer :: j, lineno, nsub, lim(2,10), error
        real(double_precision) :: h0,tstop,mass_star
        character(len=80) :: c80
        character(len=150) :: string
        !-----------------------------------------------------------------------
        ! the timestep of the simulation as writed in param.in
        open(13, file='param.in', status='old', iostat=error)
        if (error /= 0) then
            write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim('param.in')
            stop
        end if
        ! Read integration parameters
        lineno = 0
        do j = 1, 26
            ! We want the next non commented line
            do
                lineno = lineno + 1
                read (13,'(a150)') string
                if (string(1:1).ne.')') exit
            enddo

            call mio_spl (150,string,nsub,lim)
            c80(1:3) = '   '
            c80 = string(lim(1,nsub):lim(2,nsub))

            if (j.eq.3) read (c80,*) tstop
            if (j.eq.5) read (c80,*) h0
            if (j.eq.18) read (c80,*) mass_star
        enddo
        get_initial_timestep = (/tstop,abs(h0),mass_star/)
        close (13)
        !-----------------------------------------------------------------------
        return
    end function get_initial_timestep

!-------------------------------------------------------------------------------
! subroutine that write the parameters of the user_module 
! into the file 'tidesGR.out'
subroutine write_simus_properties()

    use git_infos
    implicit none
    real(double_precision), dimension(3) :: timestep
    real(double_precision) :: distance_accuracy
    integer :: j
    real(double_precision), parameter :: TWOTHIRD = 2.d0 / 3.d0
    !---------------------------------------------------------------------------

    timestep = get_initial_timestep()
    ! below this limit, with this timestep, an orbit will only contain 20 timestep or less, whiis not accurate.
    !distance_accuracy = (10. * timestep(2) / 365.25)**TWOTHIRD 
    distance_accuracy = (30. * timestep(2) *sqrt(timestep(3)*K2) / TWOPI)**TWOTHIRD 

    open(10, file='tidesGR.out')

    write(10,*) ''
    write(10,'(a)') '------------------------------------'
    write(10,'(a)') '|         Timestep stuff           |'
    write(10,'(a)') '------------------------------------'
    write(10,'(a,f12.5,a)') 'timestep = ',timestep(2), ' days'
    write(10,'(a,f12.5,a)') '  with this timestep, the simulation will not be accurate below', distance_accuracy,' AU'
    write(10,*) ''
    write(10,*) ''
    write(10,'(a)') '------------------------------------'
    write(10,'(a)') '|       Mercury Properties         |'
    write(10,'(a)') '------------------------------------'
    write(10,'(a,a)') 'branch = ', branch
    write(10,'(a,a)') 'commit = ', commit
    write(10,'(a,a)') 'tags = ', tags
    write(10,'(a)') modifs
    write(10,'(a,f12.5,a,f12.5,a)') 'With h=', timestep(2), ' days, simulation is accurate for r > ', distance_accuracy,' AU'
    write(10,*) ''
    write(10,*) ''
    write(10,'(a)') '------------------------------------'
    write(10,'(a)') '|       Special Effects            |'
    write(10,'(a)') '------------------------------------'
    write(10,*) ''
    if (tides.eq.1) write(10,'(a)') 'Tides are on'
    if (GenRel.eq.1) write(10,'(a)') 'General Relativity effects taken into account'
    if (rot_flat.eq.1) write(10,'(a)') 'Effect of rotation-induced flattening taken into account'
    write(10,*) ''
    if (brown_dwarf.eq.1) write(10,'(a)') 'The central body is an evolving Brown-dwarf'
    if (M_dwarf.eq.1) write(10,'(a)') 'The central body is an evolving M-dwarf'
    if (Sun_like_star.eq.1) write(10,'(a)') 'The central body is an evolving Sun-like star'
    if (Rscst.eq.1) write(10,'(a)') 'The central body is an non-evolving object'
    write(10,*) ''

    if (tides.eq.1) then
       write(10,'(a,i1)') 'Number of planets tidally evolving =',ntid
       do j = 2, ntid+1
          write(10,'(a,i1)') 'PLANET',j
          if ((planet_type(j-1).eq.0).or.(planet_type(j-1).eq.1)) then
             write(10,'(a,f12.5,a,f12.5)') 'k2fp =',k2fp_terr,'k2p =',k2p_terr,', rg2p =',rg2p_terr
             write(10,'(a,f12.5,a,f12.5)') 'k2pdeltap =',k2pdeltap_terr,' day, dissplan =',dissplan(j-1)
          endif
          if (planet_type(j-1).eq.2) then
             write(10,'(a,f12.5,a,f12.5)') 'k2p =',k2p_gg,', rg2p =',rg2p_gg
             write(10,'(a,f12.5,a,f12.5)') 'k2pdeltap =',k2pdeltap_gg,' day, dissplan =',dissplan(j-1)
          endif
          if (planet_type(j-1).eq.3) then
             write(10,'(a,f12.5,a,f12.5)') 'k2p =',k2tp_what(j-1),'k2p =',k2fp_what(j-1),', rg2p =',rg2p_what(j-1)
             write(10,'(a,es19.9e3,a,f12.5)') 'k2pdeltap =',k2pdeltap_what(j-1),' day, dissplan =',dissplan(j-1)
          endif
       enddo
    endif
    write(10,*) ''
    close(10)
    !---------------------------------------------------------------------------
end subroutine write_simus_properties

end module tides_constant_GR
