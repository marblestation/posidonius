!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!      CODE_FORTRAN.F90    (ErikSoft   3 May 2002)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Mercury is a general-purpose N-body integration package for problems in
! celestial mechanics.


!------------------------------------------------------------------------------


program mercury
  
    use user_module
    use physical_constant

    implicit none
  
    integer :: j,nbod,nbig

    real(double_precision) :: time,gm_tmp,dt
    real(double_precision) :: q_tmp,e_tmp,i_tmp,p_tmp,n_tmp,l_tmp
    real(double_precision) :: x_tmp,y_tmp,z_tmp,u_tmp,v_tmp,w_tmp

    real(double_precision), dimension(:), allocatable :: m,sma,ecc,inc ! (Number of bodies)
    real(double_precision), dimension(:,:), allocatable :: xh,vh ! (3,Number of bodies)
    real(double_precision), dimension(:,:), allocatable :: x,v,a ! (3,Number of bodies)
  
    !------------------------------------------------------------------------------

    ! star + planets
    nbod = 3
    ! planets
    nbig = 2

    ! Initialization
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


    ! m is here actually gm 
    m(1) = K2*1.d0
    m(2:nbod) = (/K2*3.d-6, K2*6*3.d6/)

    ! sma in AU, eccentricity and inclination
    sma(2:nbod) = (/0.01, 0.05/)
    ecc(2:nbod) = (/0.01, 0.01/)
    inc(2:nbod) = (/1.0, 1.0/)

    ! Loop on planets
    do j = 2,nbod
        gm_tmp = m(1) + m(j)
        q_tmp = sma(j)*(1.-ecc(j))
        e_tmp = ecc(j)
        i_tmp = inc(j)
        p_tmp = 0.
        n_tmp = 0.
        l_tmp = 0.
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

    enddo

    time = 0.
    dt = 0.08 ! day (re-write it in user_module)

    call mfo_user(time,nbod,nbig,m,x,v,a)

    time = time + dt
end program mercury
