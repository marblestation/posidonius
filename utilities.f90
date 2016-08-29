!******************************************************************************
! MODULE: utilities
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Modules that gather various functions about string manipulation
!! and things that are perfectly separated from mercury particuliar
!! behaviour
!
!******************************************************************************

module utilities

  use types_numeriques

  implicit none
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> John E. Chambers
!
!> @date 24 May 1997
!
! DESCRIPTION: 
!> @brief Sorts an array X, of size N, using Shell's method. Also returns an array
!! INDEX that gives the original index of every item in the sorted array X.
!
!> @note The maximum array size is 29523.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine mxx_sort (n,x,index)
  

  implicit none

  
  ! Input/Output
  integer, intent(in) :: n
  integer, intent(out) :: index(n)
  real(double_precision), intent(inout) :: x(n)
  
  ! Local
  integer :: i,j,k,l,m,inc,iy
  real(double_precision) :: y
  integer, dimension(9), parameter :: incarr = (/1,4,13,40,121,364,1093,3280,9841/)
  
  !------------------------------------------------------------------------------
  
  do i = 1, n
     index(i) = i
  end do
  
  m = 0
10 m = m + 1
  if (incarr(m).lt.n) goto 10
  m = m - 1
  
  do i = m, 1, -1
     inc = incarr(i)
     do j = 1, inc
        do k = inc, n - j, inc
           y = x(j+k)
           iy = index(j+k)
           do l = j + k - inc, j, -inc
              if (x(l).le.y) goto 20
              x(l+inc) = x(l)
              index(l+inc) = index(l)
           end do
20         x(l+inc) = y
           index(l+inc) = iy
        end do
     end do
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mxx_sort

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> John E. Chambers
!
!> @date 30 September 2000
!
! DESCRIPTION: 
!> @brief Given initial and final coordinates and velocities, the routine returns
!! the X and Y coordinates of a box bounding the motion in between the
!! end points.
!!\n\n
!! If the X or Y velocity changes sign, the routine performs a quadratic
!! interpolation to estimate the corresponding extreme value of X or Y.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine mce_box (nbod,h,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
  
  use physical_constant
  use mercury_constant

  implicit none

  
  ! Input/Output
  integer, intent(in) :: nbod !< [in] current number of bodies (1: star; 2-nbig: big bodies; nbig+1-nbod: small bodies)
  real(double_precision), intent(in) :: h !< [in] current integration timestep (days)
  real(double_precision), intent(in) :: x0(3,nbod)
  real(double_precision), intent(in) :: x1(3,nbod)
  real(double_precision), intent(in) :: v0(3,nbod)
  real(double_precision), intent(in) :: v1(3,nbod)
  real(double_precision), intent(out) :: xmin(nbod)
  real(double_precision), intent(out) :: xmax(nbod)
  real(double_precision), intent(out) :: ymin(nbod)
  real(double_precision), intent(out) :: ymax(nbod)
  
  ! Local
  integer :: j
  real(double_precision) :: temp
  
  !------------------------------------------------------------------------------
  
  do j = 2, nbod
     xmin(j) = min (x0(1,j), x1(1,j))
     xmax(j) = max (x0(1,j), x1(1,j))
     ymin(j) = min (x0(2,j), x1(2,j))
     ymax(j) = max (x0(2,j), x1(2,j))
     
     ! If velocity changes sign, do an interpolation
     if ((v0(1,j).lt.0.and.v1(1,j).gt.0).or.(v0(1,j).gt.0.and.v1(1,j).lt.0)) then
        temp = (v0(1,j)*x1(1,j) - v1(1,j)*x0(1,j)       - .5d0*h*v0(1,j)*v1(1,j)) / (v0(1,j) - v1(1,j))
        xmin(j) = min (xmin(j),temp)
        xmax(j) = max (xmax(j),temp)
     end if
     
     if ((v0(2,j).lt.0.and.v1(2,j).gt.0).or.(v0(2,j).gt.0.and.v1(2,j).lt.0)) then
        temp = (v0(2,j)*x1(2,j) - v1(2,j)*x0(2,j)       - .5d0*h*v0(2,j)*v1(2,j)) / (v0(2,j) - v1(2,j))
        ymin(j) = min (ymin(j),temp)
        ymax(j) = max (ymax(j),temp)
     end if
  end do
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mce_box

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> John E. Chambers
!
!> @date 1 December 1998
!
! DESCRIPTION: 
!> @brief Calculates minimum value of a quantity D, within an interval H, given initial
!! and final values D0, D1, and their derivatives D0T, D1T, using third-order
!! (i.e. cubic) interpolation.
!!\n\n
!! Also calculates the value of the independent variable T at which D is a
!! minimum, with respect to the epoch of D1.
!
!> @note The routine assumes that only one minimum is present in the interval H.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine mce_min (d0,d1,d0t,d1t,h,d2min,tmin)

  implicit none
  
  ! Input/Output
  real(double_precision), intent(in) :: d0
  real(double_precision), intent(in) :: d1
  real(double_precision), intent(in) :: d0t
  real(double_precision), intent(in) :: d1t
  real(double_precision), intent(in) :: h !< [in] current integration timestep (days)
  real(double_precision), intent(out) :: d2min
  real(double_precision), intent(out) :: tmin
  
  ! Local
  real(double_precision) :: a,b,c,temp,tau
  
  !------------------------------------------------------------------------------
  
  if (d0t*h.gt.0.or.d1t*h.lt.0) then
     if (d0.le.d1) then
        d2min = d0
        tmin = -h
     else
        d2min = d1
        tmin = 0.d0
     end if
  else
     temp = 6.d0*(d0 - d1)
     a = temp + 3.d0*h*(d0t + d1t)
     b = temp + 2.d0*h*(d0t + 2.d0*d1t)
     c = h * d1t
     
     temp =-.5d0*(b + sign (sqrt(max(b*b - 4.d0*a*c,0.d0)), b) )
     if (temp.eq.0) then
        tau = 0.d0
     else
        tau = c / temp
     end if
     
     ! Make sure TAU falls in the interval -1 < TAU < 0
     tau = min(tau, 0.d0)
     tau = max(tau, -1.d0)
     
     ! Calculate TMIN and D2MIN
     tmin = tau * h
     temp = 1.d0 + tau
     d2min = tau*tau*((3.d0+2.d0*tau)*d0 + temp*h*d0t)    + temp*temp*((1.d0-2.d0*tau)*d1 + tau*h*d1t)
     
     ! Make sure D2MIN is not negative
     d2min = max(d2min, 0.d0)
  end if
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mce_min

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> John E. Chambers
!
!> @date 7 July 1999
!
! DESCRIPTION: 
!> @brief Converts from Julian day number to Julian/Gregorian Calendar dates, assuming
!! the dates are those used by the English calendar.
!!\n\n
!! Algorithm taken from `Practical Astronomy with your calculator' (1988)
!! by Peter Duffett-Smith, 3rd edition, C.U.P.
!!\n\n
!! Algorithm for negative Julian day numbers (Julian calendar assumed) by
!! J. E. Chambers.
!
!> @note The output date is with respect to the Julian Calendar on or before
!! 4th October 1582, and with respect to the Gregorian Calendar on or 
!! after 15th October 1582.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine mio_jd2y (jd0,year,month,day)
  

  implicit none

  
  ! Input/Output
  integer, intent(out) :: year
  integer, intent(out) :: month
  real(double_precision), intent(in) :: jd0
  real(double_precision), intent(out) :: day
  
  ! Local
  integer :: i,a,b,c,d,e,g
  real(double_precision) :: jd,f,temp,x,y,z
  
  !------------------------------------------------------------------------------
  
  if (jd0.le.0) goto 50
  
  jd = jd0 + 0.5d0
  i = sign( dint(dabs(jd)), jd )
  f = jd - 1.d0*i
  
  ! If on or after 15th October 1582
  if (i.gt.2299160) then
     temp = (1.d0*i - 1867216.25d0) / 36524.25d0
     a = sign( dint(dabs(temp)), temp )
     temp = .25d0 * a
     b = i + 1 + a - sign( dint(dabs(temp)), temp )
  else
     b = i
  end if
  
  c = b + 1524
  temp = (1.d0*c - 122.1d0) / 365.25d0
  d = sign( dint(dabs(temp)), temp )
  temp = 365.25d0 * d
  e = sign( dint(dabs(temp)), temp )
  temp = (c-e) / 30.6001d0
  g = sign( dint(dabs(temp)), temp )
  
  temp = 30.6001d0 * g
  day = 1.d0*(c-e) + f - 1.d0*sign( dint(dabs(temp)), temp )
  
  if (g.le.13) month = g - 1
  if (g.gt.13) month = g - 13
  
  if (month.gt.2) year = d - 4716
  if (month.le.2) year = d - 4715
  
  if (day.gt.32) then
     day = day - 32
     month = month + 1
  end if
  
  if (month.gt.12) then
     month = month - 12
     year = year + 1
  end if
  return
  
50 continue
  
  ! Algorithm for negative Julian day numbers (Duffett-Smith doesn't work)
  x = jd0 - 2232101.5
  f = x - dint(x)
  if (f.lt.0) f = f + 1.d0
  y = dint(mod(x,1461.d0) + 1461.d0)
  z = dint(mod(y,365.25d0))
  month = int((z + 0.5d0) / 30.61d0)
  day = dint(z + 1.5d0 - 30.61d0*dble(month)) + f
  month = mod(month + 2, 12) + 1
  
  year = 1399 + int (x / 365.25d0)
  if (x.lt.0) year = year - 1
  if (month.lt.3) year = year + 1
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mio_jd2y

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> John E. Chambers
!
!> @date 14 November 1999
!
! DESCRIPTION: 
!> @brief Given a character string STRING, of length LEN bytes, the routine finds 
!! the beginnings and ends of NSUB substrings present in the original, and 
!! delimited by spaces. The positions of the extremes of each substring are 
!! returned in the array DELIMIT.\n
!! Substrings are those which are separated by spaces or the = symbol.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine mio_spl (length,string,nsub,delimit)

  implicit none
  
  ! Input/Output
  integer, intent(in) :: length
  integer, intent(out) :: nsub!,delimit(2,100)
  integer, dimension(:,:), intent(out) :: delimit

  character(len=1), dimension(length), intent(in) :: string
  
  ! Local
  integer :: j,k
  character(len=1) :: c
  
  ! If the first dimension of the array 'delimit' is not 2, return an error
  if (size(delimit,1) /= 2) then
    write(*,*) "mio_spl: The first dimension of 'delimit' must be of size 2"
    stop
  end if
  
  !------------------------------------------------------------------------------
  
  nsub = 0
  j = 0
  c = ' '
  delimit(1,1) = -1
  
  ! Find the start of string
10 j = j + 1
  if (j.gt.length) goto 99
  c = string(j)
  if (c.eq.' '.or.c.eq.'=') goto 10
  
  ! Find the end of string
  k = j
20 k = k + 1
  if (k.gt.length) goto 30
  c = string(k)
  if (c.ne.' '.and.c.ne.'=') goto 20
  
  ! Store details for this string
30 nsub = nsub + 1
  delimit(1,nsub) = j
  delimit(2,nsub) = k - 1
  
  if (k.lt.length) then
     j = k
     goto 10
  end if
  
99 continue
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mio_spl

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> John E. Chambers
!
!> @date 2 March 1999
!
! DESCRIPTION: 
!> @brief Calculates inverse hyperbolic cosine of an angle X (in radians).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function arcosh (x)
  

  implicit none

  
  ! Input/Output
  real(double_precision), intent(in) :: x
  real(double_precision) :: arcosh
  
  !------------------------------------------------------------------------------
  
  if (x.ge.1.d0) then
     arcosh = log (x + sqrt(x*x - 1.d0))
  else
     arcosh = 0.d0
  end if
  
  !------------------------------------------------------------------------------
  
  return
end function arcosh

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> John E. Chambers
!
!> @date 17 April 1997
!
! DESCRIPTION: 
!> @brief Calculates sin and cos of an angle X (in radians).
!
!> @warning x is modified by the routine. Without this part, outputs of mercury are different, don't ask me why.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine mco_sine (x,sx,cx)
  use physical_constant, only : PI, TWOPI

  implicit none

  
  ! Input/Output
  real(double_precision), intent(inout) :: x
  real(double_precision), intent(out) :: sx
  real(double_precision), intent(out) :: cx
    
  !------------------------------------------------------------------------------
  
  ! TODO why results of this routine are different from simple calls of intrinsec cos() and sin()
  if (x.gt.0) then
     x = mod(x,TWOPI)
  else
     x = mod(x,TWOPI) + TWOPI
  end if
  
  cx = cos(x)
  
  if (x.gt.PI) then
     sx = -sqrt(1.d0 - cx*cx)
  else
     sx =  sqrt(1.d0 - cx*cx)
  end if
  
  !------------------------------------------------------------------------------
  
  return
end subroutine mco_sine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> C. Cossou
!
!> @date 2011
!
! DESCRIPTION: 
!> @brief return the number of bodies and the number of big bodies
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine getNumberOfBodies(nb_big_bodies, nb_bodies)

  implicit none

  integer, intent(out) :: nb_bodies !< [out] total number of bodies, including the central object
  integer, intent(out) :: nb_big_bodies !< [out] number of big bodies
  
  ! Local
  integer :: j,k,lim(2,10),nsub, error
  logical test
  character(len=80) :: infile(3),filename,c80
  character(len=150) :: string
  real(double_precision), dimension(9) :: dummy ! variable to read the value and have the right position in the file without storing the values
  !------------------------------------------------------------------------------
  
  do j = 1, 80
     filename(j:j) = ' '
  end do
  do j = 1, 3
     infile(j)   = filename
  end do  
  
  ! Read in filenames and check for duplicate filenames
  inquire (file='files.in', exist=test)
  if (.not.test) write(*,*) 'Error: the file "files.in" does not exist'
  open (15, file='files.in', status='old')
  
  ! Input files
  do j = 1, 3
     read (15,'(a150)') string
     call mio_spl (150,string,nsub,lim)
     infile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
     do k = 1, j - 1
        if (infile(j).eq.infile(k)) write(*,*) 'Error: In "files.in", Some files are identical.'
     end do
  end do
  close (15)
  
  !--------------------------------------------------------

  !  READ  IN  DATA  FOR  BIG  AND  SMALL  BODIES
  
  nb_bodies = 1
  do j = 1, 2
     if (j.eq.2) nb_big_bodies = nb_bodies
     
     ! Check if the file containing data for Big bodies exists, and open it
     filename = infile(j)
     inquire (file=infile(j), exist=test)
     if (.not.test) write (*,*) filename, 'does not exist'
     open (11, file=infile(j), status='old', iostat=error)
     if (error /= 0) then
        write (*,'(/,2a)') " ERROR: Programme terminated. Unable to open ",trim(filename)
        stop
     end if
     
     ! Read data style
     do
      read (11,'(a150)') string
      if (string(1:1).ne.')') exit
    end do
     
     ! Read epoch of Big bodies
     if (j.eq.1) then
      do
        read (11,'(a150)') string
        if (string(1:1).ne.')') exit
      end do
     end if
     
     ! Read information for each object
     do
      read (11,'(a)',iostat=error) string
      if (error /= 0) exit
      
      if (string(1:1).eq.')') cycle
    
     
     nb_bodies = nb_bodies + 1
     
     
    ! we skip the line(s) that contains informations of the current planet
     do
       read (11,'(a150)') string
       if (string(1:1).ne.')') exit
     end do 
     backspace(11)
     read (11,*) dummy
     
     
     end do
    close (11)
  end do

end subroutine getNumberOfBodies

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> C. Cossou
!
!> @date 2011
!
! DESCRIPTION: 
!> @brief change cartesian coordinates into polar coordinates
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_polar_coordinates(x, y, z, radius, theta)
  use physical_constant
  
  implicit none
  
  ! Inputs
  real(double_precision), intent(in) :: x
  real(double_precision), intent(in) :: y
  real(double_precision), intent(in) :: z
  
  ! Outputs
  real(double_precision), intent(out) :: radius
  real(double_precision), intent(out) :: theta
  
  radius = sqrt(x * x + y * y + z * z)
  
  ! There will be problem if radius equal 0 because in this case, the angle can take any real value possible
  
  ! To obtain theta inside [0; 2 pi[ (source: wikipedia) -> The function atan2 do exactly this (inside [-pi ; pi[, but with a buit-in function
!~   theta = atan2(y, x)
  theta = acos(x / radius)
  
!~   write(*,*) x, y, z, radius, theta
!~   stop
  
  ! We want the output between 0 and 2*PI
!~   if (theta.lt.0.) then
!~     theta = theta + TWOPI
!~   endif
  if (y.lt.0.) then
    theta = TWOPI - theta
  endif
  
  
end subroutine get_polar_coordinates

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> C. Cossou
!
!> @date 2011
!
! DESCRIPTION: 
!> @brief function that calculate the mean value of a one dimensional array
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function get_mean(vector)

implicit none

! Input
real(double_precision), intent(in), dimension(:) :: vector

! Output
real(double_precision) :: get_mean
!--------------------------------------------------------

get_mean = sum(vector) / size(vector)

end function get_mean

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> C. Cossou
!
!> @date 2011
!
! DESCRIPTION: 
!> @brief function that calculate the standard deviation of a set of value given as a one dimensional array
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function get_stdev(vector)

implicit none

! Input
real(double_precision), intent(in), dimension(:) :: vector

! Output
real(double_precision) :: get_stdev
!--------------------------------------------------------

get_stdev = sqrt(sum((vector - get_mean(vector))**2) / size(vector))

end function get_stdev

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> C. Cossou
!
!> @date 2011
!
! DESCRIPTION: 
!> @brief Return the vectoriel product of two vectors
!! vect_product = x /\ y
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function vect_product(x, y)

  implicit none

  ! The function - output :
  real(double_precision), dimension(3) :: vect_product

  ! Inputs :
  real(double_precision), dimension(3), intent(in) :: x
  real(double_precision), dimension(3), intent(in) :: y
!--------------------------------------------------------

  vect_product(1) = x(2) * y(3) - x(3) * y(2)
  vect_product(2) = x(3) * y(1) - x(1) * y(3)
  vect_product(3) = x(1) * y(2) - x(2) * y(1)



end function vect_product

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> C. Cossou
!
!> @date 2011
!
! DESCRIPTION: 
!> @brief the subroutine get in argument a one dimension array and 
!! return an histogram of the number of bins specified in argument. By default the min and max for bins is the min and max of the 
!! data set. With optional argument it should be possible easily to accept in argument the min and max for bins (with default 
!! values) but it is not implemented yet since I have no need for this so far.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine get_histogram(datas, bin_x_values, bin_y_values) 

implicit none

! Input
real(double_precision), dimension(:), intent(in) :: datas !< [in] a one dimension array that contains data to be used for the histogram

! Output
real(double_precision), dimension(:), intent(out) :: bin_x_values !< [out] The x values for the histogram. The array must be of size 'nb_bins'
real(double_precision), dimension(:), intent(out) :: bin_y_values !< [out] The y values for the histogram. The array must be of size 'nb_bins'

! Locals
real(double_precision) :: delta_bin, max_value, min_value
integer :: index_bin, nb_points, nb_bins

integer :: i ! For loops

!--------------------------------------------------------

  ! We get the total number of points 
  nb_points = size(datas)
  
  ! We get the number of bins from the size of the arrays associated
  nb_bins = size(bin_y_values)
  
  if (.not.(size(bin_x_values).eq.nb_bins)) then
    stop 'Error in get_histogram : "bin_x_values" and "bin_y_values" do not have the same size'
  end if

  ! We initialize the values of the counting array
  bin_y_values(1:nb_bins) = 0

  ! From the list of values, we get the values for the histogram
  max_value = maxval(datas(1:nb_points))
  min_value = minval(datas(1:nb_points))
  
  delta_bin = (max_value - min_value) / float(nb_bins)
  
  if (delta_bin.eq.0.d0) then
    write(*,'(a, es6.0e2,a, es7.0e2,a, es7.0e2,a)') 'subroutine get_histogram: For ', float(nb_points), &
             ' values, the turbulent torque is between [', min_value, ' ; ', max_value, ']'
    write(*,'(a, i5, a, es7.0e2)') 'Thus, for ', nb_bins, ' bins in the histogram, the single width of a bin is : ', delta_bin
    stop 'Error in get_histogram : the min and max of the dataset are equal, resulting in a delta_bin equal to 0.'
  end if
  
  do i=1,nb_bins
    bin_x_values(i) = min_value + (i - 0.5d0) * delta_bin
  end do
    
  do i=1, nb_points
    ! if the value is exactly equal to max_value, we force the index to be nb_bins
    index_bin = min(floor((datas(i) - min_value) / delta_bin)+1, nb_bins)
    
    ! With floor, we get the immediate integer below the value. Since the index start at 1, we add 1 to the value, because the 
    ! calculation will get from 0 to the number of bins. Thus, for the max value, we will get nb_bins +1, which is not possible. 
    ! As a consequence, we take the lower value between the index and nb_bins, to ensure that for the max value, we get an index 
    ! of nb_bins.
    bin_y_values(index_bin) = bin_y_values(index_bin) + 1
  end do

  ! We normalize the histogram
  bin_y_values(1:nb_bins) = bin_y_values(1:nb_bins) / float(nb_points)

end subroutine get_histogram
  
end module utilities
