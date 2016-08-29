module spline

  !*************************************************************
  !** Modules that contains user defined modules. 
  !** Only spline will be public.
  !**
  !** Version 1.0 - june 2011
  !*************************************************************
  use types_numeriques

  implicit none

contains

  subroutine spline_b_val (ndata,tdata,ydata,tval,yval)
    !
    implicit none
    !
    integer :: ndata
    !
    real(double_precision) :: bval
    integer :: left
    integer :: right
    real(double_precision), dimension(ndata) :: tdata,ydata
    real(double_precision) :: tval,u,yval
    !
    !  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
    !
    call rvec_bracket (ndata,tdata,tval,left,right)
    !
    !  Evaluate the 5 nonzero B spline basis functions in the interval,
    !  weighted by their corresponding data values.
    !
    u = (tval-tdata(left))/(tdata(right)-tdata(left))
    yval = 0.0E+00
    !
    !  B function associated with node LEFT - 1, (or "phantom node"),
    !  evaluated in its 4th interval.
    !
    bval = ( 1.0E+00 - 3.0E+00 * u + 3.0E+00 * u**2 - u**3 ) / 6.0E+00
    if (left-1.gt.0) then
       yval = yval + ydata(left-1) * bval
    else
       yval = yval + ( 2.0E+00 * ydata(1) - ydata(2) ) * bval
    end if
    !
    !  B function associated with node LEFT,
    !  evaluated in its third interval.
    !
    bval = ( 4.0E+00 - 6.0E+00 * u**2 + 3.0E+00 * u**3 ) / 6.0E+00
    yval = yval + ydata(left) * bval
    !
    !  B function associated with node RIGHT,
    !  evaluated in its second interval.
    !
    bval = ( 1.0E+00 + 3.0E+00 * u + 3.0E+00 * u**2 - 3.0E+00 * u**3 ) / 6.0E+00
    yval = yval + ydata(right) * bval
    !
    !  B function associated with node RIGHT+1, (or "phantom node"),
    !  evaluated in its first interval.
    !
    bval = u**3 / 6.0E+00
    if ( right+1.le.ndata ) then
       yval = yval + ydata(right+1) * bval
    else
       yval = yval + ( 2.0E+00 * ydata(ndata) - ydata(ndata-1) ) * bval
    end if

    return
  end subroutine spline_b_val

  subroutine rvec_bracket (n,x,xval,left,right)
    !
    implicit none
    !
    integer :: n
    !
    integer :: i,left,right
    real(double_precision), dimension(n) :: x
    real(double_precision) :: xval
    !
    do i = 2,n-1
       if (xval.lt.x(i)) then
          left = i-1
          right = i
          return
       end if
    end do
    left = n - 1
    right = n
    return
  end subroutine rvec_bracket

end module spline
