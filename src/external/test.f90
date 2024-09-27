module angmom_lib
  implicit none
contains
  function gmosh(x, y) result(z) 
    real(8), intent(in) :: x, y
    real(8) :: z
    PRINT *, 'The value of x is:', x
    PRINT *, 'The value of y is:', y
    z = x + y
  end function gmosh
end module angmom_lib

