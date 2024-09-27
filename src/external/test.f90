module angmom_lib
  implicit none
contains
  function gmosh(x, y) result(z)
    real(8), intent(in) :: x, y
    real(8) :: z
    z = x + y
  end function gmosh
end module angmom_lib

