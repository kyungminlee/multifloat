program multifloat_test
  use, intrinsic::iso_fortran_env
  implicit none
  real(kind=16) :: var
  var = 1.0q3
  print *, var
end program
