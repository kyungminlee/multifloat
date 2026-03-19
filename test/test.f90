program multifloat_test
  use multifloat
  implicit none
  type(float64x2) :: a, b, c
  double precision :: ah, al

  ! a = 1.0 + 1.0d-10
  a = 1.0d0
  a%limbs(2) = 1.0d-10
  
  ! c = a * a = (1 + 1d-10)^2 = 1 + 2d-10 + 1d-20
  c = a * a
  
  print *, "a: ", a%limbs
  print *, "c: ", c%limbs
  print *, "c target: 1.0000000002 1.0d-20 (approx)"
  
  ! The second limb of c should be around 1.0d-20 if it works correctly.
  ! Currently, it's likely 0 or something else because mul is just limb-wise.
  ! Actually, mul is: c%limbs = a%limbs * b%limbs
  ! So c%limbs(1) = 1.0 * 1.0 = 1.0
  ! c%limbs(2) = 1.0d-10 * 1.0d-10 = 1.0d-20
  ! Wait, that's also not quite right for (a_h + a_l)^2 = a_h^2 + 2*a_h*a_l + a_l^2
  ! Current mul: c%limbs(1) = a_h * b_h, c%limbs(2) = a_l * b_l
  ! For a=b, c%limbs(1) = a_h^2, c%limbs(2) = a_l^2. Missing 2*a_h*a_l.

  b = c / a
  print *, "b = c / a: ", b%limbs
  print *, "b target: 1.0 1.0d-10"
  
end program
