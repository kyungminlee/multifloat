using MultiFloats

constants = [

    (+0x1.0000000000000p+0000, +0x1.314BACF0323FFp-0113),
    (+0x1.62E42FEFA39EFp-0001, +0x1.ABC9E3B39803Fp-0056),
    (+0x1.EBFBDFF82C58Fp-0003, -0x1.5E43A53E454F1p-0057),
    (+0x1.C6B08D704A0C0p-0005, -0x1.D3316275139AEp-0059),
    (+0x1.3B2AB6FBA4E77p-0007, +0x1.4E65DFEF67D34p-0062),
    (+0x1.5D87FE78A6731p-0010, +0x1.0717F88815ADFp-0066),
    (+0x1.430912F86C787p-0013, +0x1.BC7CDBCDC0339p-0067),
    (+0x1.FFCBFC588B0C7p-0017, -0x1.E645E286FE571p-0071),
    (+0x1.62C0223A5C863p-0020, -0x1.99EF542AA8E1Ep-0074),
    (+0x1.B5253D395E80Fp-0024, 0.0),
    (+0x1.E4CF5152FBB30p-0028, 0.0),
    (+0x1.E8CAC72F6E9E5p-0032, 0.0),
    (+0x1.C3C1919538484p-0036, 0.0),
    (+0x1.816519F74C4AFp-0040, 0.0),
   ]


num2hex(x::Float64) = string(reinterpret(UInt64, x), base=16)
num2dec(x::Float64) = string(reinterpret(UInt64, x), base=10)
for i in 14:-1:1
# for c in constants
    println("! i=$i")
    c = constants[i]
    println("! ", num2hex(c[1]), '\t', num2hex(c[2]))
    println("! ", num2dec(c[1]), '\t', num2dec(c[2]))
    println("! ", c[1], '\t', c[2])
    println("res = add(mul(res, x), float64x2((/dble(z'$(num2hex(c[1]))'), dble(z'$(num2hex(c[2]))')/)))")
end

# LOG2E
let c = (+0x1.71547652B82FEp+0000, +0x1.777D0FFDA0D24p-0056)
    println("! ", num2hex(c[1]), '\t', num2hex(c[2]))
    println("! ", num2dec(c[1]), '\t', num2dec(c[2]))
    println("! ", c[1], '\t', c[2])
    println("res = add(mul(res, x), float64x2((/dble(z'$(num2hex(c[1]))'), dble(z'$(num2hex(c[2]))')/)))")
    println("exp = ", exp2(Float64x2(c))._limbs)
end
