using Test, Nemo, Groebner, RationalUnivariateRepresentation
import AbstractAlgebra

# Univariate
# R, (x,) = polynomial_ring(QQ, ["x"])
# @test rur([x^2 - 5]) == [[-5, 0, 1], [0, 0, 2]]
# @test rur([x - 5]) == [[-5, 1]]

# Last variable is not separating, but other variables are
R, (x0,x2,x1) = polynomial_ring(AbstractAlgebra.QQ, ["x0","x2","x1"])
example = [x0^2 + 2*x1^2 + 2*x2^2 - x0, 2*x0*x1 + 2*x1*x2 - x1, x0 + 2*x1 + 2*x2 - 1]

# Check that RUR is a solution to the original system
# modulo the minimal polynomial of T.
p = 2^60 + 33
Rqq, _ = polynomial_ring(QQ, "T")
Rzp, _ = polynomial_ring(GF(p), "T")
for sys in [
	    Groebner.Examples.hexapod(),
	    Groebner.Examples.rootn(4) .^ 2,  # non-radical
	    Groebner.Examples.rootn(4) .^ 3,  # non-radical
	    Groebner.Examples.katsuran(8),
	    Groebner.Examples.noonn(5),
	    Groebner.Examples.chandran(9),    # large coefficients (20k bits)
	    example,                          # last variable not separating
	    ]
    rr, sep = zdim_parameterization(sys, get_separating_element=true)
    nv = length(gens(parent(sys[1])))
    @assert length(rr) == nv + 1 && length(sep) == nv
    sys_zp = map(f -> change_base_ring(GF(p), f), sys) 
    rr_zp = map(f -> map(c -> GF(p)(numerator(c)) // GF(p)(denominator(c)), f), rr) 
    xi_zp = [mod(Rzp(rr_zp[i]) * invmod(derivative(Rzp(rr_zp[1])), Rzp(rr_zp[1])), Rzp(rr_zp[1])) for i in 2:nv+1];
    
    # f(T) = 0, x_i = p(T) / q(T)
    # system(x_i(T)) = 0 mod f(T)
    @test all(iszero, map(eq -> mod(evaluate(eq, xi_zp), Rzp(rr_zp[1])), sys_zp))
end

# Check with trivial equations x_i = C
R, (x,y,z) = polynomial_ring(AbstractAlgebra.QQ, ["x","y","z"])
# x + 1 = 0
sys = [x^2*y^2*z^2 - 1, x - y^3 - 1, y*z - 1]
rr, sep = zdim_parameterization(sys, get_separating_element=true)
@test sep == [0, 0, 1]
@test rr == [[1//2, 0, 0, 1], [0, 0, -3], [0, 3], [0, 0, 0, 3]]

