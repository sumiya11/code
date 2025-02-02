using Test, Nemo, Groebner, RationalUnivariateRepresentation
import AbstractAlgebra

# Univariate
begin
	R, (x,) = polynomial_ring(QQ, ["x"])
	@test zdim_parameterization([x^2 - 5]) == [[-5, 0, 1], [0, 0, 2]]

	# Fails
	# @test zdim_parameterization([x - 5]) == [[-5, 1]]
end

# Access to core RUR
begin
	R, (x,y) = polynomial_ring(AbstractAlgebra.GF(2^30+3), ["x","y"])
	flag, rur = RationalUnivariateRepresentation.rur_core([x^2 + 1, y^2])
	@test !flag

	flag, rur = RationalUnivariateRepresentation.rur_core([x + 1, y^2])
	@test flag && rur == [[0x00000000, 0x00000001], [0x40000002]]

	# Fails
	# flag, rur = RationalUnivariateRepresentation.rur_core([R(0), x^2 + 1, y^2])
end

# Check that RUR is a solution to the original system
# modulo the minimal polynomial of T for both AA and Nemo
for interface in [Nemo, AbstractAlgebra]

# Last variable is not separating, but other variables are
R, (x0,x2,x1) = polynomial_ring(interface.QQ, ["x0","x2","x1"])
example = [x0^2 + 2*x1^2 + 2*x2^2 - x0, 2*x0*x1 + 2*x1*x2 - x1, x0 + 2*x1 + 2*x2 - 1]
example2 = [x0^2 - 1]

@test_throws DomainError zdim_parameterization(example2)

p = 2^60 + 33
gf_p = Nemo.GF(p)
Rqq, _ = polynomial_ring(interface.QQ, "T")
Rzp, _ = polynomial_ring(gf_p, "T")
for sys in [
	    Groebner.Examples.hexapod(np=interface),
	    Groebner.Examples.rootn(4, np=interface) .^ 2,  # non-radical
	    Groebner.Examples.rootn(4, np=interface) .^ 3,  # non-radical
	    Groebner.Examples.katsuran(8, np=interface),
	    Groebner.Examples.noonn(5, np=interface),
	    Groebner.Examples.chandran(9, np=interface),    # large coefficients (20k bits)
	    example,                          # last variable not separating
	    ]
    rr, sep = zdim_parameterization(sys, get_separating_element=true)
    nv = length(gens(parent(sys[1])))
    @assert length(rr) == nv + 1 && length(sep) == nv
    sys_zp = map(f -> change_base_ring(gf_p, f), sys)
    rr_zp = map(f -> map(c -> gf_p(numerator(c)) // gf_p(denominator(c)), f), rr)
    xi_zp = [mod(Rzp(rr_zp[i]) * invmod(derivative(Rzp(rr_zp[1])), Rzp(rr_zp[1])), Rzp(rr_zp[1])) for i in 2:nv+1];
    
    # f(T) = 0, x_i = p(T) / q(T)
    # system(x_i(T)) = 0 mod f(T)
    @test all(iszero, map(eq -> mod(evaluate(eq, xi_zp), Rzp(rr_zp[1])), sys_zp))
end
end # interface loop

# Check with trivial equations x_i = C
begin
	R, (x,y,z) = polynomial_ring(AbstractAlgebra.QQ, ["x","y","z"])
	# x + 1 = 0
	sys = [x^2*y^2*z^2 - 1, x - y^3 - 1, y*z - 1]
	rr, sep = zdim_parameterization(sys, get_separating_element=true)
	@test sep == [0, 0, 1]
	@test rr == [[1//2, 0, 0, 1], [0, 0, -3], [0, 3], [0, 0, 0, 3]]
end

# The quotient contains two elements
begin
	R, (_z__tpk2__d, _z__tpk1__d, _z__tpk3__d, _z__t29_w_t__d, _z__t29_wˍt_t__d, _z__t29_wˍtt_t__d, _z__t29_wˍttt_t__d) = polynomial_ring(AbstractAlgebra.QQ, ["_z__tpk2__d", "_z__tpk1__d", "_z__tpk3__d", "_z__t29_w_t__d", "_z__t29_wˍt_t__d", "_z__t29_wˍtt_t__d", "_z__t29_wˍttt_t__d"])
	sys = [479//200*_z__tpk2__d*_z__t29_w_t__d - 479//200*_z__tpk1__d - 837//1000, -837//1000*_z__tpk2__d*_z__t29_w_t__d + 479//200*_z__tpk2__d*_z__t29_wˍt_t__d + 837//1000*_z__tpk1__d - 651//250, -651//250*_z__tpk2__d*_z__t29_w_t__d - 837//500*_z__tpk2__d*_z__t29_wˍt_t__d + 479//200*_z__tpk2__d*_z__t29_wˍtt_t__d + 651//250*_z__tpk1__d + 43//25, 43//25*_z__tpk2__d*_z__t29_w_t__d - 7811//1000*_z__tpk2__d*_z__t29_wˍt_t__d - 251//100*_z__tpk2__d*_z__t29_wˍtt_t__d + 479//200*_z__tpk2__d*_z__t29_wˍttt_t__d - 43//25*_z__tpk1__d + 783//50, -479//200*_z__tpk2__d*_z__t29_w_t__d + _z__tpk3__d*_z__t29_w_t__d + _z__t29_wˍt_t__d, 837//1000*_z__tpk2__d*_z__t29_w_t__d - 479//200*_z__tpk2__d*_z__t29_wˍt_t__d + _z__tpk3__d*_z__t29_wˍt_t__d + _z__t29_wˍtt_t__d, 651//250*_z__tpk2__d*_z__t29_w_t__d + 837//500*_z__tpk2__d*_z__t29_wˍt_t__d - 479//200*_z__tpk2__d*_z__t29_wˍtt_t__d + _z__tpk3__d*_z__t29_wˍtt_t__d + _z__t29_wˍttt_t__d]
	rr, sep = zdim_parameterization(sys, get_separating_element=true)
	nv = 7
	R, _ = QQ["T"]
	xi = [mod(R(rr[i]) * invmod(derivative(R(rr[1])), R(rr[1])), R(rr[1])) for i in 2:nv+1]
	@test all(iszero, map(eq -> mod(evaluate(eq, xi), R(rr[1])), sys))
end

