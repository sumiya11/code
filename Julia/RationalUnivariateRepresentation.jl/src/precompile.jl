
@setup_workload begin
    # Putting some things in `setup` can reduce the size of the
    # precompile file and potentially make loading faster.
    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
	    R,(x,y,z,t) = AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ,["x","y","z","t"])
        sys = [y^2*z+2*x*y*t-2*x-z,-x^3*z+4*x*y^2*z+4*x^2*y*t+2*y^3*t+4*x^2-10*y^2+4*x*z-10*y*t+2,2*y*z*t+x*t^2-x-2*z,-x*z^3+4*y*z^2*t+4*x*z*t^2+2*y*t^3+4*x*z+4*z^2-10*y*t-10*t^2+2]
        zdim_parameterization(sys, verbose=false)
    end
end

precompile(
    zdim_parameterization,
    (Vector{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Rational{BigInt}}},)
)
