using Pkg;
pkg"activate here"
pkg"status"

###

using Primes,HostCPUFeatures,AbstractAlgebra,TimerOutputs
using BenchmarkTools
include("../../Groebner.jl/src/Groebner.jl")

include((@__DIR__)*"/../Julia/rur.jl")

nn=27;
include("../Data/Systems/noon7.jl");
sys_z = convert_sys_to_sys_z(sys);

@time rur_jam = test_param(sys_z, nn);
