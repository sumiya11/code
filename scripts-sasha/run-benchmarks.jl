using Pkg;
pkg"activate here"
pkg"status"

###

using Primes,HostCPUFeatures,AbstractAlgebra,TimerOutputs
using BenchmarkTools
include("../../Groebner.jl/src/Groebner.jl")

include((@__DIR__)*"/../Julia/rur.jl")

nn = 27

systems = Dict(
    # "chandra9" => Groebner.chandran(9, ordering=:degrevlex),
    # "chandra10" => Groebner.chandran(10, ordering=:degrevlex),
    # "henrion6" => Groebner.henrion6(ordering=:degrevlex),
)
for name in [
    # (@__DIR__)*"/../Data/Systems/kat10.jl", 
    (@__DIR__)*"/../Data/Systems/katsura12.jl", 
    # (@__DIR__)*"/../Data/Systems/noon7.jl",
    (@__DIR__)*"/../Data/Systems/eco11.jl",
    (@__DIR__)*"/../Data/Systems/eco12.jl",
    # (@__DIR__)*"/../Data/Systems/cp_d_3_n_6_p_2.jl",
    # (@__DIR__)*"/../Data/Systems/cp_d_3_n_6_p_6.jl",
    ]
    include("$name")
    systems[name] = sys
end

###############################################################
###############################################################

_gauss_reduct_jam[] = 1

@time test_param(convert_sys_to_sys_z(Groebner.eco7(ordering=:degrevlex)), nn);

@info "" _gauss_reduct_jam[]
res1 = Dict()
for (name, sys) in systems
    @info "System $name"
    sys_z = convert_sys_to_sys_z(sys)
    @time rur = test_param(sys_z, nn);
    println("Minpoly len. = $(length(rur[1])))")
    res1[name] = rur;
end;

###############################################################
###############################################################

_gauss_reduct_jam[] = 2

@time test_param(convert_sys_to_sys_z(Groebner.eco7(ordering=:degrevlex)), nn);

@info "" _gauss_reduct_jam[]
res2 = Dict()
for (name, sys) in systems
    @info "System $name"
    sys_z = convert_sys_to_sys_z(sys)
    @time rur = test_param(sys_z, nn);
    res2[name] = rur;
end;

###############################################################
###############################################################

_gauss_reduct_jam[] = 4

@time test_param(convert_sys_to_sys_z(Groebner.eco7(ordering=:degrevlex)), nn);

@info "" _gauss_reduct_jam[]
res3 = Dict()
for (name, sys) in systems
    @info "System $name"
    sys_z = convert_sys_to_sys_z(sys)
    @time rur = test_param(sys_z, nn);
    res3[name] = rur;
end;

@assert res1 == res2 == res3
