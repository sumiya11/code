using Pkg;
pkg"activate here"
pkg"status"

###

using Primes,HostCPUFeatures,AbstractAlgebra,TimerOutputs
using BenchmarkTools, CPUSummary, InteractiveUtils
include("../../Groebner.jl/src/Groebner.jl")

include((@__DIR__)*"/../Julia/rur.jl")

versioninfo()
@info "" cache_linesize() cache_size(Val(1))

nn = 27

timings = Dict()
systems = Dict(
    "chandra8" => Groebner.chandran(8, ordering=:degrevlex),
    "chandra9" => Groebner.chandran(9, ordering=:degrevlex),
    # "chandra10" => Groebner.chandran(10, ordering=:degrevlex),
    "henrion6" => Groebner.henrion6(ordering=:degrevlex),
)
for name in [
    (@__DIR__)*"/../Data/Systems/kat10.jl", 
    (@__DIR__)*"/../Data/Systems/katsura11.jl", 
    # (@__DIR__)*"/../Data/Systems/katsura12.jl", 
    (@__DIR__)*"/../Data/Systems/noon6.jl",
    # (@__DIR__)*"/../Data/Systems/noon7.jl",
    (@__DIR__)*"/../Data/Systems/eco10.jl",
    (@__DIR__)*"/../Data/Systems/eco11.jl",
    # (@__DIR__)*"/../Data/Systems/eco12.jl",
    (@__DIR__)*"/../Data/Systems/cp_d_3_n_6_p_2.jl",
    (@__DIR__)*"/../Data/Systems/cp_d_3_n_6_p_6.jl",
    ]
    include("$name")
    systems[String(last(split(name, "/")))] = sys
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
    t = @timed @time test_param(sys_z, nn);
    println("Minpoly len. = $(length(t.value[1])))")
    !haskey(timings, name) && (timings[name] = (length(t.value[1]), Dict()))
    timings[name][2][1] = t.time
    res1[name] = t.value;
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
    t = @timed @time test_param(sys_z, nn);
    println("Minpoly len. = $(length(t.value[1])))")
    !haskey(timings, name) && (timings[name] = (length(t.value[1]), Dict()))
    timings[name][2][2] = t.time
    res2[name] = t.value;
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
    t = @timed @time test_param(sys_z, nn);
    println("Minpoly len. = $(length(t.value[1])))")
    !haskey(timings, name) && (timings[name] = (length(t.value[1]), Dict()))
    timings[name][2][4] = t.time
    res3[name] = t.value;
end;

@assert res1 == res2 == res3

###############################################################
###############################################################

buf = IOBuffer()
println(buf)

sorted_keys = sort(collect(keys(timings)))
for name in sorted_keys
    data = timings[name]
    println(buf, "$name, dim=$(data[1])")
    println(buf, "  scalar\t\t$(round(Int, data[2][1])) sec")
    println(buf, "  jam2\t\t\t$(round(Int, data[2][2])) sec")
    println(buf, "  jam4\t\t\t$(round(Int, data[2][4])) sec")
    println()
end

msg = String(take!(buf))

println(msg)

io = open((@__DIR__)*"/server_timings.txt", "w")
println(io, msg)
println(io)
versioninfo(io)
close(io)
