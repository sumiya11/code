using Pkg;
pkg"activate here"
# pkg"add Nemo@0.38.3"
# pkg"add AbstractAlgebra@0.34"
# pkg"add Atomix"
# pkg"add Primes"
# pkg"add HostCPUFeatures"
# pkg"add ExprTools"
# pkg"add Combinatorics"
# pkg"add MultivariatePolynomials"
# pkg"add PrettyTables"
# pkg"add TimerOutputs"
# pkg"add JET"
# pkg"add Test"
# pkg"add LoopVectorization"
pkg"status"

###

using Primes,HostCPUFeatures,AbstractAlgebra,TimerOutputs
using BenchmarkTools
using JET
# using Traceur
include("../../Groebner.jl/src/Groebner.jl")

include((@__DIR__)*"/../Julia/rur.jl")
# include((@__DIR__)*"/rur_old.jl")

macro my_profview(ex)
    :((VSCodeServer.Profile).clear();
        VSCodeServer.Profile.init(n=10^8, delay=0.001);
    VSCodeServer.Profile.start_timer();
        $ex;
    VSCodeServer.Profile.stop_timer();
        VSCodeServer.view_profile(;))
end

###

nn=27;
tmr = TimerOutputs.TimerOutput()
sys = Groebner.eco11(k=QQ, ordering=:degrevlex);
include("../Data/Systems/eco11.jl");
sys_z = convert_sys_to_sys_z(sys);
dm,Dq,sys_T,_vars=prepare_system(sys_z, 27,R);
@time rur_jam = test_param(sys_T, nn);
@report_opt target_modules=(Main,) rur_jam = test_param(sys_T, nn);
@profview_allocs rur_jam = test_param(sys_T, nn);

d = 200
R, (x,y) = polynomial_ring(QQ, ["x","y"], ordering=:degrevlex)
sys = [sum(y^i * rand(1:1) for i in 1:d), x - sum(y^i * rand(1:1) for i in 1:div(d,2))]
#include("../Data/Systems/henrion6.jl");
gb = Groebner.groebner(sys, loglevel=0);

d = 500
R, (x,y) = polynomial_ring(QQ, ["x","y"], ordering=:degrevlex)
sys = [sum(y^i * rand(1:1) for i in 1:d), x - sum(y^i * rand(1:2) for i in 1:div(d,2))]
#include("../Data/Systems/henrion6.jl");
# gb = Groebner.groebner(sys, loglevel=-3);
sys_z = convert_sys_to_sys_z(sys);
# _, _, sys_T = prepare_system(sys_z, nn, R)
rur = test_param(sys_z, nn);
@my_profview test_param(sys_z, nn);

_gauss_reduct_jam[] = 4
# s = @report_opt target_modules=(Main,) test_param(sys_z, nn)

@time rur_jam = test_param(sys_z, nn);
@my_profview rur_jam = test_param(sys_z, nn);

begin
    reset_timer!(tmr)
    _gauss_reduct_jam[] = 4
    @time rur_jam_4 = test_param(sys_z, nn);
    tmr
end
begin
    reset_timer!(tmr)
    _gauss_reduct_jam[] = 2
    @time rur_jam_2 = test_param(sys_z, nn);
    tmr
end
begin
    reset_timer!(tmr)
    _gauss_reduct_jam[] = 1
    @time rur = test_param(sys_z, nn);
    tmr
end
rur == rur_jam_2 == rur_jam_4

###########################################
###########################################

function foo!(x,a1,v1,a2,v2,a3,v3,a4,v4)
    add_mul!(x,a1,v1)
    add_mul!(x,a2,v2)
    add_mul!(x,a3,v3)
    add_mul!(x,a4,v4)
end
function woo!(x,a1,v1,a2,v2,a3,v3,a4,v4)
    add_mul_jam4!(x,a1,v1,a2,v2,a3,v3,a4,v4)
end

n = 1_000_000
@btime foo!(x,a1,v1,a2,v2,a3,v3,a4,v4) setup=begin
    x = rand(UInt64, $n)
    a1 = rand(UInt32)
    v1 = rand(UInt32, $n)
    a2 = rand(UInt32)
    v2 = rand(UInt32, $n)
    a3 = rand(UInt32)
    v3 = rand(UInt32, $n)
    a4 = rand(UInt32)
    v4 = rand(UInt32, $n)
end
#  1.700 μs (0 allocations: 0 bytes)

@btime woo!(x,a1,v1,a2,v2,a3,v3,a4,v4) setup=begin
    x = rand(UInt64, $n)
    a1 = rand(UInt32)
    v1 = rand(UInt32, $n)
    a2 = rand(UInt32)
    v2 = rand(UInt32, $n)
    a3 = rand(UInt32)
    v3 = rand(UInt32, $n)
    a4 = rand(UInt32)
    v4 = rand(UInt32, $n)
end
#  1.280 μs (0 allocations: 0 bytes)

x = rand(UInt64, n)
a1 = rand(UInt32)
v1 = rand(UInt32, n)
a2 = rand(UInt32)
v2 = rand(UInt32, n)
a3 = rand(UInt32)
v3 = rand(UInt32, n)
a4 = rand(UInt32)
v4 = rand(UInt32, n)

io = open((@__DIR__)*"/add_mul_jam4_llvm.txt", "w")
code_llvm(io, add_mul_jam4!, map(typeof, (x,a1,v1,a2,v2,a3,v3,a4,v4)), debuginfo=:none)
close(io)

io = open((@__DIR__)*"/add_mul_jam4_native.txt", "w")
code_native(io, add_mul_jam4!, map(typeof, (x,a1,v1,a2,v2,a3,v3,a4,v4)), debuginfo=:none)
close(io)


###########################################
###########################################

@my_profview test_param(sys_z, nn);

###
  
systems = Dict(
    "chandra9" => Groebner.chandran(9, ordering=:degrevlex),
    "henrion6" => Groebner.henrion6(ordering=:degrevlex),
)
for name in ["systems/noon6", "systems/eco11"]
    include("$name.jl")
    systems[name] = sys
end

res2 = Dict()
for (name, sys) in systems
    @info "System $name"
    sys_z = convert_sys_to_sys_z(sys)
    @time rur = test_param(sys_z, nn);
    res2[name] = rur;
end;
#=
[ Info: System henrion6
[2,2](0)[4,6](0)[4,10](0)[4,14](0)[4,18](0)[4,22](0)[4,26](0)[4,30](0)[4,34](0)[4,38](0)[4,42](0)[4,46](0)[4,50](0)[4,54](0)[4,58](0)[4,62](0)[4,66](0)[4,70](0)[4,74](0)[4,78](0)[4,82](0)[8,90](0)[8,98](0)[8,106](0)[8,114](0)[8,122](0)[8,130](0)[8,138](0)[8,146](0)[8,154](0)[8,162](0)[16,178](0)[16,194](0)[16,210](0)[16,226](0)[16,242](0)[16,258](0)[16,274](0)[16,290](0)[16,306](0)[16,322](0)[32,354](0)(1) 
71.698802 seconds (9.58 M allocations: 13.167 GiB, 3.58% gc time, 0.13% compilation time)
[ Info: System systems/eco11
[2,2](0)[4,6](0)[4,10](0)[4,14](0)[4,18](0)[4,22](0)[4,26](0)[4,30](0)[4,34](0)[4,38](0)[4,42](0)[4,46](0)[4,50](0)[4,54](0)[4,58](0)[4,62](0)[4,66](0)[4,70](0)[4,74](0)[4,78](0)[4,82](0)[8,90](0)[8,98](0)[8,106](0)[8,114](0)[8,122](0)[8,130](0)(1) 
42.047049 seconds (13.54 M allocations: 7.814 GiB, 4.72% gc time)
[ Info: System chandra9
[2,2](0)[4,6](0)[4,10](0)[4,14](0)[4,18](0)[4,22](0)[4,26](0)[4,30](0)[4,34](0)[4,38](0)[4,42](0)[4,46](0)[4,50](0)[4,54](0)[4,58](0)[4,62](0)[4,66](0)[4,70](0)[4,74](0)[4,78](0)[4,82](0)[8,90](0)[8,98](0)[8,106](0)[8,114](0)[8,122](0)[8,130](0)[8,138](0)[8,146](0)[8,154](0)[8,162](0)[16,178](0)[16,194](0)[16,210](0)[16,226](0)[16,242](0)[16,258](0)[16,274](0)[16,290](0)[16,306](0)[16,322](0)[32,354](0)[32,386](0)[32,418](0)[32,450](0)[32,482](0)[32,514](0)[32,546](0)[32,578](0)[32,610](0)[32,642](0)[64,706](0)[64,770](0)[64,834](0)[64,898](0)[64,962](0)[64,1026](0)[64,1090](0)[64,1154](0)[64,1218](0)[64,1282](0)(1) 
39.240365 seconds (20.42 M allocations: 10.336 GiB, 5.64% gc time)
[ Info: System systems/noon6
[2,2](0)[4,6](0)[4,10](0)[4,14](0)[4,18](0)[4,22](0)[4,26](0)[4,30](0)[4,34](0)[4,38](0)[4,42](0)[4,46](0)[4,50](0)[4,54](0)[4,58](0)[4,62](0)[4,66](0)[4,70](0)[4,74](0)[4,78](0)[4,82](0)[8,90](0)[8,98](0)[8,106](0)[8,114](0)[8,122](0)[8,130](0)[8,138](0)[8,146](0)[8,154](0)[8,162](0)[16,178](0)[16,194](0)[16,210](0)[16,226](0)[16,242](0)[16,258](0)[16,274](0)[16,290](0)[16,306](0)[16,322](0)[32,354](0)[32,386](0)[32,418](0)[32,450](0)[32,482](0)(1)
116.138241 seconds (15.38 M allocations: 21.364 GiB, 3.48% gc time, 0.06% compilation time)

[ Info: System henrion6
[2,2](0)[4,6](0)[4,10](0)[4,14](0)[4,18](0)[4,22](0)[4,26](0)[4,30](0)[4,34](0)[4,38](0)[4,42](0)[4,46](0)[4,50](0)[4,54](0)[4,58](0)[4,62](0)[4,66](0)[4,70](0)[4,74](0)[4,78](0)[4,82](0)[8,90](0)[8,98](0)[8,106](0)[8,114](0)[8,122](0)[8,130](0)[8,138](0)[8,146](0)[8,154](0)[8,162](0)[16,178](0)[16,194](0)[16,210](0)[16,226](0)[16,242](0)[16,258](0)[16,274](0)[16,290](0)[16,306](0)[16,322](0)[32,354](0)(1) 
84.769461 seconds (9.60 M allocations: 13.167 GiB, 3.11% gc time)
[ Info: System systems/eco11
[2,2](0)[4,6](0)[4,10](0)[4,14](0)[4,18](0)[4,22](0)[4,26](0)[4,30](0)[4,34](0)[4,38](0)[4,42](0)[4,46](0)[4,50](0)[4,54](0)[4,58](0)[4,62](0)[4,66](0)[4,70](0)[4,74](0)[4,78](0)[4,82](0)[8,90](0)[8,98](0)[8,106](0)[8,114](0)[8,122](0)[8,130](0)(1) 
42.507408 seconds (13.53 M allocations: 7.813 GiB, 4.50% gc time)
[ Info: System chandra9
[2,2](0)[4,6](0)[4,10](0)[4,14](0)[4,18](0)[4,22](0)[4,26](0)[4,30](0)[4,34](0)[4,38](0)[4,42](0)[4,46](0)[4,50](0)[4,54](0)[4,58](0)[4,62](0)[4,66](0)[4,70](0)[4,74](0)[4,78](0)[4,82](0)[8,90](0)[8,98](0)[8,106](0)[8,114](0)[8,122](0)[8,130](0)[8,138](0)[8,146](0)[8,154](0)[8,162](0)[16,178](0)[16,194](0)[16,210](0)[16,226](0)[16,242](0)[16,258](0)[16,274](0)[16,290](0)[16,306](0)[16,322](0)[32,354](0)[32,386](0)[32,418](0)[32,450](0)[32,482](0)[32,514](0)[32,546](0)[32,578](0)[32,610](0)[32,642](0)[64,706](0)[64,770](0)[64,834](0)[64,898](0)[64,962](0)[64,1026](0)[64,1090](0)[64,1154](0)[64,1218](0)[64,1282](0)(1) 
45.977081 seconds (20.39 M allocations: 10.336 GiB, 6.15% gc time)
[ Info: System systems/noon6
[2,2](0)[4,6](0)[4,10](0)[4,14](0)[4,18](0)[4,22](0)[4,26](0)[4,30](0)[4,34](0)[4,38](0)[4,42](0)[4,46](0)[4,50](0)[4,54](0)[4,58](0)[4,62](0)[4,66](0)[4,70](0)[4,74](0)[4,78](0)[4,82](0)[8,90](0)[8,98](0)[8,106](0)[8,114](0)[8,122](0)[8,130](0)[8,138](0)[8,146](0)[8,154](0)[8,162](0)[16,178](0)[16,194](0)[16,210](0)[16,226](0)[16,242](0)[16,258](0)[16,274](0)[16,290](0)[16,306](0)[16,322](0)[32,354](0)[32,386](0)[32,418](0)[32,450](0)[32,482](0)(1)
131.479292 seconds (15.35 M allocations: 21.362 GiB, 3.27% gc time)
=#
