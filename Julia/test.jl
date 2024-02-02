using Pkg;
pkg"activate RUR2/scripts-sasha/here"
pkg"status"

using TimerOutputs
tmr = TimerOutputs.TimerOutput()

include("../../Groebner.jl/src/Groebner.jl")
include("rur.jl")

using AbstractAlgebra
using Test

#####################################################
# shape lemma case
#####################################################

R, (x, y) = AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ, ["x","y"], ordering=:degrevlex)

# TODO
sys = [R(1)]
sys_z = convert_sys_to_sys_z(sys)
# rur = test_param(sys_z, 27)

# TODO
sys = [x, y]
sys_z = convert_sys_to_sys_z(sys)
# rur = test_param(sys_z, 27)

sys = [x - 1, y^2 + 1]
sys_z = convert_sys_to_sys_z(sys)
rur = test_param(sys_z, 27)
@test rur == [
    [1//1, 0//1, 1//1],
    [0//1, 2//1]
]

sys = Groebner.noonn(2, ordering=:degrevlex)
sys_z = convert_sys_to_sys_z(sys)
rur = test_param(sys_z, 27)
@test rur == [
    [-11//10, 331//1100, 2//1, -11//5, -10//11, 1//1],
    [331//1210, 50//11, -6//1, -100//121, 10//11]
]

sys = Groebner.chandran(3, ordering=:degrevlex)
sys_z = convert_sys_to_sys_z(sys)
rur = test_param(sys_z, 27)
@test rur == [
    [13016198158497821231889892947963417080429844249640621965941145600//14780607353199594067361659251327125954359640739917121291481681, -4761897999097476895468565115080156576666723377371412609357578240//14780607353199594067361659251327125954359640739917121291481681, -3905565285999177228403161085968179887534211276835917014857940992//14780607353199594067361659251327125954359640739917121291481681, -41738612957090672010477264961536//3844555546900004677751591891209, 1//1],
    [3843162496864274560482853963759457628242393849363871432892169765218135421485056//852609813903724657420458088312468645897510725411277724196832567082114481175, -21334190464088875516062230172219016126022160945289137387797031114205215540117504//4263049069518623287102290441562343229487553627056388620984162835410572405875, 339298824644533344549760701685858163924876132352//1108853550823596269200473988933971789029022875, 24714223531135976//288421779135875],
    [209686554542814359838562044159012526305339401788116652656652070019310056112128//170521962780744931484091617662493729179502145082255544839366513416422896235, -399508987045470420925745310698271036481282262597663495569868502924547999137792//609007009931231898157470063080334747069650518150912660140594690772938915125, -1069444656135886602427912053268653909160117665792//1108853550823596269200473988933971789029022875, -12357111765567988//288421779135875],
]
#=
Gb in lex. is:
 H3^4 - 41738612957090672010477264961536//3844555546900004677751591891209*H3^3 - 3905565285999177228403161085968179887534211276835917014857940992//14780607353199594067361659251327125954359640739917121291481681*H3^2 - 4761897999097476895468565115080156576666723377371412609357578240//14780607353199594067361659251327125954359640739917121291481681*H3 + 13016198158497821231889892947963417080429844249640621965941145600//14780607353199594067361659251327125954359640739917121291481681
 H2 + 14780607353199594067361659251327125954359640739917121291481681//1000307844678716093896632144236723074886732225590881398292480000*H3^3 - 11876900645494429624492967440811026287247254373//74037653391698365449423501889227348466728960000*H3^2 - 21442821758858233455195955360009//12525456378617956425227304960000*H3 + 3089277941391997//2595796012222875
 H1 - 14780607353199594067361659251327125954359640739917121291481681//500153922339358046948316072118361537443366112795440699146240000*H3^3 + 11876900645494429624492967440811026287247254373//37018826695849182724711750944613674233364480000*H3^2 + 265436418719007774518919444256333//43839097325162847488295567360000*H3 - 24714223531135976//2595796012222875
=#

#####################################################
# general case
#####################################################

# Non-shape, y separates
R, (x, y) = AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ, ["x","y"], ordering=:degrevlex)
sys = [
    x^2,
    x*y,
    y^2
]
sys_z = convert_sys_to_sys_z(sys)
rur = test_param(sys_z, 27)
# @test rur == [
#     [0, 0, 1],  # y^2 = 0
#     []          # x   = 0
# ]

# This example is used for testing "biv_general" (_Z is not separating)
R, (t, _Z) = AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ, ["t", "_Z"], ordering=:degrevlex)
sys = [
    _Z^4 + _Z^2,
    t*_Z^2 + t + _Z^3 + _Z,
    t^4 + t*_Z^3 + t*_Z + _Z^2,
]
sys_z = convert_sys_to_sys_z(sys)
gb = Groebner.groebner(sys, ordering=Groebner.Lex())
rur = test_param(sys_z, 27)
# @test rur = ...

# Non-shape, _Z separates
R, (x, y, z, t, _Z) = AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ, ["x", "y", "z", "t", "_Z"], ordering=:degrevlex)
sys = [
    y^2 * z + 2 * x * y * t - 2 * x - z,
    -x^3 * z + 4 * x * y^2 * z + 4 * x^2 * y * t + 2 * y^3 * t + 4 * x^2 - 10 * y^2 + 4 * x * z - 10 * y * t + 2,
    2 * y * z * t + x * t^2 - x - 2 * z,
    -x * z^3 + 4 * y * z^2 * t + 4 * x * z * t^2 + 2 * y * t^3 + 4 * x * z + 4 * z^2 - 10 * y * t - 10 * t^2 + 2,
    x - 2y - 3z + 4t + _Z
]
sys_z = convert_sys_to_sys_z(sys)
gb = Groebner.groebner(map(f -> change_coefficient_ring(GF(2^30+3), f), sys_z), ordering=Groebner.Lex())
_, _, sys_z_sep = prepare_system(sys_z, 27, R)
rur = test_param(sys_z_sep, 27)
# @test rur = ...

