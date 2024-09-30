using Test, TestSetExtensions

using AbstractAlgebra, Groebner, SHA
using RationalUnivariateRepresentation

rur = zdim_parameterization

R, (x,) = polynomial_ring(QQ, ["x"])
@test rur([x^2 - 5]) == [[-5, 0, 1], [0, 0, 2]]
# @test rur([x - 5]) == [[-5, 1]]  # fails

R, (x,y) = polynomial_ring(QQ, ["x","y"])
@test rur([x - y^2, y^2 + 3]) == [[3, 0, 1], [0, -6], [0, 0, 2]]

c = Groebner.Examples.katsuran(8)
rr = zdim_parameterization(c)
@test sha3_256(string(rr)) == UInt8[0x95, 0xa9, 0xbb, 0x01, 0x51, 0x12, 0xd8, 0x66, 0x86, 0x16, 0x0a, 0xb3, 0x28, 0x2c, 0x79, 0x07, 0xb2, 0x99, 0xc0, 0x1e, 0x59, 0x9d, 0x6e, 0xed, 0x20, 0x89, 0x5e, 0x23, 0x4c, 0x0f, 0x68, 0xbb]

