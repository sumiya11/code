using Test, TestSetExtensions

using AbstractAlgebra
using RationalUnivariateRepresentation

rur = zdim_parameterization

R, (x,) = polynomial_ring(QQ, ["x"])
@test rur([x^2 - 5]) == [[-5, 0, 1]]
# @test rur([x - 5]) == [[-5, 1]]  # fails

R, (x,y) = polynomial_ring(QQ, ["x","y"])
@test rur([x - y^2, y^2 + 3]) == [[3, 0, 1], [0, -6]]
