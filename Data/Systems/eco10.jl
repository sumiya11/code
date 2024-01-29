R,(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)=polynomial_ring(AbstractAlgebra.QQ, [:x1,:x2,:x3,:x4,:x5,:x6,:x7,:x8,:x9,:x10], ordering=:degrevlex)
sys=[x1*x2*x10 + x2*x3*x10 + x3*x4*x10 + x4*x5*x10 + x5*x6*x10 + x6*x7*x10 + x7*x8*x10 + x8*x9*x10 + x1*x10 - 1,
x1*x3*x10 + x2*x4*x10 + x3*x5*x10 + x4*x6*x10 + x5*x7*x10 + x6*x8*x10 + x7*x9*x10 + x2*x10 - 2,
x1*x4*x10 + x2*x5*x10 + x3*x6*x10 + x4*x7*x10 + x5*x8*x10 + x6*x9*x10 + x3*x10 - 3,
x1*x5*x10 + x2*x6*x10 + x3*x7*x10 + x4*x8*x10 + x5*x9*x10 + x4*x10 - 4,
x1*x6*x10 + x2*x7*x10 + x3*x8*x10 + x4*x9*x10 + x5*x10 - 5,
x1*x7*x10 + x2*x8*x10 + x3*x9*x10 + x6*x10 - 6,
x1*x8*x10 + x2*x9*x10 + x7*x10 - 7,
x1*x9*x10 + x8*x10 - 8,
x9*x10 - 9,
x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + 1]
