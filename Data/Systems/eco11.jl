R,(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)=polynomial_ring(QQ, ["x0","x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"],ordering=:degrevlex)

sys=[x0*x1*x10+x1*x2*x10+x2*x3*x10+x3*x4*x10+x4*x5*x10+x5*x6*x10+x6*x7*x10+x7*x8*x10+x8*x9*x10+x0*x10-1,
x0*x2*x10+x1*x3*x10+x2*x4*x10+x3*x5*x10+x4*x6*x10+x5*x7*x10+x6*x8*x10+x7*x9*x10+x1*x10-2,
x0*x3*x10+x1*x4*x10+x2*x5*x10+x3*x6*x10+x4*x7*x10+x5*x8*x10+x6*x9*x10+x2*x10-3,
x0*x4*x10+x1*x5*x10+x2*x6*x10+x3*x7*x10+x4*x8*x10+x5*x9*x10+x3*x10-4,
x0*x5*x10+x1*x6*x10+x2*x7*x10+x3*x8*x10+x4*x9*x10+x4*x10-5,
x0*x6*x10+x1*x7*x10+x2*x8*x10+x3*x9*x10+x5*x10-6,
x0*x7*x10+x1*x8*x10+x2*x9*x10+x6*x10-7,
x0*x8*x10+x1*x9*x10+x7*x10-8,
x0*x9*x10+x8*x10-9,
x9*x10-10,
x0+x1+x2+x3+x4+x5+x6+x7+x8+x9+1]

