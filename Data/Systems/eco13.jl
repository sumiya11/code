R,(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)=polynomial_ring(QQ, ["x0","x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12"])

sys=[x0*x1*x12+x1*x2*x12+x2*x3*x12+x3*x4*x12+x4*x5*x12+x5*x6*x12+x6*x7*x12+x7*x8*x12+x8*x9*x12+x9*x10*x12+x10*x11*x12+x0*x12-1,
x0*x2*x12+x1*x3*x12+x2*x4*x12+x3*x5*x12+x4*x6*x12+x5*x7*x12+x6*x8*x12+x7*x9*x12+x8*x10*x12+x9*x11*x12+x1*x12-2,
x0*x3*x12+x1*x4*x12+x2*x5*x12+x3*x6*x12+x4*x7*x12+x5*x8*x12+x6*x9*x12+x7*x10*x12+x8*x11*x12+x2*x12-3,
x0*x4*x12+x1*x5*x12+x2*x6*x12+x3*x7*x12+x4*x8*x12+x5*x9*x12+x6*x10*x12+x7*x11*x12+x3*x12-4,
x0*x5*x12+x1*x6*x12+x2*x7*x12+x3*x8*x12+x4*x9*x12+x5*x10*x12+x6*x11*x12+x4*x12-5,
x0*x6*x12+x1*x7*x12+x2*x8*x12+x3*x9*x12+x4*x10*x12+x5*x11*x12+x5*x12-6,
x0*x7*x12+x1*x8*x12+x2*x9*x12+x3*x10*x12+x4*x11*x12+x6*x12-7,
x0*x8*x12+x1*x9*x12+x2*x10*x12+x3*x11*x12+x7*x12-8,
x0*x9*x12+x1*x10*x12+x2*x11*x12+x8*x12-9,
x0*x10*x12+x1*x11*x12+x9*x12-10,
x0*x11*x12+x10*x12-11,
x11*x12-12,
x0+x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+1]
