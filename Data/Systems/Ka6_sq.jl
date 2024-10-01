R,(x0,x1,x2,x3,x4,x5,x6)=polynomial_ring(QQ,["x0","x1","x2","x3","x4","x5","x6"];
sys = [-2*x0^4 + 28*x0^2*x1^2 - 8*x1^4 + 72*x0*x1^2*x2 - 8*x0^2*x2^2 + 60*x1^2*x2^2 + 152*x0*x1*x2*x3 + 104*x1*x2^2*x3 + 32*x0^2*x3^2 + 16*x1^2*x3^2 + 20*x2^2*x3^2 - 8*x3^4 + 80*x1^2*x2*x4 + 32*x0*x2^2*x4 + 216*x0*x1*x3*x4 + 72*x1*x2*x3*x4 + 72*x2*x3^2*x4 + 24*x0^2*x4^2 + 24*x1^2*x4^2 - 16*x2^2*x4^2 + 20*x3^2*x4^2 - 8*x4^4 + 112*x1*x2^2*x5 + 64*x1^2*x3*x5 + 80*x0*x2*x3*x5 + 136*x0*x1*x4*x5 + 152*x1*x2*x4*x5 + 72*x2*x3*x4*x5 + 72*x3*x4^2*x5 - 8*x0^2*x5^2 + 16*x1^2*x5^2 + 24*x2^2*x5^2 - 16*x3^2*x5^2 + 20*x4^2*x5^2 - 8*x5^4 + 32*x2^3*x6 + 144*x1*x2*x3*x6 + 80*x0*x3^2*x6 + 64*x0*x2*x4*x6 + 80*x1*x3*x4*x6 + 72*x0*x1*x5*x6 + 136*x1*x2*x5*x6 + 152*x2*x3*x5*x6 + 72*x3*x4*x5*x6 + 72*x4*x5^2*x6 - 8*x0^2*x6^2 - 16*x1^2*x6^2 + 16*x2^2*x6^2 + 24*x3^2*x6^2 - 16*x4^2*x6^2 + 20*x5^2*x6^2 - 8*x6^4 + 4*x0^3 - 28*x0*x1^2 - 36*x1^2*x2 + 8*x0*x2^2 - 76*x1*x2*x3 - 32*x0*x3^2 - 16*x2^2*x4 - 108*x1*x3*x4 - 24*x0*x4^2 - 40*x2*x3*x5 - 68*x1*x4*x5 + 8*x0*x5^2 - 40*x3^2*x6 - 32*x2*x4*x6 - 36*x1*x5*x6 + 8*x0*x6^2 - 5*x0^2 - 12*x0*x1 - 3*x1^2 - 12*x0*x2 - 24*x1*x2 - 12*x2^2 - 12*x0*x3 - 24*x1*x3 - 24*x2*x3 - 2*x3^2 - 12*x0*x4 - 24*x1*x4 - 24*x2*x4 - 24*x3*x4 - 4*x4^2 - 12*x0*x5 - 24*x1*x5 - 24*x2*x5 - 24*x3*x5 - 24*x4*x5 - 12*x5^2 - 12*x0*x6 - 24*x1*x6 - 24*x2*x6 - 24*x3*x6 - 24*x4*x6 - 24*x5*x6 - 12*x6^2 + 6*x0 + 12*x1 + 12*x2 + 12*x3 + 12*x4 + 12*x5 + 12*x6 - 3,
2*x0^4 + 8*x0^2*x1^2 + 15*x1^4 + 28*x0*x1^2*x2 + 36*x0^2*x2^2 - 8*x1^2*x2^2 + 13*x2^4 + 28*x1^3*x3 + 8*x0*x1*x2*x3 + 20*x1*x2^2*x3 - 16*x0^2*x3^2 + 64*x1^2*x3^2 - 16*x2^2*x3^2 + 8*x3^4 - 20*x1^2*x2*x4 + 76*x0*x2^2*x4 - 8*x0*x1*x3*x4 - 8*x1*x2*x3*x4 + 28*x0^2*x4^2 - 40*x1^2*x4^2 + 44*x2^2*x4^2 + 16*x3^2*x4^2 + 8*x4^4 - 28*x1*x2^2*x5 + 68*x1^2*x3*x5 - 56*x0*x2*x3*x5 + 56*x1*x3^2*x5 - 24*x0*x1*x4*x5 - 48*x1*x2*x4*x5 + 56*x2*x3*x4*x5 - 24*x0^2*x5^2 + 36*x1^2*x5^2 - 8*x2^2*x5^2 + 44*x3^2*x5^2 + 16*x4^2*x5^2 + 8*x5^4 + 20*x2^3*x6 - 72*x1*x2*x3*x6 - 48*x0*x3^2*x6 - 36*x1^2*x4*x6 + 96*x0*x2*x4*x6 + 8*x1*x3*x4*x6 + 56*x2*x4^2*x6 - 64*x0*x1*x5*x6 + 40*x1*x2*x5*x6 - 48*x2*x3*x5*x6 + 56*x3*x4*x5*x6 + 8*x0^2*x6^2 - 16*x1^2*x6^2 + 36*x2^2*x6^2 - 8*x3^2*x6^2 + 44*x4^2*x6^2 + 16*x5^2*x6^2 + 8*x6^4 - 4*x0^3 - 8*x0*x1^2 - 14*x1^2*x2 - 36*x0*x2^2 - 4*x1*x2*x3 + 16*x0*x3^2 - 38*x2^2*x4 + 4*x1*x3*x4 - 28*x0*x4^2 + 28*x2*x3*x5 + 12*x1*x4*x5 + 24*x0*x5^2 + 24*x3^2*x6 - 48*x2*x4*x6 + 32*x1*x5*x6 - 8*x0*x6^2 + 12*x0^2 + 40*x0*x1 + 40*x1^2 + 40*x0*x2 + 80*x1*x2 + 47*x2^2 + 40*x0*x3 + 80*x1*x3 + 80*x2*x3 + 34*x3^2 + 40*x0*x4 + 80*x1*x4 + 80*x2*x4 + 80*x3*x4 + 45*x4^2 + 40*x0*x5 + 80*x1*x5 + 80*x2*x5 + 80*x3*x5 + 80*x4*x5 + 32*x5^2 + 40*x0*x6 + 80*x1*x6 + 80*x2*x6 + 80*x3*x6 + 80*x4*x6 + 80*x5*x6 + 40*x6^2 - 20*x0 - 40*x1 - 40*x2 - 40*x3 - 40*x4 - 40*x5 - 40*x6 + 10,
9*x0^4 + 72*x0^2*x1^2 + 45*x1^4 + 108*x0*x1^2*x2 + 72*x0^2*x2^2 + 68*x1^2*x2^2 + 39*x2^4 + 36*x1^3*x3 + 64*x0*x1*x2*x3 + 84*x1*x2^2*x3 - 4*x0^2*x3^2 + 120*x1^2*x3^2 + 108*x2^2*x3^2 + 36*x3^4 - 44*x1^2*x2*x4 + 84*x0*x2^2*x4 + 16*x0*x1*x3*x4 + 144*x1*x2*x3*x4 + 72*x2*x3^2*x4 + 48*x0^2*x4^2 + 32*x1^2*x4^2 + 108*x2^2*x4^2 + 108*x3^2*x4^2 + 36*x4^4 - 68*x1*x2^2*x5 + 60*x1^2*x3*x5 - 8*x0*x2*x3*x5 + 72*x1*x3^2*x5 + 96*x0*x1*x4*x5 - 8*x1*x2*x4*x5 + 144*x2*x3*x4*x5 + 72*x3*x4^2*x5 + 36*x0^2*x5^2 + 84*x1^2*x5^2 + 32*x2^2*x5^2 + 108*x3^2*x5^2 + 108*x4^2*x5^2 + 36*x5^4 + 12*x2^3*x6 - 56*x1*x2*x3*x6 - 80*x0*x3^2*x6 + 36*x1^2*x4*x6 + 96*x0*x2*x4*x6 - 8*x1*x3*x4*x6 + 72*x2*x4^2*x6 + 72*x0*x1*x5*x6 + 96*x1*x2*x5*x6 - 8*x2*x3*x5*x6 + 144*x3*x4*x5*x6 + 72*x4*x5^2*x6 + 36*x0^2*x6^2 + 72*x1^2*x6^2 + 84*x2^2*x6^2 + 32*x3^2*x6^2 + 108*x4^2*x6^2 + 108*x5^2*x6^2 + 36*x6^4 - 18*x0^3 - 72*x0*x1^2 - 54*x1^2*x2 - 72*x0*x2^2 - 32*x1*x2*x3 + 4*x0*x3^2 - 42*x2^2*x4 - 8*x1*x3*x4 - 48*x0*x4^2 + 4*x2*x3*x5 - 48*x1*x4*x5 - 36*x0*x5^2 + 40*x3^2*x6 - 48*x2*x4*x6 - 36*x1*x5*x6 - 36*x0*x6^2 + 9*x0^2 + 9*x1^2 + 9*x2^2 - 10*x3^2 + 3*x4^2,
2*x0^4 + 20*x0^2*x1^2 + 16*x1^4 + 56*x0*x1^2*x2 + 40*x0^2*x2^2 + 32*x1^2*x2^2 + 18*x2^4 + 32*x1^3*x3 + 96*x0*x1*x2*x3 + 64*x1*x2^2*x3 + 12*x0^2*x3^2 + 88*x1^2*x3^2 + 28*x2^2*x3^2 + 8*x3^4 + 40*x1^2*x2*x4 + 104*x0*x2^2*x4 + 112*x0*x1*x3*x4 + 88*x1*x2*x3*x4 + 24*x2*x3^2*x4 + 48*x0^2*x4^2 + 20*x1^2*x4^2 + 48*x2^2*x4^2 + 28*x3^2*x4^2 + 8*x4^4 + 48*x1*x2^2*x5 + 112*x1^2*x3*x5 + 72*x0*x2*x3*x5 + 64*x1*x3^2*x5 + 104*x0*x1*x4*x5 + 32*x1*x2*x4*x5 + 88*x2*x3*x4*x5 + 24*x3*x4^2*x5 + 8*x0^2*x5^2 + 56*x1^2*x5^2 + 20*x2^2*x5^2 + 48*x3^2*x5^2 + 28*x4^2*x5^2 + 8*x5^4 + 40*x2^3*x6 + 88*x1*x2*x3*x6 + 8*x0*x3^2*x6 + 32*x1^2*x4*x6 + 144*x0*x2*x4*x6 + 72*x1*x3*x4*x6 + 64*x2*x4^2*x6 + 24*x0*x1*x5*x6 + 104*x1*x2*x5*x6 + 32*x2*x3*x5*x6 + 88*x3*x4*x5*x6 + 24*x4*x5^2*x6 + 8*x0^2*x6^2 + 16*x1^2*x6^2 + 56*x2^2*x6^2 + 20*x3^2*x6^2 + 48*x4^2*x6^2 + 28*x5^2*x6^2 + 8*x6^4 - 4*x0^3 - 20*x0*x1^2 - 28*x1^2*x2 - 40*x0*x2^2 - 48*x1*x2*x3 - 12*x0*x3^2 - 52*x2^2*x4 - 56*x1*x3*x4 - 48*x0*x4^2 - 36*x2*x3*x5 - 52*x1*x4*x5 - 8*x0*x5^2 - 4*x3^2*x6 - 72*x2*x4*x6 - 12*x1*x5*x6 - 8*x0*x6^2 + 5*x0^2 + 12*x0*x1 + 15*x1^2 + 12*x0*x2 + 24*x1*x2 + 20*x2^2 + 12*x0*x3 + 24*x1*x3 + 24*x2*x3 + 13*x3^2 + 12*x0*x4 + 24*x1*x4 + 24*x2*x4 + 24*x3*x4 + 22*x4^2 + 12*x0*x5 + 24*x1*x5 + 24*x2*x5 + 24*x3*x5 + 24*x4*x5 + 12*x5^2 + 12*x0*x6 + 24*x1*x6 + 24*x2*x6 + 24*x3*x6 + 24*x4*x6 + 24*x5*x6 + 12*x6^2 - 6*x0 - 12*x1 - 12*x2 - 12*x3 - 12*x4 - 12*x5 - 12*x6 + 3,
-3*x0^4 - 4*x0^2*x1^2 - 15*x1^4 + 4*x0*x1^2*x2 - 24*x0^2*x2^2 - 5*x2^4 - 12*x1^3*x3 + 24*x0*x1*x2*x3 + 44*x1*x2^2*x3 + 4*x0^2*x3^2 - 8*x1^2*x3^2 - 16*x2^2*x3^2 - 12*x3^4 + 20*x1^2*x2*x4 + 4*x0*x2^2*x4 + 104*x0*x1*x3*x4 - 8*x1*x2*x3*x4 + 16*x2*x3^2*x4 + 16*x0^2*x4^2 - 8*x1^2*x4^2 - 36*x2^2*x4^2 - 16*x3^2*x4^2 - 12*x4^4 + 60*x1*x2^2*x5 + 44*x1^2*x3*x5 + 8*x0*x2*x3*x5 - 24*x1*x3^2*x5 + 72*x0*x1*x4*x5 + 48*x1*x2*x4*x5 - 8*x2*x3*x4*x5 + 16*x3*x4^2*x5 - 12*x0^2*x5^2 + 4*x1^2*x5^2 - 8*x2^2*x5^2 - 36*x3^2*x5^2 - 16*x4^2*x5^2 - 12*x5^4 + 28*x2^3*x6 + 88*x1*x2*x3*x6 + 32*x0*x3^2*x6 - 12*x1^2*x4*x6 + 32*x0*x2*x4*x6 + 8*x1*x3*x4*x6 - 24*x2*x4^2*x6 + 16*x0*x1*x5*x6 + 72*x1*x2*x5*x6 + 48*x2*x3*x5*x6 - 8*x3*x4*x5*x6 + 16*x4*x5^2*x6 - 12*x0^2*x6^2 - 24*x1^2*x6^2 + 4*x2^2*x6^2 - 8*x3^2*x6^2 - 36*x4^2*x6^2 - 16*x5^2*x6^2 - 12*x6^4 + 6*x0^3 + 4*x0*x1^2 - 2*x1^2*x2 + 24*x0*x2^2 - 12*x1*x2*x3 - 4*x0*x3^2 - 2*x2^2*x4 - 52*x1*x3*x4 - 16*x0*x4^2 - 4*x2*x3*x5 - 36*x1*x4*x5 + 12*x0*x5^2 - 16*x3^2*x6 - 16*x2*x4*x6 - 8*x1*x5*x6 + 12*x0*x6^2 - x0^2 + 8*x0*x1 + 10*x1^2 + 8*x0*x2 + 16*x1*x2 + 5*x2^2 + 8*x0*x3 + 16*x1*x3 + 16*x2*x3 + 12*x3^2 + 8*x0*x4 + 16*x1*x4 + 16*x2*x4 + 16*x3*x4 + 15*x4^2 + 8*x0*x5 + 16*x1*x5 + 16*x2*x5 + 16*x3*x5 + 16*x4*x5 + 8*x5^2 + 8*x0*x6 + 16*x1*x6 + 16*x2*x6 + 16*x3*x6 + 16*x4*x6 + 16*x5*x6 + 8*x6^2 - 4*x0 - 8*x1 - 8*x2 - 8*x3 - 8*x4 - 8*x5 - 8*x6 + 2,
6*x0^4 + 36*x0^2*x1^2 + 25*x1^4 + 28*x0*x1^2*x2 + 28*x0^2*x2^2 + 84*x1^2*x2^2 + 31*x2^4 + 4*x1^3*x3 + 80*x0*x1*x2*x3 + 52*x1*x2^2*x3 + 48*x0^2*x3^2 + 80*x1^2*x3^2 + 68*x2^2*x3^2 + 24*x3^4 + 52*x1^2*x2*x4 + 36*x0*x2^2*x4 + 128*x0*x1*x3*x4 + 48*x1*x2*x3*x4 + 24*x2*x3^2*x4 + 52*x0^2*x4^2 + 80*x1^2*x4^2 + 52*x2^2*x4^2 + 60*x3^2*x4^2 + 24*x4^4 + 76*x1*x2^2*x5 + 60*x1^2*x3*x5 + 72*x0*x2*x3*x5 + 8*x1*x3^2*x5 + 96*x0*x1*x4*x5 + 72*x1*x2*x4*x5 + 32*x2*x3*x4*x5 + 24*x3*x4^2*x5 + 32*x0^2*x5^2 + 76*x1^2*x5^2 + 72*x2^2*x5^2 + 52*x3^2*x5^2 + 60*x4^2*x5^2 + 24*x5^4 + 28*x2^3*x6 + 120*x1*x2*x3*x6 + 48*x0*x3^2*x6 + 20*x1^2*x4*x6 + 64*x0*x2*x4*x6 + 56*x1*x3*x4*x6 + 8*x2*x4^2*x6 + 40*x0*x1*x5*x6 + 80*x1*x2*x5*x6 + 72*x2*x3*x5*x6 + 32*x3*x4*x5*x6 + 24*x4*x5^2*x6 + 24*x0^2*x6^2 + 56*x1^2*x6^2 + 76*x2^2*x6^2 + 72*x3^2*x6^2 + 52*x4^2*x6^2 + 60*x5^2*x6^2 + 24*x6^4 - 12*x0^3 - 36*x0*x1^2 - 14*x1^2*x2 - 28*x0*x2^2 - 40*x1*x2*x3 - 48*x0*x3^2 - 18*x2^2*x4 - 64*x1*x3*x4 - 52*x0*x4^2 - 36*x2*x3*x5 - 48*x1*x4*x5 - 32*x0*x5^2 - 24*x3^2*x6 - 32*x2*x4*x6 - 20*x1*x5*x6 - 24*x0*x6^2 - 4*x0^2 - 40*x0*x1 - 37*x1^2 - 40*x0*x2 - 80*x1*x2 - 39*x2^2 - 40*x0*x3 - 80*x1*x3 - 80*x2*x3 - 34*x3^2 - 40*x0*x4 - 80*x1*x4 - 80*x2*x4 - 80*x3*x4 - 33*x4^2 - 40*x0*x5 - 80*x1*x5 - 80*x2*x5 - 80*x3*x5 - 80*x4*x5 - 38*x5^2 - 40*x0*x6 - 80*x1*x6 - 80*x2*x6 - 80*x3*x6 - 80*x4*x6 - 80*x5*x6 - 40*x6^2 + 20*x0 + 40*x1 + 40*x2 + 40*x3 + 40*x4 + 40*x5 + 40*x6 - 10,
-x0^4 - 44*x0^2*x1^2 - 3*x1^4 - 76*x0*x1^2*x2 - 8*x1^2*x2^2 + 3*x2^4 + 4*x1^3*x3 + 8*x0*x1*x2*x3 - 52*x1*x2^2*x3 + 36*x0^2*x3^2 + 24*x1^2*x3^2 - 40*x2^2*x3^2 - 4*x3^4 + 84*x1^2*x2*x4 + 36*x0*x2^2*x4 + 56*x0*x1*x3*x4 - 56*x1*x2*x3*x4 - 80*x2*x3^2*x4 + 24*x0^2*x4^2 + 40*x1^2*x4^2 - 4*x2^2*x4^2 - 48*x3^2*x4^2 - 4*x4^4 + 108*x1*x2^2*x5 + 60*x1^2*x3*x5 + 104*x0*x2*x3*x5 + 8*x1*x3^2*x5 - 8*x0*x1*x4*x5 - 72*x2*x3*x4*x5 - 80*x3*x4^2*x5 + 4*x0^2*x5^2 + 20*x1^2*x5^2 + 32*x2^2*x5^2 - 4*x3^2*x5^2 - 48*x4^2*x5^2 - 4*x5^4 + 28*x2^3*x6 + 152*x1*x2*x3*x6 + 80*x0*x3^2*x6 + 20*x1^2*x4*x6 + 64*x0*x2*x4*x6 + 88*x1*x3*x4*x6 + 8*x2*x4^2*x6 - 64*x0*x1*x5*x6 - 24*x1*x2*x5*x6 - 72*x3*x4*x5*x6 - 80*x4*x5^2*x6 - 4*x0^2*x6^2 + 20*x2^2*x6^2 + 32*x3^2*x6^2 - 4*x4^2*x6^2 - 48*x5^2*x6^2 - 4*x6^4 + 2*x0^3 + 44*x0*x1^2 + 38*x1^2*x2 - 4*x1*x2*x3 - 36*x0*x3^2 - 18*x2^2*x4 - 28*x1*x3*x4 - 24*x0*x4^2 - 52*x2*x3*x5 + 4*x1*x4*x5 - 4*x0*x5^2 - 40*x3^2*x6 - 32*x2*x4*x6 + 32*x1*x5*x6 + 4*x0*x6^2 + 5*x0^2 + 24*x0*x1 + 14*x1^2 + 24*x0*x2 + 48*x1*x2 + 25*x2^2 + 24*x0*x3 + 48*x1*x3 + 48*x2*x3 + 34*x3^2 + 24*x0*x4 + 48*x1*x4 + 48*x2*x4 + 48*x3*x4 + 31*x4^2 + 24*x0*x5 + 48*x1*x5 + 48*x2*x5 + 48*x3*x5 + 48*x4*x5 + 26*x5^2 + 24*x0*x6 + 48*x1*x6 + 48*x2*x6 + 48*x3*x6 + 48*x4*x6 + 48*x5*x6 + 24*x6^2 - 12*x0 - 24*x1 - 24*x2 - 24*x3 - 24*x4 - 24*x5 - 24*x6 + 6];
