
R,(z0,z1,z2,z3,z4,z5,z6,_Z)=polynomial_ring(QQ, ["z0","z1","z2","z3","z4","z5","z6","_Z"],ordering=:degrevlex)

sys = [z0 + z1 + z2 + z3 + z4 + z5 + z6, z0*z1 + z1*z2 + z2*z3 + z3*z4 + z4*z5 + z5*z6 + z6*z0, z0*z1*z2 + z1*z2*z3 + z2*z3*z4 + z3*z4*z5 + z4*z5*z6 + z5*z6*z0 + z6*z0*z1, z0*z1*z2*z3 + z1*z2*z3*z4 + z2*z3*z4*z5 + z3*z4*z5*z6 + z4*z5*z6*z0+ z5*z6*z0*z1 + z6*z0*z1*z2, z0*z1*z2*z3*z4 + z1*z2*z3*z4*z5 + z2*z3*z4*z5*z6 + z3*z4*z5*z6*z0 + z4*z5*z6*z0*z1 + z5*z6*z0*z1*z2 + z6*z0*z1*z2*z3, z0*z1*z2*z3*z4*z5 + z1*z2*z3*z4*z5*z6 + z2*z3*z4*z5*z6*z0 + z3*z4*z5*z6*z0*z1+ z4*z5*z6*z0*z1*z2 + z5*z6*z0*z1*z2*z3 + z6*z0*z1*z2*z3*z4, z0*z1*z2*z3*z4*z5*z6 - 1,_Z - z0 - 2*z1 - 5*z2 + 10*z3 + 5*z4 - 2*z5 - 3*z6]
