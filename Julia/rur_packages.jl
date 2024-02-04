list_pkg=["AbstractAlgebra","Atomix","Combinatorics","ExprTools","HostCPUFeatures","MultivariatePolynomials","Nemo","PrettyTables","Primes","TimerOutputs","PrecompileTools","LoopVectorization","CPUTime","HostCPUFeatures","AlgebraicSolving"];

using Pkg
Pkg.gc()
for p in list_pkg
    Pkg.add(p);
end
Pkg.gc()
for p in list_pkg
    Pkg.update(p);
end

versioninfo()
using Pkg
Pkg.status()
