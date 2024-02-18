using Distributed, InteractiveUtils
versioninfo()

### Manually setup processes (only local processes from this machine) ###

# 3 local workers
addprocs(3)
@info procs()  # 1+3 in total

# can also do "@everywhere using Pkg; Pkg.develop(path=(@__DIR__)*"/../Julia/RationalUnivariateRepresentation.jl/")"
@everywhere include("/home/demin/code/Julia/RationalUnivariateRepresentation.jl/src/RationalUnivariateRepresentation.jl")

### Compute ###

using AbstractAlgebra
include("/home/demin/code/Data/Systems/katsura11.jl")

# Run only on this process
@time rur1 = RationalUnivariateRepresentation.zdim_parameterization(sys, parallelism=:multiprocessing, procs=[1])

# Run on processes [1,2,3,4]
@time rur2 = RationalUnivariateRepresentation.zdim_parameterization(sys, parallelism=:multiprocessing, procs=[1,2,3,4])

# Run on all available processes
@time rur3 = RationalUnivariateRepresentation.zdim_parameterization(sys, parallelism=:multiprocessing, procs=procs(), verbose=false);

@time rur4 = RationalUnivariateRepresentation.zdim_parameterization(sys, parallelism=:serial, verbose=false);

@assert rur1 == rur2 == rur3 == rur4
