using Distributed, InteractiveUtils
versioninfo()

### Manually setup processes (only local processes from this machine) ###

# 3 local workers
addprocs(3)
@info procs()  # 1+3 in total

# can also do "@everywhere using Pkg; Pkg.develop(path=(@__DIR__)*"/../Julia/RationalUnivariateRepresentation.jl/")"
@everywhere include((@__DIR__)*"/../Julia/RationalUnivariateRepresentation.jl/src/RationalUnivariateRepresentation.jl")

### Compute ###

using AbstractAlgebra
include((@__DIR__)*"/../Data/Systems/caprasse.jl")

# Run only on this process
@time rur1 = RationalUnivariateRepresentation.zdim_parameterization(sys, parallelism=:multiprocessing, procs=[1])

# Run on processes [1,2,3,4]
@time rur2 = RationalUnivariateRepresentation.zdim_parameterization(sys, parallelism=:multiprocessing, procs=[1,2,3,4])

# Run on all available processes
@time rur3 = RationalUnivariateRepresentation.zdim_parameterization(sys, parallelism=:multiprocessing, procs=procs(), verbose=false)

rur4 = RationalUnivariateRepresentation.zdim_parameterization(sys, parallelism=:serial)

@assert rur1 == rur2 == rur3 == rur4
