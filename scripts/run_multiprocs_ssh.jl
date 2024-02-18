using Distributed, InteractiveUtils
versioninfo()

### Manually setup processes ###

# 3 local workers
addprocs(3)

# 8 workers on l'X server via ssh :^).
# Important: the version of Julia on the server is 1.9.2, same as on the host machine.
addprocs([("demin@ron.medicis.polytechnique.fr", 8)], dir="/home/demin", exename="/home/demin/downloads/julia-1.9.2/bin/julia")

@info procs()  # 1+3+8 in total

@everywhere include((@__DIR__)*"/../Julia/RationalUnivariateRepresentation.jl/src/RationalUnivariateRepresentation.jl")

### Compute ###

using AbstractAlgebra
include((@__DIR__)*"/../Data/Systems/caprasse.jl")

# Run only on this process
RationalUnivariateRepresentation.zdim_parameterization(sys, parallelism=:multiprocessing, procs=[1])

# Run only on processes [1,2,3,4]
RationalUnivariateRepresentation.zdim_parameterization(sys, parallelism=:multiprocessing, procs=[1,2,3,4])
