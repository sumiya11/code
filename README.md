# Rational Univariate Representation in (almost) pure Julia or Maple

We propose two prototypes that implement the computation of Rational Univariate Representations for solving zero-dimensional systems of polynomial equations with rational coefficients.

This code is basically proposed to reproduce the results of an [article](Article/RUR.pdf).


If you want to get better performances with our Julia package, please aviod running several large examples simultanously on the same machine in order to avoid cache effects that might sensibly affect the computation times.

## Julia 

Install the package RationalUnivariateRepresentation.jl with

```julia
using Pkg; Pkg.develop(path="Julia/RationalUnivariateRepresentation.jl/")
```

<!-- Please first get a fresh version of the [development version of Groebner.jl package](https://github.com/sumiya11/Groebner.jl) -->

<!-- **Warning : the package does not pre-compile anything so that the first run might be slow.** -->

Then you should be able to do

```julia
using RationalUnivariateRepresentation

# Create a zero-dimensional system
include("Data/Systems/caprasse.jl")

# Find a RUR of solutions
rur = zdim_parameterization(sys)
```

## Maple
Should work with any Maple version since Maple 2018


```
read("Maple/zds.mpl");
include("Data/Systems/caprasse.mpl")

ext,coo:=zds:-rur(sys,vars)
```
