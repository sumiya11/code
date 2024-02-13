# Rational Univariate Representation in (almost) pure Julia or Maple

We propose two prototypes that implement the computation of Rational Univariate Representations for solving zero-dimensional systems of polynomial equations with rational coefficients.

Given a multivariate system S of polynomial equations depending on n variables (x1,..,xn), a Rational Univariate Representation is a linear form t=a1x1+..+anxn and a set of n+2 univariate polynomials f(T),f0(T),f1(T), ... ,fn(T), such that {(x1,..,xn) s.t. f(u)=0,xi=fi(u)/f0(u)} is in bijection with the roots of S, its inverse , from the roots of S to those of f, being defined by the linear form t.

The present source code is basically proposed to reproduce the results of an [article](Article/RUR.pdf).

## Julia 

Warning : if you want to get better performances with our Julia package, please aviod running several large examples simultanously on the same machine in order to avoid cache effects that might sensibly affect the computation times.

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
