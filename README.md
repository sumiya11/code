# Rational Univariate Representation in (almost) pure Julia or Maple

We propose two prototypes that implement the computation of Rational Univariate Representations for solving zero-dimensional systems of polynomial equations with rational coefficients.

This code is basically proposed to reproduce the results of an [article](Article/RUR.pdf).

## Julia 

Please first get a fresh version of the [development version of Groebner.jl package](https://github.com/sumiya11/Groebner.jl)

**Warning : the package does not pre-compile anything so that the first run might be slow.**

 ```
include("../Groebner.jl/src/Groebner.jl")
include("Julia/rur.jl")

QQ=AbstractAlgebra.QQ
polynomial_ring=AbstractAlgebra.polynomial_ring
include("Data/Systems/caprasse.jl")

qq=zdim_parameterization(sys);
 ```

## Maple
Should work with any Maple version since Maple 2018


 ```
read("Maple/zds.mpl");
include("Data/Systems/caprasse.mpl")

ext,coo:=zds:-rur(sys,vars)
 ```
