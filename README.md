# Rational Univariate Representation in (almost) pure Julia and Maple

We propose two prototypes that implement the computation of Rational Univariate Representation for solving zero-dimensional systems of polynomial equations with rational coefficients.

This code is basically proposed to reproduce the results of an article.

## Julia 

 ```
 import Groebner
include("Julia/rur.jl")

QQ=AbstractAlgebra.QQ
polynomial_ring=AbstractAlgebra.polynomial_ring
include("Data/Systems/caprasse.jl")

qq=zdim_parameterization(sys);
 ```

## Maple

 ```
read("Maple/zds.mpl");
include("Data/Systems/caprasse.mpl")

ext,coo:=zds:-rur(sys,vars)
 ```
