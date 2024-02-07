# Rational Univariate Representation in (almost) pure Julia and Maple

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
