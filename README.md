# Rational Univariate Representation in (almost) pure Julia or Maple

We propose two prototypes that implement the computation of Rational Univariate Representations for solving zero-dimensional systems of polynomial equations with rational coefficients.

Given a multivariate system S of polynomial equations depending on n variables (x1,..,xn), a Rational Univariate Representation is a linear form t=a1x1+..+anxn and a set of n+1 univariate polynomials f(T),f1(T), ... ,fn(T), such that {(x1,..,xn) s.t. f(u)=0,xi=fi(u)/f0(u)}, where f0 is the derivative of the squarefreepar of f, is in bijection with the roots of S, its inverse , from the roots of S to those of f, being defined by the linear form t.

The present source code is basically proposed to reproduce the results of an [article](https://arxiv.org/abs/2402.07141).

If you use this result, the right way to cite it is currently :

````
@misc{demin2024reading,
      title={Reading Rational Univariate Representations on lexicographic Groebner bases}, 
      author={Alexander Demin and Fabrice Rouillier and Joao Ruiz},
      year={2024},
      eprint={2402.07141},
      archivePrefix={arXiv},
      primaryClass={cs.SC}
}
````

## Julia 

**Warning** : if you want to get better performances with our Julia package, please aviod running several large examples simultanously on the same machine in order to avoid cache effects that might sensibly affect the computation times.

Our package is based on the Gröbner engine [Groebner.jl](https://github.com/sumiya11/Groebner.jl), you might need to update it if it is already installed on your system.

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
Note that the polynomials of the RUR are given as vectors of coefficients to allow conversions in your favorite package.

## Maple

The following should work with any Maple version since Maple 2018

Note that the output is a little bit more sophisticated than in Julia and is of the form f=0,[xi=fi/f0]

```
read("Maple/zds.mpl");
include("Data/Systems/caprasse.mpl")

ext,coo:=zds:-rur(sys,vars)
```

You can also use the "isolate" function that will isolate the real roots of the system directly.

**IMPORTANT** : for historical reasons, the default `RootFinding[Isolate]` package from Maple is using a small limit (8096 bits) for the hybrid (fast) computations when using the 20-years old algorithm 'RS'. To get better results avoiding this switch, please  set this limit to 2^30 using 
```
fgbrs:-rs_set_hybrid_limit(2^30);
```
and call it with a squarefree polynomial (the internal computation of the squafreepart is old and slow).
