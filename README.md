# Rational Univariate Representation in (almost) pure Julia or Maple

We propose two prototypes that implement the computation of Rational Univariate Representations for solving zero-dimensional systems of polynomial equations with rational coefficients.

Given a multivariate system $S$ of polynomial equations depending on $n$ variables $(x_1,..,x_n)$, a Rational Univariate Representation is a linear form $t=a_1x_1+..+a_nx_n$ and a set of $n+1$ univariate polynomials $\{f(T),f_1(T), ... ,f_n(T)\}$, such that $\{(x_1,..,x_n)\; s.t.\; f(u)=0,x_i=\frac{f_i(u)}{f_0(u)}\}$, where $f_0$ is the derivative of the squarefreepart of $f$, is in bijection with the roots of $S$, its inverse , from the roots of $S$ to those of $f$, being defined by the linear form $t$.

The present source codes are basically proposed to reproduce the results of an [article](https://arxiv.org/abs/2402.07141).

If you use this result, the right way to cite it is currently :

````
@misc{demin2025readingrationalunivariaterepresentations,
      title={Reading Rational Univariate Representations on lexicographic Groebner bases}, 
      author={Alexander Demin and Fabrice Rouillier and Joao Ruiz},
      year={2025},
      eprint={2402.07141},
      archivePrefix={arXiv},
      primaryClass={cs.SC},
      url={https://arxiv.org/abs/2402.07141}, 
}
````

## Julia 

See https://newrur.gitlabpages.inria.fr/RationalUnivariateRepresentation.jl/ for installation instructions, source code and examples of use.

The package can be combined with RS.jl (see https://pace.gitlabpages.inria.fr/rs.jl/) which exports 
some few functions of RS C-Library used by the Maple function Isolate (see section Maple below) or can be access 
from PACE.jl (https://pace.gitlabpages.inria.fr/pace.jl/) in a general framework for solving systems of algebraic equations.

## Maple

The following should work with any Maple version since Maple 2018

Note that the output is a little bit more sophisticated than in Julia and is of the form $f=0,[x_i=f_i/f_0]$

```
read("Maple/zds.mpl");
include("Data/Systems/caprasse.mpl")

ext,coo:=zds:-rur(sys)
```

You can also use the "isolate" function that will isolate the real roots of the system directly.
```
iso:=zds:-isolate(sys,vars):
```
If you want to get the same output as with `RootFinding[Isolate]` then you might apply the following transform 
```
iso:=map(v->map(u->lhs(u)=evalf((rhs(u)[1]+rhs(u)[2])/2,Digits),[op(v)]),iso);
```

**IMPORTANT for fans of Benchmarks** : for those who like performing benchmarks to try to illustrate theoretical informations, it must be precised that for historical reasons, mainly linked to the fact that the algorithm is able to work with approximate coefficients :  
- the default `RootFinding[Isolate]` package from Maple, when run with the 20-year old algorithm 'RS', is using a small limit (8096 bits) for the hybrid  computations. On large examples, you will get better results by avoiding the switch to the slow algorithm, by modifying this limit to a greater value, either by setting the environment variable USPMAXPREC before launching Maple (the saftest way) or using the function `fgbrs:-rs_set_hybrid_limit(<any integer < 2^31>);`.
- call the isolate function with a squarefree polynomial (the internal computation of the squafreepart is old and slow, the one of Maple is way faster), otherwise you will essentially measure the computation of the gcd.

