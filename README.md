

# Computer-assisted proofs of solitary waves, branches of solutions, eigencouples and spectral stability in the Whitham and capillary-gravity Whitham equations.



Table of contents:


* [Introduction](#introduction)
* [The capillary-gravity Whitham equation](#the-capillary-gravity-whitham-equation)
   * [Proof of existence of solitary waves and stability](#Proof-of-existence-of-solitary-waves-and-stability)
       * [Proof of the first 2 eigencouples](#Proof-of-the-first-2-eigencouples)
       * [Proof of spectral stability](#Proof-of-spectral-stability)
   * [Proof of a branch of solutions](#Proof-of-a-branch-of-solutions)
* [Utilisation and References](#utilisation-and-references)
* [License and Citation](#license-and-citation)
* [Contact](#contact)



# Introduction

This Julia code is a complement to the article 

#### [[1]](https://arxiv.org/abs/2403.18718) : "Constructive proofs of existence and stability of solitary waves in the Whitham and capillary-gravity Whitham equations", Matthieu Cadiot [ArXiv Link](https://arxiv.org/abs/2403.18718)

as it provides the necessary rigorous computations along the paper. Such computations are performed using the package [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl). The mathematical objects (spaces, sequences, operators,...) are built using the package [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl). 


# The Whitham and capillary-gravity Whitham equations

The capillary-gravity Whitham equtation is a nonlocal equation modelling surface water waves. The equation possesses travelling waves solutions, which can be studied throughout the following equation
$$\mathbb{M}_Tu - cu + u^2 = 0,$$
where $c$ is the is the velocity of the wave and $\mathbb{M}_T$ is a Fourier multiplier operator given by its symbol
$$\mathcal{F}(\mathbb{M}_Tu)(\xi)  = \sqrt{\frac{\tanh(2\pi\xi)(1+T(2\pi\xi)^2)}{2\pi\xi}}\hat{u}(\xi).$$
$T \geq 0$ is the Bond number accounting for the surface tension. In particular, if $T=0$, then the equation is purely gravitational and is often refered as the so-called Whitham equation in the literature.

Both the case $T>0$ and $T=0$ are known to possess solitary waves, that is travelling waves on the real line vanishing at infinity. These are the central topic of study of this code.

## Proof of existence of solitary waves and stability

The code "proof_full_whitham.jl" provides the computer-assisted component for the constructive proofs of existence of solitary waves using the analysis of [[1]](https://arxiv.org/abs/2403.18718). Four approximate solutions  (cf. the .jld2 files) have already been computed and the main code provides the algorithmic details for a proof of existence in a neighborhood of such approximate solutions (i.e. we prove Theorems 4.7 and 4.8). The file W_025_08_50.jld2 corresponds to an approximate solution for $T=0.25$ and $c=0.8$, W_05_08_40.jld2 corresponds to $T=0.5$ and $c=0.8$,  W_0_11_50.jld2 corresponds to  $T=0$ and $c=1.1$ and W_3_08_100_500 to $T=3$ and $c=0.8$.

The code will compute rigorously the needed bounds of Section 3 of [[1]](https://arxiv.org/abs/2403.18718) and validate (or not) the computer-assisted proof. If the computer-assisted proof succeeds, the radius for the smallest and biggest ball of contraction are displayed, as well as the bounds of Theorem 3.1.

### Proof of the first 2 eigencouples

If the proof of the soliton is achieved, the code will then compute approximations for the first 2 eigencouples of the linearization around the proved soliton. Then, the needed bounds for the proof of eigencouples are computed, following the analysis of Section 5 of [[1]](https://arxiv.org/abs/2403.18718) (cf. Theorem 5.6). If the proof of the eigencouple succeeds, we obtain that the corresponding eigenvalue is simple and the interval in which it is unique is displayed.
 
 
 ### Proof of spectral stability

If the proof of the first 2 eigencouples is achieved, then the code will try to verify that these 2 eigencouples indeed correspond to the first two eigenvalues $\lambda_1, \lambda_2$ of the linearization. To do so, we use the analysis of Section 5.2 and prove non-existence of negative eigenvalues away from $\lambda_1$ and $\lambda_2$. In particular, we verify that zero is in the interval of uniqueness of $\lambda_2$, proving that $\lambda_2=0$ (cf. Remark 5.1). Moreover, we verify that $\lambda_1$ is negative. Finally, we verify that the Vakhitov-Kolokolov quantity introduced in Lemma 5.1 is negative. If that is the case, using Lemma 5.1 combined with Theorem 5.3, we obtain the spectral stability of the solitary wave.



## Proof of a branch of solutions


The code "proof_continuation.jl" performs the verification of an existence proof for a branch of solitary waves in the case $T=0$. Following Section 4.6, we look for a branch of solutions parameterized by the velocity $c$. We obtain an approximated branch using a Chebyshev expansion in $c$ and compute the bounds introduced in Theorem 4.9. The proof is achieved in five steps : one first branch on $[1.05, 1.17]$, then $[1.17,1.18]$, $[1.18,1.19]$, $[1.19,1.20$ and finally $[1.20,1.21]$. For the first part, the solutions have a small amplitude and a large range of values for $c$ can be covered. Then, as $c$ gets bigger, the amplitude increases and the solution starts peaking. When reaching the range $[1.20,1.21]$ the number of Fourier coefficients has to be increased to 2000 in order to accurately approximate the cusping behavior.   


 
 # Utilisation and References

 The interested user can directly go to line 630 of the main code in order to use the code. The code is general and provides rigorous computations of the required bounds for any values of $c$ and $T$ satisfying Assumption 1 in [[1]](https://arxiv.org/abs/2403.18718). Illustrations are provided for both files W_0_11_50.jld2 and W_05_08_40.jld2.
  The code is built using the following packages :
 - [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl) 
 - [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl)
 - [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
 - [JLD2](https://github.com/JuliaIO/JLD2.jl)
 - [PrettyTables](https://ronisbr.github.io/PrettyTables.jl/stable/).
 
 
 # License and Citation
 
  This code is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).
  
If you wish to use this code in your publication, research, teaching, or other activities, please cite it using the following BibTeX template:

```
@software{WhithamSoliton.jl,
  author = {Matthieu Cadiot},
  title  = {WhithamSoliton.jl},
  url    = {https://github.com/matthieucadiot/WhithamSoliton.jl},
  note = {\url{ https://github.com/matthieucadiot/WhithamSoliton.jl},
  year   = {2024}
}
```

# Contact

You can contact me at :

matthieu.cadiot@mail.mcgill.ca
