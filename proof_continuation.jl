using MAT, JLD2, RadiiPolynomial, IntervalArithmetic, LinearAlgebra 

include("list_of_functions.jl")


####################################################################################################
setprecision(128)    ; Id = UniformScaling(interval(1))
N = 1000
d= 50          # size of the domain = half period of the functions
T = 0  # Bond number
c0 = 1.05
fourierC = CosFourier(N, π/d)   ; fourierC_ = CosFourier(N, interval(π)/interval(d))
a0b = interval(big(1.5)) ; db = interval(big(d))

 W = load("W_0_105_40_500.jld2","W")       # we load the approximate solution obtained at c=1.1
 W = Sequence(fourierC, coefficients(project(W,CosFourier(N, frequency(W)))))

 U0 = newton_whitham(W, c0, T, fourierC)    # we compute a better approximated solution 


 # Now we compute an approximate branch

 Nc = 8   # Chebyshev order
 c1 = 1.07
 c = Sequence(CosFourier(0, π/d)⊗ Chebyshev(1), [0.5*(c0 + c1), 0.25*(-c0 + c1)])  # function c mapsto c in Chebyshev polynomials 
 
 W = Sequence(fourierC ⊗ Chebyshev(0), coefficients(U0))
 W = project(W, fourierC ⊗ Chebyshev(Nc))   # Initial guess for newton
 
  W, _ = newton(W; maxiter = 10) do W
     return F_W(W, c,T), DF_W(W, c,T)  # Newton method to compute an approximate branch with a Cheb expansion in c
  end
 W_first = W
 c0 = c1
 
 
 Nc =  8  # Chebyshev order
 c1 = 1.11
 c = Sequence(CosFourier(0, π/d)⊗ Chebyshev(1), [0.5*(c0 + c1), 0.25*(-c0 + c1)])   # function c mapsto c in Chebyshev polynomials 
 
 W = Sequence(fourierC ⊗ Chebyshev(0), coefficients(W(nothing,1)))
 W = project(W, fourierC ⊗ Chebyshev(Nc))   # Initial guess for newton
 
  W, _ = newton(W; maxiter = 10) do W
     return F_W(W, c,T), DF_W(W, c,T)  # Newton method to compute an approximate branch with a Cheb expansion in c
  end
 W_second = W
 c0 = c1
 
 
  Nc = 8   # Chebyshev order
 c1 = 1.14
 c = Sequence(CosFourier(0, π/d)⊗ Chebyshev(1), [0.5*(c0 + c1), 0.25*(-c0 + c1)])   # function c mapsto c in Chebyshev polynomials 
 
 W = Sequence(fourierC ⊗ Chebyshev(0), coefficients(W(nothing,1)))
 W = project(W, fourierC ⊗ Chebyshev(Nc))   # Initial guess for newton
 
  W, _ = newton(W; maxiter = 10) do W
     return F_W(W, c,T), DF_W(W, c,T)  # Newton method to compute an approximate branch with a Cheb expansion in c
  end
 W_third = W
 c0 = c1
 
 # we do the same thing a second time
  
 c1 = 1.17
 c = Sequence(CosFourier(0, π/d)⊗ Chebyshev(1), [0.5*(c0 + c1), 0.25*(-c0 + c1)]) 
 
 W = Sequence(fourierC ⊗ Chebyshev(0), coefficients(W(nothing,1)))
 W = project(W, fourierC ⊗ Chebyshev(Nc))
  W, _ = newton(W; maxiter = 10) do W
     return F_W(W, c,T), DF_W(W, c,T)
  end
 W_fourth = W
 
 

 ###################################################################################################################################################
 
 c0 = 1.05 ; c2 = 1.07 ; c3 = 1.11 ; c4 = 1.14
 c1 = 1.17
 Nc = 10
 N_ = 3Nc
 N_fft_ = nextpow(2, 2N_ + 1)
 npts_ = N_fft_ ÷ 2 + 1
 
 c = [0.5 * (c0 + c1) + 0.5cospi(2k / N_fft_)*(c1 - c0) for k ∈ 0:npts_-1]
 V = Vector{Sequence{CosFourier,Vector{Float64}}}(undef, npts_)
 
 for n=1:npts_
     if c[n]<= c2
         V[n] = Sequence(fourierC_,W_first(nothing,2/(c2-c0)*c[n] + 1-2*c2/(c2-c0))[(:,0)])
     elseif (c[n]<= c3)&&(c[n]>c2)
        V[n] = Sequence(fourierC_,W_second(nothing,2/(c3-c2)*c[n] + 1-2*c3/(c3-c2))[(:,0)])
     elseif (c[n]<= c4)&&(c[n]>c3)
        V[n] = Sequence(fourierC_,W_third(nothing,2/(c4-c3)*c[n] + 1-2*c4/(c4-c3))[(:,0)])  
     else
        V[n] = Sequence(fourierC_,W_fourth(nothing,2/(c1-c4)*c[n] + 1-2*c1/(c1-c4))[(:,0)])         
     end
 end
 
 # we verified that a, σ0 and σ1 satisfy (39) in Proposition 4.1 using the function check_constants. 
 a = interval(1) ; σ0 = interval(0.0228) ; σ1 = interval(0.1) ; d= interval(d)
  ν_ = interval(4)/interval(π)^2 ; ν=ν_
 
 
 D² = project(Derivative(2), fourierC_, fourierC_,Interval{Float64})
 d2 = sqrt.(-diag(coefficients(D²))) ; d2[1] = interval(1) ; d2 = (tanh.(d2)).*(interval(1).+interval(T)*d2.^2)./d2 ; d2[1] = interval(1) ; d2 =  sqrt.(d2)
 M0 = LinearOperator(fourierC_,fourierC_,Diagonal(d2)) 
 dL = diag(Matrix(coefficients(M0)))
 
 Lν = diag(coefficients(Id - ν_*D²))
 D1 = convert(Vector{Interval{Float64}},interval.(exp2cos(N))) ; D2 = interval.(ones((N+1)))./D1
 
 # we do a proof from c0 = 1.05 to c1 = 1.17. Here we can take a wide interval since solutions are smooth and with "small amplitude". At least the decay of the Fourier coefficients is strong.
 V0_ini, V0_end, r0_min, r0_max = proof_continuation_whitham(V,c0,c1,N,Nc,d,a,σ0,σ1)  
 
 


########### Second branch ####################
c0 = c1
c1 = 1.18
W = W(nothing,1)   # we take the end point of the previous branch
W = Sequence(fourierC, coefficients(project(W,CosFourier(N, frequency(W)))))

U0 = newton_whitham(W, c0, T, fourierC)    # we compute a better approximated solution 

Nc = 3   # Chebyshev order
c = Sequence(CosFourier(0, mid(π/d))⊗ Chebyshev(1), [0.5*(c0 + c1), 0.25*(-c0 + c1)])  # function c mapsto c in Chebyshev polynomials 

W = Sequence(fourierC ⊗ Chebyshev(0), coefficients(U0))
W = project(W, fourierC ⊗ Chebyshev(Nc))   # Initial guess for newton

 W, _ = newton(W; maxiter = 10) do W
    return F_W(W, c,T), DF_W(W, c,T)  # Newton method to compute an approximate branch with a Cheb expansion in c
 end


###################################################################################################################################################

Nc = 4
N_ = 3Nc
N_fft_ = nextpow(2, 2N_ + 1)
npts_ = N_fft_ ÷ 2 + 1

c = [0.5 * (c0 + c1) + 0.5cospi(2k / N_fft_)*(c1 - c0) for k ∈ 0:npts_-1]
V = Vector{Sequence{CosFourier,Vector{Float64}}}(undef, npts_)

for n=1:npts_
        V[n] = Sequence(fourierC_,W(nothing,2/(c1-c0)*c[n] + 1-2*c1/(c1-c0))[(:,0)])
end

# we verified that a, σ0 and σ1 satisfy (39) in Proposition 4.1 using the function check_constants. 
a = interval(1) ; σ0 = interval(0.0228) ; σ1 = interval(0.1) ; d= interval(d)
 ν_ = interval(4)/interval(π)^2 ; ν=ν_


D² = project(Derivative(2), fourierC_, fourierC_,Interval{Float64})
d2 = sqrt.(-diag(coefficients(D²))) ; d2[1] = interval(1) ; d2 = (tanh.(d2)).*(interval(1).+interval(T)*d2.^2)./d2 ; d2[1] = interval(1) ; d2 =  sqrt.(d2)
M0 = LinearOperator(fourierC_,fourierC_,Diagonal(d2)) 
dL = diag(Matrix(coefficients(M0)))

Lν = diag(coefficients(Id - ν_*D²))
D1 = convert(Vector{Interval{Float64}},interval.(exp2cos(N))) ; D2 = interval.(ones((N+1)))./D1

V1_ini, V1_end, r1_min, r1_max  = proof_continuation_whitham(V,c0,c1,N,Nc,d,a,σ0,σ1)   # proof of the branch between c1 = 1.17 and c2 = 1.18




########### Third branch ####################
c0 = c1
c1 = 1.19
W = W(nothing,1)   # we take the end point of the previous branch
N = 1200
fourierC = CosFourier(N, mid(π/d))   ; fourierC_ = CosFourier(N, interval(π)/interval(d))
W = Sequence(fourierC, coefficients(project(W,CosFourier(N, frequency(W)))))

U0 = newton_whitham(W, c0, T, fourierC)    # we compute a better approximated solution 


Nc = 3   # Chebyshev order
c = Sequence(CosFourier(0, mid(π/d))⊗ Chebyshev(1), [0.5*(c0 + c1), 0.25*(-c0 + c1)])  # function c mapsto c in Chebyshev polynomials 

W = Sequence(fourierC ⊗ Chebyshev(0), coefficients(U0))
W = project(W, fourierC ⊗ Chebyshev(Nc))   # Initial guess for newton

 W, _ = newton(W; maxiter = 10) do W
    return F_W(W, c,T), DF_W(W, c,T)  # Newton method to compute an approximate branch with a Cheb expansion in c
 end


###################################################################################################################################################

Nc = 4
N_ = 3Nc
N_fft_ = nextpow(2, 2N_ + 1)
npts_ = N_fft_ ÷ 2 + 1

c = [0.5 * (c0 + c1) + 0.5cospi(2k / N_fft_)*(c1 - c0) for k ∈ 0:npts_-1]
V = Vector{Sequence{CosFourier,Vector{Float64}}}(undef, npts_)

for n=1:npts_
        V[n] = Sequence(fourierC_,W(nothing,2/(c1-c0)*c[n] + 1-2*c1/(c1-c0))[(:,0)])
end

# we verified that a, σ0 and σ1 satisfy (39) in Proposition 4.1 using the function check_constants. 
a = interval(1) ; σ0 = interval(0.0228) ; σ1 = interval(0.1) ; d= interval(d)
 ν_ = interval(4)/interval(π)^2 ; ν=ν_


D² = project(Derivative(2), fourierC_, fourierC_,Interval{Float64})
d2 = sqrt.(-diag(coefficients(D²))) ; d2[1] = interval(1) ; d2 = (tanh.(d2)).*(interval(1).+interval(T)*d2.^2)./d2 ; d2[1] = interval(1) ; d2 =  sqrt.(d2)
M0 = LinearOperator(fourierC_,fourierC_,Diagonal(d2)) 
dL = diag(Matrix(coefficients(M0)))

Lν = diag(coefficients(Id - ν_*D²))
D1 = convert(Vector{Interval{Float64}},interval.(exp2cos(N))) ; D2 = interval.(ones((N+1)))./D1

V2_ini, V2_end, r2_min, r2_max  = proof_continuation_whitham(V,c0,c1,N,Nc,d,a,σ0,σ1)   # proof of the branch between c2 = 1.18 and c3 = 1.19





########### Fourth branch ####################
c0 = c1
c1 = 1.20
W = W(nothing,1)   # we take the end point of the previous branch
N = 1800
fourierC = CosFourier(N, mid(π/d))   ; fourierC_ = CosFourier(N, interval(π)/interval(d))
W = Sequence(fourierC, coefficients(project(W,CosFourier(N, frequency(W)))))

U0 = newton_whitham(W, c0, T, fourierC)    # we compute a better approximated solution 


Nc = 3   # Chebyshev order
c = Sequence(CosFourier(0, mid(π/d))⊗ Chebyshev(1), [0.5*(c0 + c1), 0.25*(-c0 + c1)])  # function c mapsto c in Chebyshev polynomials 

W = Sequence(fourierC ⊗ Chebyshev(0), coefficients(U0))
W = project(W, fourierC ⊗ Chebyshev(Nc))   # Initial guess for newton

 W, _ = newton(W; maxiter = 10) do W
    return F_W(W, c,T), DF_W(W, c,T)  # Newton method to compute an approximate branch with a Cheb expansion in c
 end


###################################################################################################################################################

Nc = 4
N_ = 3Nc
N_fft_ = nextpow(2, 2N_ + 1)
npts_ = N_fft_ ÷ 2 + 1

c = [0.5 * (c0 + c1) + 0.5cospi(2k / N_fft_)*(c1 - c0) for k ∈ 0:npts_-1]
V = Vector{Sequence{CosFourier,Vector{Float64}}}(undef, npts_)

for n=1:npts_
        V[n] = Sequence(fourierC_,W(nothing,2/(c1-c0)*c[n] + 1-2*c1/(c1-c0))[(:,0)])
end

# we verified that a, σ0 and σ1 satisfy (39) in Proposition 4.1 using the function check_constants. 
a = interval(1) ; σ0 = interval(0.0228) ; σ1 = interval(0.1) ; d= interval(d)
 ν_ = interval(4)/interval(π)^2 ; ν=ν_


D² = project(Derivative(2), fourierC_, fourierC_,Interval{Float64})
d2 = sqrt.(-diag(coefficients(D²))) ; d2[1] = interval(1) ; d2 = (tanh.(d2)).*(interval(1).+interval(T)*d2.^2)./d2 ; d2[1] = interval(1) ; d2 =  sqrt.(d2)
M0 = LinearOperator(fourierC_,fourierC_,Diagonal(d2)) 
dL = diag(Matrix(coefficients(M0)))

Lν = diag(coefficients(Id - ν_*D²))
D1 = convert(Vector{Interval{Float64}},interval.(exp2cos(N))) ; D2 = interval.(ones((N+1)))./D1

V3_ini, V3_end, r3_min, r3_max  = proof_continuation_whitham(V,c0,c1,N,Nc,d,a,σ0,σ1)   # proof of the branch between c3 = 1.19 and c4 = 1.20






########### Fifth branch ####################
c0 = c1
c1 = 1.21
W = W(nothing,1)   # we take the end point of the previous branch
N = 2000
fourierC = CosFourier(N, mid(π/d))   ; fourierC_ = CosFourier(N, interval(π)/interval(d))
W = Sequence(fourierC, coefficients(project(W,CosFourier(N, frequency(W)))))

U0 = newton_whitham(W, c0, T, fourierC)    # we compute a better approximated solution 


Nc = 3   # Chebyshev order
c = Sequence(CosFourier(0, mid(π/d))⊗ Chebyshev(1), [0.5*(c0 + c1), 0.25*(-c0 + c1)])  # function c mapsto c in Chebyshev polynomials 

W = Sequence(fourierC ⊗ Chebyshev(0), coefficients(U0))
W = project(W, fourierC ⊗ Chebyshev(Nc))   # Initial guess for newton

 W, _ = newton(W; maxiter = 10) do W
    return F_W(W, c,T), DF_W(W, c,T)  # Newton method to compute an approximate branch with a Cheb expansion in c
 end


###################################################################################################################################################

Nc = 4
N_ = 3Nc
N_fft_ = nextpow(2, 2N_ + 1)
npts_ = N_fft_ ÷ 2 + 1

c = [0.5 * (c0 + c1) + 0.5cospi(2k / N_fft_)*(c1 - c0) for k ∈ 0:npts_-1]
V = Vector{Sequence{CosFourier,Vector{Float64}}}(undef, npts_)

for n=1:npts_
        V[n] = Sequence(fourierC_,W(nothing,2/(c1-c0)*c[n] + 1-2*c1/(c1-c0))[(:,0)])
end

# we verified that a, σ0 and σ1 satisfy (39) in Proposition 4.1 using the function check_constants. 
a = interval(1) ; σ0 = interval(0.0228) ; σ1 = interval(0.1) ; d= interval(d)
 ν_ = interval(4)/interval(π)^2 ; ν=ν_


D² = project(Derivative(2), fourierC_, fourierC_,Interval{Float64})
d2 = sqrt.(-diag(coefficients(D²))) ; d2[1] = interval(1) ; d2 = (tanh.(d2)).*(interval(1).+interval(T)*d2.^2)./d2 ; d2[1] = interval(1) ; d2 =  sqrt.(d2)
M0 = LinearOperator(fourierC_,fourierC_,Diagonal(d2)) 
dL = diag(Matrix(coefficients(M0)))

Lν = diag(coefficients(Id - ν_*D²))
D1 = convert(Vector{Interval{Float64}},interval.(exp2cos(N))) ; D2 = interval.(ones((N+1)))./D1

V4_ini, V4_end, r4_min, r4_max  = proof_continuation_whitham(V,c0,c1,N,Nc,d,a,σ0,σ1)   # proof of the branch between c4 = 1.20 and c5 = 1.21


# Now we verify that the branch is continuous. In particular, we use that \|u\|_{Hcon} <= (c0+c1)/2*\|Lν u\|_2

if sup(sqrt(2d)*(c0+c1)/2*norm(Lν.*(V0_end-V1_ini),2) + (c0+c1)/2*r0_min) < inf(r1_max)
   display("continuity of branch 0 and 1") 
else 
   display("discontinuity of branch 0 and 1")
end

if sup(sqrt(2d)*(c0+c1)/2*norm(Lν.*(V1_end-V2_ini),2) + (c0+c1)/2*r1_min) < inf(r2_max)
   display("continuity of branch 1 and 2") 
else 
   display("discontinuity of branch 1 and 2")
end

if sup(sqrt(2d)*(c0+c1)/2*norm(Lν.*(V2_end-V3_ini),2) + (c0+c1)/2*r2_min) < inf(r3_max)
   display("continuity of branch 2 and 3") 
else 
   display("discontinuity of branch 2 and 3")
end

if sup(sqrt(2d)*(c0+c1)/2*norm(Lν.*(V3_end-V4_ini),2) + (c0+c1)/2*r3_min) < inf(r4_max)
   display("continuity of branch 3 and 4") 
else 
   display("discontinuity of branch 3 and 4")
end
