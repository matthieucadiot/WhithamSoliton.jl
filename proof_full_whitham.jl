using  JLD2, RadiiPolynomial, IntervalArithmetic, LinearAlgebra 


include("list_of_functions.jl")


######################## MAIN CODE ############################################################
#############################################################################################################################

setprecision(100)


################### Loading the candidate ##############################


### Candidate for T = 0, c = 1.1, N = 800 and d = 50
#   W = load("W_0_11_50.jld2","W")
#   N = 800 ; N0 = 800
#   dn= 50 ; db = interval(big(dn)); d = interval(dn)
# c = interval(1.1)
# Tn = 0; T= abs(interval(Tn)) ; Tb = abs(interval(big(Tn))) 
# #### candidate values for a, σ0 and σ1. These values were obtained studying the graph of the function |mT-c| 
# #### values for T=0 and c=1.1
# a = 1;  ab = interval(big(a)) ; a = interval(a)
# σ0 = interval(0.03); σ1 = interval(0.1);




# ### Candidate for T = 0.5, c = 0.8, N= 800 and d = 40
#     W = load("W_05_08_40.jld2","W")
#     N = 800 ; N0 = 800
#   dn= 40 ; db = interval(big(dn)); d = interval(dn)
# c = interval(0.8)
# Tn = 0.5; T= abs(interval(Tn)) ; Tb = abs(interval(big(Tn))) 
# ### candidate values for a, σ0 and σ1. These values were obtained studying the graph of the function |mT-c| 
# ####values for T=0.5 and c=0.8
# a = 1;  ab = interval(big(a)) ; a = interval(a)
# σ1 = interval(0.1)
# σ0 = interval(0.08)





W = load("W_3_08_100_500.jld2","W")
N = 500 ; N0 = 500
dn= 100 ; db = interval(big(dn)); d = interval(dn)
c = interval(0.8)
Tn = 3; T= abs(interval(Tn)) ; Tb = abs(interval(big(Tn))) 
### candidate values for a, σ0 and σ1. These values were obtained studying the graph of the function |mT-c| 
####values for T=0.5 and c=0.8
a = 0.5;  ab = interval(big(a)) ; a = interval(a)
σ1 = interval(0.1) ; σ0 = interval(0.2)
xmin = interval(0.1) ; ymin = interval(0.05) ; ν = 1/sqrt(T) ; νb = 1/sqrt(Tb) 





# W = load("W_025_08_50.jld2","W")
# N = 1000 ; N0 = 1000
# dn= 50 ; db = interval(big(dn)); d = interval(dn)
# c = interval(0.8)
# Tn = 0.25; T= interval(Tn) ; Tb = interval(big(Tn))
# ### candidate values for a, σ0 and σ1. These values were obtained studying the graph of the function |mT-c| 
# ####values for T=0.5 and c=0.8
# a = 0.5;  ab = interval(big(a)) ; a = interval(a)
# σ1 = interval(0.1) ; σ0 = interval(0.1)
# xmin = interval(0.1) ; ymin = interval(0.05) ; ν = 1/sqrt(T) ; νb = 1/sqrt(Tb) 


# verify that the chosen constants verify their related conditions
P = check_constants(a,σ0,σ1,xmin,ymin,ν)

if P==1
  display("the constants a, σ0 and σ1 satisfy the required conditions")
else
  display("the constants a, σ0 and σ1 do not satisfy the required conditions")
  Breakdown = BreakNow
end


## we construct multiple operators and objects which are necessary along the proof
fourier = CosFourier(N,π/d); fourier0 = fourier
ν = interval(4)/interval(π)^2
νb = interval(big(4))/interval(big(π))^2

fourierE = CosFourier(2*N, interval(π)/d)
  # Fourier coefficients of cosh(2ax)
Eb = Sequence(fourierE ,interval.(big.(zeros(2*N+1))))
for n = 1:2N  
    Eb[n] = real(interval(big((-1)^n))*( (interval(big(1))-exp(-interval(big(4))*ab*db))/(interval(big(2))*ab-im*n*π/db) ))*(interval(big(1))-exp(-interval(big(4))*ab*db))
end
Eb[0] = interval(big(1))/(interval(big(2))*ab)*(interval(big(1))-exp(-interval(big(4))*ab*db)); Eb = Eb/(interval(big(2))*db) ;

E = Sequence(fourierE ,interval.((zeros(2*N+1))))
for n = 1:2N  
    E[n] = real(interval(((-1)^n))*( (interval((1))-exp(-interval((4))*a*d))/(interval((2))*a-im*n*π/d) ))*(interval((1))-exp(-interval((4))*a*d))
end
E[0] = interval((1))/(interval((2))*a)*(interval((1))-exp(-interval((4))*a*d)); E = E/(interval(2)*d) ;


#################################################

# Construction of the operator Lν
Id = UniformScaling(interval(1))
D² = project(Derivative(2), fourier0, fourier0,Interval{Float64})
Lν = diag(coefficients(Id - ν*D²))

# Construction of the operator L
d2 = diag(coefficients(D²))  ;   dd = sqrt.(-d2) ; d0=dd;  d0[1] = interval(1);
dL = tanh.(dd).*(ones(N0+1)+T*dd.^2)./d0 ;  dL[1] = interval(1) ; dL = sqrt.(dL)
MT = LinearOperator(fourier0,fourier0,Diagonal(dL))

D1 = convert(Vector{Interval{Float64}},interval.(exp2cos(N0)))
D2 = interval.(ones((N0+1)))./D1

D² = project(Derivative(2), CosFourier(N, π/db), CosFourier(N, π/db),Interval{BigFloat})
Lνb = I - νb*D²
Lνb = diag(coefficients(Lνb)).*2

# computer-assisted proof of a solitary wave
rmin, Yu, norm_DF_inv = proof_soliton(W,N,d,c,T,a,σ0,σ1)



# ###########################        PROOF OF STABILITY       ##################################################################

# construction of a projection on the trace zero functions
setprecision(100)
D = interval.(big.((1:N+1).^4))
S = trace(N); C = S' ; S =  interval.(big.(S)) ;  C =  interval.(big.(C))
        
V = interval.(big.(coefficients(W)))
W0 = (D.*C)*solve_linear(S*(D.*C),vec(S*V),100)   # the function solve_linear solves a linear system rigorously
U0b = Sequence(fourier,interval(big.(coefficients(W))) - vec(W0))
U0 = interval.(Float64.(inf.(U0b),RoundDown),Float64.(sup.(U0b),RoundUp) )

# Passage from cosine series to exponential series
fourier_f = Fourier(N,π/d)
U0f = coefficients(U0) ; U0f = [reverse(U0f[2:N+1]); U0f[1] ; U0f[2:N+1]] ; U0f = Sequence(fourier_f,U0f)
U0fb = coefficients(U0b) ; U0fb = [reverse(U0fb[2:N+1]); U0fb[1] ; U0fb[2:N+1]] ; U0fb = Sequence(fourier_f,U0fb)


# Construction of the operator L
D² = real.(project(Derivative(2), fourier_f, fourier_f,Complex{Interval{Float64}}))
dd = diag(coefficients(D²))  ;   dd = sqrt.(-dd) ; d0=dd;  d0[N+1] = interval(1);
dL = tanh.(dd).*(ones(2*N+1)+T*dd.^2)./d0 ;  dL[N+1] = interval(1) ; dL = sqrt.(dL)
MT = LinearOperator(fourier_f,fourier_f,Diagonal(dL))

if inf(T)>0
    L = MT - c*I
    Linv = LinearOperator(fourier_f,fourier_f,Diagonal(ones(2*N+1)./diag(coefficients(L))))
    DG = interval(2)*real.(project(Multiplication(U0f),fourier_f, fourier_f,Complex{Interval{Float64}}))
else
  L = c*I - MT
  Linv = LinearOperator(fourier_f,fourier_f,Diagonal(ones(2*N+1)./diag(coefficients(L))))
  DG = -interval(2)*real.(project(Multiplication(U0f),fourier_f, fourier_f,Complex{Interval{Float64}}))
end

E = eigen(Hermitian(coefficients(mid.(L) + mid.(DG))),1:2)

pspace = ParameterSpace()×fourier_f

# Construction of the eigencouples
λ1 = E.values[1] ; λ2 = E.values[2] # Theoretically we expect the second eigenvalue to be exactly zero with eigenvector u' 
V1 = E.vectors[:,1] ; V2 = E.vectors[:,2]
U1 = Sequence(pspace,vec([λ1;V1])) ; U2 = Sequence(pspace,vec([λ2;V2]))


# Proof of the eigencouples 
include("C:/Users/Matthieu/Desktop/code/julia/Whitham/list_of_functions.jl")
P1 = proof_eigen(U1,U0f,c,T,ν,d,rmin)
P2 = proof_eigen(U2,U0f,c,T,ν,d,rmin)

P = no_eigen_outside(P1[2],P2[2],U1,U2,U0f,rmin,c,T,ν)

#verify that 0 is the second eigenvalue (we know that zero is an eigenvalue, we check that it is the second one)

if P==1
  if (sup(component(real(U2),1)[1]-P2[2])<0)&&(inf(component(real(U2),1)[1]+P2[2])>0)
    display("zero is the second eigenvalue")
    display("There is no eigenvalue in between 0 and λi")
    display("Moreover, there is no eigenvalue smaller than λi")
  else
    display("there exists a second negative eigenvalue")
  end
else
  display("The inverse becomes singular for some value of λ")
end


# we verify that the Vakhitov-Kolokolov quantity is negative
P0, τ = VK_quantity(U0,rmin,σ0, norm_DF_inv)


if P0==1
  display("The Vakhitov-Kolokolov quantity is negative")
  display("The solution is stable")
else
  display("We could not prove that the Vakhitov-Kolokolov quantity is negative")
end
