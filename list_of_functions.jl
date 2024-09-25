




function continuation_whitham(u1,u2,c1,c2,T,d,N,N0,σ0,σ1,a)

    fourier = CosFourier(N,π/d) 
    fourier0 = CosFourier(N0,π/d) 

    u1 = Sequence(fourier,coefficients(u1)[1:N+1]) ; u2 = Sequence(fourier,coefficients(u2)[1:N+1])
  
  c = c1

  L = Lν.*(MT - c*Id)
  Linv = LinearOperator(fourier0,fourier0,Diagonal(ones(N0+1)./diag(coefficients(L))))
  
  #Construction of the fourier series of cosh(2ax) on Ω
 
        
    V = interval.(big.(coefficients(u1)))
        
    W0 = (D.*C)*solve_linear(S*(D.*C),vec(S*V),100)   # the function solve_linear solves a linear system rigorously
    display("U0 created")
        #
    U0b = interval.(big.(u1)) - Sequence(fourier,vec(W0))
    U0 = interval.(Float64.(inf.(U0b),RoundDown),Float64.(sup.(U0b),RoundUp) )
  
  
  V2 = interval(2)*U0
  V1 = -interval(4)*ν*π/d*(0:N).*U0
  V0 = -interval(2)*ν*π^2/d^2*((0:N).^2).*U0
  
  
  
  ############ Computation of A (Section 3.2) ########################################
  
  DG = Lν.*project(Multiplication(V2),fourier0, fourier0,Interval{Float64})
  B = interval.(inv(I + mid.(DG)*mid.(Linv)))
  Badj = LinearOperator(fourier0,fourier0,coefficients(B)')
  
  # Construction of WT
  e0 = Sequence(fourier0 ,interval.(zeros(N0+1)))
  e0[0] = interval(1)
  if inf(T)>0
      WT = e0
  else
    e02 = Sequence(CosFourier(2*N0, π/d) ,interval.(zeros(2*N0+1)))
    e02[0] = interval(1)
    WT = interval.(mid.(project(Multiplication(mid.(e0 - 1/mid(c)*V2)),fourier0, CosFourier(2*N0, π/d),Interval{Float64}))\mid.(e02))
  end
  
  # computation of the norm of B
  MWT = project(Multiplication(WT*WT),fourier0, fourier0,Interval{Float64}) - project(Multiplication(WT),fourier0, fourier0,Interval{Float64})^2
  
  
  n_WT = sqrt(opnorm(LinearOperator(coefficients(D1.*(MWT).*D2')),2))   # computation of \|π_NW_Tπ^N\|_2
  n_BN = (opnorm(LinearOperator(coefficients(D1.*((B*Badj)^8).*D2')),2))^(interval(1)/interval(16))   #computation of \|B^N_T\|_2
  norm_B = maximum([interval(1) maximum([norm(WT,1) n_BN])+n_WT])  #computation of the norm of B
  
  display("norm of B")
  display(norm_B)
  
  MWT = Nothing
  
  ################# Computation of Z1 ####################################
  
  W0 = WT*V0 ;              
  W1 = WT*V1 ;  
  W2 = WT*V2 ;          
  
  if inf(T)>0
    L0 = Linv ;             l0 = (interval(1)+ν*(interval(N0+1)*π/d)^2)*(maximum(vec((abs.(dL.-c)))))
    L1 = Linv*sqrt.(-D²) ;  l1 = (interval(1)+ν*(interval(N0+1)*π/d)^2)*(maximum(vec((abs.(dL.-c)))))/(interval(N0+1)*π/d)
    L2 = Linv.*Lν' ;          l2 = maximum(vec((abs.(dL.-c))))
  else
    L0 = Linv ;             l0 = (interval(1)+ν*(interval(N0+1)*π/d)^2)*(maximum(vec((abs.(dL.-c)))))
    L1 = Linv*sqrt.(-D²) ;  l1 = (interval(1)+ν*(interval(N0+1)*π/d)^2)*(maximum(vec((abs.(dL.-c)))))/(interval(N0+1)*π/d)
    L2 = Linv.*Lν' + interval(1)/c*Id ;  l2 = maximum(vec(abs.((c*(dL.-c))./dL)))
  end
  
  #MWi corresponds to the operator πᴺ(Vi*WT)(π^{3N}-πᴺ)(Vi*WT)πᴺ = πᴺ(Vi*WT*Vi*WT)πᴺ - (πᴺ(Vi*WT)πᴺ)^2 (recall that Wi = Vi*WT)
  MW0 = project(Multiplication(W0*W0),fourier0, fourier0,Interval{Float64}) - project(Multiplication(W0),fourier0, fourier0,Interval{Float64})^2
  MW1 = project(Multiplication(W1*W1),fourier0, fourier0,Interval{Float64}) - project(Multiplication(W1),fourier0, fourier0,Interval{Float64})^2
  MW2 = project(Multiplication(W2*W2),fourier0, fourier0,Interval{Float64}) - project(Multiplication(W2),fourier0, fourier0,Interval{Float64})^2
  
  #MVi corresponds to the operator πᴺVi(π^{2N}-πᴺ)Viπᴺ = πᴺ(Vi*Vi)πᴺ - (πᴺViπᴺ)^2
  MV0 = project(Multiplication(V0*V0),fourier0, fourier0,Interval{Float64}) - project(Multiplication(V0),fourier0, fourier0,Interval{Float64})^2
  MV1 = project(Multiplication(V1*V1),fourier0, fourier0,Interval{Float64}) - project(Multiplication(V1),fourier0, fourier0,Interval{Float64})^2
  MV2 = project(Multiplication(V2*V2),fourier0, fourier0,Interval{Float64}) - project(Multiplication(V2),fourier0, fourier0,Interval{Float64})^2
  
  # computation of each component of Z1
  Z11 =  opnorm(LinearOperator(coefficients(D1.*(Id - B*(Id+DG*Linv)).*D2')),2)^2
  Z11 = Z11 + ( sqrt(opnorm(LinearOperator(coefficients(D1.*(L0*MW0*L0).*D2')),2)) + sqrt(opnorm(LinearOperator(coefficients(D1.*(L1*MW1*L1).*D2')),2)) + sqrt(opnorm(LinearOperator(coefficients(D1.*(L2*MW2*L2).*D2')),2)) )^2
  Z11 = sqrt(Z11)
  
  display("value of Z11")
  display(Z11)
  
  Z12 = sqrt(opnorm(LinearOperator(coefficients(D1.*B*MV0*Badj.*D2')),2))/l0 + sqrt(opnorm(LinearOperator(coefficients(D1.*B*MV1*Badj.*D2')),2))/l1 + sqrt(opnorm(LinearOperator(coefficients(D1.*B*MV2*Badj.*D2')),2))/l2
  
  display("value of Z12")
  display(Z12)
  
  Z13 = norm(W0,1)/l0 + norm(W1,1)/l1 + norm(W2,1)/l2
  
  if inf(T)>0
    Z14 = 0
  else
    Z14 = norm(e0 - WT*(e0-1/c*V2),1)
  end
  
  Z1 = sqrt(Z11^2 + Z12^2 + Z13^2 + Z14^2)
  
  display("value of Z1")
  display(Z1)
  
  
  # To gain computational memory, we decided to cut truncate the number of coefficients of the multiplication operator associated to U0.
  # Indeed, we use a truncation of size N0 instead of N. This introduces an extra bound in the computation of 𝒵1, which we denote Z1inf
  # and accounts for such a truncation. Specifically it provides an upper bound for \|BΛν(U0 - π^N0 U0)L^{-1}Λν^{-1}\|_2  
  
  Z1inf  = (opnorm(D1.*B*project(Multiplication((V0-project(V0,fourier0))^2),fourier0, fourier0,Interval{Float64})^2*Badj.*D2',2) + norm(V0-project(V0,fourier0),1)^2)^(1/2) 
  + (opnorm(D1.*B*project(Multiplication((V1-project(V1,fourier0))^2),fourier0, fourier0,Interval{Float64})^2*Badj.*D2',2)+norm(V0-project(V0,fourier0),1)^2)^(1/2) 
  +  (opnorm(D1.*B*project(Multiplication((V2-project(V2,fourier0))^2),fourier0, fourier0,Interval{Float64})^2*Badj.*D2',2)+norm(V0-project(V0,fourier0),1)^2)^(1/2)
  
  Z1inf = 1/σ0*Z1inf
  
  display("value of Z1inf")
  display(Z1inf)
  
  MV0 = Nothing ; MV1 = Nothing ; MV2 = Nothing; MW0 = Nothing; MW1 = Nothing; MW2 = Nothing
  
  ######################## Computation of Zu ##############################
  
  
    ##### Construction of the needed constants
  # Construction of ξ0 
  Ca = (interval(1) + abs(cos(interval(2)*a)))/(interval(1)-abs(cos(interval(2)*a)))
  
  if inf(T)>0
    ξ0 = maximum([interval(1) interval(1)/sqrt(T) c^2*Ca*interval(4)/T])
    ξ0 = maximum([interval(3)*T/(interval(2)*tanh(ξ0)) interval(2)*c^2/(T*tanh(ξ0))])
  else
    ξ0 = interval(1)
  end
  
  # Constants for the bound Y0
  if inf(T)>0
    a0 = minimum([interval(1)/sqrt(T) interval(1)/interval(2)*π])
    a0b = minimum([interval(1)/sqrt(Tb) interval(big(1))/interval(big(2))*π])
  else
    a0 = interval(1.5)
    a0b = interval(big(1.5))
  end
  
  
  # Notice that max_{s>0}min{√2 + √s, √2/(1-exp(-πs)) } ≤ max{1+√2, √2/(1-exp(-π))}, so we can compute CY0 as follows (cf. Lemma 4.1)
  pi_i = interval(π)
  if T==0
    CY0 = interval(2.28)*interval(0.5)*maximum([interval(1)+sqrt(interval(2)) sqrt(interval(2))/(interval(1)-exp(-pi_i))])
  else
    CY0 =maximum([interval(1)+sqrt(interval(2)) sqrt(interval(2))/(interval(1)-exp(-pi_i))])*interval(1)/(sqrt(pi_i)*T^(1/4))*(interval(2)/sqrt(interval(1)+a0*sqrt(T)) + interval(1)/sqrt(interval(2)))
  end
  
  # Other constants of Lemma 4.1
  C0 = interval(1)/(π*σ0*(interval(1)-ν*a^2)) + interval(1)/(π*ν*σ0)
  
  if inf(T)>0
    C1 = interval(1)/(interval(2)*π)*( interval(2)*(interval(1)+a)/(σ0*(interval(1)-ν*a^2)) + 4*(interval(1)+a)/(σ1*sqrt(T)*ν) )
    K1 = interval(2)*ξ0/(π*σ0) + interval(2)*sqrt(ξ0)*(interval(1)+abs(c))/(π*σ0*sqrt(T)) + interval(2)*(interval(1)/(3*T) + abs(c)/(4*T^(3/2)) + interval(2)*c^2/T )/(π*sqrt(tanh(ξ0)*T*ξ0)) + abs(c)/(T*π)*(interval(2)+3*log(ξ0)) + interval(1)/sqrt(interval(2)*π*T) 
    K2 = Ca*abs(ξ0+im*a)/( interval(2)*π*σ0^2*(interval(1)-T*a^2)^2)*( interval(1)/a^2 + T + Ca/a + Ca*T*abs(ξ0+im*a) ) + interval(2)*Ca^2/π*( interval(2)*(interval(1)+T)/sqrt(T*ξ0) + (interval(2)+a)/interval(2)*exp(-interval(2)ξ0))
    C2 = maximum([K2 K1*exp(a)])  
  else
    C1 = interval(1)/(interval(2)*abs(c)*ν) + Ca^(interval(1)/interval(4))*(interval(1)+sqrt(a))/(π*abs(c)*σ0)*( interval(1)/(interval(1)-ν*a^2) + interval(2)/ν )
    K1 = (interval(2)+4*exp(-interval(2)*ξ0))/(minimum([interval(1) abs(c)^3])*π*sqrt(ξ0)*σ0) + interval(1)/(π*σ0*abs(c)) + interval(2)/(π*c^2) + interval(1)/(π*abs(c)^3)*(interval(2)+3*log(3*ξ0)) + interval(1)/(c^2*sqrt(interval(2)*π))
    K2 = Ca^(interval(1)/4)/π*( interval(1)/(interval(2)*σ0^2)*( interval(2) + a^(-3/2)) + interval(1)/(4*σ0^2*(interval(1)-abs(cos(interval(2)a)))^2*sqrt(a)) )
    C2 = maximum([K2 K1*exp(a)]) 
  end
  
  
  
  
  
  
  ###################### Computation of 𝒵u #####################################
  
  #### Computation of the constants C(d) and C1(d) in Lemmas 4.2 and 4.6
  
  C1d = exp(-interval(2)*a*d)*( interval(2)*sqrt(π)/(sqrt(interval(4)*a*d)*(interval(1)-exp(-interval(2)*a*d))) + interval(4)/(interval(1)-exp(-interval(2)*a*d)) )
  Cd = exp(-interval(2)*a*d)*( interval(4)*d + interval(4)*exp(-a*d)/(a*(interval(1)-exp(-interval(3)*a*d/interval(2)))) + interval(2)/(a*(interval(1)-exp(-a*d*interval(2)))))
  
  ### Inner product (V0,E*V0)
  D1N = convert(Vector{Interval{Float64}},interval.(exp2cos(N)))
  D2N = interval.(ones((N+1)))./D1N
  PS0 = abs( (coefficients(D1N.*project(V0*E,fourier))')*coefficients(D1N.*V0) )
  
  ### Inner product (V1,Efull*V1)
  k = (-N:N)*π/db
  fourier_f = Fourier(N,π/db)
  V1f = interval(4)*νb*coefficients(U0b) ; V1f = [reverse(V1f[2:N+1]); V1f[1] ; V1f[2:N+1]].*k ; V1f = Sequence(fourier_f,V1f)
  Ef = coefficients(Eb); Ef = [reverse(Ef[2:2*N+1]); Ef[1] ; Ef[2:2*N+1]] ; Ef = Sequence(Fourier(2*N,π/db),Ef)
  PS1 = abs( (coefficients(project(V1f*Ef,fourier_f))')*coefficients(V1f) )
  
  ### Inner product (V2,E*V2)
  PS2 = interval(8)*abs( (coefficients(project(U0b*Eb,fourier))')*coefficients(U0b) )
  
  # Computation of ∫|v2'|^2 from d-1 to d
  char = Sequence(Fourier(2*N,π/db) ,interval.(big.(zeros(4*N+1))))
  for n = [-2N:-1;1:2N]  
      char[n] = real(interval(big((-1)^n))*(exp(im*n*π/db)-interval(big(1)))/(interval(big(2))*im*n*π))
  end
  char[0] = interval(big(1))/(interval(big(2))*db)
  Elog =  abs( (coefficients(project(V1f*char,fourier_f))')*coefficients(V1f) )
  
  # Computation of the Zu_i
  Zu0 = interval(4)*d*C0^2/a*PS0 + interval(2)*d*Cd*C0^2*PS0
  Zu1 = interval(4)*d*C1^2/a*PS1 + interval(2)*d*Cd*C1^2*PS1
  Zu2 = interval(2)*C2^2*( interval(2)*d*PS2/a + interval(4)*log(interval(2))*Elog ) + interval(2)*d*C1d*C2^2*PS2
  
  𝒵u = norm_B*sqrt( Zu0 + Zu1 + Zu2 )
  
  display("value of 𝒵u")
  display(𝒵u)
  
  ######################## Computation of ℨ1 ##############################
  
  #in addition we add the continuous term
  
  𝒵1 = Z1 + 𝒵u + Z1inf + abs(c2-c1)*norm_B/2
  display("value of 𝒵1")
  display(𝒵1)
  
  
  ################### Computation of Z2 ######################
  if (T==0)||(inf(T)>1/3)
    κT = interval(1)/(abs(interval(1)-c)*two_norm_inverse_l(d,ν,c,T))
  else
    κT = interval(1)/(sqrt(ν)*σ0^2)
  end
  𝒵2 = interval(2)*κT*norm_B
  
  display("value of 𝒵2")
  display(𝒵2)
  
  
  #################### Computation of Y0 ##################
  
  Cd = exp(-interval(2)*a0*d)*( interval(4)*d + interval(4)*exp(-a0*d)/(a0*(interval(1)-exp(-interval(3)*a0*d/interval(2)))) + interval(2)/(a0*(interval(1)-exp(-a0*d*interval(2)))))
  
  D1 = convert(Vector{Interval{BigFloat}},interval.(big.(exp2cos(N))))
  D2 = interval.(big.(ones((N+1))))./D1

  
  
  Eb = Sequence(fourierE ,interval.(big.(zeros(2*N+1))))
  for n = 1:2N  
      Eb[n] = real(interval(big((-1)^n))*(interval(big(1))/(im*interval(big(n))/db*interval(big(pi)) + interval(big(2))a0b) - interval(big(1))/(im*interval(big(n))/db*interval(big(pi))-interval(big(2))*a0b)))
  end
  Eb[0] = interval(big(1))/a0b
  Eb = Eb*(interval(big(1))-exp(-interval(big(4))*a0b*db))/(interval(big(4))*db)
  
  
  ### Inner product part 
  PSY0 = abs( (coefficients(D1.*project(((Lνb.*U0b)*Eb),fourier))')*coefficients(D1.*(Lνb.*U0b)) )
  
  Yu1 = sqrt(interval(4)*d^2*CY0*PSY0) ;
  
  display("bound Yu1")
  display(Yu1)
  
  D² = project(Derivative(2), CosFourier(2N, π/db), CosFourier(2N, π/db),Interval{Float64})
  Lν2 = Id - ν*D²
  Lν2 = diag(coefficients(Lν2))
  
  D² = project(Derivative(2), fourier, fourier,Interval{Float64})
  # Construction of the operator L
  d2 = diag(coefficients(D²))  ;   dd = sqrt.(-d2) ; d0=dd;  d0[1] = interval(1);
  dL = tanh.(dd).*(ones(N+1)+T*dd.^2)./d0 ;  dL[1] = interval(1) ; dL = sqrt.(dL)
  MT = LinearOperator(fourier,fourier,Diagonal(dL))
  L = Lν2[1:N+1].*diag(coefficients((MT - c*Id)))
  
  FU0 = L.*U0 + Lν2.*(U0*U0)
  Y0 = sqrt(interval(2)*d)*sqrt( norm(B*project(FU0,fourier),2)^2 + norm(WT*FU0-project(WT*FU0,fourier),2)^2)
  u = (u2-u1)/(c2-c1)
  
  Yc = sqrt(interval(2)*d)*sqrt( norm(B*(L*u - Lν.*u1 + 2*Lν.*(u1*u)),2)^2 + norm(WT*2*Lν.*(u1*u)-project(WT*2*Lν.*(u1*u),fourier),2)^2 )*abs(c2-c1) + sqrt(interval(2)*d)*sqrt(norm(B*Lν.*(u^2-u),2)^2 + norm(WT*Lν.*(u^2-u) - project(WT*Lν.*(u^2-u),fourier),2)^2)*abs(c2-c1)^2
  
  display("value of Yc")
  display(Yc)
  
  𝒴0 = Y0 + Yu1*sqrt(interval(1) + norm_B^2*(interval(1)+Cd)) + Yc
  display("value of 𝒴0")
  display(𝒴0)
  
  α = interval( sup( interval(1)/interval(2) + 𝒴0*𝒵2/(interval(1)-𝒵1)^2 ) )
  r = interval(2)*α*𝒴0/(interval(1)-𝒵1)
  
  if inf(1- 𝒵1)>0
    if sup(𝒵2*r+𝒵1) < 1
        if inf(1/2*𝒵2*r^2 - (1-𝒵1)*r + 𝒴0) < 0
            display("The computer-assisted proof was successful")
            display("Radius of contraction :")
            display(r)
        else
        display("failure: discriminant is negative")
        end
    else
        display("r is too big")
    end
  else
    display("failure: 1-z > 0")
  end
  
  end
  
  
  
  
  
  
  
  










  
function trace(N)
    w = 0:N;
    Y = w.^2;
    X = 2*(-1).^w; X[1] = 1;
  
    f = [X';
         (X.*Y)']
    return f
  end
  
  
  
  function ini(x,c,T,b)
    c= c-1;
    @. 2*c/cosh(b*sqrt(abs(3*c/(2*(1-3*T))))*x)^2;
  end
  
  
  
  
  
  
  function compute_boundary(V,N,l,a)
  S = interval(big(0));
  p = interval(big(pi))
  for n = 1:N
    n1 = interval(big(n))
    for m = 1:N
  
                m1 = interval(big(m))
               b = 4*a*(1/((p/l*(n1-m1))^2+4*a^2) + 1/((p/l*(n1+m1))^2+4*a^2));
               u= V[n];
               v= V[m];
               S = S+ (-1)^(n-m)*u*v*b;
    end
               b = 4*a/((p/l*(n1))^2+4*a^2) ;
               u= V[n];
               v= V[0];
               S = S+(-1)^(n)*u*v*b;
  end
  
       for m = 1:N
           m1 = interval(big(m))
               b = 4*a/((p/l*(m1))^2+4*a^2)  ;
               u= V[0];
               v= V[m];
               S = S+ (-1)^(m)*u*v*b;
       end
        b =  1/(2*a) ;
        u= V[0];
        v= V[0];
        S = S+ u*v*b;
  
  f = (sqrt(abs(S)));
  end
  
  
  
  
  
  
  
  function compute_boundary_1(V,N,l,a)
  S = interval(big(0));
  p = interval(big(pi))
  for n = 1:N
    for m = 1:N
                n1 = interval(big(n))
                m1 = interval(big(m))
               b = 4*a*(1/((p/l*(n1-m1))^2+4*a^2) - 1/((p/l*(n1+m1))^2+4*a^2)) ;
               u= n1*p/l*V[n];
               v= m1*p/l*V[m];
               S = S+ (-1)^(n-m)*u*v*b;
    end
  
  end
  
  f = (sqrt(abs(S)));
  end
  
  
  
  
  
  
  
  
  
  
  function compute_boundary_total_new(V,N,d,a,C1,C0,ep,diff_norm)
  
  
  ######## computation of the constants A_0(w)
    A0 =   compute_boundary(V,N,d,a)
    A1 =  compute_boundary_1(V,N,d,a)
    A2 =  compute_boundary_1(V,N,d,a/2)
  
  ####### computation of C1 and C3
  C_0 = C0*(2*a*d+3)/(a^2*ep)*A1    #Lw (L1 norm)
  C_1 = 2*C0*(2*a*d+1)/(a*ep)*A1   #Lw(d) and Lw(-d)
  C_2 = sqrt(2/a)*C0/ep*A2         #Lw (l2 norm)
  
  ######### tilde part
  
  C_3 = (4*a+2)*C0/a*A0           #Lw (l1 norm)
  C_4 = 2*C1/a*A0                 #dxLw
  C_5 = C0/sqrt(a)*A0             #Lw (l2 norm)
  C_6 = (4*a+2)*C0/a*A1           #Lw' (l1 norm)
  
  ######## computation of C
  
  
  N = 16
  
  
  I1 = sqrt((2*N+1)/(2*d)*C_0^2 + 2*d/(pi^2*N)*(C_6 + C_1)^2 )
  I2 = sqrt((2*N+1)/(2*d)*C_3^2 + 2*d/(pi^2*N)*C_4^2 )
  
  C = I1 + I2 + C_5 + C_2 + 2*diff_norm
  
    return C
  end
  
  
  
  
  
  
  function compute_boundary_Y(V,W,N,d,a,C1,diff,ep)
  
  
  ######## computation of the constants A_0(w)
    A0 =   compute_boundary(V,N,d,a)
    A1 =  compute_boundary_1(V,N,d,a)
    A2 =   compute_boundary(W,N,d,a)
  
  ####### computation of C1 and C3
  
  C_0 = 2*C1*(sqrt(2*d)*(exp(-a)+1/(ep*(1+ep)))*(A0+A1) + diff/a*ep^(5/2) )
  
  C_1 = 2*C1*(sqrt(2*d)*(exp(-a)+1/(ep*(1+ep)))*(A1+A2) + 3*diff/a*ep^(3/2) )
  
  C_2 = 2*C1*sqrt((ep*exp(-2a)/3)*(A0^2+A1^2) + diff^2/a^2*ep^(5) )
  C_2 = C_2 + 2*C1*sqrt((  ((1+ep)^3-1)/(ep^3*(1+ep)^3) +exp(-2a)/3 )*(A0^2+A1^2))
  
  ######## computation of C
  
  N = 16
  
  C = sqrt( (2*N+1)/(2*d)*C_0^2 + 2*d/(pi^2*N)*C_1^2 ) + C_2
  
    return C
  end
  
  
  
  function eval_f_cos(V,N,d,dx,ep)
  
    di = interval(big(d))
    p = interval(big(pi))
    S = V[0];
    M= 0;
  for i=d-ep:dx:d
    x = interval(big(i),big(i+dx))
    for m = 1:N
        S = S + 2*V[m]*cos(big(m)*p/di*x)
    end
    S = abs(S)
    M = max(sup(S),M)
    S = V[0];
  end
  display(M)
  return M
  end
  
  
  
  
  
  
  
  function exp2cos(N)
  
    d = 2*(ones((N+1)))
  
    d[1] =1;
  
    return d
  end
  
  
  
  
  function test_constants(T,ν, c,a,σ0,σ1,x,di,dj,xmin,ymin)
  
    setprecision(100)
  
    if inf(T)>0
        # first requirement on a
        if inf(a) >= sup(minimum([π/2 1/sqrt(ν)]))
            return 0,0
        end
        # we start by verifying σ0 and σ1, for these we can fix the imaginary part y= a
        X =  -sup(x):dj:sup(x)-dj 
        Z = interval.(X,X.+dj*ones(size(X))).+im*a
  
        m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z).*( 1 .+ T*Z.*Z)).-c).-σ0
        m2 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z).*( 1 .+ T*Z.*Z)).-c).-σ1*sqrt.(T*abs.(X))
         
        if (minimum(inf.(m1))<=0)||(minimum(inf.(m2))<=0)
          display(minimum(inf.(m1)))
            return 0,inf.(m1)
        end
  
        
        # Verify that |m_T-c| >0 on S_x and |y| ≥ ymin
        ddi = 0.01
        ddj = 0.01
        y1 = -sup(a) ; y2 = -sup(ymin)
        x1 = -sup(x) ; x2 = sup(x) 
        k=0
        while y1<y2
            while x1<x2
                Z = interval(x1,x1+ddi)+im*interval(y1,y1+ddj)
                m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z).*( 1 .+ T*(Z.*Z))).-c)
                while (minimum(inf.(m1))<=0)&&(k<6) 
                    ddi = ddi/10
                    ddj=ddj/10
                    Z = interval(x1,x1+ddi)+im*interval(y1,y1+ddj)
                    m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z).*( 1 .+ T*(Z.*Z))).-c)
                    k=k+1
                end
                if k==6
                    display("value for a is too big test 1")
                    return 0,0
                end
                k=0
                x1 = x1+ddi
                y1 = y1+ddj
                ddi = 0.01
                ddj = 0.01
            end
        end
  
  
        
        # Verify that |m_T-c| >0 on S_x and |y| ≥ ymin
        ddi = 0.01
        ddj = 0.01
        y1 = -sup(ymin) ; y2 = sup(ymin)
        x1 = -sup(x) ; x2 = -sup(xmin) 
        k=0
        while y1<y2
            while x1<x2
                Z = interval(x1,x1+ddi)+im*interval(y1,y1+ddj)
                m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z).*( 1 .+ T*(Z.*Z))).-c)
                while (minimum(inf.(m1))<=0)&&(k<4) 
                    ddi = ddi/10
                    ddj=ddj/10
                    Z = interval(x1,x1+ddi)+im*interval(y1,y1+ddj)
                    Z2 = real(Z)^2-imag(Z)^2 + 2*im*real(Z)*imag(Z)
                    m1 = abs.(sqrt.((exp.(interval(2)*Z).-interval(1.0))./((exp.(interval(2)*Z).+interval(1.0)).*Z).*( interval(1) .+ T*Z*Z)).-c)
                  
                    k=k+1
                end
                if k==4
                    display("value for a is too big test 2")
                  
                    return 0,0
                end
                k=0
                x1 = x1+ddi
                y1 = y1+ddj
                ddi = 0.01
                ddj = 0.01
            end
        end
  
  
  
        
        # Verify that |m_T-c| >0 on S_x and |y| ≥ ymin
        ddi = 0.01
        ddj = 0.01
        y1 = -sup(ymin) ; y2 = sup(ymin)
        x1 = inf(xmin) ; x2 = sup(x) 
        k=0
        while y1<y2
            while x1<x2
                Z = interval(x1,x1+ddi)+im*interval(y1,y1+ddj)
                m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z).*( 1 .+ T*Z.*Z)).-c)
                while (minimum(inf.(m1))<=0)&&(k<4) 
                    ddi = ddi/10
                    ddj=ddj/10
                    Z = interval(x1,x1+ddi)+im*interval(y1,y1+ddj)
                    m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z).*( 1 .+ T*Z*Z)).-c)
                    k=k+1
                end
                if k==4
                    display("value for a is too big test 3")
                    return 0,0
                end
                k=0
                x1 = x1+ddi
                y1 = y1+ddj
                ddi = 0.01
                ddj = 0.01
            end
        end
          
        # # Verify that |m_T-c| >0 on S_x and |y| ≤ ymin, |x| ≥ xmin
        # Y= -sup(ymin):di:inf(ymin)
        # X= [-sup(x):dj:-sup(xmin) ; inf(xmin):dj:sup(x)-dj]
        # NX = length(X) ;  NY = length(Y) ;
        # X = X'.*ones(NY); Y = ones(NX)'.*Y ;
        # Z = interval.(X,X.+dj*ones(size(X)))+im*interval.(Y,Y.+di*ones(size(Y)))
  
        # # m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z).*( 1 .+ T*Z.^2)).-c)
        # m1 = abs.(sqrt.(tanh.(Z).*(1 .+ T*Z.^2)./Z).-c)
        
        # if (minimum(inf.(m1))<=0)
        #     display("second test")
        #     return 0,inf.(m1)
        # end
             
  
    else
        if sup(a) >= π/2
            return 0,0
        end
  
         # we start by verifying σ0 and σ1, for these we can fix the imaginary part y= a
        X =  -sup(x):dj/4:sup(x)-dj/4 ; 
        Z = interval.(X,X.+dj/4*ones(size(X))).+im*a
  
        m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z)).-c).-σ0
        if minimum(inf.(m1))<=0
            display("fails for σ0")
            return 0,inf.(m1)
        end
  
        # Verify that |m_T-c| >0 on S_x and |y| ≥ ymin
        ddi = 0.01
        ddj = 0.01
        y1 = -sup(a) ; y2 = -sup(ymin)
        x1 = -sup(x) ; x2 = sup(x) 
        k=0
        while y1<y2
            while x1<x2
                Z = interval(x1,x1+ddi)+im*interval(y1,y1+ddj)
                m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z)).-c)
                while (minimum(inf.(m1))<=0)&&(k<4) 
                    ddi = ddi/10
                    ddj=ddj/10
                    Z = interval(x1,x1+ddi)+im*interval(y1,y1+ddj)
                    m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z)).-c)
                    k=k+1
                end
                if k==4
                    display("value for a is too big 4")
                    return 0,0
                end
                k=0
                x1 = x1+ddi
                y1 = y1+ddj
                ddi = 0.01
                ddj = 0.01
            end
        end
          
        # Verify that |m_T-c| >0 on S_x and |y| ≤ ymin, |x| ≥ xmin
        Y= -sup(ymin):di:inf(ymin)
        X= [-sup(x):dj:-sup(xmin) ; inf(xmin):dj:sup(x)-dj]
        NX = length(X) ;  NY = length(Y) ;
        X = X'.*ones(NY); Y = ones(NX)'.*Y ;
        Z = interval.(X,X.+dj*ones(size(X)))+im*interval.(Y,Y.+di*ones(size(Y)))
  
        m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z)).-c)
        if (minimum(inf.(m1))<=0)
            display("second test")
            return 0,inf.(m1)
        end
            
    end
    return 1,1
  end
  
  
  
  
  using MATLAB
  
  function PlotCoeffs1D(U0,a,b)
    #U0 is a sequence in 2D
    #a,b,c,d are the endpoints of the interval [a,b] × [c,d]
    y1 = a:0.01:b
    m=length(y1)
    U = zeros(m)
    for b₁ = 1:m
            U[b₁] = real(U0(y1[b₁]))
    end
    mat"
    h = plot($y1,$U)"
  end


  function PlotCoeffs2D(U0,a1,a2,b1,b2)
    #U0 is a sequence in 2D
    #a,b,c,d are the endpoints of the interval [a,b] × [c,d]
    y1 = LinRange(a1,a2,200)
    y2 = LinRange(b1,b2,3000)
    m1=length(y1) ; m2=length(y2)
    U = zeros(m1,m2)
    for k1 = 1:m1
      for k2 = 1:m2
            U[k1,k2] = real(U0(2/(a2-a1)*y1[k1] + 1-2*a2/(a2-a1),y2[k2]))
      end
    end
    mat"
    h = surf($y2,$y1,$U)
     set(h,'LineStyle','none')"
  end
  
  

  function Plot_branch(U0,a1,a2,b1,b2)
    #U0 is a sequence in 2D
    #a,b,c,d are the endpoints of the interval [a,b] × [c,d]
    y1 = LinRange(a1,a2,5)
    y2 = LinRange(b1,b2,3000)
    m1=length(y1) ; m2=length(y2)
    U = zeros(m1,m2)
    for k1 = 1:m1
      for k2 = 1:m2
            U[k1,k2] = real(U0(2/(a2-a1)*y1[k1] + 1-2*a2/(a2-a1),y2[k2]))
      end
    end
    mat" 
    hold on
    for i = 1:5
      plot($y2,$U(i,:))
    end
     "
  end
  
  
  
  function two_norm_inverse_l(d,ν,c,T)
  
      S = interval(0)
      N = floor(mid(100*d/π))
      for n=1:N
        mT = sqrt(tanh(n*π/d)*(1+T*(n*π/d)^2)/(n*π/d))
        S = S + 1/((mT-c)^2*(1+ν*(n*π/d)^2)^2)
      end
      S = S + 1/((1-c)^2) + d^4/(3*π^4*ν^2*N^3*sqrt(tanh((N+1)*π/d)*(1+T*((N+1)*π/d)^2)/(N*π/d))^2)
      S = S/d 
      return sqrt(S)
  
    end
  
  
  
  
  function proof_eigen(U,U0,c,T,ν,d,rmin)
  
    setprecision(100)
  
    pspace = space(U)
    fourier = space(U0)
    N = order(U0)
  
    V = component(U,2)  #(V correspond to Ψ in Section 5)
    λ = interval(real(component(U,1)[1]))

    Vb = interval.(big.(real.(V))) + im*interval.(big.(imag.(V)))
  
    #construction of the linear operator MT
    D² = real.(project(Derivative(2), fourier, fourier,Complex{Interval{Float64}}))
    dd = diag(coefficients(D²))  ;   dd = sqrt.(-dd) ; d0=dd;  d0[N+1] = interval(1);
    dL = tanh.(dd).*(ones(2*N+1)+T*dd.^2)./d0 ;  dL[N+1] = interval(1) ; dL = sqrt.(dL)
    MT = LinearOperator(fourier,fourier,Diagonal(dL))
    #construction of L
    if inf(T)>0
        L = MT - c*Id - λ*Id
        Linv = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(L))))
        DG = 2*real.(project(Multiplication(U0),fourier, fourier,Complex{Interval{Float64}}))
    else
        L = c*Id - MT - λ*Id
        Linv = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(L))))
        DG = -2*real.(project(Multiplication(U0),fourier, fourier,Complex{Interval{Float64}}))
    end
  
    # projection of V on trace zero (V correspond to Ψ in Section 5)
    S = trace_full(N,d); Sᵀ = S' ;
    ddd= ones(2N+1)./(1:2*N+1)
    Dtrace = Matrix(interval.(big.(Diagonal(ddd)^3)))
  
    Vb = Vb - Sequence(fourier,Dtrace*Sᵀ*solve_linear(S*Dtrace*Sᵀ,S*coefficients(Vb),100))
    V = interval.(Float64.(inf.(real.(Vb)),RoundDown),Float64.(sup.(real.(Vb)),RoundUp) ) + im*interval.(Float64.(inf.(imag.(Vb)),RoundDown),Float64.(sup.(imag.(Vb)),RoundUp) ) 
  
    U = Sequence(pspace,vec([λ; coefficients(V)]))

    
  
    # 
    M = LinearOperator(pspace,pspace, [[0   -2*d*coefficients(L*(L*V))'] ; [-coefficients(V) coefficients(L+DG)]] )
    D = [[interval(1) zeros(1,2*N+1)] ; [zeros(2*N+1,1) coefficients(L)]]
    Di = [[interval(1) zeros(1,2*N+1)] ; [zeros(2*N+1,1) coefficients(Linv)]]
    D = LinearOperator(pspace,pspace,D)
    Di = LinearOperator(pspace,pspace,Di)
   
  
  # # # # ################ Z BOUND ######################################################
  
    B = inv(mid.(M)*mid.(Di)) ; B = interval.(real.(B)) + im*interval.(imag.(B))
    Badj = LinearOperator(pspace,pspace,coefficients(B)')

   
  
    e0 = Sequence(fourier,interval.(zeros(2*N+1)))
    e0[0] = 1
    if inf(T)>0
        WT = e0
    else
        e02 = Sequence(Fourier(2*N, π/d),interval.(zeros(4*N+1)))
        e02[0] = 1
        WT = interval.(mid.(project(Multiplication(mid.(e0 - 2/mid(c-λ)*U0)),fourier,Fourier(2*N, π/d),Interval{Float64}))\mid.(e02))
    end


  
  # computation of the norm of B
    MWT = project(Multiplication(WT*WT),fourier, fourier,Interval{Float64}) - project(Multiplication(WT),fourier, fourier,Interval{Float64})*project(Multiplication(WT),fourier, fourier,Interval{Float64})
    n_WT = sqrt(opnorm(LinearOperator(coefficients(MWT)),2))   # computation of \|π_NW_Tπ^N\|_2
    n_BN = sqrt(opnorm(LinearOperator(coefficients((B*Badj))),2))   #computation of \|B^N_T\|_2
    norm_B = maximum([1 maximum([norm(WT,1) n_BN])+n_WT])  #computation of the norm of B

  
  
    if inf(T)>0
      L2 = Linv ;          l2 = maximum(vec((abs.(dL.-(c+λ)))))
    else
      L2 = Linv + 1/(c-λ)*Id ;  l2 = maximum(vec(abs.(((c-λ)*(dL.-(c-λ)))./dL)))
    end
  
      W2 = 2*WT*U0
      #MWi corresponds to the operator πᴺ(Vi*WT)(π^{3N}-πᴺ)(Vi*WT)πᴺ = πᴺ(Vi*WT*Vi*WT)πᴺ - (πᴺ(Vi*WT)πᴺ)^2 (recall that Wi = Vi*WT)
      MW2 = project(Multiplication(W2*W2),fourier, fourier,Interval{Float64}) - project(Multiplication(W2),fourier, fourier,Interval{Float64})^2
  #   MVi corresponds to the operator πᴺVi(π^{2N}-πᴺ)Viπᴺ = πᴺ(Vi*Vi)πᴺ - (πᴺViπᴺ)^2
      MV2 = 4*project(Multiplication(U0*U0),fourier, fourier,Interval{Float64}) - 4*project(Multiplication(U0),fourier, fourier,Interval{Float64})^2
      MV2 =  LinearOperator(pspace,pspace, [[0   zeros(1,2*N+1)] ; [zeros(2*N+1,1) coefficients(MV2)]] )
  
    Z11 =  opnorm(LinearOperator(coefficients((Id - B*M*Di))),2)^2
    Z11 = Z11 + sqrt(opnorm(LinearOperator(coefficients(L2*MW2*L2)),2))
    Z11 = sqrt(Z11)

    Z12 = sqrt(opnorm(LinearOperator(coefficients(B*MV2*Badj)),2))/l2 
  
    Z13 =  norm(W2,1)/l2
  
    if inf(T)>0
      Z14 = 0
    else
      Z14 = norm(e0 - WT*(e0-2/(c-λ)*U0),1)
    end
  
    Z1 = sqrt(Z11^2 + Z12^2 + Z13^2 + Z14^2)
  
  
  
  
  #computation of the supremum of Zu for all λ ∈ (0,λmin)
  # We first compute the supremum of Cλ
  λmin = abs(1-c)-2*norm(U0,1)-rmin/(4*sqrt(ν)*abs(1-c)) # λmin is a lower bound for the eigenvalues
  if inf(T)>0
    Cλ = computation_Cλ(maximum([c c+λmin]),a,σ0,T)
  else
    Cλ = computation_Cλ(c-λmin,a,σ0,T)
  end
  
  # coefficients of the characteristic function
  char = Sequence(Fourier(2*N,π/d) ,interval.(big.(zeros(4*N+1))))
  for n = [-2N:-1;1:2N]  
    char[n] = real((-1)^n*(exp(im*n*π/db)-1)/(2*im*n*π))
  end
  char[0] = 1/(2*db)

  a0b = interval(big(1.5))
  
  # Coefficients of cosh(2ax)
  Eb = Sequence(Fourier(2*N,π/d),interval.(big.(zeros(4*N+1))))
  for n = -2N:2N  
    Eb[n] = a0b/db*(-1)^n*(1-exp(-4*a0b*db))/(4*a0b^2 + (n*interval(big(π)/db))^2 )
  end
  
  
  # For all λ ≤0, a-λ ≥ a and so the maximal 𝒵λ is achieved for λ=0
  C1d = exp(-2*a*d)*( 2*sqrt(π)/(sqrt(4*a*d)*(1-exp(-2*a*d))) + 4/(1-exp(-2*a*d)) )
  DU0 = U0fb.*( (-N:N)*π/db )
  Elog =  abs( (coefficients(project(DU0*char,fourier))')*coefficients(DU0) )
  𝒵u0 = abs( (coefficients(project(U0fb*Eb,fourier))')*coefficients(U0fb) )
  
  #Computation of the supremum of 𝒵λ for all λ ≤ 0
  𝒵u = 8*Cλ^2*(d*𝒵u0*(2/a+C1d) + 4*log(interval(2))*Elog )
  𝒵u = sqrt(𝒵u)
  
  𝒵u = interval.(Float64(inf(𝒵u),RoundDown),Float64(sup(𝒵u),RoundUp) ) # go back to standard precision
  
  
    ## Computation of 𝒵1
    𝒵1 = sup(𝒵u*norm_B + Z1 + norm_B*rmin/(4*sqrt(ν)*σ0*(σ0-λ)))
  
  
  ################# Y BOUND ######################################################
  a0 = interval(1.5)
  Cd = exp(-2*a0*d)*( 4*d + 4*exp(-a0*d)/(a0*(1-exp(-3*a0*d/2))) + 2/(a0*(1-exp(-a0*d*2))))
  
  D² = real.(project(Derivative(2), Fourier(N, π/db), Fourier(N, π/db),Complex{Interval{BigFloat}}))
  Lνb = νb*D²
  Lνb = diag(interval(big(1)).-coefficients(Lνb)).*2
  
  Eb = Sequence(Fourier(2*N, π/d) ,interval.((zeros(4*N+1))))
  for n = -2N:2N  
      Eb[n] = real((-1)^n*(1/(im*n/db*π + 2a0b) - 1/(im*n/db*π-2*a0b)))
  end
  Eb[0] = 1/a0b
  Eb = Eb*(1-exp(-4*a0b*db))/(4*db)
  
  ### Inner product part 
  PSY0 = abs( (coefficients(project(((Lνb.*Vb)*Eb),fourier))')*coefficients((Lνb.*Vb)) )
  pi_i = interval(π)
  if T==0
    CY0 = interval(2.28)*interval(0.5)*maximum([interval(1)+sqrt(interval(2)) sqrt(interval(2))/(interval(1)-exp(-pi_i))])
  else
    CY0 =maximum([interval(1)+sqrt(interval(2)) sqrt(interval(2))/(interval(1)-exp(-pi_i))])*interval(1)/(sqrt(pi_i)*T^(1/4))*(interval(2)/sqrt(interval(1)+a0*sqrt(T)) + interval(1)/sqrt(interval(2)))
  end
  Yu = sqrt((4*d^2*CY0*PSY0)*(1+norm_B^2*(Cd+1)))
 
  B22 = component(B,2,2)
  B12 = component(B,1,2)
  Y0 = sqrt(2*interval(d))*( sqrt(2*norm(B22*(L*V+DG*V),2)^2 + 2*norm(2*U0*V - project(2*U0*V,fourier))^2)) + 2*sqrt(2*interval(d))*norm(L*V+ DG*V,2)*norm(L*V,2)

  
  𝒴0 = Yu + Y0 + 2*rmin/σ0*maximum([norm(B12*V,2)/(sqrt(2*d))  opnorm(LinearOperator(coefficients(B22)),2)*norm(V,1) ])


  ################## Z2 BOUND ######################################################
    Z2 = norm_B
 
  ################## Verification ###################################################
  # # # r_star = 1.0;
  # # # interval_proof = interval_of_existence(Y,Interval(0),Z2,r_star,C²Condition())
  # # # display(interval_proof)
  # # #
  x = 1;
  if inf(1- 𝒵1)>0
    if inf(1-sup(4*𝒴0*Z2)) > 0
      rmin=sup((1-𝒵1 - sqrt((1-𝒵1 )^2-4*𝒴0*Z2))/(2*Z2));
      rmax=inf((1-𝒵1  + sqrt((1-𝒵1 )^2-4*𝒴0*Z2))/(2*Z2));
      if rmin<rmax
        display("The computer-assisted proof was successful for the eigenvalue ν with value")
        display(Float64.(mid(λ)))
        # display("Value of the minimal radius for the contraction")
        # display(Float64(rmin))
        # display("Value of the maximal radius for the contraction")
        # display(Float64(rmax))
        
      else
        display("rmin>=rmax")
        x=0
      end
    else
      display("failure: discriminant is negative")
      x=0
    end
  else
      display("failure: 1-z > 0")
      x=0
  end
  
  return [x==1,rmax,rmin]
  
  end
  
  



function VK_quantity(U,rmin,σ0, norm_DF_inv)

  #upper bound for \|DF^{-1}\|
    𝒞 = σ0*norm_DF_inv/(1-2*rmin/(4*sqrt(ν)*σ0))
  
    D² = project(Derivative(2), fourier0, fourier0,Interval{Float64})
    Lν = Id - ν*D²
    
    # Construction of the operator L
    d2 = diag(coefficients(D²))  ;   dd = sqrt.(-d2) ; d0=dd;  d0[1] = interval(1);
    dL = tanh.(dd).*(ones(N0+1)+T*dd.^2)./d0 ;  dL[1] = interval(1) ; dL = sqrt.(dL)
    MT = LinearOperator(fourier0,fourier0,Diagonal(dL))
    
    L =MT - c*Id  
    DG = interval(2)*project(Multiplication(U),fourier0, fourier0,Interval{Float64})
  
    A = interval.(inv(mid.(L + DG)))

    V0 = mid.(A)*mid.(U) 
    V = interval.(big.(coefficients(V0)))
    D = interval.(big.((1:N+1).^4))
    S = trace(N); C = S' ; S =  interval.(big.(S)) ;  C =  interval.(big.(C))
    V_p = (D.*C)*solve_linear(S*(D.*C),vec(S*V),100)   # the function solve_linear solves a linear system rigorously
    Vb = Sequence(fourier,V - vec(V_p))
    V = interval.(Float64.(inf.(Vb),RoundDown),Float64.(sup.(Vb),RoundUp) )



    a0b = interval(big(1.5)) ; a0= interval(1.5)
  
    pi_i = interval(π)
    if T==0
      CY0 = interval(2.28)*interval(0.5)*maximum([interval(1)+sqrt(interval(2)) sqrt(interval(2))/(interval(1)-exp(-pi_i))])
    else
      CY0 =maximum([interval(1)+sqrt(interval(2)) sqrt(interval(2))/(interval(1)-exp(-pi_i))])*interval(1)/(sqrt(pi_i)*T^(1/4))*(interval(2)/sqrt(interval(1)+a0*sqrt(T)) + interval(1)/sqrt(interval(2)))
    end
    D1 = convert(Vector{Interval{BigFloat}},interval.(big.(exp2cos(N))))
    D2 = interval.(big.(ones((N+1))))./D1
    D² = project(Derivative(2), CosFourier(N, π/db), CosFourier(N, π/db),Interval{BigFloat})
    Lνb =  νb*D²
    Lνb = diag(interval(big(1)).-coefficients(Lνb)).*2
    
    fourierE = CosFourier(2*N,interval(π)/d)
    Eb = Sequence(fourierE ,interval.(big.(zeros(2*N+1))))
    for n = 1:2N  
        Eb[n] = real(interval(big((-1)^n))*(interval(big(1))/(im*interval(big(n))/db*interval(big(pi)) + interval(big(2))a0b) - interval(big(1))/(im*interval(big(n))/db*interval(big(pi))-interval(big(2))*a0b)))
    end
    Eb[0] = interval(big(1))/a0b
    Eb = Eb*(interval(big(1))-exp(-interval(big(4))*a0b*db))/(interval(big(4))*db)
    
    
    ### Inner product part 
    PSY0 = abs( (coefficients(D1.*project(((Lνb.*Vb)*Eb),fourier))')*coefficients(D1.*(Lνb.*Vb)) )
    
    Yu1 = sqrt(interval(4)*d^2*CY0*PSY0) ;


   ϵ = 𝒞*(Yu1 + norm(U0 - L*V - 2*V*U0,2)) + rmin/(interval(2)*sqrt(ν)*σ0)*sqrt(interval(2)*d)*𝒞*norm(V,2)

   display((coefficients(U)')*coefficients(V))

   τ = interval(2)*d*(coefficients(U)')*coefficients(V) + ϵ*sqrt(interval(2)*d)*norm(U,2) + interval(2)*𝒞*(sqrt(interval(2)*d)*norm(U,2) + rmin)*rmin
    
   if sup(τ)<0
    return 1, τ
   else
    return 0, τ
   end
end




  
  
  function PlotmT(c,T)
    #U0 is a sequence in 2D
    #a,b,c,d are the endpoints of the interval [a,b] × [c,d]
    y1 = 0.01:0.01:10
    y2 = 0.01:0.01:10
    m1=length(y1)
    m2=length(y2)
    U = complex(zeros(m1,m2))
    for b₁ = 1:m1
        for b₂ = 1:m2
          z = y1[b₁]+im*y2[b₂]
            U[b₁,b₂] = abs(sqrt(tanh(z)*(1+T*z^2)/z)-c)
        end
    end
    U = real.(U)
    mat"
    h = surf($y1,$y2,$U)
    set(h,'LineStyle','none')"
  end
  
  
  
  
  
  function no_eigen_outside(r1,r2,U1,U2,U0,rmin,c,T,ν)
  
    setprecision(100)
  
    fourier = space(U0)
    N = order(U0)
  
    λ1 = interval(component(real(U1),1)[1])
    λ2 = interval(component(real(U2),1)[1])

    r1 = interval.(Float64(inf(r1),RoundDown),Float64(sup(r1),RoundUp) )
    r2 = interval.(Float64(inf(r2),RoundDown),Float64(sup(r2),RoundUp) )
    λ1 = interval.(Float64(inf(λ1),RoundDown),Float64(sup(λ1),RoundUp) )
    λ2 = interval.(Float64(inf(λ2),RoundDown),Float64(sup(λ2),RoundUp) )
  
    # Construction of the quantities that do not change with λ
  
    D² = real.(project(Derivative(2), fourier_f, fourier_f,Complex{Interval{Float64}}))
    dd = diag(coefficients(D²))  ;   dd = sqrt.(-dd) ; d0=dd;  d0[N+1] = interval(1);
    dL = tanh.(dd).*(ones(2*N+1)+T*dd.^2)./d0 ;  dL[N+1] = interval(1) ; dL = sqrt.(dL)
    MT = LinearOperator(fourier_f,fourier_f,Diagonal(dL))
  
    if inf(T)>0
        L = MT - c*Id
        Linv = LinearOperator(fourier_f,fourier_f,Diagonal(ones(2*N+1)./diag(coefficients(L))))
        DG = 2*real.(project(Multiplication(U0),fourier_f, fourier_f,Complex{Interval{Float64}}))
    else
        L = c*Id - MT
        Linv = LinearOperator(fourier_f,fourier_f,Diagonal(ones(2*N+1)./diag(coefficients(L))))
        DG = -2*real.(project(Multiplication(U0),fourier_f, fourier_f,Complex{Interval{Float64}}))
    end
    
    # Construction of WT
    e0 = Sequence(fourier ,interval.(zeros(2*N+1)))
    e0[0] = interval(1)
   
  
  # Now that the above quantities are fixed, we prove that there is no eigenvalue between λ1 and λ2 and that there is no eigenvalue between λ2 and λmin. This implies that λ1 and λ2 are the smallest eigenvalues of DF(u)
  
  λmin = abs(1-c)-2*norm(U0,1)-rmin/(4*sqrt(ν)*abs(1-c))
  λ = interval(inf(λ2-r2))
  P = 1
  k=1
  
  #computation of the supremum of Zλ for all λ ∈ (0,λmin)
  # We first compute the supremum of Cλ
  if inf(T)>0
    Cλ = computation_Cλ(maximum([c c+λmin]),a,σ0,T)
  else
    Cλ = computation_Cλ(c-λmin,a,σ0,T)
  end
 
  
  # coefficients of the characteristic function
  char = Sequence(Fourier(2*N,π/d) ,interval.(big.(zeros(4*N+1))))
  for n = [-2N:-1;1:2N]  
    char[n] = real((-1)^n*(exp(im*n*π/db)-1)/(2*im*n*π))
  end
  char[0] = 1/(2*db)
  
  # Coefficients of cosh(2ax)
  Eb = Sequence(Fourier(2*N,π/d),interval.(big.(zeros(4*N+1))))
  for n = -2N:2N  
    Eb[n] = real((-1)^n*(1/(im*n/db*pi + 2ab) - 1/(im*n/db*pi-2*ab)))
  end
  Eb[0] = 1/ab
  Eb = Eb*(1-exp(-4*ab*db))/(4*db)

  
  
  # For all λ ≤0, a-λ ≥ a and so the maximal 𝒵λ is achieved for λ=0
  C1d = exp(-2*a*d)*( 2*sqrt(π)/(sqrt(4*a*d)*(1-exp(-2*a*d))) + 4/(1-exp(-2*a*d)) )
  DU0 = U0.*( (-N:N)*π/db )
  Elog =  abs( (coefficients(project(DU0*char,fourier))')*coefficients(DU0) )

 
  𝒵u0 = abs( (coefficients(project(U0fb*Eb,fourier))')*coefficients(U0fb) )

 
  
  #Computation of the supremum of 𝒵λ for all λ ≤ 0
  𝒵λ = 8*Cλ^2*(d*𝒵u0*(2/a+C1d) + 4*log(interval(2))*Elog )
  𝒵λ = sqrt(𝒵λ)
  
  𝒵λ = interval.(Float64(inf(𝒵λ),RoundDown),Float64(sup(𝒵λ),RoundUp) )

  #MV2 corresponds to the operator πᴺVi(π^{2N}-πᴺ)Viπᴺ = πᴺ(Vi*Vi)πᴺ - (πᴺViπᴺ)^2
  MV2 = interval(4)*project(Multiplication(U0*U0),fourier, fourier,Interval{Float64}) - interval(4)*project(Multiplication(U0),fourier, fourier,Interval{Float64})^2
  
 
  

    while  (P==1)&&(sup(λ)>=inf(λ1+r1))
        Lλ = L - interval(λ)*Id
        Li = ones(2*N+1)./diag(coefficients(Lλ))
  
        Df = Lλ + DG
        B = interval.(inv(mid.(Df).*mid.(Li)'))
        Badj = LinearOperator(fourier,fourier,coefficients(B)')
        n_BN = sqrt(opnorm( LinearOperator(coefficients(B*Badj)),2 ))   #computation of \|B^N_T\|_2
        
        if inf(T)>0
            WT = e0
            n_WT = 0
        else
            e02 = Sequence(CosFourier(2*N, π/d) ,interval.(zeros(2*N+1)))
            e02[0] = 1
            WT = interval.(mid.(project(Multiplication(mid.(e0 - 2/mid(c-λ)*U0)),fourier,Fourier(2*N, π/d),Interval{Float64}))\mid.(e02))
            MWT = project(Multiplication(WT*WT),fourier, fourier,Interval{Float64}) - project(Multiplication(WT),fourier, fourier,Interval{Float64})*project(Multiplication(WT),fourier, fourier,Interval{Float64})
    
            n_WT = sqrt(opnorm(LinearOperator(coefficients(D1.*(MWT).*D2')),2))   # computation of \|π_NW_Tπ^N\|_2
            MWT = Nothing
        end
    
        # computation of the norm of B
        norm_B = maximum([1 maximum([norm(WT,1) n_BN])+n_WT]) # Computation for an upper bound of the norm of B
  
        if inf(T)>0
            lT = maximum(vec((abs.(dL.-c-λ))))
        else
            lT = maximum(vec(abs.(((c-λ)*(dL.-(c-λ))./dL))))
        end 

  
        W2 = interval(2)*WT*U0 ;          
  
  #MW2 corresponds to the operator πᴺ(Vi*WT)(π^{3N}-πᴺ)(Vi*WT)πᴺ = πᴺ(Vi*WT*Vi*WT)πᴺ - (πᴺ(Vi*WT)πᴺ)^2 (recall that Wi = Vi*WT)
        MW2 = project(Multiplication(W2*W2),fourier, fourier,Interval{Float64}) - project(Multiplication(W2),fourier, fourier,Interval{Float64})^2
  
  
        if inf(T)>0
            Z14 = 0
        else
            Z14 = norm(e0 - WT*(e0-interval(2)/(c-λ)*U0),1)
        end
  
       
        # computation of each component of Z1
        Z11 =  opnorm(LinearOperator(coefficients(Id - B*(Df.*Li'))),2)^2
        Z11 = Z11 + opnorm(LinearOperator(coefficients((Li.*MW2.*Li'))),2)
        Z11 = sqrt(Z11)
  
        Z12 = sqrt(opnorm(LinearOperator(coefficients(B*MV2*Badj)),2))/lT
  
        
        Z13 = norm(W2,1)/lT
  
        Z1 = sqrt(Z11^2 + Z12^2 + Z13^2 + Z14^2)
  
        𝒞 = (1-norm_B*𝒵λ-Z1)/norm_B*(abs(1-c) -λ) - rmin/(4*sqrt(ν)*abs(1-c)) 
        if inf(𝒞) <= 0
            P=0
            return P==1
        end
  
        λ = λ - interval.(Float64(inf(𝒞),RoundDown),Float64(sup(𝒞),RoundUp) ) 
      
    end
    
    λ = λ1-r1
    k=2
    
    while  (P==1)&&(sup(λ)>=inf(λmin))
        
        Lλ = L - λ*Id
        Li = interval.(ones(2*N+1))./diag(coefficients(Lλ))
  
        Df = Lλ + DG
        B = interval.(inv(mid.(Df).*mid.(Li)'))
        Badj = LinearOperator(fourier,fourier,coefficients(B)')
        n_BN = sqrt(opnorm( LinearOperator(coefficients(B*Badj)),2 ))   #computation of \|B^N_T\|_2

        if inf(T)>0
          WT = e0
          n_WT = 0
      else
          e02 = Sequence(CosFourier(2*N, π/d) ,interval.(zeros(2*N+1)))
          e02[0] =interval(1)
          WT = interval.(mid.(project(Multiplication(mid.(e0 - interval(2)/mid(c-λ)*U0)),fourier,Fourier(2*N, π/d),Interval{Float64}))\mid.(e02))
          MWT = project(Multiplication(WT*WT),fourier, fourier,Interval{Float64}) - project(Multiplication(WT),fourier, fourier,Interval{Float64})*project(Multiplication(WT),fourier, fourier,Interval{Float64})
  
          n_WT = sqrt(opnorm(LinearOperator(coefficients(D1.*(MWT).*D2')),2))   # computation of \|π_NW_Tπ^N\|_2
          MWT = Nothing
      end
      
        norm_B = maximum([interval(1) maximum([norm(WT,1) n_BN])+n_WT]) # Computation for an upper bound of the norm of B
       
        if inf(T)>0
            lT = maximum(vec((abs.(dL.-(c+λ)))))
        else
            lT = maximum(vec(abs.(((c-λ)*(dL.-(c-λ))./dL))))
        end 
  
        W2 = interval(2)*WT*U0 ; 
  #MW2 corresponds to the operator πᴺ(Vi*WT)(π^{3N}-πᴺ)(Vi*WT)πᴺ = πᴺ(Vi*WT*Vi*WT)πᴺ - (πᴺ(Vi*WT)πᴺ)^2 (recall that Wi = Vi*WT)
  MW2 = project(Multiplication(W2*W2),fourier, fourier,Interval{Float64}) - project(Multiplication(W2),fourier, fourier,Interval{Float64})^2
  
  
        
        if inf(T)>0
            Z14 = 0
        else
            Z14 = norm(e0 - WT*(e0-interval(2)/(c-λ)*U0),1)
        end
  
       
        # computation of each component of Z1
        Z11 =  opnorm(LinearOperator(coefficients(Id - B*(Df.*Li'))),2)^2
        Z11 = Z11 + opnorm(LinearOperator(coefficients((Li.*MW2.*Li'))),2)
        Z11 = sqrt(Z11)
  
       
        Z12 = sqrt(opnorm(LinearOperator(coefficients(B*MV2*Badj)),2))/lT
  
        
        Z13 = norm(W2,1)/lT
  
        Z1 = sqrt(Z11^2 + Z12^2 + Z13^2 + Z14^2)
  
        𝒞 = (1-norm_B*𝒵λ-Z1)/norm_B*(abs(1-c) -λ) - rmin/(4*sqrt(ν)*abs(1-c)) 
        if inf(𝒞) <= 0
            P=0
            return P==1
        end
  
        λ = λ - interval.(Float64(inf(𝒞),RoundDown),Float64(sup(𝒞),RoundUp) ) 
      
    end
  
    return P==1
  
  
  end
  
  
  
  function computation_Cλ(c,a,σ0,T)
  
    setprecision(100)
  
  ##### Construction of the needed constants
  # Construction of ξ0 
    Ca = (1 + abs(cos(2*a)))/(1-abs(cos(2*a)))
  
  if inf(T)>0
  ξ0 = maximum([interval(1) 1/sqrt(T) c^2*Ca*4/T])
  ξ0 = maximum([3*T/(2*tanh(ξ0)) 2*c^2/(T*tanh(ξ0))])
  else
  ξ0 = interval(1)
  end
  
  if inf(T)>0
  K1 = 2*ξ0/(π*σ0) + 2*sqrt(ξ0)*(1+abs(c))/(π*σ0*sqrt(T)) + 2*(1/(3*T) + abs(c)/(4*T^(3/2)) + 2*c^2/T )/(π*sqrt(tanh(ξ0)*T*ξ0)) + abs(c)/(T*π)*(2+3*log(ξ0)) + 1/sqrt(2π*T) 
  K2 = Ca*abs(ξ0+im*a)/( 2π*σ0^2*(1-T*a^2)^2)*( 1/a^2 + T + Ca/a + Ca*T*abs(ξ0+im*a) ) + 2*Ca^2/π*( 2*(1+T)/sqrt(T*ξ0) + (2+a)/2*exp(-2ξ0))
  C2 = maximum([K2 K1*exp(a)])  
  else
  K1 = (2+4*exp(-2*ξ0))/(minimum([1 abs(c)^3])*π*sqrt(ξ0)*σ0) + 1/(π*σ0*abs(c)) + 2/(π*c^2) + 1/(π*abs(c)^3)*(2+3*log(3*ξ0)) + 1/(c^2*sqrt(2*π))
  K2 = Ca^(1/4)/π*( 1/(2*σ0^2)*( 2 + a^(-3/2)) + 1/(4*σ0^2*(1-abs(cos(2a)))^2*sqrt(a)) )
  C2 = maximum([K2 K1*exp(a)]) 
  end
  
  return C2
  
  end
  
  
  
  function trace_full(N,d)
    w = -N:N;
    p = interval(pi)
    d = interval(d)
    Y2 = (p/d)*interval.(w);
    Y1 = (p^2/d^2)*interval.(w).^2;
    X = interval.((-1).^w);
  
    f = [X';
        (X.*Y1)';
         (X.*Y2)']
    return f
  end
  
  
  function solve_linear(M,b,precis)
    setprecision(precis)
  
    x = interval.(big.(Float64.(mid.(M))\Float64.(mid.(b)))) 
    Minv = interval.(big.(inv(Float64.(mid.(M)))))
    N = size(b)[1]
    Id = interval.(1.0*I[1:N,1:N])
  
    Minv2 = interval.(Float64.(inf.(Minv),RoundDown),Float64.(sup.(Minv),RoundUp) )
    M2 = interval.(Float64.(inf.(M),RoundDown),Float64.(sup.(M),RoundUp) )
    Z1 = opnorm(LinearOperator(Id - Minv2*M2),Inf)
  
    if inf(interval(1)-Z1) >0
        Y0 = norm(Minv*(M*x-b),Inf)
        rmin = sup(Y0/inf(1-Z1))
        return interval.(inf.(x).-rmin,sup.(x).+rmin)
    else
        return Nan
    end
  end

  




  function newton_whitham(U, c, T, fourier)

    D² = project(Derivative(2), fourier, fourier,Float64)
    d2 = diag(coefficients(D²))
    d2 = sqrt.(-d2)
    d2[1] = 1.0
    d2 = (tanh.(d2)).*((1.0).+T*d2.^2)./d2
    d2[1] = (1.0)
    d2 =  sqrt.(d2)
    L = LinearOperator(fourier,fourier,Diagonal(d2)) - c*I
    l = diag(coefficients(L)) ;
    F = l.*U + project(U*U,fourier)

    n_F = norm(F,2)
    display(n_F)
    err = 1.0e-15
    k=0
    while (n_F > err)&&(k<10)
        DG = project(Multiplication(2*U),fourier, fourier,Float64)
        DF = L + DG
        U = U - DF\F
        
        F = l.*U + project(U*U,fourier)
        @show n_F = norm(F,2) # macro pour montrer le résultat du calcul
        k = k+1
    end
    
    return U
end









function check_constants(a,σ0,σ1,xmin,ymin,ν)

  
    K = 10
    f = interval(2)
    if T==0 #if T=0 and c>1, we have c>m_T so we want f such that |m0| ≤ f and so |c-mT| ≥ c-|mT| ≥ c-f
      for k=1:K
        f = f + interval(2)^(k+1)*abs(xmin+im*ymin)^k/(interval(factorial(k+1)))
      end
      f = sqrt((f+  interval(1)/interval(factorial(K+1))*interval(2)^(K+1)*abs(xmin+im*ymin)^K)*(interval(1)+T*abs(xmin+im*ymin)^2)/(interval(1)+cos(interval(2)*xmin)))
    
      if inf(c-f)>0
        display("the constants xmin and ymin satisfy the required conditions")
      else  
        display("the constants xmin and ymin do not satisfy the required conditions")
        Breakdown = BreakNow
        # if the proof fails, we can decrease xmin and ymin for it to work out
      end
    
    else #if T>0 and c<1, we have c<m_T so we want f such that |mT| ≥ f and so |c-mT| ≥ |mT|-c ≥ f-c
      for k=1:K
        f = f - interval(2)^(k+1)*abs(xmin+im*ymin)^k/(interval(factorial(k+1)))
      end
      f = sqrt((f- interval(1)/interval(factorial(K+1))*interval(2)^(K+1)*abs(xmin+im*ymin)^K)*(interval(1)+T*abs(xmin+im*ymin)^2)/(interval(2)+sin(interval(2)*xmin)))
    
      if inf(f-c)>0
        display("the constants xmin and ymin satisfy the required conditions")
      else
        display("the constants xmin and ymin do not satisfy the required conditions")
        Breakdown = BreakNow
        # if the proof fails, we can decrease xmin and ymin for it to work out
      end
    end
    
    x = a
    P1 = -interval(1.0)
    P2 = -interval(1.0)
    while (inf(P1)<=0)&&(inf(P2)<=0)
      if inf(T)>0
        x = x+interval(0.1)
        P1 = x^2 - sqrt((x^2+a^2)/(T^2*(cosh(2x)-interval(1))/(interval(1)+cosh(interval(2)*x))))*(c+σ0)^2
        P2 = (cosh(interval(2)x)-interval(1))/(cosh(interval(2)x)+interval(1))*T^2*x^4/(x^2+a^2)-interval(2)^4*c^4
      else
        x = x+interval(0.1)
        P1 = x^2 - interval(1)+abs(cos(interval(2)a))/(cosh(interval(2)x)*(c+σ0)^4*(interval(1)-abs(cos(interval(2)a))/cosh(interval(2)x)))
        P2 = interval(1)
      end
    end
    
    ######## We verify that the constants satisfy the desired inequalities
    di = 0.0001
    dj = 0.0001
    P,m1 = test_constants(T,ν,c,a,σ0,σ1,x,di,dj,xmin,ymin)
    
    
    
      return P
    end  

    




    


    function proof_soliton(W,N,d,c,T,a,σ0,σ1)



        fourier = CosFourier(N, π/d)
        fourier0 = CosFourier(N0,pi/d)
        W = Sequence(fourier,coefficients(W))
        
        # Choose ν = T if T>0 and ν = 1 if T=0
        if inf(T)>0
            ν = T
            νb = Tb
        else
            ν = interval(4)/interval(π)^2
            νb = interval(big(4))/interval(big(π))^2
        end
        
        
        ##################### VERIFICATION OF THE CONSTANTS a, σ0 and σ1 AS DESCRIBED IN APPENDIX 7.2 #########################################
        
        # We start by some analysis around the point 0. We prove that |m_T(ξ)-c|≥ σ0 for all ξ = x+iy and |x|≤xmin, |y|≤ymin
        # Then we use that m_T(0)=1 and the series extension of |tanh(ξ)/ξ| = |1-exp(-2ξ)/(ξ+ξexp(-2ξ))| at zero to find that |m_T(ξ)-c| \geq ||m_T|-|c|| ≥ |1-c|- ϵ and ϵ is given by f below. In particular, if xmin and ymin are small enough, then |m_T(ξ)-c|≥ σ0 will be satisfied and we display a verified message. Otherwise the code breaks.
        
        # xmin = interval(0.45)
        # ymin = interval(0.05)
        
        xmin = interval(0.1)
        ymin = interval(0.05)
        
        ######################################################################
        
        #Construction of the fourier series of cosh(2ax) on Ω
        fourierE = CosFourier(2*N, interval(π)/d)
        
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
        D² = project(Derivative(2), fourier0, fourier0,Interval{Float64})
        Lν = Id - ν*D²
        
        # Construction of the operator L
        d2 = diag(coefficients(D²))  ;   dd = sqrt.(-d2) ; d0=dd;  d0[1] = interval(1);
        dL = tanh.(dd).*(ones(N0+1)+T*dd.^2)./d0 ;  dL[1] = interval(1) ; dL = sqrt.(dL)
        MT = LinearOperator(fourier0,fourier0,Diagonal(dL))
        
        L = Lν*(MT - c*Id)
        Linv = LinearOperator(fourier0,fourier0,Diagonal(interval.(ones(N0+1))./diag(coefficients(L))))
        
        D1 = convert(Vector{Interval{Float64}},interval.(exp2cos(N0)))
        D2 = interval.(ones((N0+1)))./D1
        
        ########### Projection on zero trace functions (Section 3.1) #############################
        setprecision(100)
        D = interval.(big.((1:N+1).^4))
        S = trace(N); C = S' ;
        S =  interval.(big.(S)) ;  C =  interval.(big.(C))
        
        V = interval.(big.(coefficients(W)))
        
        W0 = (D.*C)*solve_linear(S*(D.*C),vec(S*V),100)   # the function solve_linear solves a linear system rigorously
        display("U0 created")
        #
        U0b = interval.(big.(W)) - Sequence(fourier,vec(W0))
        U0 = interval.(Float64.(inf.(U0b),RoundDown),Float64.(sup.(U0b),RoundUp) )
        
        V2 = interval(2)*U0
        V1 = -interval(4)*ν*π/d*(0:N).*U0
        V0 = -interval(2)*ν*π^2/d^2*((0:N).^2).*U0
        
        
        ### Computation of epsilon (Assumption 2)
        if T==0
          ϵ = c/interval(2) - norm(U0,1)
          display("value of ϵ")
          display(inf(ϵ))
          if inf(ϵ)<=0
            display("epsilon is negative, the proof will fail")
            Breakdown = BreakNow
          end
        end
        
        
        ############ Computation of A (Section 3.2) ########################################
        
        DG = Lν*project(Multiplication(V2),fourier0, fourier0,Interval{Float64})
        B = interval.(inv(I + mid.(DG)*mid.(Linv)))
        Badj = LinearOperator(fourier0,fourier0,coefficients(B)')
        
        # Construction of WT
        e0 = Sequence(fourier0 ,interval.(zeros(N0+1)))
        e0[0] = interval(1)
        if inf(T)>0
            WT = e0
        else
          e02 = Sequence(CosFourier(2*N0, π/d) ,interval.(zeros(2*N0+1)))
          e02[0] = interval(1)
          WT = interval.(mid.(project(Multiplication(mid.(e0 - 1/mid(c)*V2)),fourier0, CosFourier(2*N0, π/d),Interval{Float64}))\mid.(e02))
        end
        
        # computation of the norm of B
        MWT = project(Multiplication(WT*WT),fourier0, fourier0,Interval{Float64}) - project(Multiplication(WT),fourier0, fourier0,Interval{Float64})^2
        
        
        n_WT = sqrt(opnorm(LinearOperator(coefficients(D1.*(MWT).*D2')),2))   # computation of \|π_NW_Tπ^N\|_2
        n_BN = (opnorm(LinearOperator(coefficients(D1.*((B*Badj)^8).*D2')),2))^(interval(1)/interval(16))   #computation of \|B^N_T\|_2
        norm_B = maximum([interval(1) maximum([norm(WT,1) n_BN])+n_WT])  #computation of the norm of B
        
        display("norm of B")
        display(norm_B)
        
        MWT = Nothing
        
        ################# Computation of Z1 ####################################
        
        W0 = WT*V0 ;              
        W1 = WT*V1 ;  
        W2 = WT*V2 ;          
        
        if inf(T)>0
          L0 = Linv ;             l0 = (interval(1)+ν*(interval(N0+1)*π/d)^2)*(maximum(vec((abs.(dL.-c)))))
          L1 = Linv*sqrt.(-D²) ;  l1 = (interval(1)+ν*(interval(N0+1)*π/d)^2)*(maximum(vec((abs.(dL.-c)))))/(interval(N0+1)*π/d)
          L2 = Linv*Lν ;          l2 = maximum(vec((abs.(dL.-c))))
        else
          L0 = Linv ;             l0 = (interval(1)+ν*(interval(N0+1)*π/d)^2)*(maximum(vec((abs.(dL.-c)))))
          L1 = Linv*sqrt.(-D²) ;  l1 = (interval(1)+ν*(interval(N0+1)*π/d)^2)*(maximum(vec((abs.(dL.-c)))))/(interval(N0+1)*π/d)
          L2 = Linv*Lν + interval(1)/c*Id ;  l2 = maximum(vec(abs.((c*(dL.-c))./dL)))
        end
        
        #MWi corresponds to the operator πᴺ(Vi*WT)(π^{3N}-πᴺ)(Vi*WT)πᴺ = πᴺ(Vi*WT*Vi*WT)πᴺ - (πᴺ(Vi*WT)πᴺ)^2 (recall that Wi = Vi*WT)
        MW0 = project(Multiplication(W0*W0),fourier0, fourier0,Interval{Float64}) - project(Multiplication(W0),fourier0, fourier0,Interval{Float64})^2
        MW1 = project(Multiplication(W1*W1),fourier0, fourier0,Interval{Float64}) - project(Multiplication(W1),fourier0, fourier0,Interval{Float64})^2
        MW2 = project(Multiplication(W2*W2),fourier0, fourier0,Interval{Float64}) - project(Multiplication(W2),fourier0, fourier0,Interval{Float64})^2
        
        #MVi corresponds to the operator πᴺVi(π^{2N}-πᴺ)Viπᴺ = πᴺ(Vi*Vi)πᴺ - (πᴺViπᴺ)^2
        MV0 = project(Multiplication(V0*V0),fourier0, fourier0,Interval{Float64}) - project(Multiplication(V0),fourier0, fourier0,Interval{Float64})^2
        MV1 = project(Multiplication(V1*V1),fourier0, fourier0,Interval{Float64}) - project(Multiplication(V1),fourier0, fourier0,Interval{Float64})^2
        MV2 = project(Multiplication(V2*V2),fourier0, fourier0,Interval{Float64}) - project(Multiplication(V2),fourier0, fourier0,Interval{Float64})^2
        
        # computation of each component of Z1
        Z11 =  opnorm(LinearOperator(coefficients(D1.*(Id - B*(Id+DG*Linv)).*D2')),2)^2
        Z11 = Z11 + ( sqrt(opnorm(LinearOperator(coefficients(D1.*(L0*MW0*L0).*D2')),2)) + sqrt(opnorm(LinearOperator(coefficients(D1.*(L1*MW1*L1).*D2')),2)) + sqrt(opnorm(LinearOperator(coefficients(D1.*(L2*MW2*L2).*D2')),2)) )^2
        Z11 = sqrt(Z11)
        
        display("value of Z11")
        display(Z11)
        
        Z12 = sqrt(opnorm(LinearOperator(coefficients(D1.*B*MV0*Badj.*D2')),2))/l0 + sqrt(opnorm(LinearOperator(coefficients(D1.*B*MV1*Badj.*D2')),2))/l1 + sqrt(opnorm(LinearOperator(coefficients(D1.*B*MV2*Badj.*D2')),2))/l2
        
        display("value of Z12")
        display(Z12)
        
        Z13 = norm(W0,1)/l0 + norm(W1,1)/l1 + norm(W2,1)/l2
        
        if inf(T)>0
          Z14 = 0
        else
          Z14 = norm(e0 - WT*(e0-1/c*V2),1)
        end
        
        Z1 = sqrt(Z11^2 + Z12^2 + Z13^2 + Z14^2)
        
        display("value of Z1")
        display(Z1)
        
        
        # To gain computational memory, we decided to cut truncate the number of coefficients of the multiplication operator associated to U0.
        # Indeed, we use a truncation of size N0 instead of N. This introduces an extra bound in the computation of 𝒵1, which we denote Z1inf
        # and accounts for such a truncation. Specifically it provides an upper bound for \|BΛν(U0 - π^N0 U0)L^{-1}Λν^{-1}\|_2  
        
        Z1inf  = (opnorm(D1.*B*project(Multiplication((V0-project(V0,fourier0))^2),fourier0, fourier0,Interval{Float64})^2*Badj.*D2',2) + norm(V0-project(V0,fourier0),1)^2)^(1/2) 
        + (opnorm(D1.*B*project(Multiplication((V1-project(V1,fourier0))^2),fourier0, fourier0,Interval{Float64})^2*Badj.*D2',2)+norm(V0-project(V0,fourier0),1)^2)^(1/2) 
        +  (opnorm(D1.*B*project(Multiplication((V2-project(V2,fourier0))^2),fourier0, fourier0,Interval{Float64})^2*Badj.*D2',2)+norm(V0-project(V0,fourier0),1)^2)^(1/2)
        
        Z1inf = 1/σ0*Z1inf
        
        
        MV0 = Nothing ; MV1 = Nothing ; MV2 = Nothing; MW0 = Nothing; MW1 = Nothing; MW2 = Nothing
        
        ######################## Computation of Zu ##############################
        
        
          ##### Construction of the needed constants
        # Construction of ξ0 
        Ca = (interval(1) + abs(cos(interval(2)*a)))/(interval(1)-abs(cos(interval(2)*a)))
        
        if inf(T)>0
          ξ0 = maximum([interval(1) interval(1)/sqrt(T) c^2*Ca*interval(4)/T])
          ξ0 = maximum([interval(3)*T/(interval(2)*tanh(ξ0)) interval(2)*c^2/(T*tanh(ξ0))])
        else
          ξ0 = interval(1)
        end
        
        # Constants for the bound Y0
        if inf(T)>0
          a0 = minimum([interval(1)/sqrt(T) interval(1)/interval(2)*π])
          a0b = minimum([interval(1)/sqrt(Tb) interval(big(1))/interval(big(2))*π])
        else
          a0 = interval(1.5)
          a0b = interval(big(1.5))
        end
        
        
        # Notice that max_{s>0}min{√2 + √s, √2/(1-exp(-πs)) } ≤ max{1+√2, √2/(1-exp(-π))}, so we can compute CY0 as follows (cf. Lemma 4.1)
        pi_i = interval(π)
        if T==0
          CY0 = interval(2.28)*interval(0.5)*maximum([interval(1)+sqrt(interval(2)) sqrt(interval(2))/(interval(1)-exp(-pi_i))])
        else
          CY0 =maximum([interval(1)+sqrt(interval(2)) sqrt(interval(2))/(interval(1)-exp(-pi_i))])*interval(1)/(sqrt(pi_i)*T^(1/4))*(interval(2)/sqrt(interval(1)+a0*sqrt(T)) + interval(1)/sqrt(interval(2)))
        end
        
        # Other constants of Lemma 4.1
        C0 = interval(1)/(π*σ0*(interval(1)-ν*a^2)) + interval(1)/(π*ν*σ0)
        
        if inf(T)>0
          C1 = interval(1)/(interval(2)*π)*( interval(2)*(interval(1)+a)/(σ0*(interval(1)-ν*a^2)) + 4*(interval(1)+a)/(σ1*sqrt(T)*ν) )
          K1 = interval(2)*ξ0/(π*σ0) + interval(2)*sqrt(ξ0)*(interval(1)+abs(c))/(π*σ0*sqrt(T)) + interval(2)*(interval(1)/(3*T) + abs(c)/(4*T^(3/2)) + interval(2)*c^2/T )/(π*sqrt(tanh(ξ0)*T*ξ0)) + abs(c)/(T*π)*(interval(2)+3*log(ξ0)) + interval(1)/sqrt(interval(2)*π*T) 
          K2 = Ca*abs(ξ0+im*a)/( interval(2)*π*σ0^2*(interval(1)-T*a^2)^2)*( interval(1)/a^2 + T + Ca/a + Ca*T*abs(ξ0+im*a) ) + interval(2)*Ca^2/π*( interval(2)*(interval(1)+T)/sqrt(T*ξ0) + (interval(2)+a)/interval(2)*exp(-interval(2)ξ0))
          C2 = maximum([K2 K1*exp(a)])  
        else
          C1 = interval(1)/(interval(2)*abs(c)*ν) + Ca^(interval(1)/interval(4))*(interval(1)+sqrt(a))/(π*abs(c)*σ0)*( interval(1)/(interval(1)-ν*a^2) + interval(2)/ν )
          K1 = (interval(2)+4*exp(-interval(2)*ξ0))/(minimum([interval(1) abs(c)^3])*π*sqrt(ξ0)*σ0) + interval(1)/(π*σ0*abs(c)) + interval(2)/(π*c^2) + interval(1)/(π*abs(c)^3)*(interval(2)+3*log(3*ξ0)) + interval(1)/(c^2*sqrt(interval(2)*π))
          K2 = Ca^(interval(1)/4)/π*( interval(1)/(interval(2)*σ0^2)*( interval(2) + a^(-3/2)) + interval(1)/(4*σ0^2*(interval(1)-abs(cos(interval(2)a)))^2*sqrt(a)) )
          C2 = maximum([K2 K1*exp(a)]) 
        end
        
        
        
        
        
        
        ###################### Computation of 𝒵u #####################################
        
        #### Computation of the constants C(d) and C1(d) in Lemmas 4.2 and 4.6
        
        C1d = exp(-interval(2)*a*d)*( interval(2)*sqrt(π)/(sqrt(interval(4)*a*d)*(interval(1)-exp(-interval(2)*a*d))) + interval(4)/(interval(1)-exp(-interval(2)*a*d)) )
        Cd = exp(-interval(2)*a*d)*( interval(4)*d + interval(4)*exp(-a*d)/(a*(interval(1)-exp(-interval(3)*a*d/interval(2)))) + interval(2)/(a*(interval(1)-exp(-a*d*interval(2)))))
        
        ### Inner product (V0,E*V0)
        D1N = convert(Vector{Interval{Float64}},interval.(exp2cos(N)))
        D2N = interval.(ones((N+1)))./D1N
        PS0 = abs( (coefficients(D1N.*project(V0*E,fourier))')*coefficients(D1N.*V0) )
        
        ### Inner product (V1,Efull*V1)
        k = (-N:N)*π/db
        fourier_f = Fourier(N,π/db)
        V1f = interval(4)*νb*coefficients(U0b) ; V1f = [reverse(V1f[2:N+1]); V1f[1] ; V1f[2:N+1]].*k ; V1f = Sequence(fourier_f,V1f)
        Ef = coefficients(Eb); Ef = [reverse(Ef[2:2*N+1]); Ef[1] ; Ef[2:2*N+1]] ; Ef = Sequence(Fourier(2*N,π/db),Ef)
        PS1 = abs( (coefficients(project(V1f*Ef,fourier_f))')*coefficients(V1f) )
        
        ### Inner product (V2,E*V2)
        PS2 = interval(8)*abs( (coefficients(project(U0b*Eb,fourier))')*coefficients(U0b) )
        
        # Computation of ∫|v2'|^2 from d-1 to d
        char = Sequence(Fourier(2*N,π/db) ,interval.(big.(zeros(4*N+1))))
        for n = [-2N:-1;1:2N]  
            char[n] = real(interval(big((-1)^n))*(exp(im*n*π/db)-interval(big(1)))/(interval(big(2))*im*n*π))
        end
        char[0] = interval(big(1))/(interval(big(2))*db)
        Elog =  abs( (coefficients(project(V1f*char,fourier_f))')*coefficients(V1f) )
        
        # Computation of the Zu_i
        Zu0 = interval(4)*d*C0^2/a*PS0 + interval(2)*d*Cd*C0^2*PS0
        Zu1 = interval(4)*d*C1^2/a*PS1 + interval(2)*d*Cd*C1^2*PS1
        Zu2 = interval(2)*C2^2*( interval(2)*d*PS2/a + interval(4)*log(interval(2))*Elog ) + interval(2)*d*C1d*C2^2*PS2
        
        𝒵u = norm_B*sqrt( Zu0 + Zu1 + Zu2 )
        
        display("value of 𝒵u")
        display(𝒵u)
        
        ######################## Computation of ℨ1 ##############################
        
        
        𝒵1 = Z1 + 𝒵u + Z1inf
        display("value of 𝒵1")
        display(𝒵1)
        
        
        ################### Computation of Z2 ######################
        if (T==0)||(inf(T)>1/3)
          κT = interval(1)/(abs(interval(1)-c)*two_norm_inverse_l(d,ν,c,T))
        else
          κT = interval(1)/(sqrt(ν)*σ0^2)
        end
        𝒵2 = interval(2)*κT*norm_B
        
        display("value of 𝒵2")
        display(𝒵2)
        
        
        #################### Computation of Y0 ##################
        
        Cd = exp(-interval(2)*a0*d)*( interval(4)*d + interval(4)*exp(-a0*d)/(a0*(interval(1)-exp(-interval(3)*a0*d/interval(2)))) + interval(2)/(a0*(interval(1)-exp(-a0*d*interval(2)))))
        
        D1 = convert(Vector{Interval{BigFloat}},interval.(big.(exp2cos(N))))
        D2 = interval.(big.(ones((N+1))))./D1
        D² = project(Derivative(2), CosFourier(N, π/db), CosFourier(N, π/db),Interval{BigFloat})
        Lνb =  νb*D²
        Lνb = diag(interval(big(1)).-coefficients(Lνb)).*2
        
        
        Eb = Sequence(fourierE ,interval.(big.(zeros(2*N+1))))
        for n = 1:2N  
            Eb[n] = real(interval(big((-1)^n))*(interval(big(1))/(im*interval(big(n))/db*interval(big(pi)) + interval(big(2))a0b) - interval(big(1))/(im*interval(big(n))/db*interval(big(pi))-interval(big(2))*a0b)))
        end
        Eb[0] = interval(big(1))/a0b
        Eb = Eb*(interval(big(1))-exp(-interval(big(4))*a0b*db))/(interval(big(4))*db)
        
        
        ### Inner product part 
        PSY0 = abs( (coefficients(D1.*project(((Lνb.*U0b)*Eb),fourier))')*coefficients(D1.*(Lνb.*U0b)) )
        
        Yu1 = sqrt(interval(4)*d^2*CY0*PSY0) ;
        
        display("bound Yu1")
        display(Yu1)
        
        D² = project(Derivative(2), CosFourier(2N, π/db), CosFourier(2N, π/db),Interval{Float64})
        Lν2 = Id - ν*D²
        Lν2 = diag(coefficients(Lν2))
        
        D² = project(Derivative(2), fourier, fourier,Interval{Float64})
        # Construction of the operator L
        d2 = diag(coefficients(D²))  ;   dd = sqrt.(-d2) ; d0=dd;  d0[1] = interval(1);
        dL = tanh.(dd).*(ones(N+1)+T*dd.^2)./d0 ;  dL[1] = interval(1) ; dL = sqrt.(dL)
        MT = LinearOperator(fourier,fourier,Diagonal(dL))
        L = Lν2[1:N+1].*diag(coefficients((MT - c*Id)))
        
        FU0 = L.*U0 + Lν2.*(U0*U0)
        Y0 = sqrt(interval(2)*d)*sqrt( norm(B*project(FU0,fourier),2)^2 + norm(WT*FU0-project(WT*FU0,fourier),2)^2)
        
        𝒴0 = Y0 + Yu1*sqrt(interval(1) + norm_B^2*(interval(1)+Cd))
        display("value of 𝒴0")
        display(𝒴0)
        
        α = interval( sup( interval(1)/interval(2) + 𝒴0*𝒵2/(interval(1)-𝒵1)^2 ) )
        r = interval(2)*α*𝒴0/(interval(1)-𝒵1)
        
        if inf(1- 𝒵1)>0
          if sup(𝒵2*r+𝒵1) < 1
              if inf(1/2*𝒵2*r^2 - (1-𝒵1)*r + 𝒴0) < 0
                  display("The computer-assisted proof was successful")
                  display("Radius of contraction :")
                  display(r)
              else
              display("failure: discriminant is negative")
              end
          else
              display("r is too big")
          end
        else
          display("failure: 1-z > 0")
        end
        
        
        ######### Regularity for the solution (using Proposition 4.1)
        
        if T==0
          if sup(r) <= inf(interval(4)*ϵ*sqrt(ν)*minimum([abs(c-interval(1)),abs(c)]))
            display("the solution is infinitely differentiable")
          else
            display("the solution is H^(3/2) (continuous at least), but could not prove more regularity")
          end
        end
        
        return r, Yu1, norm_B/(1- 𝒵1)
        
        end    
  


        function F_W(U,c,T)
          fourier = space(U)
          D² = project(Derivative(2,0), fourier, fourier,Float64)
      
          for i=0:order(U)[2]
              D²[(0,i),(0,i)] = -1
          end
      
          d2 = diag(coefficients(D²))
          d2 = sqrt.(-d2)
          d2 = (tanh.(d2)).*((1.0).+T*d2.^2)./d2
          d2 =  sqrt.(d2)
      
          L = LinearOperator(fourier,fourier,Diagonal(d2))
      
          for i=0:order(U)[2]
              L[(0,i),(0,i)] = 1
          end
          
          F = L*U-project(c*U,fourier) + project(U*U,fourier)
          return F
      end
      
      
      
      
      function F_W(U::Vector{Vector{Interval{Float64}}},c,T)
          
          fourier = fourierC_
          D² = project(Derivative(2), fourier, fourier,Interval{Float64})
      
          ν_ = interval(4)/interval(π)^2
          D² = project(Derivative(2), fourierC_, fourierC_,Interval{Float64})
          Lν = diag(coefficients(Id - ν_*D²))
          d2 = diag(coefficients(D²))
          d2 = sqrt.(-d2)
          d2[1] = interval(1)
          d2 = (tanh.(d2)).*(interval(1).+T*d2.^2)./d2
          d2 =  sqrt.(d2)
          d2[1] = interval(1)
      
          L = LinearOperator(fourier,fourier,Diagonal(d2))
      
          F = [coefficients(Lν.*((L-c[i]*Id)*Sequence(fourier,U[i])+project(Sequence(fourier,U[i])^2,fourier))) for i=1:length(U)]
          
          return F
      end
      
      
      
      function F_W(U::Vector{Vector{Float64}},c,T)
          
          fourier = fourierC_
           ν_ = interval(4)/interval(π)^2
          D² = project(Derivative(2), fourierC_, fourierC_,Interval{Float64})
          Lν = diag(coefficients(Id - ν_*D²))
      
          d2 = diag(coefficients(D²))
          d2 = sqrt.(-d2)
          d2[1] = interval(1)
          d2 = (tanh.(d2)).*(interval(1).+T*d2.^2)./d2
          d2 =  sqrt.(d2)
          d2[1] = interval(1)
      
          L = LinearOperator(fourier,fourier,Diagonal(d2))
      
          F = [coefficients(Lν.*((L-c[i]*Id)*Sequence(fourier,U[i])+project(Sequence(fourier,U[i])^2,fourier))) for i=1:length(U)]
          
          return F
      end
      
      
      
      function F_Whitham(U,c,T)
          fourier = space(U)
          D² = project(Derivative(2), fourier, fourier,Float64)
      
      
          d2 = diag(coefficients(D²))
          d2 = sqrt.(-d2)
          d2[1] = 1
          d2 = (tanh.(d2)).*((1.0).+T*d2.^2)./d2
          d2 =  sqrt.(d2)
          d2[1] = 1
      
          L = LinearOperator(fourier,fourier,Diagonal(d2))
          
          F = L*U-c*U + project(U*U,fourier)
          return F
      end
      
      
      
      function DF_W(U,c::Sequence{TensorSpace{Tuple{CosFourier{Float64},Chebyshev}},Vector{Float64}},T)
          fourier = space(U)
          D² = project(Derivative(2,0), fourier, fourier,Float64)
      
          for i=0:order(U)[2]
              D²[(0,i),(0,i)] = -1
          end
      
          d2 = diag(coefficients(D²))
          d2 = sqrt.(-d2)
          d2 = (tanh.(d2)).*((1.0).+T*d2.^2)./d2
          d2 =  sqrt.(d2)
      
          L = LinearOperator(fourier,fourier,Diagonal(d2))
      
          for i=0:order(U)[2]
              L[(0,i),(0,i)] = 1
          end
          
          L = L - project(Multiplication(c),fourier, fourier,Float64)
          DF = L + 2*project(Multiplication(U),fourier, fourier,Float64)
          return DF
      end
      
      
      function DF_W(U,c::Float64,T)
          fourier = space(U)
          D² = project(Derivative(2), fourier, fourier,Float64)
          d2 = diag(coefficients(D²))
          d2 = sqrt.(-d2)
          d2[1] = 1.0
          d2 = (tanh.(d2)).*((1.0).+T*d2.^2)./d2
          d2[1] = (1.0)
          d2 =  sqrt.(d2)
          L = LinearOperator(fourier,fourier,Diagonal(d2)) - c*I
          DF = L + 2*project(Multiplication(U),fourier, fourier,Float64)
          return DF
      end
      
      function DF_W(U,c,T,fourier)
         
          D² = project(Derivative(2), fourier, fourier,Interval{Float64})
          d2 = diag(coefficients(D²))
          d2 = sqrt.(-d2)
          d2[1] = 1.0
          d2 = (tanh.(d2)).*((1.0).+T*d2.^2)./d2
          d2[1] = (1.0)
          d2 =  sqrt.(d2)
          L = LinearOperator(fourier,fourier,Diagonal(d2)) - c*Id
          DF = L + 2*project(Multiplication(Sequence(fourier,coefficients(U))),fourier, fourier,Interval{Float64})
          return DF
      end
      
      
      function norm_cheb(U::VecOrMat{<:Sequence},fourier)
          
          F = interval(0)
          for i=1:length(U)
              F = F + norm(U[i],2)
          end
          return F
      end
      
      function norm_cheb(U::VecOrMat{<:Vector},fourier)
          
          F = interval(0)
          for i=1:length(U)
              F = F + norm(Sequence(fourier,U[i]),2)
          end
          return F
      end
      
      
      
      
      function project_trace(U)
          N = order(U)
          U[0] = interval(2)*sum(coefficients(U)[3:end].*(interval(2:N).^2 .-interval(1) ).*interval.((-1).^(2:N)))
          U[1] = sum(coefficients(U)[3:end].*(interval(2:N).^2 ).*interval.((-1).^(2:N)))
          return U
      end    
      
      
      
      
      
      function newton_whitham(U, c, T, fourier,kmax)
      
          D² = project(Derivative(2), fourier, fourier,Float64)
          d2 = diag(coefficients(D²))
          d2 = sqrt.(-d2)
          d2[1] = 1.0
          d2 = (tanh.(d2)).*((1.0).+T*d2.^2)./d2
          d2[1] = (1.0)
          d2 =  sqrt.(d2)
          L = LinearOperator(fourier,fourier,Diagonal(d2)) - c*I
          l = diag(coefficients(L)) ;
          F = l.*U + project(U*U,fourier)
      
          n_F = norm(F,2)
          display(n_F)
          err = 1.0e-15
          k=0
          while (n_F > err)&&(k<kmax)
              DG = project(Multiplication(2*U),fourier, fourier,Float64)
              DF = L + DG
              U = U - DF\F
              
              F = l.*U + project(U*U,fourier)
              @show n_F = norm(F,2) # macro pour montrer le résultat du calcul
              k = k+1
          end
          
          return U
      end
      
      function cheb2grid(x::VecOrMat{<:Sequence}, N_fft)
          vals = fft.(x, N_fft)
          return [real.(getindex.(vals, i)) for i ∈ eachindex(vals[1])]
      end
      
      function grid2cheb(x_fft,Nc)
          return [rifft!(complex.([x_fft[j][n] for j = eachindex(x_fft)]), Chebyshev(Nc)) for n=eachindex(x_fft[1])]
      end    

      function grid2cheb0(B_grid,Nc)
        F = [rifft!(complex.([B_grid[k][i,j] for k = 1:length(B_grid)]), Chebyshev(Nc)) for i = axes(B_grid,1), j = axes(B_grid,2)]
        return [real.(getindex.(F, i)) for i ∈ eachindex(coefficients(F[1])).-1] 
      end    
    
      
      
      
      
      
      function computation_Z1(B,WT,U,c0,c1,c,N)

        fourier0 = fourierC_
          ν = interval(4)/interval(π)^2
          dn = interval(2)*interval.(ones(N_+1)) ; dn[1] = interval(1)
          V2 = interval(2)*U
          V1 = [-interval(4)*ν*π/d*(0:N).*U[i] for i = eachindex(U)]
          V0 = [-interval(2)*ν*π^2/d^2*((0:N).^2).*U[i] for i = eachindex(U)]

          D² = project(Derivative(2), fourierC_, fourierC_,Interval{Float64})
          Lν = diag(coefficients(Id - ν*D²))
      
          Badj =  [LinearOperator(fourier0,fourier0,coefficients(B[i])') for i = eachindex(B)]
          L = Lν.*(M0 - c0*Id)
          L0 = Lν.*(M0 - interval(0.5)*(c0+c1)*Id)
          Linv = interval.(ones(N+1))./diag(coefficients(L))
          L0inv = interval.(ones(N+1))./diag(coefficients(L0))
          
          W0 = [WT[i]*V0[i] for i = eachindex(WT)]               
          W1 = [WT[i]*V1[i] for i = eachindex(WT)] 
          W2 = [WT[i]*V2[i] for i = eachindex(WT)]  
      
            L0 = Linv ;             l0 = (interval(1)+ν*(interval(N+1)*π/d)^2)*(maximum(vec((abs.(dL.-c0)))))
            L1 = Linv.*diag(coefficients(sqrt.(-D²))) ;  l1 = (interval(1)+ν*(interval(N+1)*π/d)^2)*(maximum(vec((abs.(dL.-c0)))))/(interval(N+1)*π/d)
            L2 = Linv.*Lν + interval(1)/c0*interval.(ones(N+1)) ;  l2 = maximum(vec(abs.((c0*(dL.-c0))./dL)))
          
            DG = [Lν.*project(Multiplication(V2[i]),fourier0, fourier0,Interval{Float64}) for i=eachindex(V2)]
          #MWi corresponds to the operator πᴺ(Vi*WT)(π^{3N}-πᴺ)(Vi*WT)πᴺ = πᴺ(Vi*WT*Vi*WT)πᴺ - (πᴺ(Vi*WT)πᴺ)^2 (recall that Wi = Vi*WT)
          MW0 = [project(Multiplication(W0[i]*W0[i]),fourier0, fourier0,Interval{Float64}) - project(Multiplication(W0[i]),fourier0, fourier0,Interval{Float64})^2 for i = eachindex(W0)]
          MW1 = [project(Multiplication(W1[i]*W1[i]),fourier0, fourier0,Interval{Float64}) - project(Multiplication(W1[i]),fourier0, fourier0,Interval{Float64})^2 for i = eachindex(W0)]
          MW2 = [project(Multiplication(W2[i]*W2[i]),fourier0, fourier0,Interval{Float64}) - project(Multiplication(W2[i]),fourier0, fourier0,Interval{Float64})^2 for i = eachindex(W0)]
          
          
          # computation of each component of Z1

          Z11 = [coefficients(D1.*(L0inv.*(M0-c[i]*Id) - B[i]*(L0inv.*(M0-c[i]*Id)+DG[i].*L0inv')).*D2') for i = eachindex(c)]
          Z11 = grid2cheb0(Z11,N_) ; Z11 = norm(dn.*[opnorm(LinearOperator(Z11[i]),2) for i = eachindex(Z11)],1)
          
          ZW0 = [coefficients(D1.*(L0.*MW0[i].*L0').*D2') for i=eachindex(MW0)] ; ZW0 = grid2cheb0(ZW0,N_) 
          ZW1 = [coefficients(D1.*(L1.*MW1[i].*L1').*D2') for i=eachindex(MW1)] ; ZW1 = grid2cheb0(ZW1,N_) 
          ZW2 = [coefficients(D1.*(L2.*MW2[i].*L2').*D2') for i=eachindex(MW2)] ; ZW2 = grid2cheb0(ZW2,N_)

          Z11 = Z11 + norm(dn.*[sqrt(opnorm(LinearOperator(ZW0[i]),2)) for i = eachindex(ZW0)],1) 
          Z11 = Z11 + norm(dn.*[sqrt(opnorm(LinearOperator(ZW1[i]),2)) for i = eachindex(ZW1)],1) 
          Z11 = Z11 + norm(dn.*[sqrt(opnorm(LinearOperator(ZW2[i]),2)) for i = eachindex(ZW2)],1) 
          
          display("value of Z11")
          display(Z11) ; MW0 = Nothing; MW1 = Nothing; MW2 = Nothing

        #MVi corresponds to the operator πᴺVi(π^{2N}-πᴺ)Viπᴺ = πᴺ(Vi*Vi)πᴺ - (πᴺViπᴺ)^2
        MV0 = [project(Multiplication(V0[i]*V0[i]),fourier0, fourier0,Interval{Float64}) - project(Multiplication(V0[i]),fourier0, fourier0,Interval{Float64})^2 for i = eachindex(V0)]
        MV1 = [project(Multiplication(V1[i]*V1[i]),fourier0, fourier0,Interval{Float64}) - project(Multiplication(V1[i]),fourier0, fourier0,Interval{Float64})^2 for i = eachindex(V1)]
        MV2 = [project(Multiplication(V2[i]*V2[i]),fourier0, fourier0,Interval{Float64}) - project(Multiplication(V2[i]),fourier0, fourier0,Interval{Float64})^2 for i = eachindex(V2)]

        ZV0 = [coefficients(D1.*B[i]*MV0[i]*Badj[i].*D2') for i = eachindex(MV0)] ; ZV0 = grid2cheb0(ZV0,N_) 
        ZV1 = [coefficients(D1.*B[i]*MV1[i]*Badj[i].*D2') for i = eachindex(MV1)] ; ZV1 = grid2cheb0(ZV1,N_) 
        ZV2 = [coefficients(D1.*B[i]*MV2[i]*Badj[i].*D2') for i = eachindex(MV2)] ; ZV2 = grid2cheb0(ZV2,N_) 

        Z12 =       norm(dn.*[sqrt(opnorm(LinearOperator(ZV0[i]),2)) for i = eachindex(ZV0)],1)/l0 
        Z12 = Z12 + norm(dn.*[sqrt(opnorm(LinearOperator(ZV1[i]),2)) for i = eachindex(ZV1)],1)/l1 
        Z12 = Z12 + norm(dn.*[sqrt(opnorm(LinearOperator(ZV2[i]),2)) for i = eachindex(ZV2)],1)/l2 
          
          display("value of Z12")
          display(Z12)
          
          Z13 =       norm(dn.*[norm(grid2cheb([coefficients(W0[i]) for i = eachindex(W0)],N_)[j],1) for j = eachindex(ZV0)],1)/l0 
          Z13 = Z13 + norm(dn.*[norm(grid2cheb([coefficients(W1[i]) for i = eachindex(W1)],N_)[j],1) for j = eachindex(ZV0)],1)/l1  
          Z13 = Z13 + norm(dn.*[norm(grid2cheb([coefficients(W2[i]) for i = eachindex(W2)],N_)[j],1) for j = eachindex(ZV0)],1)/l2  

          display("value of Z13")
          display(Z13)
          e0 = Sequence(fourierC_ ,interval.(zeros(N+1)))
          e0[0] = interval(1)
          Z14 = [coefficients(c[i]*e0 - WT[i]*(c[i]*e0-V2[i])) for i ∈ eachindex(c)] ; Z14 = grid2cheb(Z14,N_)
          Z14 = interval(4)*norm([norm(Z14[i],1)/c0 for i = eachindex(Z14)],1)

          display("value of Z14")
          display(Z14)
      
          return sqrt(Z11^2 + Z12^2 + Z13^2 + Z14^2)
      end
      
      
      
      function computation_norm_B(B,W,N,i)
        D1 = convert(Vector{Interval{Float64}},interval.(exp2cos(N)))
          D2 = interval.(ones((N+1)))./D1
        if i <= Nc
              MWT = project(Multiplication(W*W),fourierC_, fourierC_,Interval{Float64}) - project(Multiplication(W),fourierC_, fourierC_,Interval{Float64})^2
              n_W = sqrt(opnorm(LinearOperator(coefficients(D1.*(MWT).*D2')),2))   # computation of \|π_NW_Tπ^N\|_2
              n_BN = (opnorm(LinearOperator(D1.*((B*B')).*D2'),2))^(interval(1)/interval(2))   #computation of \|B^N_T\|_2
              return  maximum([norm(W,1) n_BN])+n_W 

        else
          return  norm(W,1) + (opnorm(LinearOperator(D1.*((B*B')).*D2'),2))^(interval(1)/interval(2))

        end
      end    
      
      
      
      function computation_Yu(U,a0b,db,N)
      
          fourierE = CosFourier(2*N, interval(π)/d)
          CY0 = interval(2.28)*interval(0.5)*maximum([interval(1)+sqrt(interval(2)) sqrt(interval(2))/(interval(1)-exp(-interval(π)))])
          D1 = convert(Vector{Interval{BigFloat}},interval.(big.(exp2cos(N))))
          νb = interval(big(4))/interval(big(π))^2
          D² = project(Derivative(2), CosFourier(N,interval(big(π))/interval(big(d))), CosFourier(N,interval(big(π))/interval(big(d))),Interval{BigFloat})
          Lνb = diag(interval(big(1)).-coefficients(νb*D²))
          Eb = Sequence(fourierE ,interval.(big.(zeros(2*N+1))))
          for n = 1:2N  
              Eb[n] = real(interval(big((-1)^n))*(interval(big(1))/(im*interval(big(n))/db*interval(big(pi)) + interval(big(2))a0b) - interval(big(1))/(im*interval(big(n))/db*interval(big(pi))-interval(big(2))*a0b)))
          end
          Eb[0] = interval(big(1))/a0b
          Eb = Eb*(interval(big(1))-exp(-interval(big(4))*a0b*db))/(interval(big(4))*db)
          
          
          ### Inner product part 
          PSY0 = abs( (coefficients(D1.*project(((Lνb.*U)*Eb),fourierC_))')*coefficients(D1.*(Lνb.*U)) )
          return sqrt(interval(4)*d^2*CY0*PSY0)
      end
      

 function computation_Zu(U,a,d,c0) 

        d=interval(d)
        fourier = fourierC_
        U = Sequence(fourierC_,vec(U))
      
        V0 = -interval(2)*ν*π^2/d^2*((0:N).^2).*U
          ##### Construction of the needed constants
        # Construction of ξ0 
        Ca = (interval(1) + abs(cos(interval(2)*a)))/(interval(1)-abs(cos(interval(2)*a)))
        ξ0 = interval(1)
        
        # Other constants of Lemma 4.1
        C0 = interval(1)/(π*σ0*(interval(1)-ν*a^2)) + interval(1)/(π*ν*σ0)

          C1 = interval(1)/(interval(2)*abs(c0)*ν) + Ca^(interval(1)/interval(4))*(interval(1)+sqrt(a))/(π*abs(c0)*σ0)*( interval(1)/(interval(1)-ν*a^2) + interval(2)/ν )
          K1 = (interval(2)+4*exp(-interval(2)*ξ0))/(minimum([interval(1) abs(c0)^3])*π*sqrt(ξ0)*σ0) + interval(1)/(π*σ0*abs(c0)) + interval(2)/(π*c0^2) + interval(1)/(π*abs(c0)^3)*(interval(2)+3*log(3*ξ0)) + interval(1)/(c0^2*sqrt(interval(2)*π))
          K2 = Ca^(interval(1)/4)/π*( interval(1)/(interval(2)*σ0^2)*( interval(2) + a^(-3/2)) + interval(1)/(4*σ0^2*(interval(1)-abs(cos(interval(2)a)))^2*sqrt(a)) )
          C2 = maximum([K2 K1*exp(a)]) 
        
        
        ###################### Computation of 𝒵u #####################################
        
        #### Computation of the constants C(d) and C1(d) in Lemmas 4.2 and 4.6
        
        C1d = exp(-interval(2)*a*d)*( interval(2)*sqrt(π)/(sqrt(interval(4)*a*d)*(interval(1)-exp(-interval(2)*a*d))) + interval(4)/(interval(1)-exp(-interval(2)*a*d)) )
        Cd = exp(-interval(2)*a*d)*( interval(4)*d + interval(4)*exp(-a*d)/(a*(interval(1)-exp(-interval(3)*a*d/interval(2)))) + interval(2)/(a*(interval(1)-exp(-a*d*interval(2)))))
        
        ### Inner product (V0,E*V0)

        fourierE = CosFourier(2*N, interval(π)/interval(d)) ; E = Sequence(fourierE,interval.(zeros(2*N+1)))
          for n = 1:2N  
              E[n] = real(interval(((-1)^n))*(interval((1))/(im*interval(big(n))/d*interval(big(pi)) + interval(big(2))a) - interval(big(1))/(im*interval((n))/d*interval((pi))-interval((2))*a)))
          end
          E[0] = interval((1))/a
          E = E*(interval((1))-exp(-interval((4))*a*d))/(interval((4))*d)

        D1N = convert(Vector{Interval{Float64}},interval.(exp2cos(N)))
        D2N = interval.(ones((N+1)))./D1N
        PS0 = abs( (coefficients(D1N.*project(V0*E,fourier))')*coefficients(D1N.*V0) )
        
        ### Inner product (V1,Efull*V1)
        k = (-N:N)*π/db
        fourier_f = Fourier(N,π/d)
        V1f = interval(4)*ν_*coefficients(U) ; V1f = [reverse(V1f[2:N+1]); V1f[1] ; V1f[2:N+1]].*k ; V1f = Sequence(fourier_f,V1f)
        Ef = coefficients(E); Ef = [reverse(Ef[2:2*N+1]); Ef[1] ; Ef[2:2*N+1]] ; Ef = Sequence(Fourier(2*N,π/d),Ef)
        PS1 = abs( (coefficients(project(V1f*Ef,fourier_f))')*coefficients(V1f) )
        
        ### Inner product (V2,E*V2)
        PS2 = interval(8)*abs( (coefficients(project(U*E,fourier))')*coefficients(U) )
        
        # Computation of ∫|v2'|^2 from d-1 to d
        char = Sequence(Fourier(2*N,π/d) ,interval.(big.(zeros(4*N+1))))
        for n = [-2N:-1;1:2N]  
            char[n] = real(interval(big((-1)^n))*(exp(im*n*π/db)-interval(big(1)))/(interval(big(2))*im*n*π))
        end
        char[0] = interval(big(1))/(interval(big(2))*db)
        Elog =  abs( (coefficients(project(V1f*char,fourier_f))')*coefficients(V1f) )
        
        # Computation of the Zu_i
        Zu0 = interval(4)*d*C0^2/a*PS0 + interval(2)*d*Cd*C0^2*PS0
        Zu1 = interval(4)*d*C1^2/a*PS1 + interval(2)*d*Cd*C1^2*PS1
        Zu2 = interval(2)*C2^2*( interval(2)*d*PS2/a + interval(4)*log(interval(2))*Elog ) + interval(2)*d*C1d*C2^2*PS2
        
        𝒵u = sqrt(Zu0 + Zu1 + Zu2) 
        

        return 𝒵u
end












function proof_continuation_whitham(V,c0,c1,N,Nc,d,a,σ0,σ1)

  a0b = interval(big(1.5)) ; db = interval(big(d))
  fourierC = CosFourier(N, π/d)   ; fourierC_ = CosFourier(N, interval(π)/interval(d))
N_ = 3Nc
N_fft_ = nextpow(2, 2N_ + 1)
npts_ = N_fft_ ÷ 2 + 1



c = [0.5 * (c0 + c1) + 0.5cospi(2k / N_fft_)*(c1 - c0) for k ∈ 0:npts_-1]
c_ = [interval(0.5) * (interval(c0) + interval(c1)) + interval(0.5)*cospi(interval(2)*interval(k) / interval(N_fft_)) * (interval(c1) - interval(c0)) for k ∈ 0:npts_-1]


V_fft= [V ; reverse(V)[2:npts_-1]] ; c_fft = [c ; reverse(c)[2:npts_-1]] ; c_fft_ = [c_ ; reverse(c_)[2:npts_-1]]
V_cheb = Vector{Sequence{Chebyshev,Vector{Interval{Float64}}}}(undef, N+1)

for n=0:N
    V_cheb[n+1] = rifft!(complex.([V_fft[j][n] for j = 1:2*(npts_-1)]), Chebyshev(Nc))  #computation of the Chebyshev coefficients
end
V_cheb = [real.(getindex.(V_cheb, i)) for i ∈ eachindex(coefficients(V_cheb[1])).-1] # reorganization of the coefficients
 V_cheb0_b = [project_trace(Sequence(fourierC_,interval.(big.(mid.(V_cheb[i]))))) for i = eachindex(V_cheb)]  # Rigorous projection of the Chebyshev coefficients into the trace zero functions in big float
 V_cheb0 = [project_trace(Sequence(fourierC_,interval.((mid.(V_cheb[i]))))) for i = eachindex(V_cheb)]  # Rigorous projection of the Chebyshev coefficients into the trace zero functions in float64
V_cheb = [Sequence(Chebyshev(Nc),[V_cheb0[j][i] for j = 1:length(V_cheb0)]) for i = 0:N]
V_cheb_b = [Sequence(Chebyshev(Nc),[V_cheb0_b[j][i] for j = 1:length(V_cheb0_b)]) for i = 0:N]
V_fft_ = cheb2grid(V_cheb, N_fft_)

V_cheb0 = [real.(getindex.(V_cheb, i)) for i ∈ eachindex(coefficients(V_cheb[1])).-1] 
V_cheb0_b = [real.(getindex.(V_cheb_b, i)) for i ∈ eachindex(coefficients(V_cheb_b[1])).-1] 

# # V_cheb_seq = [coefficients(V_cheb[i])[j] for i = eachindex(V_cheb), j = 1:N+1] 

# V_cheb_seq = [(V_cheb0[i])[j] for i = eachindex(V_cheb0), j = 1:N+1] 
#  V_cheb_seq = Sequence(Chebyshev(Nc)⊗fourierC,mid.(vec(V_cheb_seq)))
# # V_fft_ = V_cheb_seq.(c_fft_,nothing) ; V_fft_ = [coefficients(V_fft_[i]) for i= eachindex(V_fft_)]
############# construction of the coefficients Bn ###################
 
B_grid = complex.(inv.([mid.(Lν).*mid.(DF_W(mid.(V_fft[i]),c_fft[i],T,fourierC_))*inv((mid.(M0)-c_fft[i]*I).*mid.(Lν)') for i=1:2(npts_-1)]))
B_cheb = [interval.(mid.(rifft!(complex.([B_grid[k][i,j] for k = 1:length(B_grid)]), Chebyshev(Nc)))) for i = indices(codomain(B_grid[1])), j = indices(domain(B_grid[1]))]
B_cheb0 = [real.(getindex.(B_cheb, i)) for i ∈ eachindex(coefficients(B_cheb[1])).-1] 


BF_fft_ = real.(cheb2grid(B_cheb, N_fft_)) .* F_W(V_fft_, c_fft_,T)
BF_cheb = grid2cheb(complex.(BF_fft_), N_)
BF_cheb = [real.(getindex.(BF_cheb, i)) for i ∈ eachindex(coefficients(BF_cheb[1])).-1]

e0 = Sequence(fourierC_ ,interval.(zeros(N+1)))
e0[0] = interval(1)
e02 = Sequence(CosFourier(2*N, interval(π)/interval(d)) ,interval.(zeros(2*N+1)))
e02[0] = interval(1)
W_fft = [interval.(mid.(project(Multiplication(e0 - 2/mid(c_fft[i])*Sequence(fourierC_,V_fft_[i])),fourierC_, CosFourier(2*N, interval(π)/interval(d)),Interval{Float64}))\mid.(e02)) for i=1:length(c_fft)]
W_fft0 = [coefficients(W_fft[i]) for i=1:length(c_fft)]
W_cheb = grid2cheb(complex.(W_fft0), Nc) 
W_cheb = [real.(getindex.(W_cheb, i)) for i ∈ eachindex(coefficients(W_cheb[1])).-1]
W_cheb0 = [Sequence(fourierC_,W_cheb[i]) for  i∈eachindex(W_cheb)]


## Computation of the norm of the components of B(c)
 norm_B = [computation_norm_B(B_cheb0[i],W_cheb0[i],N,i) for i=1:length(B_cheb0)]
 norm_B[1] = maximum([interval(1) norm_B[1]])

 dn = interval(2)*interval.(ones(Nc+1)) ; dn[1] = interval(1)
 norm_B_sup = norm(dn.*norm_B,1)
 display("norm of B")
 display(norm_B_sup)


################################# Computation of the Y0 bound #############################



D²_2 = project(Derivative(2), CosFourier(2N, interval(π)/interval(d)), CosFourier(2N, interval(π)/interval(d)),Interval{Float64})
Lν2 = Id - ν_*D²_2
Lν2 = diag(coefficients(Lν2))
WF_fft = [coefficients(W_fft[i]*(Lν2.*(Sequence(fourierC_,V_fft_[i])^2 - project(Sequence(fourierC_,V_fft_[i])^2,fourierC_)))) for i=eachindex(W_fft)]

WF_cheb = [coefficients(rifft!(complex.([(WF_fft[k])[i] for k = 1:length(WF_fft)]), Chebyshev(N_))) for i=eachindex((WF_fft[1]))]
WF_cheb = [real.(getindex.(WF_cheb, i)) for i ∈ eachindex(WF_cheb[1])]

Y0 = (norm_cheb(BF_cheb,CosFourier(N,interval(π)/interval(d)))^2 + norm_cheb(WF_cheb,CosFourier(3N,interval(π)/interval(d)))^2)^0.5

display("value of Y0")
display(Y0)

a0 = interval(1.5)
Cd = exp(-interval(2)*a0*d)*( interval(4)*d + interval(4)*exp(-a0*d)/(a0*(interval(1)-exp(-interval(3)*a0*d/interval(2)))) + interval(2)/(a0*(interval(1)-exp(-a0*d*interval(2)))))
Yu1 = [computation_Yu(Sequence(fourierC_,V_cheb0_b[i]),a0b,db,N) for i=eachindex(V_cheb0_b)]
dn = interval(2)*interval.(ones(Nc+1)) ; dn[1] = interval(1)
Yu = norm(dn.*Yu1*(interval(1)+norm_B_sup^2*(interval(1)+Cd))^(0.5),1)

display("value of Yu")
display(Yu)

𝒴0 = Yu + Y0


##################### Computation of the Z1 bound ##############################


Z1 = computation_Z1(B_grid,W_fft,[Sequence(fourierC_,V_fft_[i]) for i∈eachindex(V_fft_)],interval(c0),interval(c1),c_fft_,N)

display("value of Z1")
display(Z1)

Zu = norm(dn.*[computation_Zu(V_cheb0[i],a,d,interval(c0)) for i ∈ eachindex(V_cheb0)].*norm_B,1)

display("value of Zu")
display(Zu)

𝒵1 = Z1 + Zu

display("value of 𝒵1")
display(𝒵1)

################### Computation of the Z2 bound #####################################

# construction of the function 1-c as a Chebyshev polynomial 
f_c = Sequence(Chebyshev(Nc),interval.(zeros(Nc+1)))
f_c[0] = interval(1) + interval(2)*(c1-c0)/interval(4)*(interval(1)-interval(2)*c1/(c1-c0))
f_c[1] = -(c1-c0)/interval(4)

# function 1
e0 = Sequence(Chebyshev(2*Nc),interval.(zeros(2*Nc+1))) ; e0[0] = interval(1)
# approximate inverse of f_c
c_inv =  interval.(mid.(project(Multiplication(mid.(f_c)),Chebyshev(Nc), Chebyshev(2*Nc),Interval{Float64}))\mid.(e0))

n_c = interval(1)/(interval(1)-norm(e0-f_c*c_inv,1))   # rigorous diffect in the Neumann series
# computation of B(c)/(1-c)
Bf_c = real.(cheb2grid(B_cheb, N_fft_)) ; c_inv_grid =  real.(fft(c_inv,N_fft_)) ;  Bf_c = [c_inv_grid[i]*Bf_c[i] for i = eachindex(c_inv_grid)]
Bf_c = [rifft!(complex.([Bf_c[k][i,j] for k = 1:length(Bf_c)]), Chebyshev(2*Nc)) for i = axes(Bf_c[1],1), j = axes(Bf_c[1],2)]
Bf_c = [real.(getindex.(Bf_c, i)) for i ∈ eachindex(coefficients(Bf_c[1])).-1]

W_c = grid2cheb(complex.([c_inv_grid[i]*W_fft0[i] for i = eachindex(c_inv_grid)]), 2*Nc)  
W_c = [real.(getindex.(W_c, i)) for i ∈ eachindex(coefficients(W_c[1])).-1]
W_c = [Sequence(fourierC_,W_c[i]) for  i∈eachindex(W_c)]

n_Bf_c = [computation_norm_B(Bf_c[i],W_c[i],N,i) for i=1:length(Bf_c)]

dn2 = interval(2)*interval.(ones(2*Nc+1)) ; dn2[1] = interval(1)
𝒵2 = interval(2)*two_norm_inverse_l(d,ν,c0,T)*maximum([1/(c0-1) norm(dn2.*n_Bf_c,1)])*n_c

display("value of 𝒵2")
display(𝒵2)

display(interval(2)*two_norm_inverse_l(d,ν,c0,T)*1/(c0-1)*norm_B_sup)


α = interval( sup( interval(1)/interval(2) + 𝒴0*𝒵2/(interval(1)-𝒵1)^2 ) )
r = interval(2)*α*𝒴0/(interval(1)-𝒵1)

if inf(1- 𝒵1)>0
  if sup(𝒵2*r+𝒵1) < 1
      if inf(1/2*𝒵2*r^2 - (1-𝒵1)*r + 𝒴0) < 0
          display("The computer-assisted proof was successful")
          display("Radius of contraction :")
          display(r)
      else
      display("failure: discriminant is negative")
      end
  else
      display("r is too big")
  end
else
  display("failure: 1-z > 0")
end


######### Regularity for the solution (using Proposition 4.1)

    dn = interval(2)*interval.(ones(Nc+1)) ; dn[1] = interval(1)
    ϵ = interval(0.5)*interval(c0) - (interval(2)/interval(π)*log(interval(Nc+1))+interval(1))*maximum([norm(V_fft_[i],1) for i = eachindex(V_fft_)])   # Lebesgue constant to compute an upper bound for the infinity norm (cf. Trefethen Approximation Theory and Approximation Practice)
    ϵ0 = interval(0.5)*interval(c0) - norm(dn.*[norm(V_cheb0[i],1) for i= eachindex(V_cheb0)],1) ; ϵ = maximum([ϵ,ϵ0])
    display(inf(ϵ))
    if inf(ϵ)<=0
      display("epsilon is negative, the proof fails")
    end


  if sup(r) <= inf(interval(4)*ϵ*sqrt(ν)*σ0)
    display("the solution is infinitely differentiable")
  else
    display("the solution is H^(3/2) (continuous at least), but could not prove more regularity")
  end

    return V_fft_[npts_], V_fft_[1], r , (1-𝒵1)/𝒵2


end