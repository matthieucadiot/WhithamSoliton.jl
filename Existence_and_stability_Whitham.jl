using  JLD2, PrettyTables, RadiiPolynomial, IntervalArithmetic, LinearAlgebra 


function trace(N)
  w = 0:N;
  Y = w.^2;
  X = 2*(-1).^w; X[1] = 1;

  f = [X';
       (X.*Y)']
  return f
end



function exp2cos(N)

  d = 2*(ones((N+1)))

  d[1] =1;

  return d
end


function test_constants(T,ν, c,a,σ0,σ1,x,di,dj,xmin,ymin)

  setprecision(Precis)

  if T>0
      # first requirement on a
      if a >= minimum([π/2 1/sqrt(ν)])
          return 0,0
      end
      # we start by verifying σ0 and σ1, for these we can fix the imaginary part y= a
      X =  -sup(x):dj:sup(x)-dj 
      Z = interval.(X,X.+dj*ones(size(X))).+im*a

      m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z).*( 1 .+ T*Z.^2)).-c).-σ0
      m2 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z).*( 1 .+ T*Z.^2)).-c).-σ1*sqrt.(T*abs.(X))
       
      if (minimum(inf.(m1))<=0)||(minimum(inf.(m2))<=0)
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
              m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z).*( 1 .+ T*Z.^2)).-c)
              while (minimum(inf.(m1))<=0)&&(k<4) 
                  ddi = ddi/10
                  ddj=ddj/10
                  Z = interval(x1,x1+ddi)+im*interval(y1,y1+ddj)
                  m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z).*( 1 .+ T*Z.^2)).-c)
                  k=k+1
              end
              if k==4
                  display("value for a is too big")
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

      m1 = abs.(sqrt.((exp.(2Z).-1.0)./((exp.(2Z).+1.0).*Z).*( 1 .+ T*Z.^2)).-c)
      if (minimum(inf.(m1))<=0)
          display("second test")
          return 0,inf.(m1)
      end
           

  else
      if a >= π/2
          return 0,0
      end

       # we start by verifying σ0 and σ1, for these we can fix the imaginary part y= a
      X =  -sup(x):dj:sup(x)-dj ; 
      Z = interval.(X,X.+dj*ones(size(X))).+im*a

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
                  display("value for a is too big")
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




function proof_eigen(U,Vb,U0,c,T,ν,d,r0)

  setprecision(Precis)

  fourier = space(U0)
  N = order(U0)

  V = component(U,2)  #(V correspond to Ψ in Section 5)
  λ = real(component(U,1)[1])
  
  #construction of the linear operator MT
  D² = real.(project(Derivative(2), fourier, fourier,Complex{Interval{Float64}}))
  dd = diag(coefficients(D²))  ;   dd = sqrt.(-dd) ; d0=dd;  d0[N+1] = interval(1);
  dL = tanh.(dd).*(ones(2*N+1)+T*dd.^2)./d0 ;  dL[N+1] = interval(1) ; dL = sqrt.(dL)
  MT = LinearOperator(fourier,fourier,Diagonal(dL))
  #construction of L
  if T>0
      L = MT - c*I - λ*I
      Linv = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(L))))
      DG = 2*real.(project(Multiplication(U0),fourier, fourier,Complex{Interval{Float64}}))
  else
      L = c*I - MT - λ*I
      Linv = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(L))))
      DG = -2*real.(project(Multiplication(U0),fourier, fourier,Complex{Interval{Float64}}))
  end

  pspace = ParameterSpace()×fourier
 
  M = LinearOperator(pspace,pspace, [[0   -2*d*coefficients(V)'] ; [-coefficients(V) coefficients(L+DG)]] )
  D = [[interval(1) zeros(1,2*N+1)] ; [zeros(2*N+1,1) coefficients(L)]]
  Di = [[interval(1) zeros(1,2*N+1)] ; [zeros(2*N+1,1) coefficients(Linv)]]
  P1 = [[interval(1) zeros(1,2*N+1)] ; [zeros(2*N+1,1) sqrt(2d)*Diagonal(ones(2*N+1))]]
  P2 = [[interval(1) zeros(1,2*N+1)] ; [zeros(2*N+1,1) 1/sqrt(2d)*Diagonal(ones(2*N+1))]]
  D = LinearOperator(pspace,pspace,D)
  Di = LinearOperator(pspace,pspace,Di)
  P1 = LinearOperator(pspace,pspace,P1)
  P2 = LinearOperator(pspace,pspace,P2)


# # # # ################ Z BOUND ######################################################

  B = inv(mid.(M)*mid.(Di)) ; B = interval.(real.(B)) + im*interval.(imag.(B))
  Badj = LinearOperator(pspace,pspace,coefficients(B)')

  e0 = Sequence(fourier,interval.(zeros(2*N+1)))
  e0[0] = 1
  if T>0
      WT = e0
  else
      e02 = Sequence(Fourier(2*N, π/d),interval.(zeros(4*N+1)))
      e02[0] = 1
      WT = interval.(mid.(project(Multiplication(mid.(e0 - 2/mid(c-λ)*U0)),fourier,Fourier(2*N, π/d),Interval{Float64}))\mid.(e02))
  end

# computation of the norm of B
  MWT = project(Multiplication(WT*WT),fourier, fourier,Interval{Float64}) - project(Multiplication(WT),fourier, fourier,Interval{Float64})*project(Multiplication(WT),fourier, fourier,Interval{Float64})
  n_WT = sqrt(opnorm(LinearOperator(coefficients(MWT)),2))   # computation of \|π_NW_Tπ^N\|_2
  n_BN = sqrt(opnorm(LinearOperator(coefficients((P1*B*P2^2*Badj*P1))),2))   #computation of \|B^N_T\|_2
  norm_B = maximum([1 maximum([norm(WT,1) n_BN])+n_WT])  #computation of the norm of B

 
  if T>0
    L2 = Linv ;          l2 = maximum(vec((abs.(dL.-(c+λ)))))
  else
    L2 = Linv + 1/(c-λ)*I ;  l2 = maximum(vec(abs.(((c-λ)*(dL.-(c-λ)))./dL)))
  end

    W2 = 2*WT*U0
    #MWi corresponds to the operator πᴺ(Vi*WT)(π^{3N}-πᴺ)(Vi*WT)πᴺ = πᴺ(Vi*WT*Vi*WT)πᴺ - (πᴺ(Vi*WT)πᴺ)^2 (recall that Wi = Vi*WT)
    MW2 = project(Multiplication(W2*W2),fourier, fourier,Interval{Float64}) - project(Multiplication(W2),fourier, fourier,Interval{Float64})^2
#   MVi corresponds to the operator πᴺVi(π^{2N}-πᴺ)Viπᴺ = πᴺ(Vi*Vi)πᴺ - (πᴺViπᴺ)^2
    MV2 = 4*project(Multiplication(U0*U0),fourier, fourier,Interval{Float64}) - 4*project(Multiplication(U0),fourier, fourier,Interval{Float64})^2
    MV2 =  LinearOperator(pspace,pspace, [[0   zeros(1,2*N+1)] ; [zeros(2*N+1,1) coefficients(MV2)]] )

  Z11 =  opnorm(LinearOperator(coefficients(P1*(I - B*M*Di)*P2)),2)^2
  Z11 = Z11 + sqrt(opnorm(LinearOperator(coefficients(L2*MW2*L2)),2))
  Z11 = sqrt(Z11)

  Z12 = sqrt(opnorm(LinearOperator(coefficients(P1*B*MV2*Badj*P2)),2))/l2 

  Z13 =  norm(W2,1)/l2

  if T>0
    Z14 = 0
  else
    Z14 = norm(e0 - WT*(e0-2/(c-λ)*U0),1)
  end

  Z1 = sqrt(Z11^2 + Z12^2 + Z13^2 + Z14^2)


#computation of the supremum of Zu for all λ ∈ (0,λmin)
# We first compute the supremum of Cλ
λmin = abs(1-c)-2*norm(U0,1)-r0/(4*sqrt(ν)*abs(1-c)) # λmin is a lower bound for the eigenvalues
if T>0
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
  Eb[n] = ab/db*(-1)^n*(1-exp(-4*ab*db))/(4*ab^2 + (n*interval(big(π)/db))^2 )
end


# For all λ ≤0, a-λ ≥ a and so the maximal 𝒵λ is achieved for λ=0
C1d = exp(-2*a*d)*( 2*sqrt(π)/(sqrt(4*a*d)*(1-exp(-2*a*d))) + 4/(1-exp(-2*a*d)) )
DU0 = U0fb.*( (-N:N)*π/db )
Elog0 =  abs( (coefficients(project(DU0*char,fourier))')*coefficients(DU0) )
DV = Vb.*( (-N:N)*π/db )
Elogphi =  abs( (coefficients(project(DV*char,fourier))')*coefficients(DV) )
𝒵u0 = abs( (coefficients(project(U0fb*Eb,fourier))')*coefficients(U0fb) )
𝒵uphi = abs( (coefficients(project(Vb*Eb,fourier))')*coefficients(Vb) )

#Computation of the supremum of 𝒵λ for all λ ≤ 0
𝒵u = maximum([2*sqrt(2*Cλ^2*(d*𝒵u0*(2/a+C1d) + 8*log(interval(2))*Elog0 ))  sqrt(2*d)*sqrt(2*Cλ^2*(d*𝒵uphi*(2/a+C1d) + 8*log(interval(2)) *Elogphi)) ])
𝒵u = Interval.(Float64(inf(𝒵u),RoundDown),Float64(sup(𝒵u),RoundUp) ) # go back to standard precision

  ## Computation of 𝒵1
  𝒵1 = sup(𝒵u*norm_B + Z1 + norm_B*r0/(4*sqrt(ν)*σ0*(σ0-λ)))

################# Y BOUND ######################################################

Cd = exp(-2*a0*d)*( 4*d + 4*exp(-a0*d)/(a0*(1-exp(-3*a0*d/2))) + 2/(a0*(1-exp(-a0*d*2))))

D² = real.(project(Derivative(2), Fourier(N, π/db), Fourier(N, π/db),Complex{Interval{BigFloat}}))
Lνb = I - νb*D²
Lνb = diag(coefficients(Lνb)).*2

Eb = Sequence(Fourier(2*N, π/d) ,interval.((zeros(4*N+1))))
for n = -2N:2N  
    Eb[n] = a0b/db*(-1)^n*(1-exp(-4*a0b*db))/(4*a0b^2 + (n*interval(big(π)/db))^2 )
end

### Inner product part 
PSY0 = abs( (coefficients(project(((Lνb.*Vb)*Eb),fourier))')*coefficients((Lνb.*Vb)) )

Yu = sqrt((4*d^2*CY0*PSY0)*(1+norm_B^2*(Cd+1)))


B22 = component(B,2,2)
B12 = component(B,1,2)
Y0 = sqrt(2*interval(d))*( sqrt(2*norm(B22*(L*V+DG*V),2)^2 + 2*norm(2*U0*V - project(2*U0*V,fourier))^2)) + 2*sqrt(2*interval(d))*norm(L*V+ DG*V,2)*norm(L*V,2)

𝒴0 = Yu + Y0 + 2*r0/σ0*maximum([norm(B12*V,2)/(sqrt(2*d))  opnorm(LinearOperator(coefficients(B22)),2)*norm(V,1) ])

################## Z2 BOUND ######################################################
  Z2 = norm_B
  
################## Verification ###################################################
# # # r_star = 1.0;
# # # interval_proof = interval_of_existence(Y,Interval(0),Z2,r_star,C²Condition())
# # # display(interval_proof)
# # #
x = 1;
if T>0
  λmax = σ0
else
  λmax = Interval.(Float64(inf(minimum([σ0 c-2*(norm(U0,1)+r0/(4*sqrt(ν)*σ0))])),RoundDown))
end

if 1- 𝒵1>0
  if 1-sup(4*𝒴0*Z2) > 0
    rmin=Float64(sup((1-𝒵1 - sqrt((1-𝒵1 )^2-4*𝒴0*Z2))/(2*Z2)),RoundUp);
    if rmin<rmax
      display("The computer-assisted proof was successful for the eigenvalue λ with value")
      display(Float64.(mid(λ)))
      display("Value of the minimal radius for the contraction")
      display(Float64(rmin))
      if (1-𝒵1-rmin*norm_B/((σ0 - λ)^2)>0)&&(Float64(sup(λ+rmin),RoundUp)<λmax)
        R = Interval.(Float64(inf((σ0 - λ)*(1-𝒵1-rmin*norm_B/((σ0 - λ)^2))/(norm_B)),RoundDown))
        display("λ is the only eigenvalue in")
        display([sup(λ-R) inf(minimum([λ+R λmax]))])
      else
        R = None
      end
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
min_rad_constraction = minimum([R λmax-λ])
return [x==1,min_rad_constraction,rmin]

end




function no_eigen_outside(r1,r2,U1,U2,U0,rmin,c,T,ν)

  setprecision(Precis)

  fourier = space(U0)
  N = order(U0)

  λ1 = real(component(U1,1)[1])
  λ2 = real(component(U2,1)[1])

  # Construction of the quantities that do not change with λ

  D² = real.(project(Derivative(2), fourier, fourier,Complex{Interval{Float64}}))
  dd = diag(coefficients(D²))  ;   dd = sqrt.(-dd) ; d0=dd;  d0[N+1] = interval(1);
  dL = tanh.(dd).*(ones(2*N+1)+T*dd.^2)./d0 ;  dL[N+1] = interval(1) ; dL = sqrt.(dL)
  MT = LinearOperator(fourier,fourier,Diagonal(dL))

  if T>0
      L = MT - c*I
      Linv = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(L))))
      DG = 2*real.(project(Multiplication(U0),fourier, fourier,Complex{Interval{Float64}}))
  else
      L = c*I - MT
      Linv = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(L))))
      DG = -2*real.(project(Multiplication(U0),fourier, fourier,Complex{Interval{Float64}}))
  end

# Now that the above quantities are fixed, we prove that there is no eigenvalue between λ1 and λ2 and that there is no eigenvalue between λ2 and λmin. This implies that λ1 and λ2 are the smallest eigenvalues of DF(u)

λmin = abs(1-c)-2*norm(U0,1)-rmin/(4*sqrt(ν)*abs(1-c))
λ = interval(Float64(sup(λ2-r2),RoundUp))
P = 1
k=1

#computation of the supremum of Zλ for all λ ∈ (0,λmin)
# We first compute the supremum of Cλ
if T>0
  Cλ = computation_Cλ(maximum([c c+λmin]),a,σ0,T)
else
  Cλ = computation_Cλ(c-λmin,a,σ0,T)
end

# coefficients of the characteristic function
char = Sequence(Fourier(2*N,π/d) ,interval.(big.(zeros(4*N+1))))
for n = [-2N:-1;1:2N]  
  char[n] = real((-1)^n*(exp(im*n*interval(big(π))/db)-1)/(2*im*n*interval(big(π))))
end
char[0] = 1/(2*db)

# Coefficients of cosh(2ax)
Eb = Sequence(Fourier(2*N,π/d),interval.(big.(zeros(4*N+1))))
for n = -2N:2N  
  Eb[n] = real((-1)^n*(1/(im*n/db*interval(big(pi)) + 2ab) - 1/(im*n/db*interval(big(pi))-2*ab)))
end
Eb[0] = 1/ab
Eb = Eb*(1-exp(-4*ab*db))/(4*db)

# For all λ ≤0, a-λ ≥ a and so the maximal 𝒵λ is achieved for λ=0
C1d = exp(-2*a*d)*( 2*sqrt(interval(big(π)))/(sqrt(4*a*d)*(1-exp(-2*a*d))) + 4/(1-exp(-2*a*d)) )
DU0 = U0.*( (-N:N)*π/db )
Elog =  abs( (coefficients(project(DU0*char,fourier))')*coefficients(DU0) )
𝒵u0 = abs( (coefficients(project(U0fb*Eb,fourier))')*coefficients(U0fb) )

#Computation of the supremum of 𝒵λ for all λ ≤ 0
𝒵λ = 8*Cλ^2*(d*𝒵u0*(2/a+C1d) + 4*log(interval(2))*Elog )
𝒵λ = sqrt(𝒵λ)

𝒵λ = Interval.(Float64(inf(𝒵λ),RoundDown),Float64(sup(𝒵λ),RoundUp) )
 # Construction of the unit vector e0
 e0 = Sequence(fourier ,interval.(zeros(2*N+1)))
 e0[0] = 1


  while  (P==1)&&(sup(λ)>=inf(λ1+r1))
      Lλ = L - λ*I
      Li = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(Lλ))))

      Df = Lλ + DG
      B = interval.(inv(mid.(Df*Li)))
      Badj = LinearOperator(fourier,fourier,coefficients(B)')
      n_BN = sqrt(opnorm( LinearOperator(coefficients(B*Badj)),2 ))   #computation of \|B^N_T\|_2

      if T>0
          WT = e0
          n_WT = 0
      else
          e02 = Sequence(Fourier(2*N, π/d) ,interval.(zeros(4*N+1)))
          e02[0] = 1
          WT = interval.(mid.(project(Multiplication(mid.(e0 - 2/mid(c-λ)*U0)),fourier,Fourier(2*N, π/d),Interval{Float64}))\mid.(e02))
          MWT = project(Multiplication(WT*WT),fourier, fourier,Interval{Float64}) - project(Multiplication(WT),fourier, fourier,Interval{Float64})*project(Multiplication(WT),fourier, fourier,Interval{Float64})
  
          n_WT = sqrt(opnorm(LinearOperator(coefficients(MWT)),2))   # computation of \|π_NW_Tπ^N\|_2
          MWT = Nothing
      end
  
      # computation of the norm of B
      norm_B = maximum([1 maximum([norm(WT,1) n_BN])+n_WT]) # Computation for an upper bound of the norm of B

      if T>0
          lT = maximum(vec((abs.(dL.-(c+λ)))))
      else
          lT = maximum(vec(abs.(((c-λ)*(dL.-(c-λ))./dL))))
      end 
      V2 = 2*U0
      W2 = WT*V2 ;          

#MW2 corresponds to the operator πᴺ(Vi*WT)(π^{3N}-πᴺ)(Vi*WT)πᴺ = πᴺ(Vi*WT*Vi*WT)πᴺ - (πᴺ(Vi*WT)πᴺ)^2 (recall that Wi = Vi*WT)
      MW2 = project(Multiplication(W2*W2),fourier, fourier,Interval{Float64}) - project(Multiplication(W2),fourier, fourier,Interval{Float64})^2

#MV2 corresponds to the operator πᴺVi(π^{2N}-πᴺ)Viπᴺ = πᴺ(Vi*Vi)πᴺ - (πᴺViπᴺ)^2
      MV2 = project(Multiplication(V2*V2),fourier, fourier,Interval{Float64}) - project(Multiplication(V2),fourier, fourier,Interval{Float64})^2

      if T>0
          Z14 = 0
      else
          Z14 = norm(e0 - WT*(e0-1/(c-λ)*V2),1)
      end

      # computation of each component of Z1
      Z11 =  opnorm(LinearOperator(coefficients(I - B*(Df*Li))),2)^2
      Z11 = Z11 + opnorm(LinearOperator(coefficients((Li*MW2*Li))),2)
      Z11 = sqrt(Z11)

      Z12 = sqrt(opnorm(LinearOperator(coefficients(B*MV2*Badj)),2))/lT

      Z13 = norm(W2,1)/lT

      Z1 = sqrt(Z11^2 + Z12^2 + Z13^2 + Z14^2)
      𝒞 = (1-norm_B*𝒵λ-Z1)/norm_B*(abs(1-c) -λ) - rmin/(4*sqrt(ν)*abs(1-c)) 
      if inf(𝒞) <= 0
          P=0
          return P==1
      end

      λ = interval(Float64(sup(λ - 𝒞),RoundUp))
  end

  λ = λ1-r1
  k=2
  
  while  (P==1)&&(sup(λ)>=inf(λmin))
      
    Lλ = L - λ*I
    Li = LinearOperator(fourier,fourier,Diagonal(ones(2*N+1)./diag(coefficients(Lλ))))

    Df = Lλ + DG
    B = interval.(inv(mid.(Df*Li)))
    Badj = LinearOperator(fourier,fourier,coefficients(B)')
    n_BN = sqrt(opnorm( LinearOperator(coefficients(B*Badj)),2 ))   #computation of \|B^N_T\|_2

    if T>0
        WT = e0
        n_WT = 0
    else
        e02 = Sequence(Fourier(2*N, π/d) ,interval.(zeros(4*N+1)))
        e02[0] = 1
        WT = interval.(mid.(project(Multiplication(mid.(e0 - 2/mid(c-λ)*U0)),fourier,Fourier(2*N, π/d),Interval{Float64}))\mid.(e02))
        MWT = project(Multiplication(WT*WT),fourier, fourier,Interval{Float64}) - project(Multiplication(WT),fourier, fourier,Interval{Float64})*project(Multiplication(WT),fourier, fourier,Interval{Float64})

        n_WT = sqrt(opnorm(LinearOperator(coefficients(MWT)),2))   # computation of \|π_NW_Tπ^N\|_2
        MWT = Nothing
    end

    # computation of the norm of B
    norm_B = maximum([1 maximum([norm(WT,1) n_BN])+n_WT]) # Computation for an upper bound of the norm of B

    if T>0
        lT = maximum(vec((abs.(dL.-(c+λ)))))
    else
        lT = maximum(vec(abs.(((c-λ)*(dL.-(c-λ))./dL))))
    end 
    V2 = 2*U0
    W2 = WT*V2 ;          

#MW2 corresponds to the operator πᴺ(Vi*WT)(π^{3N}-πᴺ)(Vi*WT)πᴺ = πᴺ(Vi*WT*Vi*WT)πᴺ - (πᴺ(Vi*WT)πᴺ)^2 (recall that Wi = Vi*WT)
    MW2 = project(Multiplication(W2*W2),fourier, fourier,Interval{Float64}) - project(Multiplication(W2),fourier, fourier,Interval{Float64})^2

#MV2 corresponds to the operator πᴺVi(π^{2N}-πᴺ)Viπᴺ = πᴺ(Vi*Vi)πᴺ - (πᴺViπᴺ)^2
    MV2 = project(Multiplication(V2*V2),fourier, fourier,Interval{Float64}) - project(Multiplication(V2),fourier, fourier,Interval{Float64})^2

    if T>0
        Z14 = 0
    else
        Z14 = norm(e0 - WT*(e0-1/(c-λ)*V2),1)
    end

    # computation of each component of Z1
    Z11 =  opnorm(LinearOperator(coefficients(I - B*(Df*Li))),2)^2
    Z11 = Z11 + opnorm(LinearOperator(coefficients((Li*MW2*Li))),2)
    Z11 = sqrt(Z11)

    Z12 = sqrt(opnorm(LinearOperator(coefficients(B*MV2*Badj)),2))/lT

    Z13 = norm(W2,1)/lT

    Z1 = sqrt(Z11^2 + Z12^2 + Z13^2 + Z14^2)

    𝒞 = (1-norm_B*𝒵λ-Z1)/norm_B*(abs(1-c) -λ) - rmin/(4*sqrt(ν)*abs(1-c)) 
    if inf(𝒞) <= 0
        P=0
        return P==1
    end

    λ = interval(Float64(sup(λ - 𝒞),RoundUp))
  end

  return P==1


end



function computation_Cλ(c,a,σ0,T)

  setprecision(Precis)

##### Construction of the needed constants
# Construction of ξ0 
  Ca = (1 + abs(cos(2*a)))/(1-abs(cos(2*a)))

if T>0
ξ0 = maximum([interval(1) 1/sqrt(T) c^2*Ca*4/T])
ξ0 = maximum([3*T/(2*tanh(ξ0)) 2*c^2/(T*tanh(ξ0))])
else
ξ0 = interval(1)
end

if T>0
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
  p = interval(big(π))
  Y2 = (p/d)*interval.(big.(w));
  Y1 = (p^2/d^2)*interval.(big.(w)).^2;
  Y3 = (p^3/d^3)*interval.(big.(w)).^3;
  X = interval.(big.((-1).^w));

  f = [X';
      (X.*Y1)';
       (X.*Y2)';
       (X.*Y3)']
  return f
end









######################## MAIN CODE ############################################################
#############################################################################################################################
Base.abs2(x::Interval) = (y=abs(x);return y*y)

Precis = 100 
setprecision(Precis)  #set up the precision for the bigfloat computations


################### Loading the candidate ##############################


### Candidate for T = 0, c = 1.1, N = 800 and d = 50
#   W = load("W_0_11_50.jld2","W")
#   N = 800 ; 
#   dn= 50 ; db = interval(big(dn)); d = interval(dn)
# c = interval(1.1)
# Tn = 0; T= abs(interval(Tn)) ; Tb = abs(interval(big(Tn))) 
# #### candidate values for a, σ0 and σ1. These values were obtained studying the graph of the function |mT-c| 
# #### values for T=0 and c=1.1
# a = 1;  ab = interval(big(a)) ; a = interval(a)
# σ0 = interval(0.03); σ1 = interval(0.1);

### Candidate for T = 0.5, c = 0.8, N= 800 and d = 40
    W = load("W_05_08_40.jld2","W")
    N = 800 ; 
  dn= 40 ; db = interval(big(dn)); d = interval(dn)
c = interval(0.8)
Tn = 0.5; T= abs(interval(Tn)) ; Tb = abs(interval(big(Tn))) 
### candidate values for a, σ0 and σ1. These values were obtained studying the graph of the function |mT-c| 
####values for T=0.5 and c=0.8
a = 1;  ab = interval(big(a)) ; a = interval(a)
σ1 = interval(0.1)
σ0 = interval(0.08)

fourier = CosFourier(N, π/d)
W = Sequence(fourier,coefficients(W))

# Choose ν = T if T>0 and ν = 1 if T=0
if T>0
    ν = T
    νb = Tb
else
    ν = interval(4/π^2)
    νb = interval(big(4/π^2))
end


##################### VERIFICATION OF THE CONSTANTS a, σ0 and σ1 AS DESCRIBED IN APPENDIX 7.2 #########################################

# We start by some analysis around the point 0. We prove that |m_T(ξ)-c|≥ σ0 for all ξ = x+iy and |x|≤xmin, |y|≤ymin
# Then we use that m_T(0)=1 and the series extension of |tanh(ξ)/ξ| = |1-exp(-2ξ)/(ξ+ξexp(-2ξ))| at zero to find that |m_T(ξ)-c| \geq ||m_T|-|c|| ≥ |1-c|- ϵ and ϵ is given by f below. In particular, if xmin and ymin are small enough, then |m_T(ξ)-c|≥ σ0 will be satisfied and we display a verified message. Otherwise the code breaks.

xmin = interval(0.12)
ymin = interval(0.12)

K = 10
f = interval(2)
if T==0 #if T=0 and c>1, we have c>m_T so we want f such that |m0| ≤ f and so |c-mT| ≥ c-|mT| ≥ c-f
  for k=1:K
    global f
    f = f + 2^(k+1)*abs(xmin+im*ymin)^k/(interval(factorial(k+1)))
  end
  f = sqrt((f+ 1/interval(factorial(K+1))*2^(K+1)*abs(xmin+im*ymin)^K)*(1+T*abs(xmin+im*ymin)^2)/(1+cos(2*xmin)))

  if c-f<0
    display("the constants xmin and ymin do not satisfy the required conditions")
    Breakdown = BreakNow
    # if the proof fails, we can decrease xmin and ymin for it to work out
  end

else #if T>0 and c<1, we have c<m_T so we want f such that |mT| ≥ f and so |c-mT| ≥ |mT|-c ≥ f-c
  for k=1:K
    global f
    f = f - 2^(k+1)*abs(xmin+im*ymin)^k/(interval(factorial(k+1)))
  end
  f = sqrt((f- 1/interval(factorial(K+1))*2^(K+1)*abs(xmin+im*ymin)^K)*(1+T*abs(xmin+im*ymin)^2)/(2+sin(2*xmin)))

  if f-c>0
    display("the constants xmin and ymin satisfy the required conditions")
  else
    display("the constants xmin and ymin do not satisfy the required conditions")
    Breakdown = BreakNow
    # if the proof fails, we can decrease xmin and ymin for it to work out
  end
end

x = a
P1 = -1.0
P2 = -1.0
while (inf(P1)<=0)&&(inf(P2)<=0)
  if T>0
    global x, P1, P2
    x = x+interval(1.0)
    P1 = x^2 - sqrt((x^2+a^2)/(T^2*(cosh(2x)-1)/(1+cosh(2x))))*(c+σ0)^2
    P2 = (cosh(2x)-1)/(cosh(2x)+1)*T^2*x^4/(x^2+a^2)-2^4*c^4
  else
    global x, P1, P2
    x = x+interval(1.0)
    P1 = x^2 - 1+abs(cos(2a))/(cosh(2x)*(c+σ0)^4*(1-abs(cos(2a))/cosh(2x)))
    P2 = 1
  end
end

######## We verify that the constants satisfy the desired inequalities
di = 0.001
dj = 0.001
P,m1 = test_constants(T,ν,c,a,σ0,σ1,x,di,dj,xmin,ymin)

if P==1
  display("the constants a, σ0 and σ1 satisfy the required conditions")
else
  display("the constants a, σ0 and σ1 do not satisfy the required conditions")
  Breakdown = BreakNow
end


######################################################################

#Construction of the fourier series of cosh(2ax) on Ω
fourierE = CosFourier(2*N, π/d)

Eb = Sequence(fourierE ,interval.(big.(zeros(2*N+1))))
for n = 0:2N  
    Eb[n] = ab/db*(-1)^n*(1-exp(-4*ab*db))/(4*ab^2 + (n*pi/db)^2) 
end

E = Sequence(fourierE ,interval.((zeros(2*N+1))))
for n = 0:2N  
    E[n] = a/d*(-1)^n*(1-exp(-4*a*d))/(4*a^2 + (n*pi/d)^2)
end

#################################################

# Construction of the operator Lν
D² = project(Derivative(2), fourier, fourier,Interval{Float64})
Lν = I - ν*D²

# Construction of the operator L
dd = diag(coefficients(D²))  ;   dd = sqrt.(-dd) ; d0=dd;  d0[1] = interval(1);
dL = tanh.(dd).*(ones(N+1)+T*dd.^2)./d0 ;  dL[1] = interval(1) ; dL = sqrt.(dL)
MT = LinearOperator(fourier,fourier,Diagonal(dL))

L = Lν*(MT - c*I)
Linv = LinearOperator(fourier,fourier,Diagonal(ones(N+1)./diag(coefficients(L))))

D1 = convert(Vector{Interval{Float64}},interval.(exp2cos(N)))
D2 = ones((N+1))./D1

########### Projection on zero trace functions (Section 3.1) #############################
setprecision(Precis)
D = interval.(big.(diag(coefficients(mid.(Linv))).^10))
S = trace(N); C = S' ;
S =  interval.(big.(S)) ;  C =  interval.(big.(C))

V = interval.(big.(coefficients(W)))

W0 = (D.*C)*solve(S*(D.*C),vec(S*V))
#
U0b = interval.(big.(W)) - Sequence(fourier,vec(W0))
U0 = Interval.(Float64.(inf.(U0b),RoundDown),Float64.(sup.(U0b),RoundUp) )

V2 = 2*U0
V1 = -4*ν*sqrt.(-D²)*U0
V0 = -2*ν*D²*U0


### Computation of epsilon (Assumption 2)
if T==0
  ϵ = c/2 - norm(U0,1)
  display("value of ϵ")
  display(inf(ϵ))
  if inf(ϵ)<=0
    display("epsilon is negative, the proof will fail")
    Breakdown = BreakNow
  end
end


############ Computation of A (Section 3.2) ########################################

DG = Lν*project(Multiplication(V2),fourier, fourier,Interval{Float64})
B = interval.(inv(I + mid.(DG)*mid.(Linv)))
Badj = LinearOperator(fourier,fourier,coefficients(B)')

# Construction of WT
e0 = Sequence(fourier ,interval.(zeros(N+1)))
e0[0] = 1
if T>0
    WT = e0
else
  e02 = Sequence(CosFourier(2*N, π/d) ,interval.(zeros(2*N+1)))
  e02[0] = 1
  WT = interval.(mid.(project(Multiplication(mid.(e0 - 1/mid(c)*V2)),fourier, CosFourier(2*N, π/d),Interval{Float64}))\mid.(e02))
end

# computation of the norm of B
MWT = project(Multiplication(WT*WT),fourier, fourier,Interval{Float64}) - project(Multiplication(WT),fourier, fourier,Interval{Float64})*project(Multiplication(WT),fourier, fourier,Interval{Float64})


n_WT = sqrt(opnorm(LinearOperator(coefficients(D1.*(MWT).*D2')),2))   # computation of \|π_NW_Tπ^N\|_2
n_BN = sqrt(opnorm(LinearOperator(coefficients(D1.*(B*Badj).*D2')),2))   #computation of \|B^N_T\|_2
norm_B = maximum([1 maximum([norm(WT,1) n_BN])+n_WT])  #computation of the norm of B

display("norm of B")
display(norm_B)

MWT = Nothing

################# Computation of Z1 ####################################

W0 = WT*V0 ;              
W1 = WT*V1 ;  
W2 = WT*V2 ;          

if T>0
  L0 = Linv ;             l0 = (1+ν*((N+1)*π/d)^2)*(maximum(vec((abs.(dL.-c)))))
  L1 = Linv*sqrt.(-D²) ;  l1 = (1+ν*((N+1)*π/d)^2)*(maximum(vec((abs.(dL.-c)))))/((N+1)*π/d)
  L2 = Linv*Lν ;          l2 = maximum(vec((abs.(dL.-c))))
else
  L0 = Linv ;             l0 = (1+ν*((N+1)*π/d)^2)*(maximum(vec((abs.(dL.-c)))))
  L1 = Linv*sqrt.(-D²) ;  l1 = (1+ν*((N+1)*π/d)^2)*(maximum(vec((abs.(dL.-c)))))/((N+1)*π/d)
  L2 = Linv*Lν + 1/c*I ;  l2 = maximum(vec(abs.((c*(dL.-c))./dL)))
end

#MWi corresponds to the operator πᴺ(Vi*WT)(π^{3N}-πᴺ)(Vi*WT)πᴺ = πᴺ(Vi*WT*Vi*WT)πᴺ - (πᴺ(Vi*WT)πᴺ)^2 (recall that Wi = Vi*WT)
MW0 = project(Multiplication(W0*W0),fourier, fourier,Interval{Float64}) - project(Multiplication(W0),fourier, fourier,Interval{Float64})^2
MW1 = project(Multiplication(W1*W1),fourier, fourier,Interval{Float64}) - project(Multiplication(W1),fourier, fourier,Interval{Float64})^2
MW2 = project(Multiplication(W2*W2),fourier, fourier,Interval{Float64}) - project(Multiplication(W2),fourier, fourier,Interval{Float64})^2

#MVi corresponds to the operator πᴺVi(π^{2N}-πᴺ)Viπᴺ = πᴺ(Vi*Vi)πᴺ - (πᴺViπᴺ)^2
MV0 = project(Multiplication(V0*V0),fourier, fourier,Interval{Float64}) - project(Multiplication(V0),fourier, fourier,Interval{Float64})^2
MV1 = project(Multiplication(V1*V1),fourier, fourier,Interval{Float64}) - project(Multiplication(V1),fourier, fourier,Interval{Float64})^2
MV2 = project(Multiplication(V2*V2),fourier, fourier,Interval{Float64}) - project(Multiplication(V2),fourier, fourier,Interval{Float64})^2

# computation of each component of Z1
Z11 =  opnorm(LinearOperator(coefficients(D1.*(I - B*(I+DG*Linv)).*D2')),2)^2
Z11 = Z11 + ( sqrt(opnorm(LinearOperator(coefficients(D1.*(L0*MW0*L0).*D2')),2)) + sqrt(opnorm(LinearOperator(coefficients(D1.*(L1*MW1*L1).*D2')),2)) + sqrt(opnorm(LinearOperator(coefficients(D1.*(L2*MW2*L2).*D2')),2)) )^2
Z11 = sqrt(Z11)

Z12 = sqrt(opnorm(LinearOperator(coefficients(D1.*B*MV0*Badj.*D2')),2))/l0 + sqrt(opnorm(LinearOperator(coefficients(D1.*B*MV1*Badj.*D2')),2))/l1 + sqrt(opnorm(LinearOperator(coefficients(D1.*B*MV2*Badj.*D2')),2))/l2

Z13 = norm(W0,1)/l0 + norm(W1,1)/l1 + norm(W2,1)/l2

if T>0
  Z14 = 0
else
  Z14 = norm(e0 - WT*(e0-1/c*V2),1)
end

Z1 = sqrt(Z11^2 + Z12^2 + Z13^2 + Z14^2)

display("value of Z1")
display(Z1)


MV0 = Nothing ; MV1 = Nothing ; MV2 = Nothing; MW0 = Nothing; MW1 = Nothing; MW2 = Nothing

######################## Computation of Zu ##############################


  ##### Construction of the needed constants
# Construction of ξ0 
Ca = (1 + abs(cos(2*a)))/(1-abs(cos(2*a)))

if T>0
  ξ0 = maximum([interval(1) 1/sqrt(T) c^2*Ca*4/T])
  ξ0 = maximum([3*T/(2*tanh(ξ0)) 2*c^2/(T*tanh(ξ0))])
else
  ξ0 = interval(1)
end

# Constants for the bound Y0
if T>0
  a0 = minimum([1/sqrt(T) interval(1/2)*pi])
  a0b = minimum([1/sqrt(Tb) interval(big(1/2))*pi])
else
  a0 = interval(1)
  a0b = interval(big(1))
end


# Notice that max_{s>0}min{√2 + √s, √2/(1-exp(-πs)) } ≤ max{1+√2, √2/(1-exp(-π))}, so we can compute CY0 as follows (cf. Lemma 4.1)
pi_i = interval(π)
if T==0
  CY0 = maximum([1+sqrt(interval(2)) sqrt(interval(2))/(1-exp(-pi_i))])
else
  CY0 =maximum([1+sqrt(interval(2)) sqrt(interval(2))/(1-exp(-pi_i))])*1/(sqrt(pi_i)*T^(1/4))*(2/sqrt(1+a0*sqrt(T)) + 1/sqrt(interval(2)))
end

# Other constants of Lemma 4.1
C0 = 1/(π*σ0*(1-ν*a^2)) + 1/(π*ν*σ0)

if T>0
  C1 = 1/(2*π)*( 2*(1+a)/(σ0*(1-ν*a^2)) + 4*(1+a)/(σ1*sqrt(T)*ν) )
  K1 = 2*ξ0/(π*σ0) + 2*sqrt(ξ0)*(1+abs(c))/(π*σ0*sqrt(T)) + 2*(1/(3*T) + abs(c)/(4*T^(3/2)) + 2*c^2/T )/(π*sqrt(tanh(ξ0)*T*ξ0)) + abs(c)/(T*π)*(2+3*log(ξ0)) + 1/sqrt(2π*T) 
  K2 = Ca*abs(ξ0+im*a)/( 2π*σ0^2*(1-T*a^2)^2)*( 1/a^2 + T + Ca/a + Ca*T*abs(ξ0+im*a) ) + 2*Ca^2/π*( 2*(1+T)/sqrt(T*ξ0) + (2+a)/2*exp(-2ξ0))
  C2 = maximum([K2 K1*exp(a)])  
else
  C1 = 1/(2*abs(c)*ν) + Ca^(1/4)*(1+sqrt(a))/(π*abs(c)*σ0)*( 1/(1-ν*a^2) + 2/ν )
  K1 = (2+4*exp(-2*ξ0))/(minimum([1 abs(c)^3])*π*sqrt(ξ0)*σ0) + 1/(π*σ0*abs(c)) + 2/(π*c^2) + 1/(π*abs(c)^3)*(2+3*log(3*ξ0)) + 1/(c^2*sqrt(2*π))
  K2 = Ca^(1/4)/π*( 1/(2*σ0^2)*( 2 + a^(-3/2)) + 1/(4*σ0^2*(1-abs(cos(2a)))^2*sqrt(a)) )
  C2 = maximum([K2 K1*exp(a)]) 
end






###################### Computation of 𝒵u #####################################

#### Computation of the constants C(d) and C1(d) in Lemmas 4.2 and 4.6

C1d = exp(-2*a*d)*( 2*sqrt(π)/(sqrt(4*a*d)*(1-exp(-2*a*d))) + 4/(1-exp(-2*a*d)) )
Cd = exp(-2*a*d)*( 4*d + 4*exp(-a*d)/(a*(1-exp(-3*a*d/2))) + 2/(a*(1-exp(-a*d*2))))

### Inner product (V0,E*V0)
PS0 = abs( (coefficients(D1.*project(V0*E,fourier))')*coefficients(D1.*V0) )

### Inner product (V1,Efull*V1)
k = (-N:N)*π/db
fourier_f = Fourier(N,π/db)
V1f = 4*νb*coefficients(U0b) ; V1f = [reverse(V1f[2:N+1]); V1f[1] ; V1f[2:N+1]].*k ; V1f = Sequence(fourier_f,V1f)
Ef = coefficients(Eb); Ef = [reverse(Ef[2:2*N+1]); Ef[1] ; Ef[2:2*N+1]] ; Ef = Sequence(Fourier(2*N,π/db),Ef)
PS1 = abs( (coefficients(project(V1f*Ef,fourier_f))')*coefficients(V1f) )

### Inner product (V2,E*V2)
PS2 = 8*abs( (coefficients(D1.*project(U0b*Eb,fourier))')*coefficients(D1.*U0b) )

# Computation of ∫|v2'|^2 from d-1 to d
char = Sequence(Fourier(2*N,π/db) ,interval.(big.(zeros(4*N+1))))
for n = [-2N:-1;1:2N]  
    char[n] = real((-1)^n*(exp(im*n*π/db)-1)/(2*im*n*π))
end
char[0] = 1/(2*db)
Elog =  abs( (coefficients(project(V1f*char,fourier_f))')*coefficients(V1f) )

# Computation of the Zu_i
Zu0 = 4d*C0^2/a*PS0 + 2*d*Cd*C0^2*PS0
Zu1 = 4d*C1^2/a*PS1 + 2*d*Cd*C1^2*PS1
Zu2 = 2*C2^2*( 2d*PS2/a + 4*log(interval(2))*Elog ) + 2*d*C1d*C2^2*PS2

𝒵u = norm_B*sqrt( Zu0 + Zu1 + Zu2 )

display("value of 𝒵u")
display(𝒵u)

######################## Computation of ℨ1 ##############################


𝒵1 = Z1 + 𝒵u
display("value of 𝒵1")
display(𝒵1)


################### Computation of Z2 ######################
if (T==0)||(T>1/3)
  κT = 1/(sqrt(ν)*abs(1-c)^2)
else
  κT = 1/(sqrt(ν)*σ0^2)
end
𝒵2 = 2*κT*norm_B

display("value of 𝒵2")
display(𝒵2)


#################### Computation of Y0 ##################
setprecision(Precis)
Cd = exp(-2*a0*d)*( 4*d + 4*exp(-a0*d)/(a0*(1-exp(-3*a0*d/2))) + 2/(a0*(1-exp(-a0*d*2))))

D1 = convert(Vector{Interval{BigFloat}},interval.(big.(exp2cos(N))))
D2 = ones((N+1))./D1
D² = project(Derivative(2), CosFourier(N, π/db), CosFourier(N, π/db),Interval{BigFloat})
Lνb = I - νb*D²
Lνb = diag(coefficients(Lνb)).*2


Eb = Sequence(fourierE ,interval.(big.(zeros(2*N+1))))
for n = 1:2N  
    Eb[n] = real((-1)^n*(1/(im*interval(big(n))/db*interval(big(pi)) + 2a0b) - 1/(im*interval(big(n))/db*interval(big(pi))-2*a0b)))
end
Eb[0] = 1/a0b
Eb = Eb*(1-exp(-4*a0b*db))/(4*db)

### Inner product part 
PSY0 = abs( (coefficients(D1.*project(((Lνb.*U0b)*Eb),fourier))')*coefficients(D1.*(Lνb.*U0b)) )

Yu1 = sqrt(8*d^2*CY0*PSY0) ; Yu2 = sqrt(4*d^2*CY0*PSY0*Cd)

FU0 = L*U0 + Lν*(U0*U0)
Y0 = sqrt(2d)*sqrt( norm(B*project(FU0,fourier),2)^2 + norm(WT*FU0-project(WT*FU0,fourier),2)^2)

𝒴0 = Y0 + sqrt(Yu1^2 + norm(B)^2*Yu2^2)
display("value of 𝒴0")
display(𝒴0)


if 1- 𝒵1>0
  if (1-𝒵1)^2-4*𝒴0*𝒵2 > 0
    rmin=abs((1- 𝒵1 - sqrt((1- 𝒵1)^2-4*𝒴0*𝒵2))/(2*𝒵2))
    rmax=abs((1- 𝒵1 + sqrt((1- 𝒵1)^2-4*𝒴0*𝒵2))/(2*𝒵2))
    if sup(rmin)<inf(rmax)
      display("The proof was successful, minimal radius of the contraction :")
      display(Float64(sup(rmin)))
      display("Maximal radius of the contraction :")
      display(Float64(inf(rmax)))
    else
      display("rmin>=rmax")
    end
  else
    display("failure: discriminant is negative")
  end
else
    display("failure: linear term is positive")
end

######### Regularity for the solution (using Proposition 4.1)

if T==0
  if rmin <= 4*ϵ*sqrt(ν)*minimum([abs(c-1),abs(c)])
    display("the solution is infinitely differentiable")
  else
    display("the solution is H^(3/2) (continuous at least), but could not prove more regularity")
  end
end



###############################################################################################################################
###########################        PROOF OF STABILITY       ##################################################################

# Passage from cosine series to exponential series
fourier_f = Fourier(N,π/d)
U0f = coefficients(U0) ; U0f = [reverse(U0f[2:N+1]); U0f[1] ; U0f[2:N+1]] ; U0f = Sequence(fourier_f,U0f)
U0fb = coefficients(U0b) ; U0fb = [reverse(U0fb[2:N+1]); U0fb[1] ; U0fb[2:N+1]] ; U0fb = Sequence(fourier_f,U0fb)


# Construction of the operator L
D² = real.(project(Derivative(2), fourier_f, fourier_f,Complex{Interval{Float64}}))
dd = diag(coefficients(D²))  ;   dd = sqrt.(-dd) ; d0=dd;  d0[N+1] = interval(1);
dL = tanh.(dd).*(ones(2*N+1)+T*dd.^2)./d0 ;  dL[N+1] = interval(1) ; dL = sqrt.(dL)
MT = LinearOperator(fourier_f,fourier_f,Diagonal(dL))

if T>0
    L = MT - c*I
    Linv = LinearOperator(fourier_f,fourier_f,Diagonal(ones(2*N+1)./diag(coefficients(L))))
    DG = 2*real.(project(Multiplication(U0f),fourier_f, fourier_f,Complex{Interval{Float64}}))
else
  L = c*I - MT
  Linv = LinearOperator(fourier_f,fourier_f,Diagonal(ones(2*N+1)./diag(coefficients(L))))
  DG = -2*real.(project(Multiplication(U0f),fourier_f, fourier_f,Complex{Interval{Float64}}))
end

E = eigen(Hermitian(coefficients(mid.(L) + mid.(DG))),1:2)

pspace = ParameterSpace()×fourier_f

# Construction of the eigencouples
λ1 = E.values[1] ; λ2 = E.values[2] # Theoretically we expect the second eigenvalue to be exactly zero with eigenvector u' 
V1 = E.vectors[:,1] ; V2 = E.vectors[:,2] 

  # projection of V on trace zero (V correspond to Ψ in Section 5)
  S = trace_full(N,db); Sᵀ = S' ;
  ddd= mid(d/π)*ones(2N+1)./(1:2*N+1)
  Dtrace = Matrix(interval.(big.(Diagonal(ddd)^2))) ; ddd = Nothing

  V1b = Sequence(fourier_f,interval.(big.(V1))-Dtrace*Sᵀ*inv(S*Dtrace*Sᵀ)*S*interval.(big.(V1)))
  V1 = Interval.(Float64.(inf.(real.(V1b)),RoundDown),Float64.(sup.(real.(V1b)),RoundUp) ) + im*Interval.(Float64.(inf.(imag.(V1b)),RoundDown),Float64.(sup.(imag.(V1b)),RoundUp) )   ; 

  V2b = Sequence(fourier_f,interval.(big.(V2))-Dtrace*Sᵀ*inv(S*Dtrace*Sᵀ)*S*interval.(big.(V2)))
  V2 = Interval.(Float64.(inf.(real.(V2b)),RoundDown),Float64.(sup.(real.(V2b)),RoundUp) ) + im*Interval.(Float64.(inf.(imag.(V2b)),RoundDown),Float64.(sup.(imag.(V2b)),RoundUp) )   ; Dtrace = Nothing ; S= Nothing; Sᵀ = Nothing

  U1 = Sequence(pspace,vec([interval(λ1);coefficients(V1)])) ; U2 = Sequence(pspace,vec([interval(λ2);coefficients(V2)]))
# Proof of the eigencouples 
P1 = proof_eigen(U1,V1b,U0f,c,T,ν,d,rmin)
P2 = proof_eigen(U2,V2b,U0f,c,T,ν,d,rmin)

P = no_eigen_outside(P1[2],P2[2],U1,U2,U0f,rmin,c,T,ν)

#verify that 0 is the second eigenvalue (we know that zero is an eigenvalue, we check that it is the second one)

if P==1
  if (sup(component(real(U2),1)[1]-P2[2])<0)&&(inf(component(real(U2),1)[1]+P2[2])>0)
    display("zero is the second eigenvalue")
    display("There is no eigenvalue in between 0 and λi")
    display("Moreover, there is no eigenvalue smaller than λi")
    display("Finally, the soliton is stable !")
  else
    display("there exists a second negative eigenvalue")
  end
else
  display("The inverse becomes singular for some value of λ")
end



display("Values in the computer-assisted proof of the soliton")

header1 = (["||DF(u0)^{-1}||_{2,l}", "Z1", "Zu ", "Y0 ", "Minimal radius of contraction", "Maximal value of contraction"])
data = [Float64(sup(norm_B/(1-𝒵1))) Float64(sup(𝒵1)) Float64(sup(𝒵u)) Float64(sup(𝒴0)) Float64(sup(rmin)) Float64(inf(rmax)) ]

pretty_table(data;
             header = header1
             )


display("Values in the computer-assisted proof of the eigencouples")

header = (["Eigenvalue", "Minimal radius of contraction", "Radius of uniqueness"])
data = [mid(λ1)  sup(P1[3]) inf(P1[2]);
        mid(λ2)  sup(P2[3]) inf(P2[2])]
             
pretty_table(data;
            header = header )
