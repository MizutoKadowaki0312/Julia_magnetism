using Plots
using LaTeXStrings

"""
χ_free の Plot
"""

χ_free(gj , μB , J , kB , T) = ((gj * μB)^2 * J * (J + 1))/(3 * kB * T)
gJ(S , L , J) = 3/2 + (S*(S+1) - L*(L+1))/(2*J*(J+1))

"""
変数の設定
"""

S = 1/2
L = 3

if S <= 7/2
    J = L - S
else
    J = L + S
end

println("S = $S")
println("L = $L")
println("J = $J")

gj = gJ(S , L , J)
println("gj = $gj")

μB = 9.284764 * 10^(-24)
println("μB = $μB")

kB = 1.380649 * 10^(-23)
println("kB = $kB")

N_A = 6 * BigInt(10)^(23)
println("N_A = $N_A")

T = range(0,300,length = 30001)

χ_free_mol(N_A , gj , μB , J , kB , T) = N_A * χ_free(gj , μB , J , kB , T)

f(T) = 1/χ_free_mol.(N_A , gj , μB , J , kB , T)

plot(T , f , xlabel = L"T(K)" , ylabel = L"1/\chi" , label = L"1/\chi_{free}" , legend =:topleft)

"""
Γ7 Ground State
"""

Γ_seven(Δ,T) = (5 + 26*exp(-Δ/(kB*T)) + 32*kB*T*(1-exp(-Δ/(kB*T)))/Δ) / (21 * (1 + 2*exp(-Δ/(kB*T))))

χ_Γ7(gj , μB , J , kB , N_A , T , Δ) = χ_free_mol(N_A , gj , μB , J , kB , T) * Γ_seven(Δ,T)

Δ = 200

g(T) = 1/χ_Γ7.(gj , μB , J , kB , N_A , T , Δ)

plot!(T,g,label = L"1/\chi_{\chi_{\Gamma_7}}",title="Δ = $Δ")