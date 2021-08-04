using Plots
gr()
using Plots.PlotMeasures

"""
1つのイオンの磁化率
"""

χ_free(gj , μB , J , kB , T) = ((gj * μB)^2 * J * (J + 1))/(3 * kB * T)

"""
Landeのg因子
"""
gJ(S , L , J) = 3/2 + (S*(S+1) - L*(L+1))/(2*J*(J+1))

"""
変数の設定
"""

S = 1
L = 5

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

T = range(0,300,length = 301)

Δ = 73
println("Δ = $Δ")

E = kB * Δ
println("E = $E")

"""
モル磁化率
"""
χ_free_mol(N_A , gj , μB , J , kB , T) = N_A * χ_free(gj , μB , J , kB , T)
display(χ_free_mol.(N_A , gj , μB , J , kB , T))


plot(T , inv.(χ_free_mol.(N_A , gj , μB , J , kB , T)) , title = "χ_free_mol" , label="1/χ_free_mol" , legend =:outerright , ls=:dash)

#savefig("output-file/χ_free_mol_Pr3+.pdf")
#savefig("output-file/χ_free_mol_Pr3+.png")

"""
分配関数
"""
Z(kB , E , T) = 1 + 3*exp(-E/(kB*T))
"""
display(Z.(kB , E , T))
plot(T , Z.(kB , E , T))
"""

"""
磁化率
"""

Co_one(kB , E , T) = (3/40 * exp(-E/(kB*T)) + (2*kB*T/E)*(1 - exp(-E/(kB*T))))/Z(kB , E , T)


χ(T) = Co_one.(kB , E , T) * χ_free_mol.(N_A , gj , μB , J , kB , T)

plot(T , inv.(χ_free_mol.(N_A , gj , μB , J , kB , T)) , xlabel = "T(K)" , ylabel = "1/χ" , title = "Γ1 ground state" , label="1/χ_free" , legend =:topleft , ls=:dash , color =:black)

plot!(T , inv.(χ.(T)) , label="1/χ_Γ1" , color =:red)

plot!(twinx() , χ_free_mol.(N_A , gj , μB , J , kB , T) , label = "χ_free" , legend=:bottomright , ylabel = "χ(emu/mol)" , color =:blue)

plot!(T , χ.(T) , label = "χ_Γ1" , color =:green)
plot!(right_margin = 10mm)

#savefig("output-file/myplot_Gamma1-GS_Pr3+.pdf")
#savefig("output-file/myplot_Gamma1-GS_Pr3+.png")

plot(T , χ.(T) , label = "χ_Γ1" , color =:green)
#savefig("output-file/myplot_Gamma1-only_Pr3+.pdf")
#savefig("output-file/myplot_Gamma1-only_Pr3+.png")


"""
比熱
"""
R = 8.3
C(R , Z , E , kB , T) = R*(E/(kB*T))^2/(Z(kB , E , T))^2 * (3*Z(kB , E , T)*exp(-E/(kB*T)) - 9* exp(-2*(E/(kB*T))))

display(C.(R , Z , E , kB , T))
plot(T , C.(R , Z , E , kB , T) , xlabel = "T(K)" , ylabel = "C(T)")

savefig("output-file/myplot_heat-capacity_Pr3+.pdf")
savefig("output-file/myplot_heat-capacity_Pr3+.png")