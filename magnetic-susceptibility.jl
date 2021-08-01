using Plots
using Plots.PlotMeasures
using LaTeXStrings

"""
Ce-ion(3+) について，Γ7基底状態とΓ8基底状態それぞれについて
磁化率をPlotする．
"""

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

T = range(0,300,length = 301)

Δ = 200
println("Δ = $Δ")

print("χ_free = "); flush(stdout); display(χ_free.(gj , μB , J , kB , T))

"""
モル磁化率
"""
χ_free_mol(N_A , gj , μB , J , kB , T) = N_A * χ_free(gj , μB , J , kB , T)

print("χ_free_mol = "); flush(stdout); display(χ_free_mol.(N_A , gj , μB , J , kB , T))


plot(T , inv.(χ_free_mol.(N_A , gj , μB , J , kB , T)) , title = "χ_free_mol" , label="1/χ_free_mol" , legend =:outerright , ls=:dash)

savefig("χ_free_mol.pdf")
savefig("χ_free_mol.png")


"""
Γ7 基底状態
"""
Co_seven(Δ , kB , T) = (5 + 26 * exp(-Δ/(kB*T)) + 32 * kB * T / Δ * (1 - exp(-Δ/(kB*T)))) / (21 * (1 + 2*exp(-Δ/(kB*T))))
χ(T) = Co_seven.(Δ , kB , T) * χ_free_mol.(N_A , gj , μB , J , kB , T)
scatter(T , Co_seven.(Δ , kB , T))
scatter!(T , Co_seven.(Δ , 1 , T))
p = plot()
#plot!(right_margin = 20mm)
for n in [kB , 0.1 , 0.3 , 0.5 , 0.7 , 1.0]
    plot!(T , Co_seven.(Δ , n , T) , label = "kB = $n" , legend =:topleft)
end
p

χ(T) = Co_seven.(Δ , 1 , T) * χ_free_mol.(N_A , gj , μB , J , kB , T)

plot(T , inv.(χ_free_mol.(N_A , gj , μB , J , kB , T)) , xlabel = "T(K)" , ylabel = "1/χ" , title = "Γ7 ground state" , label="1/χ_free" , legend =:topleft , ls=:dash , color =:black)
plot!(T , inv.(χ.(T)) , label="1/χ_Γ7" , color =:red)
plot!(twinx() , χ_free_mol.(N_A , gj , μB , J , kB , T) , label = "χ_free" , legend=:bottomright , ylabel = "χ(emu/mol)" , color =:blue)
plot!(T , χ.(T) , label = "χ_Γ7" , color =:green)
plot!(right_margin = 10mm)

savefig("myplot_Gamma7.pdf")
savefig("myplot_Gamma7.png")

p = plot()
plot!(T , inv.(χ_free_mol.(N_A , gj , μB , J , kB , T)) , title = "Γ7 ground state" , label="1/χ_free" , legend =:outerright , ls=:dash)

for n in [kB , 0.1 , 0.3 , 0.5 , 0.7 , 1.0]
    χ(T) = Co_seven.(Δ , n , T) * χ_free_mol.(N_A , gj , μB , J , kB , T)
    plot!(T , inv.(χ.(T)) , ylims = (0,60) , label = "kB = $n" , legend =:bottomright)
end
p
savefig("myplot_Gamma7_diff_kB.pdf")
savefig("myplot_Gamma7_diff_kB.png")


"""
Γ8 基底状態
"""
Co_eight(Δ , kB , T) = (26 + 5 * exp(-Δ/(kB*T)) + 32 * kB * T / Δ * (1 - exp(-Δ/(kB*T)))) / (21 * (2 + exp(-Δ/(kB*T))))
χ(T) = Co_eight.(Δ , kB , T) * χ_free_mol.(N_A , gj , μB , J , kB , T)

scatter(T , Co_eight.(Δ , kB , T))
scatter!(T , Co_eight.(Δ , 1 , T))


p = plot()
for n in [kB , 0.1 , 0.3 , 0.5 , 0.7 , 1.0]
    plot!(T , Co_eight.(Δ , n , T) , label = "kB = $n" , legend =:topleft)
end
p


χ(T) = Co_eight.(Δ , 1 , T) * χ_free_mol.(N_A , gj , μB , J , kB , T)
display(χ.(T))
display(inv.(χ.(T)))

plot(T , inv.(χ_free_mol.(N_A , gj , μB , J , kB , T)) , xlabel = "T(K)" , ylabel = "1/χ" , title = "Γ8 ground state" , label="1/χ_free" , legend =:topleft , ls=:dash , color =:black)
plot!(T , inv.(χ.(T)) , label="1/χ_Γ8" , color =:red)
plot!(twinx() , χ_free_mol.(N_A , gj , μB , J , kB , T) , label = "χ_free" , legend=:bottomright , ylabel = "χ(emu/mol)" , color =:blue)
plot!(T , χ.(T) , label = "χ_Γ8" , color =:green)
plot!(right_margin = 10mm)

savefig("myplot_Gamma8.pdf")
savefig("myplot_Gamma8.png")

p = plot()
plot!(T , inv.(χ_free_mol.(N_A , gj , μB , J , kB , T)) , title = "Γ8 ground state" , label="1/χ_free" , legend =:outerright , ls=:dash)

for n in [kB , 0.1 , 0.3 , 0.5 , 0.7 , 1.0]
    χ(T) = Co_eight.(Δ , n , T) * χ_free_mol.(N_A , gj , μB , J , kB , T)
    plot!(T , inv.(χ.(T)) , ylims = (0,60) , label = "kB = $n" , legend =:bottomright)
end
p

savefig("myplot_Gamma8_diff_kB.pdf")
savefig("myplot_Gamma8_diff_kB.png")


"""
1つのグラフに出力
"""

χ(T) = Co_seven.(Δ , 1 , T) * χ_free_mol.(N_A , gj , μB , J , kB , T)
χ.(T)
inv.(χ.(T))

plot(T , inv.(χ_free_mol.(N_A , gj , μB , J , kB , T)) , 
    xlabel = "T(K)" ,ylabel = "1/χ" , title = "magnetic susceptibility of Ce^3+" , label="1/χ_free" ,legend =:topleft , ls=:dash)
plot!(T , inv.(χ.(T)) , label="1/χ_Γ7")
plot!(twinx() , χ_free_mol.(N_A , gj , μB , J , kB , T) , label = "χ_free" , legend=:bottomright , ylabel = "χ(emu/mol)")
plot!(T , χ.(T) , label = "χ_Γ7")
plot!(right_margin = 10mm)

χ(T) = Co_eight.(Δ , 1 , T) * χ_free_mol.(N_A , gj , μB , J , kB , T)

plot!(T , inv.(χ.(T)) , label="1/χ_Γ8")
plot!(T , χ.(T) , label = "χ_Γ8")
plot!(right_margin = 10mm)

savefig("myplot_Gamma7-and-Gamma8.pdf")
savefig("myplot_Gamma7-and-Gamma8.png")