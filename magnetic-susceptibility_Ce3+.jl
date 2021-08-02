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

E = kB * Δ
println("E = $E")

print("χ_free = "); flush(stdout); display(χ_free.(gj , μB , J , kB , T))

"""
モル磁化率
"""
χ_free_mol(N_A , gj , μB , J , kB , T) = N_A * χ_free(gj , μB , J , kB , T)

print("χ_free_mol = "); flush(stdout); display(χ_free_mol.(N_A , gj , μB , J , kB , T))


plot(T , inv.(χ_free_mol.(N_A , gj , μB , J , kB , T)) , title = "χ_free_mol" , label="1/χ_free_mol" , legend =:outerright , ls=:dash)

savefig("output-file/χ_free_mol.pdf")
savefig("output-file/χ_free_mol.png")


"""
Γ7 基底状態
"""
Co_seven(E , T) = (5 + 26 * exp(-E/(kB*T)) + 32 * kB * T / Δ * (1 - exp(-E/(kB*T)))) / (21 * (1 + 2*exp(-E/(kB*T))))
χ(T) = Co_seven.(E , T) * χ_free_mol.(N_A , gj , μB , J , kB , T)

plot(T , inv.(χ_free_mol.(N_A , gj , μB , J , kB , T)) , xlabel = "T(K)" , ylabel = "1/χ" , title = "Γ7 ground state" , label="1/χ_free" , legend =:topleft , ls=:dash , color =:black)

plot!(T , inv.(χ.(T)) , label="1/χ_Γ7" , color =:red)

plot!(twinx() , χ_free_mol.(N_A , gj , μB , J , kB , T) , label = "χ_free" , legend=:bottomright , ylabel = "χ(emu/mol)" , color =:blue)

plot!(T , χ.(T) , label = "χ_Γ7" , color =:green)
plot!(right_margin = 10mm)

savefig("output-file/myplot_Gamma7-GS.pdf")
savefig("output-file/myplot_Gamma7-GS.png")

"""
Γ8 基底状態
"""
Co_eight(E , T) = (26 + 5 * exp(-E/(kB*T)) + 32 * kB * T / Δ * (1 - exp(-E/(kB*T)))) / (21 * (2 + exp(-E/(kB*T))))
χ(T) = Co_eight.(E , T) * χ_free_mol.(N_A , gj , μB , J , kB , T)

plot(T , inv.(χ_free_mol.(N_A , gj , μB , J , kB , T)) , xlabel = "T(K)" , ylabel = "1/χ" , title = "Γ8 ground state" , label="1/χ_free" , legend =:topleft , ls=:dash , color =:black)

plot!(T , inv.(χ.(T)) , label="1/χ_Γ8" , color =:red)

plot!(twinx() , χ_free_mol.(N_A , gj , μB , J , kB , T) , label = "χ_free" , legend=:bottomright , ylabel = "χ(emu/mol)" , color =:blue)

plot!(T , χ.(T) , label = "χ_Γ8" , color =:green)

plot!(right_margin = 10mm)

savefig("output-file/myplot_Gamma8-GS.pdf")
savefig("output-file/myplot_Gamma8-GS.png")


"""
1つのグラフに出力
"""

χ(T) = Co_seven.(E , T) * χ_free_mol.(N_A , gj , μB , J , kB , T)

plot(T , inv.(χ_free_mol.(N_A , gj , μB , J , kB , T)) , 
    xlabel = "T(K)" ,ylabel = "1/χ" , title = L"\texttt{magnetic~susceptibility~of}~Ce{^3+}" , label="1/χ_free" ,legend =:topleft , ls=:dash)
plot!(T , inv.(χ.(T)) , label="1/χ_Γ7")
plot!(twinx() , χ_free_mol.(N_A , gj , μB , J , kB , T) , label = "χ_free" , legend=:bottomright , ylabel = "χ")
plot!(T , χ.(T) , label = "χ_Γ7")
plot!(right_margin = 10mm)

χ(T) = Co_eight.(E , T) * χ_free_mol.(N_A , gj , μB , J , kB , T)

plot!(T , inv.(χ.(T)) , label="1/χ_Γ8")
plot!(T , χ.(T) , label = "χ_Γ8")
plot!(right_margin = 10mm)

savefig("output-file/myplot_both.pdf")
savefig("output-file/myplot_both.png")