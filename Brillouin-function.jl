using Plots
using LaTeXStrings

B(x) = (2*J+1)/(2*J)*coth((2*J+1)/(2*J)*x) - 1/(2*J)*coth(x/(2*J))
x = -3:0.1:3

p= plot()

for n in 1:20
    for J in n/2
        B(x) = (2*J+1)/(2*J)*coth((2*J+1)/(2*J)*x) - 1/(2*J)*coth(x/(2*J))
        plot!(x,B,xlabel=L"x",ylabel=L"B_J(x)",label="J=$n/2",legend=:outertopright)
    end
end

p
#savefig("Brillouin-function.png")