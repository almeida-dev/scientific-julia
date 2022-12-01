include("gauss_seidel.jl") #address of file that contains external functions
using LinearAlgebra, Plots, DataFrames

# heat parameters
k = 190.0
h = 5000.0
q₀ = 1e5
Tₑₙᵥ = 15.0
tol = 1e-5

# optimized relaxation parameter by Optim.jl
λ = 1.798374

# iteration counting
iterations = 0

# creating 6x5 zero matrix as initial guess
T = ones(6, 5)

for i ∈ eachindex(T[:1])

    for j ∈ eachindex(T[1])

        T[i, j] = 26.25 #optimized temperature guess

    end

end

# iterating gauss_seidel! function for matrix T
while true

    xₙ = copy(T)

    gauss_seidel!(T, k, h, q₀, Tₑₙᵥ) #gauss seidel iteration
    
    global iterations += 1
    
    norm(T - xₙ, Inf) > tol || break #infinite norm and tolerance
    
    xₙ₊₁ = relax(T, xₙ, λ) # relaxation
    
    global T = xₙ₊₁ # updating iteration

end

# printing results as txt file
result = """
Numerical parameters:
k = $k W/mK
h = $h W/m²K
q"₀ = $q₀ W/m²
Tₑₙᵥ = $Tₑₙᵥ °C
tolerance = $tol
iterations = $iterations """

write("ibmresults.txt", result)

Tc = copy(T)
df = DataFrame(T, :auto)
allowmissing!(df)

df[3, 1] = missing #missing values for unnecessary data
df[3, 2] = missing

iostr = sprint(show, df)
write("Temperatures.txt", iostr)

# ploting annotated heatmap
T[3, 1] = 10.0
T[3, 2] = 10.0
S = reverse(T, dims=1)
texts = [(j, i, text(round(S[i, j], digits=2), 8, :black, :center))
         for i in 1:6 for j in 1:5]

heatmap(S, aspect_ratio=1.6, color=:thermal, xlim=(0, 6))
annotate!(texts, linecolor=:blue)