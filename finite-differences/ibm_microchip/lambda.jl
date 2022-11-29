include("gauss_seidel.jl")
using LinearAlgebra

# heat parameters
k = 190.0
h = 5000.0
q₀ = 1e5
Tₑₙᵥ = 15.0
tol = 1e-5

function lambda(λ)
    @assert 1.0 <= λ <= 2.0 "overrelaxation λ must be a value between 1.0 and 2.0"
    
    # iteration counting
    iterations = 0

    # creating 6x5 zero matrix as initial guess
    T = zeros(6, 5)
    for i in 1:6
        for j in 1:5
            T[i, j] = 26.25
        end
    end

    # iterating gauss_seidel! function for matrix T
    while true

        xₙ = copy(T)
        gauss_seidel!(T, k, h, q₀, Tₑₙᵥ)
        iterations += 1
        norm(T - xₙ, Inf) > tol || break
        xₙ₊₁ = relax(T, xₙ, λ)
        T = xₙ₊₁

    end

    return iterations

end