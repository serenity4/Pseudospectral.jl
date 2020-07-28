function diff_matrix!(M::AbstractArray{T} where T, n)
    leg_coefs = zeros(n + 1)
    leg_coefs[end] = 1
    leg = Fun(Jacobi(0, 0), leg_coefs)
    M[1, 1] = n * (n + 1) / 4
    M[end, end] = -M[1, 1]
    τ, _ = gausslobatto(n + 1)
    for (i, τᵢ) in enumerate(τ), (j, τⱼ) in enumerate(τ)
        if τᵢ != τⱼ
            M[i, j] = leg(τᵢ) / (leg(τⱼ) * (τᵢ - τⱼ))
        end
    end
end