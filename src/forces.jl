function lj_forces(positions::AbstractMatrix, box::CubicBox;
                   ϵ::Real=1.0, σ::Real=1.0, rcut::Real=Inf,
                   shift::Bool=false, return_potential::Bool=false)
    N, D = size(positions)
    F = zeros(Float64, N, D)
    U = 0.0

    σ2 = float(σ)^2
    rcut2 = float(rcut)^2

    Uc = 0.0
    if shift && isfinite(rcut)
        s2c = σ2 / rcut2
        s6c = s2c^3
        Uc  = 4*float(ϵ)*(s6c^2 - s6c)
    end

    Δ = zeros(Float64, D)

    @inbounds for i in 1:N-1
        ri = @view positions[i, :]
        for j in i+1:N
            rj = @view positions[j, :]
            Δ .= ri .- rj
            minimum_image!(Δ, box)
            r2 = dot(Δ, Δ)
            if r2 == 0.0
                continue
            end
            if r2 <= rcut2
                invr2 = 1 / r2
                s2 = σ2 * invr2
                s6 = s2^3
                fr_over_r = 24*float(ϵ)*(2*s6^2 - s6) * invr2
                for k in 1:D
                    f = fr_over_r * Δ[k]
                    F[i,k] += f
                    F[j,k] -= f
                end
                if return_potential
                    U += 4*float(ϵ)*(s6^2 - s6) - Uc
                end
            end
        end
    end
    return return_potential ? (F, U) : F
end
