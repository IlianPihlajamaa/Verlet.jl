using Test
using Verlet

using LinearAlgebra

@testset "ParticleSystem construction" begin
    ps = ParticleSystem([0.0 0.0], [1.0 0.0], [1.0])
    @test size(ps.positions) == (1,2)
    @test size(ps.velocities) == (1,2)
    @test ps.masses == [1.0]
end

@testset "Velocity Verlet free particle" begin
    forces(r) = zeros(size(r))
    ps = ParticleSystem([0.0 0.0], [1.0 0.0], [1.0])
    velocity_verlet!(ps, forces, 0.1)
    @test isapprox(ps.positions[1,1], 0.1; atol=1e-12)
    @test isapprox(ps.velocities[1,1], 1.0; atol=1e-12)
end

@testset "Velocity Verlet harmonic oscillator" begin
    forces(r) = -r
    ps = ParticleSystem([1.0 0.0], [0.0 0.0], [1.0])
    velocity_verlet!(ps, forces, 0.1)
    @test ps.positions[1,1] < 1.0   # should move left
end

@testset "kinetic energy" begin
    # v = (3,4) => v^2 = 25; KE = 0.5 * m * v^2
    ps = ParticleSystem([0.0 0.0], [3.0 4.0], [2.0])
    @test isapprox(kinetic_energy(ps), 0.5 * 2.0 * 25.0; atol=1e-12)
end

@testset "potential energy conventions" begin
    # A force that supports keyword return_potential
    function f_kw(r; return_potential::Bool=false)
        F = -r
        if return_potential
            U = 0.5 * sum(r.^2)  # harmonic potential: 0.5 * Σ r^2
            return (F, U)
        else
            return F
        end
    end
    ps = ParticleSystem([1.0 0.0], [0.0 0.0], [1.0])
    @test isapprox(potential_energy(ps, f_kw), 0.5; atol=1e-12)

    # A force that does NOT provide potential should error
    f_plain(r) = -r
    @test_throws ErrorException potential_energy(ps, f_plain)
end



@testset "CubicBox + minimum_image!" begin
    box = CubicBox(10.0)
    Δ = [ 6.0, -6.0,  0.1]  # in 3D
    minimum_image!(Δ, box)
    @test Δ[1] ==  6.0 - 10.0    # -> -4.0
    @test Δ[2] == -6.0 + 10.0    # ->  4.0
    @test isapprox(Δ[3], 0.1; atol=1e-12)
end

@testset "LJ two-particle (no PBC, analytic check)" begin
    # Place two particles distance r along x
    r  = 1.5
    σ  = 1.0
    ϵ  = 2.0
    R  = [0.0 0.0 0.0;
          r   0.0 0.0]
    box = CubicBox(100.0) # effectively no wrapping
    F, U = lj_forces(R, box; ϵ=ϵ, σ=σ, rcut=Inf, return_potential=true)

    # The component on particle 1 (at x=0) should be **positive** for r > r_min (attraction).
    # Negate the coefficient so `fmag` is the expected +x component on particle 1.
    s    = σ/r
    s6   = s^6
    fmag = -24*ϵ*(2*s6^2 - s6)/r
    @test isapprox(F[1,1],  fmag; atol=1e-10)
    @test isapprox(F[2,1], -fmag; atol=1e-10)
    @test isapprox(F[1,2], 0.0; atol=1e-12)
    @test isapprox(F[2,3], 0.0; atol=1e-12)

    # Analytic potential
    Uref = 4*ϵ*(s6^2 - s6)
    @test isapprox(U, Uref; atol=1e-10)
end

@testset "LJ minimum-image under PBC" begin
    # Put particles near opposite faces of a small box; distance should wrap
    L = 5.0
    box = CubicBox(L)
    R = [ -2.4   0.0  0.0;   # ~ -L/2 + 0.1
           2.4   0.0  0.0]   # ~  L/2 - 0.1
    F = lj_forces(R, box; ϵ=1.0, σ=1.0, rcut=Inf, return_potential=false)
    # Displacement should be 0.2 along x after wrapping; force pushes apart
    @test F[1,1] > 0.0
    @test F[2,1] < 0.0
    @test isapprox(F[1,1], -F[2,1]; atol=1e-12)
end

@testset "LJ cutoff and energy shift" begin
    box = CubicBox(50.0)
    R   = [0.0 0.0; 1.2 0.0]  # 2D for variety
    rcut = 1.25
    # With shift
    F1, U1 = lj_forces(R, box; rcut=rcut, shift=true, return_potential=true)
    # Move just beyond cutoff: forces zero, potential ~0 with shift=true
    R2  = [0.0 0.0; 1.26 0.0]
    F2, U2 = lj_forces(R2, box; rcut=rcut, shift=true, return_potential=true)
    @test isapprox(norm(F2), 0.0; atol=1e-12)
    @test U2 ≈ 0.0 atol=1e-10
    @test U1 < 0.0  # attractive well inside cutoff
end

@testset "Compatibility with velocity_verlet! API" begin
    # Use lj_forces as the `forces` function with potential return
    box = CubicBox(20.0)
    forces(r; return_potential=false) = return_potential ?
        lj_forces(r, box; return_potential=true) :
        lj_forces(r, box; return_potential=false)

    ps = ParticleSystem([0.9 0.0 0.0; 0.0 0.0 0.0],
                        zeros(2,3),
                        ones(2))
    E0 = kinetic_energy(ps) + potential_energy(ps, forces)
    velocity_verlet!(ps, forces, 1e-3)
    E1 = kinetic_energy(ps) + potential_energy(ps, forces)
    @test isfinite(E0) && isfinite(E1)
end

    # --- Neighbor list tests ---
    include("test_neighborlist.jl")