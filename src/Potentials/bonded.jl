# Bonded potentials
export HarmonicBond, Bond
export HarmonicAngle, Angle
export PeriodicDihedral, Dihedral
import ..Core: compute_forces!, System


# Concrete implementations of bonded potentials
struct HarmonicBond{T} <: AbstractBondPotential
    k::T    # Spring constant
    r0::T   # Equilibrium distance
end

struct HarmonicAngle{T} <: AbstractAnglePotential
    k::T    # Spring constant
    θ0::T   # Equilibrium angle
end

struct PeriodicDihedral{T} <: AbstractDihedralPotential
    k::T    # Barrier height
    n::Int  # Periodicity
    ϕ0::T   # Phase offset
end

# Structs to hold the particle indices for each interaction
# These are what will be stored in the System's specific_potentials list
struct Bond{T<:AbstractBondPotential}
    i::Int
    j::Int
    potential::T
end

struct Angle{T<:AbstractAnglePotential}
    i::Int
    j::Int
    k::Int
    potential::T
end

struct Dihedral{T<: AbstractDihedralPotential}
    i::Int
    j::Int
    k::Int
    l::Int
    potential::T
end

# Potential energy and force implementations for bonded potentials

# Harmonic Bond
function potential_energy(bond::HarmonicBond, r)
    dr = r - bond.r0
    return 0.5 * bond.k * dr^2
end

function compute_forces!(bond::Bond, system)
    vec_ij = system.positions[bond.j] - system.positions[bond.i]
    r = norm(vec_ij)

    dr = r - bond.potential.r0
    f_mag = -bond.potential.k * dr
    f_vec = f_mag * (vec_ij / r)

    system.forces[bond.i] -= f_vec
    system.forces[bond.j] += f_vec
end

# Harmonic Angle
function potential_energy(angle::HarmonicAngle, θ)
    dθ = θ - angle.θ0
    return 0.5 * angle.k * dθ^2
end

function compute_forces!(angle::Angle, system)
    r_ji = system.positions[angle.i] - system.positions[angle.j]
    r_jk = system.positions[angle.k] - system.positions[angle.j]

    r_ji_norm = norm(r_ji)
    r_jk_norm = norm(r_jk)

    cos_theta = dot(r_ji, r_jk) / (r_ji_norm * r_jk_norm)
    # Ensure cos_theta is within the valid domain for acos
    cos_theta = clamp(cos_theta, -1.0, 1.0)
    theta = acos(cos_theta)

    d_theta = theta - angle.potential.θ0
    f_mag = -angle.potential.k * d_theta

    # Forces on i and k using formula from http://docs.openmm.org/latest/userguide/theory.html#harmonic-angle-force
    f_i = (f_mag / (r_ji_norm * sin(theta))) * (cos_theta * (r_ji / r_ji_norm) - (r_jk / r_jk_norm))
    f_k = (f_mag / (r_jk_norm * sin(theta))) * (cos_theta * (r_jk / r_jk_norm) - (r_ji / r_ji_norm))

    system.forces[angle.i] += f_i
    system.forces[angle.k] += f_k
    system.forces[angle.j] -= (f_i + f_k)
end


# Periodic Dihedral
# See http://docs.openmm.org/latest/userguide/theory.html#periodic-torsion-force
function potential_energy(dihedral::PeriodicDihedral, ϕ)
    return dihedral.k * (1 + cos(dihedral.n * ϕ - dihedral.ϕ0))
end

function compute_forces!(dihedral::Dihedral, system::System{T,IT,3}) where {T,IT}
    r_ij = system.positions[dihedral.j] - system.positions[dihedral.i]
    r_jk = system.positions[dihedral.k] - system.positions[dihedral.j]
    r_kl = system.positions[dihedral.l] - system.positions[dihedral.k]

    # Normals to the planes defined by r_ij, r_jk and r_jk, r_kl
    m = cross(r_ij, r_jk)
    n = cross(r_jk, r_kl)
    m_norm = norm(m)
    n_norm = norm(n)

    # Angle between the planes
    cos_phi = dot(m, n) / (m_norm * n_norm)
    cos_phi = clamp(cos_phi, -1.0, 1.0)

    # Determine the sign of the angle
    sign_val = dot(r_ij, n) <= 0 ? 1.0 : -1.0
    phi = acos(cos_phi) * sign_val

    # Force calculation
    d_phi = dihedral.potential.n * phi - dihedral.potential.ϕ0
    f_mag = -dihedral.potential.k * dihedral.potential.n * sin(d_phi)

    f_i = f_mag * (norm(r_jk) / (m_norm^2)) * m
    f_l = -f_mag * (norm(r_jk) / (n_norm^2)) * n

    r_jk_norm_sq = dot(r_jk, r_jk)
    f_j = ((dot(r_ij, r_jk) / r_jk_norm_sq) - 1) * f_i - (dot(r_kl, r_jk) / r_jk_norm_sq) * f_l
    f_k = - (dot(r_ij, r_jk) / r_jk_norm_sq) * f_i + ((dot(r_kl, r_jk) / r_jk_norm_sq) - 1) * f_l

    system.forces[dihedral.i] += f_i
    system.forces[dihedral.j] += f_j
    system.forces[dihedral.k] += f_k
    system.forces[dihedral.l] += f_l
end

function compute_forces!(dihedral::Dihedral, system::System{T,IT,D}) where {T,IT,D}
    throw(DomainError(D, "Periodic dihedral forces require 3 spatial dimensions; got D=$D."))
end
