"""
    System(T<:AbstractFloat, IT<:Integer, Dims)

A flexible and type-stable container for particle-based simulations.

# Type Parameters
- `T<:AbstractFloat`: Floating-point type for particle properties (e.g., `Float32`, `Float64`).
- `IT<:Integer`: Integer type for particle indices and types (e.g., `Int32`, `Int64`).
- `Dims`: The number of spatial dimensions.

# Fields
- `positions::Vector{SVector{Dims, T}}`: Particle positions.
- `velocities::Vector{SVector{Dims, T}}`: Particle velocities.
- `forces::Vector{SVector{Dims, T}}`: Forces acting on particles.
- `masses::Vector{T}`: Particle masses.
- `box::AbstractBox{T}`: Simulation box defining the periodic boundary conditions.
- `types::Vector{IT}`: Particle type identifiers (integers).
- `type_names::Dict{IT, Symbol}`: Mapping from type identifiers to descriptive names (e.g., `1 => :H`).
- `natoms::IT`: Total number of atoms.
- `specific_potentials::Tuple`: Optional bonded interactions.
- `forcefield`: Object providing `compute_forces!`; typically a `ForceField` or `nothing`.
"""
struct System{T<:AbstractFloat,IT<:Integer,Dims, BOX<:AbstractBox{T}, FF}
    positions::Vector{SVector{Dims,T}}
    velocities::Vector{SVector{Dims,T}}
    forces::Vector{SVector{Dims,T}}
    masses::Vector{T}
    box::BOX
    types::Vector{IT}
    type_names::Dict{IT,Symbol}
    natoms::IT
    specific_potentials::Tuple
    forcefield::FF

    function System(
        positions::Vector{SVector{Dims,T}},
        velocities::Vector{SVector{Dims,T}},
        forces::Vector{SVector{Dims,T}},
        masses::Vector{T},
        box::AbstractBox{T},
        types::Vector{IT},
        type_names::Dict{IT,Symbol};
        specific_potentials::Tuple=(),
        forcefield=nothing
    ) where {T<:AbstractFloat,IT<:Integer,Dims}
        natoms = length(positions)
        @assert length(velocities) == natoms "velocities must be same size as positions"
        @assert length(forces) == natoms "forces must be same size as positions"
        @assert length(masses) == natoms "masses must be same size as positions"
        @assert length(types) == natoms "types must be same size as positions"
        new{T,IT,Dims, typeof(box), typeof(forcefield)}(positions, velocities, forces, masses, box, types, type_names, natoms, specific_potentials, forcefield)
    end
end

"""
    natoms(sys::System) -> Integer

Get the number of atoms in the system.
"""
natoms(sys::System) = sys.natoms

"""
    natomtypes(sys::System) -> Int

Get the number of unique atom types in the system.
"""
natomtypes(sys::System) = length(sys.type_names)

function system_volume(sys::System{T,IT,Dims}) where {T,IT,Dims}
    return box_volume(sys.box, Dims)
end

volume(sys::System) = system_volume(sys)

"""
    kinetic_energy(sys::System) -> T

Total kinetic energy: `∑ ½ mᵢ ‖vᵢ‖²`.
"""
function kinetic_energy(sys::System{T,IT,Dims}) where {T,IT,Dims}
    Ekin = zero(T)
    for i in 1:natoms(sys)
        Ekin += T(0.5) * sys.masses[i] * dot(sys.velocities[i], sys.velocities[i])
    end
    return Ekin
end
