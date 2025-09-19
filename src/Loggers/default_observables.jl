struct DensityObservable <: Observable end
struct TemperatureObservable <: Observable end
struct KineticEnergyObservable <: Observable end
struct PotentialEnergyObservable <: Observable end
struct TotalEnergyObservable <: Observable end
struct VelocityObservable <: Observable end
struct ForceObservable <: Observable end
struct VolumeObservable <: Observable end
struct PressureObservable <: Observable end

observe!(::DensityObservable, system::System, step::Integer) = density(system)
density(system::System) = Float64(natoms(system)) / volume(system)
out_type(::DensityObservable) = Float64
observed_quantity(::DensityObservable) = :density

observe!(::TemperatureObservable, system::System, step::Integer) = Float64(instantaneous_temperature(system))
out_type(::TemperatureObservable) = Float64
observed_quantity(::TemperatureObservable) = :temperature

observe!(::KineticEnergyObservable, system::System, step::Integer) = Float64(kinetic_energy(system))
out_type(::KineticEnergyObservable) = Float64
observed_quantity(::KineticEnergyObservable) = :kinetic_energy

function observe!(::PotentialEnergyObservable, system::System, step::Integer)
    return Float64(compute_potential_energy(system))
end
out_type(::PotentialEnergyObservable) = Float64
observed_quantity(::PotentialEnergyObservable) = :potential_energy

function observe!(::TotalEnergyObservable, system::System, step::Integer)
    ekin = observe!(KineticEnergyObservable(), system, step)
    epot = observe!(PotentialEnergyObservable(), system, step)
    return ekin + epot
end
out_type(::TotalEnergyObservable) = Float64
observed_quantity(::TotalEnergyObservable) = :total_energy

function _center_of_mass_velocity(system::System)
    total_mass = sum(system.masses)
    total_mass == 0 && throw(ArgumentError("total mass is zero"))
    momentum = zero(system.velocities[1])
    for i in 1:natoms(system)
        momentum += system.masses[i] * system.velocities[i]
    end
    com_velocity = momentum / total_mass
    return Float64.(collect(com_velocity))
end

observe!(::VelocityObservable, system::System, step::Integer) = _center_of_mass_velocity(system)
out_type(::VelocityObservable) = Vector{Float64}
observed_quantity(::VelocityObservable) = :velocity

function _net_force(system::System)
    total_force = zero(system.forces[1])
    for i in 1:natoms(system)
        total_force += system.forces[i]
    end
    return Float64.(collect(total_force))
end

observe!(::ForceObservable, system::System, step::Integer) = _net_force(system)
out_type(::ForceObservable) = Vector{Float64}
observed_quantity(::ForceObservable) = :force

observe!(::VolumeObservable, system::System, step::Integer) = Float64(volume(system))
out_type(::VolumeObservable) = Float64
observed_quantity(::VolumeObservable) = :volume

function observe!(::PressureObservable, system::System{T,IT,DimsVal,BOX,FF}, step::Integer) where {T<:AbstractFloat,IT<:Integer,DimsVal,BOX<:AbstractBox{T},FF}
    vol = volume(system)
    vol == 0 && throw(ArgumentError("volume is zero"))
    ekin = kinetic_energy(system)
    prefactor = 2 / (DimsVal * vol)
    return Float64(prefactor * ekin)
end
out_type(::PressureObservable) = Float64
observed_quantity(::PressureObservable) = :pressure
