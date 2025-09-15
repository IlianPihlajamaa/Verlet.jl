using Verlet
using StaticArrays
using Test

@testset "Bonded Interactions" begin

    function potential_energy_total(system, interaction)
        if interaction isa Bond
            vec_ij = system.positions[interaction.j] - system.positions[interaction.i]
            r = norm(vec_ij)
            return Verlet.Potentials.potential_energy(interaction.potential, r)
        elseif interaction isa Angle
            r_ji = system.positions[interaction.i] - system.positions[interaction.j]
            r_jk = system.positions[interaction.k] - system.positions[interaction.j]
            cos_theta = dot(r_ji, r_jk) / (norm(r_ji) * norm(r_jk))
            theta = acos(clamp(cos_theta, -1.0, 1.0))
            return Verlet.Potentials.potential_energy(interaction.potential, theta)
        elseif interaction isa Dihedral
            r_ij = system.positions[interaction.j] - system.positions[interaction.i]
            r_jk = system.positions[interaction.k] - system.positions[interaction.j]
            r_kl = system.positions[interaction.l] - system.positions[interaction.k]
            m = cross(r_ij, r_jk)
            n = cross(r_jk, r_kl)
            cos_phi = dot(m, n) / (norm(m) * norm(n))
            sign_val = dot(r_ij, n) <= 0 ? 1.0 : -1.0
            phi = acos(clamp(cos_phi, -1.0, 1.0)) * sign_val
            return Verlet.Potentials.potential_energy(interaction.potential, phi)
        end
        return 0.0
    end

    @testset "HarmonicBond Finite Difference" begin
        p1 = SVector(0.0, 0.0, 0.0)
        p2 = SVector(1.1, 0.2, -0.3)
        positions = [p1, p2]

        bond_potential = HarmonicBond(100.0, 1.0)
        bond = Bond(1, 2, bond_potential)

        system = System(positions, [SVector(0.0,0.0,0.0) for _ in 1:2], [SVector(0.0,0.0,0.0) for _ in 1:2], [1.0, 1.0], CubicBox(10.0), [1, 1], Dict(1 => :A); specific_potentials=(bond,))

        # Calculate analytical forces
        ff = Verlet.Neighbors.ForceField(())
        Verlet.Neighbors.compute_all_forces!(system, ff)
        f_analytical = copy(system.forces)

        # Calculate numerical forces using finite differences
        f_numerical = [SVector(0.0, 0.0, 0.0) for _ in 1:2]
        ε = 1e-6
        for i in 1:2
            for j in 1:3
                pos_plus = deepcopy(positions)
                pos_minus = deepcopy(positions)

                pos_plus[i] = pos_plus[i] + SVector(j==1 ? ε : 0, j==2 ? ε : 0, j==3 ? ε : 0)
                pos_minus[i] = pos_minus[i] - SVector(j==1 ? ε : 0, j==2 ? ε : 0, j==3 ? ε : 0)

                sys_plus = System(pos_plus, system.velocities, system.forces, system.masses, system.box, system.types, system.type_names, specific_potentials=(bond,))
                sys_minus = System(pos_minus, system.velocities, system.forces, system.masses, system.box, system.types, system.type_names, specific_potentials=(bond,))

                U_plus = potential_energy_total(sys_plus, bond)
                U_minus = potential_energy_total(sys_minus, bond)

                f_numerical[i] = f_numerical[i] + SVector(j==1 ? -(U_plus - U_minus) / (2ε) : 0, j==2 ? -(U_plus - U_minus) / (2ε) : 0, j==3 ? -(U_plus - U_minus) / (2ε) : 0)
            end
        end

        @test f_analytical[1] ≈ f_numerical[1] atol=1e-5
        @test f_analytical[2] ≈ f_numerical[2] atol=1e-5
    end

    @testset "HarmonicAngle Finite Difference" begin
        positions = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(1.0, 1.0, 0.0)]
        angle_potential = HarmonicAngle(100.0, π / 2)
        angle = Angle(1, 2, 3, angle_potential)
        system = System(positions, [SVector(0.0,0.0,0.0) for _ in 1:3], [SVector(0.0,0.0,0.0) for _ in 1:3], ones(3), CubicBox(10.0), ones(Int, 3), Dict(1 => :A); specific_potentials=(angle,))

        ff = Verlet.Neighbors.ForceField(())
        Verlet.Neighbors.compute_all_forces!(system, ff)
        f_analytical = copy(system.forces)

        f_numerical = [SVector(0.0, 0.0, 0.0) for _ in 1:3]
        ε = 1e-6
        for i in 1:3
            for j in 1:3
                pos_plus = deepcopy(positions)
                pos_minus = deepcopy(positions)
                pos_plus[i] = pos_plus[i] + SVector(j==1 ? ε : 0, j==2 ? ε : 0, j==3 ? ε : 0)
                pos_minus[i] = pos_minus[i] - SVector(j==1 ? ε : 0, j==2 ? ε : 0, j==3 ? ε : 0)
                sys_plus = System(pos_plus, system.velocities, system.forces, system.masses, system.box, system.types, system.type_names, specific_potentials=(angle,))
                sys_minus = System(pos_minus, system.velocities, system.forces, system.masses, system.box, system.types, system.type_names, specific_potentials=(angle,))
                U_plus = potential_energy_total(sys_plus, angle)
                U_minus = potential_energy_total(sys_minus, angle)
                f_numerical[i] = f_numerical[i] + SVector(j==1 ? -(U_plus - U_minus) / (2ε) : 0, j==2 ? -(U_plus - U_minus) / (2ε) : 0, j==3 ? -(U_plus - U_minus) / (2ε) : 0)
            end
        end

        for i in 1:3
            @test f_analytical[i] ≈ f_numerical[i] atol=1e-5
        end
    end

    @testset "PeriodicDihedral Finite Difference" begin
        positions = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(1.0, 1.0, 0.0), SVector(2.0, 1.0, 0.0)]
        dihedral_potential = PeriodicDihedral(10.0, 1, 0.0)
        dihedral = Dihedral(1, 2, 3, 4, dihedral_potential)
        system = System(positions, [SVector(0.0,0.0,0.0) for _ in 1:4], [SVector(0.0,0.0,0.0) for _ in 1:4], ones(4), CubicBox(10.0), ones(Int, 4), Dict(1 => :A); specific_potentials=(dihedral,))

        ff = Verlet.Neighbors.ForceField(())
        Verlet.Neighbors.compute_all_forces!(system, ff)
        f_analytical = copy(system.forces)

        f_numerical = [SVector(0.0, 0.0, 0.0) for _ in 1:4]
        ε = 1e-6
        for i in 1:4
            for j in 1:3
                pos_plus = deepcopy(positions)
                pos_minus = deepcopy(positions)
                pos_plus[i] = pos_plus[i] + SVector(j==1 ? ε : 0, j==2 ? ε : 0, j==3 ? ε : 0)
                pos_minus[i] = pos_minus[i] - SVector(j==1 ? ε : 0, j==2 ? ε : 0, j==3 ? ε : 0)
                sys_plus = System(pos_plus, system.velocities, system.forces, system.masses, system.box, system.types, system.type_names, specific_potentials=(dihedral,))
                sys_minus = System(pos_minus, system.velocities, system.forces, system.masses, system.box, system.types, system.type_names, specific_potentials=(dihedral,))
                U_plus = potential_energy_total(sys_plus, dihedral)
                U_minus = potential_energy_total(sys_minus, dihedral)
                f_numerical[i] = f_numerical[i] + SVector(j==1 ? -(U_plus - U_minus) / (2ε) : 0, j==2 ? -(U_plus - U_minus) / (2ε) : 0, j==3 ? -(U_plus - U_minus) / (2ε) : 0)
            end
        end

        for i in 1:4
            @test f_analytical[i] ≈ f_numerical[i] atol=1e-5
        end
    end
end
