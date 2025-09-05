# Forces & Potentials

A **force function** takes positions `r::Matrix{Float64}` of shape `(N×D)` and returns a force matrix of the same shape.

```@example forces
# Free particle (no forces)
forces_free(r) = zeros(size(r))

# Linear spring to the origin (Hooke's law, k = 1)
forces_ho(r) = -r

# With potential-energy support:
function forces_ho_with_U(r; return_potential=false)
    F = -r
    return return_potential ? (F, 0.5 * sum(abs2, r)) : F
end
nothing
```

To make [`potential_energy`](@ref) work, support the keyword protocol:

```@example forces
# Called internally as: forces(r; return_potential=true) => (F, U)
nothing
```

If your force does **not** provide `return_potential=true`, `potential_energy` will throw with a clear message.