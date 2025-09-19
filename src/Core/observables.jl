abstract type Observable end

function observe!(obs::Observable, system::System, step::Integer)
    throw(MethodError(observe!, (obs, system, step)))
end

observed_quantity(::Observable) = Symbol()

out_type(::Observable) = Any
