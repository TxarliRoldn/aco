using Printf
using Statistics
using TOML

mutable struct Particle
    id::Int64
    m::Float64
    x::Float64
    v::Float64
    rho::Float64
    e::Float64
    h::Float64
    dv::Float64
    drho::Float64
    de::Float64
end

dkernel_near(r, invh) = 0.6666666666666666*invh^2*(-3.0*(r*invh) + 2.25*(r*invh)^2)
dkernel_far(r, invh)  = -0.5*invh^2*(2.0 - r*invh)^2

function create_system(ms, xs, vs, rhos, es, NPART)
    length(ms)   == 1 ? ms   = ms   .* ones(NPART) : nothing
    length(xs)   == 1 ? xs   = xs   .* ones(NPART) : nothing
    length(vs)   == 1 ? vs   = vs   .* ones(NPART) : nothing
    length(rhos) == 1 ? rhos = rhos .* ones(NPART) : nothing
    length(es)   == 1 ? es   = es   .* ones(NPART) : nothing

    system = Array{Particle}(undef, NPART)
    Threads.@threads for i = 1:NPART
        system[i] = Particle(i, ms[i], xs[i], vs[i], rhos[i], es[i], 0.0, 0.0, 0.0, 0.0)

    end
    return system
end

function update_particle(system, p, dt)
    system[p].x   += dt * system[p].v
    system[p].v   += dt * system[p].dv
    system[p].rho += dt * system[p].drho
    system[p].e   += dt * system[p].de
    return
end

function get_deriv(system, p, NPART, NSPH)
    signs     = Array{Float64}(undef, NPART)
    distances = Array{Float64}(undef, NPART)
    for i in 1:NPART
        distances[i] = system[p].x - system[i].x
        signs[i]     = sign(distances[i])
        distances[i] = abs(distances[i])
    end
    indexs      = sortperm(distances)
    system[p].h = distances[indexs[NSPH + 1]]
    system[p].dv, system[p].drho, system[p].de = 0.0, 0.0, 0.0

    td_sph = 2.0*system[p].h
    id_sph = 1.0/system[p].h
    temp0  = system[p].e/system[p].rho
    for i in 2:(NSPH+1)
        temp1 = system[indexs[i]].e/system[indexs[i]].rho + temp0
        temp2 = system[indexs[i]].v - system[p].v
        temp3 = signs[indexs[i]]*dkernel_near(distances[indexs[i]], id_sph)

        system[p].dv   -= system[indexs[i]].m*temp1*temp3
        system[p].drho -= system[indexs[i]].m*temp2*temp3
        system[p].de   -= system[indexs[i]].m*temp1*temp2*temp3
    end
    for i in (NSPH+2):NPART
        if distances[indexs[i]] > td_sph; break; end
        temp1 = system[indexs[i]].e/system[indexs[i]].rho + temp0
        temp2 = system[indexs[i]].v - system[p].v
        temp3 = signs[indexs[i]]*dkernel_far(distances[indexs[i]], id_sph)

        system[p].dv   -= system[indexs[i]].m*temp1*temp3
        system[p].drho -= system[indexs[i]].m*temp2*temp3
        system[p].de   -= system[indexs[i]].m*temp1*temp2*temp3
    end
    system[p].dv *= 0.4
    system[p].de *= 0.2
    return
end

function integrate_euler(system, NPART, NSPH, NSTEPS, dt)
    for _ in 1:NSTEPS
        Threads.@threads for i = 1:NPART
            get_deriv(system, i, NPART, NSPH)
        end
        Threads.@threads for i = 1:NPART
            update_particle(system, i, dt)
        end
    end
    return
end

function main()
    params = TOML.parsefile("input.toml")
    NPART  = params["NPART"]
    NSTEPS = params["NSTEPS"]
    NSPH   = params["NSPH"]
    TSTART = params["TSTART"]
    TEND   = params["TEND"]
    NRUNS  = params["NRUNS"]

    dx             = 1.0/(NPART - 1.0)
    dt             = (TEND - TSTART)/(NSTEPS - 1.0)
    xs             = range(0.0, stop=1.0, length=NPART)
    es             = 1e-5 * ones(NPART)
    es[NPART >> 1] = 1.0
    system         = create_system(1.0, xs, 0.0, 1.0/dx, es, NPART)

    integrate_euler(system, NPART, NSPH, NSTEPS, dt) # Compile

    t_eu = Array{Float64}(undef, NRUNS)
    for i in 1:NRUNS
        t_s = time_ns()
        system = create_system(1.0, xs, 0.0, 1.0/dx, es, NPART)
        integrate_euler(system, NPART, NSPH, NSTEPS, dt)
        t_e = time_ns()
        t_eu[i] = (t_e - t_s)*1.0e-9
    end

    @printf("%ld\t%.8e\t%.8e\n",NPART, mean(t_eu), sqrt(var(t_eu)/NRUNS))
end

main()