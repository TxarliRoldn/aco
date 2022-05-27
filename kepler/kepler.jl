using Printf
using StaticArrays
using TOML
using Statistics

struct Phsp <: FieldVector{4, Float64}
    x::Float64
    y::Float64
    vx::Float64
    vy::Float64
end

const mGRAV = -39.48009132

function grav_position(phsp)
    return Phsp(phsp.vx, phsp.vy, 0.0, 0.0)
end

function grav_momenta(phsp)
    grmt = mGRAV/sqrt(phsp.x^2 + phsp.y^2)^3

    return Phsp(0.0, 0.0, phsp.x*grmt, phsp.y*grmt)
end

function rk4(fs, traj, i, dt)
    hdt = 0.5*dt
    sdt = 0.16666666666666666*dt

    m1, k1 = fs[1](traj[i]),          fs[2](traj[i])
    m2, k2 = fs[1](traj[i] + hdt*k1), fs[2](traj[i] + hdt*m1)
    m3, k3 = fs[1](traj[i] + hdt*k2), fs[2](traj[i] + hdt*m2)
    m4, k4 = fs[1](traj[i] +  dt*k3), fs[2](traj[i] +  dt*m3)

    traj[i+1] = traj[i] + sdt*(k1 + 2*k2 + 2*k3 + k4 + m1 + 2*m2 + 2*m3 + m4)
    return
end

function leapfrog(fs, traj, i, dt)
    hdt = 0.5*dt

    half      = traj[i] + hdt*fs[1](traj[i])
    half      = half    + dt*fs[2](half)
    traj[i+1] = half    + hdt*fs[1](half)
    return
end

function integrate(fs, phsp0, integrator, ns, dt)
    traj    = Array{Phsp}(undef, ns)
    traj[1] = phsp0
    for i in 1:(ns-1)
        integrator(fs, traj, i, dt)
    end
    return traj
end

function main()
    params = TOML.parsefile("input.toml")

    phsp0 = Phsp(params["X0"], params["Y0"], params["VX0"], params["VY0"])
    ns    = params["NSTEPS"]
    dt    = (params["TIME1"] - params["TIME0"])/(ns - 1)
    NRUNS = params["NRUNS"]

    t_rk = Array{Float64}(undef, NRUNS)
    t_lp = Array{Float64}(undef, NRUNS)

    # Compile
    integrate((grav_momenta, grav_position), phsp0, rk4, ns, dt)
    integrate((grav_momenta, grav_position), phsp0, leapfrog, ns, dt)
    
    for i = 1:NRUNS
        t_s  = time_ns()
        traj_rk = integrate((grav_momenta, grav_position), phsp0, rk4, ns, dt)
        t_e  = time_ns()
        t_rk[i] = (t_e - t_s)*1.0e-9

        t_s  = time_ns()
        traj_lp = integrate((grav_momenta, grav_position), phsp0, leapfrog, ns, dt)
        t_e  = time_ns()
        t_lp[i] = (t_e - t_s)*1.0e-9
    end
    
    @printf("%ld\t%.8e\t%.8e\t%.8e\t%.8e\n", ns, mean(t_rk), mean(t_lp), sqrt(var(t_rk)/NRUNS), sqrt(var(t_lp)/NRUNS))
    return
end

main()