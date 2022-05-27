import time
import tomli
import numpy as np

MGRAV = -39.48009132

def grav_position(phsp):
    return np.array([phsp[2], phsp[3], 0.0, 0.0])

def grav_momenta(phsp):
    grmt = MGRAV/np.sqrt(phsp[0]**2 + phsp[1]**2)**3
    return np.array([0.0, 0.0, phsp[0]*grmt, phsp[1]*grmt])

def rk4(traj, i, dt):
    hdt = 0.5*dt
    sdt = 0.16666666666666666*dt

    m1, k1 = grav_momenta(traj[i,:]),          grav_position(traj[i,:])
    m2, k2 = grav_momenta(traj[i,:] + hdt*k1), grav_position(traj[i,:] + hdt*m1)
    m3, k3 = grav_momenta(traj[i,:] + hdt*k2), grav_position(traj[i,:] + hdt*m2)
    m4, k4 = grav_momenta(traj[i,:] +  dt*k3), grav_position(traj[i,:] +  dt*m3)

    traj[i+1,:] = traj[i,:] + sdt*(k1 + 2*k2 + 2*k3 + k4 + m1 + 2*m2 + 2*m3 + m4)


def leapfrog(traj, i, dt):
    hdt = 0.5*dt

    half        = traj[i,:] + hdt*grav_momenta(traj[i,:])
    half        = half      +  dt*grav_position(half)
    traj[i+1,:] = half      + hdt*grav_momenta(half)

def integrate(phsp0, integrator, ns, dt):
    traj      = np.zeros((ns, 4))
    traj[0,:] = phsp0
    for i in range(ns-1):
        integrator(traj, i, dt)
    return traj

if __name__ == '__main__':
    with open("input.toml", "rb") as f:
        params = tomli.load(f)

    phsp0 = np.array([params["X0"], params["Y0"], params["VX0"], params["VY0"]])
    ns    = params["NSTEPS"]
    dt    = (params["TIME1"] - params["TIME0"])/(ns - 1)
    NRUNS = params["NRUNS"]

    t_rk = np.empty(NRUNS)
    t_lp = np.empty(NRUNS)

    # Compile
    integrate(phsp0, rk4, ns, dt)
    integrate(phsp0, leapfrog, ns, dt)

    for i in range(NRUNS):
        t_s  = time.time_ns()
        traj_rk = integrate(phsp0, rk4, ns, dt)
        t_e  = time.time_ns()
        t_rk[i] = (t_e - t_s)*1.0e-9

        t_s  = time.time_ns()
        traj_lp = integrate(phsp0, leapfrog, ns, dt)
        t_e  = time.time_ns()
        t_lp[i] = (t_e - t_s)*1.0e-9
    
    print("%ld\t%.8e\t%.8e\t%.8e\t%.8e" % (ns, np.mean(t_rk), np.mean(t_lp), np.sqrt(np.var(t_rk, ddof=1)/NRUNS), np.sqrt(np.var(t_lp, ddof=1)/NRUNS)))