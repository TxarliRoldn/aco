import time
import tomli
import numpy as np

def dkernel_near(r, invh):
    return 0.6666666666666666*invh**2*(-3.0*(r*invh) + 2.25*(r*invh)**2)

def dkernel_far(r, invh):
    return -0.5*invh**2*(2.0 - r*invh)**2

def create_system(ms, xs, vs, rhos, es, NPART):
    ms =   ms   if isinstance(ms, np.ndarray)   else ms * np.ones(NPART)
    xs =   xs   if isinstance(xs, np.ndarray)   else xs * np.ones(NPART)
    vs =   vs   if isinstance(vs, np.ndarray)   else vs * np.ones(NPART)
    rhos = rhos if isinstance(rhos, np.ndarray) else rhos * np.ones(NPART)
    es =   es   if isinstance(es, np.ndarray)   else es * np.ones(NPART)

    system = np.empty((NPART, 9))
    for i in range(NPART):
        system[i, :] = np.array([ms[i], xs[i], vs[i], rhos[i], es[i], 0.0, 0.0, 0.0, 0.0])
    return system

def update_particle(system, p, dt):
    system[p][1] += dt * system[p][2]
    system[p][2] += dt * system[p][6]
    system[p][3] += dt * system[p][7]
    system[p][4] += dt * system[p][8]

def get_deriv(system, p, NPART, NSPH):
    signs     = np.empty(NPART)
    distances = np.empty(NPART)
    for i in range(NPART):
        distances[i] = system[p][1] - system[i][1]
        signs[i]     = np.sign(distances[i])
        distances[i] = np.absolute(distances[i])
    indexs       = np.argsort(distances)
    system[p][5] = distances[indexs[NSPH]]
    system[p][6], system[p][7], system[p][8] = 0.0, 0.0, 0.0

    td_sph = 2.0*system[p][5]
    id_sph = 1.0/system[p][5]
    temp0  = system[p][4]/system[p][3]
    for i in range(1, NSPH+1):
        temp1 = system[indexs[i]][4]/system[indexs[i]][3] + temp0
        temp2 = system[indexs[i]][2] - system[p][2]
        temp3 = signs[indexs[i]]*dkernel_near(distances[indexs[i]], id_sph)

        system[p][6] -= system[indexs[i]][0]*temp1*temp3
        system[p][7] -= system[indexs[i]][0]*temp2*temp3
        system[p][8] -= system[indexs[i]][0]*temp1*temp2*temp3
    for i in range(NSPH+1, NPART):
        if distances[indexs[i]] > td_sph: 
            break
        temp1 = system[indexs[i]][4]/system[indexs[i]][3] + temp0
        temp2 = system[indexs[i]][2] - system[p][2]
        temp3 = signs[indexs[i]]*dkernel_far(distances[indexs[i]], id_sph)

        system[p][6] -= system[indexs[i]][0]*temp1*temp3
        system[p][7] -= system[indexs[i]][0]*temp2*temp3
        system[p][8] -= system[indexs[i]][0]*temp1*temp2*temp3
    system[p][6] *= 0.4
    system[p][8] *= 0.2

def integrate_euler(system, NPART, NSPH, NSTEPS, dt):
    for _ in range(NSTEPS):
        for i in range(NPART):
            get_deriv(system, i, NPART, NSPH)
        for i in range(NPART):
            update_particle(system, i, dt)

if __name__ == '__main__':
    with open("input.toml", "rb") as f:
        params = tomli.load(f)

    NPART  = params["NPART"]
    NSTEPS = params["NSTEPS"]
    NSPH   = params["NSPH"]
    TSTART = params["TSTART"]
    TEND   = params["TEND"]
    NRUNS  = params["NRUNS"]

    dx             = 1.0/(NPART - 1.0)
    dt             = (TEND - TSTART)/(NSTEPS - 1.0)
    xs             = np.linspace(0.0, 1.0, num=NPART)
    es             = 1e-5 * np.ones(NPART)
    es[NPART >> 1] = 1.0
    system         = create_system(1.0, xs, 0.0, 1.0/dx, es, NPART)

    integrate_euler(system, NPART, NSPH, NSTEPS, dt) # Compile

    t_eu = np.empty(NRUNS)
    for i in range(NRUNS):
        t_s = time.time_ns()
        system = create_system(1.0, xs, 0.0, 1.0/dx, es, NPART)
        integrate_euler(system, NPART, NSPH, NSTEPS, dt)
        t_e = time.time_ns()
        t_eu[i] = (t_e - t_s)*1.0e-9

    print("%ld\t%.8e\t%.8e" % (NPART, np.mean(t_eu), np.sqrt(np.var(t_eu, ddof=1)/NRUNS)))