import time
import tomli
import numpy as np

def isconvergent_complex(c_real, c_imag, MAXITER):
    c = c_real + c_imag*1j
    z = 0.0    + 0.0*1j

    for i in range(MAXITER):
        z = z**2 + c
        if z.real**2 + z.imag**2 > 4:
            break

    return i

def isconvergent_real(c_real, c_imag, MAXITER):
    zr = 0.0
    zi = 0.0

    for i in range(MAXITER):
        zr = zr**2 - zi**2 + c_real
        zi = 2*zr*zi + c_imag
        if zr**2 + zi**2 > 4:
            break

    return i

def compute(converged, c_real, c_imag, method, NPTS_DIM, MAXITER):
    for i in range(NPTS_DIM):
        for j in range(NPTS_DIM):
            converged[i,j] = method(c_real[i], c_imag[j], MAXITER)

if __name__ == '__main__':
    with open("input.toml", "rb") as f:
        params = tomli.load(f)

    NPTS_DIM   = params["NPTS_DIM"]
    MAXITER    = params["MAXITER"]
    c_real_min = params["REAL_MIN"]
    c_imag_min = params["IMAG_MIN"]
    c_real_max = params["REAL_MAX"]
    c_imag_max = params["IMAG_MAX"]
    NRUNS      = params["NRUNS"]

    dc_real = (c_real_max - c_real_min)/(NPTS_DIM - 1.0)
    dc_imag = (c_imag_max - c_imag_min)/(NPTS_DIM - 1.0)

    c_real    = np.array([c_real_min + i*dc_real for i in range(NPTS_DIM)])
    c_imag    = np.array([c_imag_min + i*dc_imag for i in range(NPTS_DIM)])

    converged = np.empty((NPTS_DIM, NPTS_DIM), dtype=np.int64)
    t_real    =  np.empty(NRUNS)
    t_complex =  np.empty(NRUNS)

    # Compile
    compute(converged, c_real, c_imag, isconvergent_real, NPTS_DIM, MAXITER)
    compute(converged, c_real, c_imag, isconvergent_complex, NPTS_DIM, MAXITER)

    for i in range(NRUNS):
        t_s = time.time_ns()
        compute(converged, c_real, c_imag, isconvergent_real, NPTS_DIM, MAXITER)
        t_e = time.time_ns()
        t_real[i] = (t_e - t_s)*1.0e-9

        t_s = time.time_ns()
        compute(converged, c_real, c_imag, isconvergent_complex, NPTS_DIM, MAXITER)
        t_e = time.time_ns()
        t_complex[i] = (t_e - t_s)*1.0e-9
    
    print("%ld\t%.8e\t%.8e\t%.8e\t%.8e" % (NPTS_DIM, np.mean(t_real), np.mean(t_complex), np.sqrt(np.var(t_real, ddof=1)/NRUNS), np.sqrt(np.var(t_complex, ddof=1)/NRUNS)))