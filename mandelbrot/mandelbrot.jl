using Printf
using Statistics
using TOML

function isconvergent_complex(c_real, c_imag, MAXITER)
    c = c_real + c_imag*1im
    z = 0.0    + 0.0*1im

    for i in 0:(MAXITER-1)
        z = z^2 + c
        if abs2(z) > 4; return i; end
    end

    return MAXITER-1
end

function isconvergent_real(c_real, c_imag, MAXITER)
    zr = 0.0
    zi = 0.0

    for i in 0:(MAXITER-1)
        zr = zr^2 - zi^2 + c_real
        zi = 2*zr*zi + c_imag
        if zr^2 + zi^2 > 4; return i; end
    end

    return MAXITER-1
end


function compute!(converged, c_real, c_imag, method, NPTS_DIM, MAXITER)
    Threads.@threads for j in 1:NPTS_DIM
        for i in 1:NPTS_DIM
            converged[i,j] = method(c_real[j], c_imag[i], MAXITER)
        end
    end
end

function main()
    params = TOML.parsefile("input.toml")

    NPTS_DIM   = params["NPTS_DIM"]
    MAXITER    = params["MAXITER"]
    c_real_min = params["REAL_MIN"]
	c_imag_min = params["IMAG_MIN"]
	c_real_max = params["REAL_MAX"]
	c_imag_max = params["IMAG_MAX"]
    NRUNS      = params["NRUNS"]

	dc_real = (c_real_max - c_real_min)/(NPTS_DIM - 1.0)
	dc_imag = (c_imag_max - c_imag_min)/(NPTS_DIM - 1.0)

    c_real    = [c_real_min + i*dc_real for i in 0:(NPTS_DIM-1)]
    c_imag    = [c_imag_min + i*dc_imag for i in 0:(NPTS_DIM-1)]

    converged = Array{Int64}(undef, NPTS_DIM, NPTS_DIM)
    t_real    = Array{Float64}(undef, NRUNS)
    t_complex = Array{Float64}(undef, NRUNS)

    # Compile
    compute!(converged, c_real, c_imag, isconvergent_real, NPTS_DIM, MAXITER)
    compute!(converged, c_real, c_imag, isconvergent_complex, NPTS_DIM, MAXITER)

    for i = 1:NRUNS
        t_s = time_ns()
	    compute!(converged, c_real, c_imag, isconvergent_real, NPTS_DIM, MAXITER)
        t_e = time_ns()
        t_real[i] = (t_e - t_s)*1.0e-9

        t_s = time_ns()
	    compute!(converged, c_real, c_imag, isconvergent_complex, NPTS_DIM, MAXITER)
        t_e = time_ns()
        t_complex[i] = (t_e - t_s)*1.0e-9
    end
    
    @printf("%ld\t%.8e\t%.8e\t%.8e\t%.8e\n",NPTS_DIM, mean(t_real), mean(t_complex), sqrt(var(t_real)/NRUNS), sqrt(var(t_complex)/NRUNS))

end

main()