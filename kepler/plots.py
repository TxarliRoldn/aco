import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

# mpl.use("pgf")

mpl.rc("figure", figsize = (11,9))
# mpl.rc("text", usetex = True)
mpl.rc("text", usetex = False)
# mpl.rc("pgf", texsystem = "lualatex")
mpl.rc("font", family = "Computer Modern Roman", weight = "normal")
mpl.rc("axes", labelsize = "24", titlesize = "28")
mpl.rc("xtick", labelsize = "20")
mpl.rc("ytick", labelsize = "20")
mpl.rc("legend", title_fontsize = "24", fontsize = "20")

cc = pd.read_csv("ccfile", sep='\t')
co = pd.read_csv("cofile", sep='\t')
jl = pd.read_csv("jlfile", sep='\t')
py = pd.read_csv("pyfile", sep='\t')
pp = pd.read_csv("ppfile", sep='\t')

fig0, ax0 = plt.subplots()
ax0.errorbar(cc.NSTEPS, cc.t_rk4, yerr=cc.err_t_rk4, color="k", marker="o", markersize=10, capsize=5, label=r"gcc")
ax0.errorbar(co.NSTEPS, co.t_rk4, yerr=co.err_t_rk4, color="b", marker="o", markersize=10, capsize=5, label=r"gcc -O3")
ax0.errorbar(jl.NSTEPS, jl.t_rk4, yerr=jl.err_t_rk4, color="r", marker="o", markersize=10, capsize=5, label=r"julia")
ax0.errorbar(py.NSTEPS, py.t_rk4, yerr=py.err_t_rk4, color="y", marker="o", markersize=10, capsize=5, label=r"python + numpy + numba")
ax0.errorbar(pp.NSTEPS, pp.t_rk4, yerr=pp.err_t_rk4, color="c", marker="o", markersize=10, capsize=5, label=r"python + numpy")
ax0.set(xlim=[1000, 100000], xlabel=r"Number of steps", ylabel=r"Time (s)", title=r"Runge-Kutta 4", xscale="log", yscale="log")
ax0.legend()

fig1, ax1 = plt.subplots()
ax1.errorbar(cc.NSTEPS, cc.t_leapfrog, yerr=cc.err_t_leapfrog, color="k", marker="o", markersize=10, capsize=5, label=r"gcc")
ax1.errorbar(co.NSTEPS, co.t_leapfrog, yerr=co.err_t_leapfrog, color="b", marker="o", markersize=10, capsize=5, label=r"gcc -O3")
ax1.errorbar(jl.NSTEPS, jl.t_leapfrog, yerr=jl.err_t_leapfrog, color="r", marker="o", markersize=10, capsize=5, label=r"julia")
ax1.errorbar(py.NSTEPS, py.t_leapfrog, yerr=py.err_t_leapfrog, color="y", marker="o", markersize=10, capsize=5, label=r"python + numpy + numba")
ax1.errorbar(pp.NSTEPS, pp.t_leapfrog, yerr=pp.err_t_leapfrog, color="c", marker="o", markersize=10, capsize=5, label=r"python + numpy")
ax1.set(xlim=[1000, 100000], xlabel=r"Number of steps", ylabel=r"Time (s)", title=r"Leapfrog", xscale="log", yscale="log")
ax1.legend()

fig0.tight_layout()
fig1.tight_layout()
# fig0.savefig("rk4.pgf", backend="pgf")
# fig1.savefig("leapfrog.pgf", backend="pgf")
fig0.savefig("kepler_rk4.png")
fig1.savefig("kepler_leapfrog.png")
