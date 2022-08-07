import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

mpl.use("pgf")

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
ax0.errorbar(cc.NPART, cc.t_eu, yerr=cc.err_t_eu, color="k", marker="o", markersize=10, capsize=5, label=r"gcc")
ax0.errorbar(co.NPART, co.t_eu, yerr=co.err_t_eu, color="b", marker="o", markersize=10, capsize=5, label=r"gcc -O3")
ax0.errorbar(jl.NPART, jl.t_eu, yerr=jl.err_t_eu, color="r", marker="o", markersize=10, capsize=5, label=r"julia")
ax0.errorbar(py.NPART, py.t_eu, yerr=py.err_t_eu, color="y", marker="o", markersize=10, capsize=5, label=r"python + numpy + numba")
ax0.errorbar(pp.NPART, pp.t_eu, yerr=pp.err_t_eu, color="c", marker="o", markersize=10, capsize=5, label=r"python + numpy")
ax0.set(xlim=[50, 2000], xlabel=r"Number of particles", ylabel=r"Time (s)", title=r"SPH Euler Integrator", xscale="log", yscale="log")
ax0.legend()

fig0.tight_layout()
# fig0.savefig("euler.pgf", backend="pgf")
fig0.savefig("sph_euler.png")
