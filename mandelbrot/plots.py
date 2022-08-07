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
ax0.errorbar(cc.NPTS_DIM, cc.t_real, yerr=cc.err_t_real, color="k", marker="o", markersize=10, capsize=5, label=r"gcc")
ax0.errorbar(co.NPTS_DIM, co.t_real, yerr=co.err_t_real, color="b", marker="o", markersize=10, capsize=5, label=r"gcc -O3")
ax0.errorbar(jl.NPTS_DIM, jl.t_real, yerr=jl.err_t_real, color="r", marker="o", markersize=10, capsize=5, label=r"julia")
ax0.errorbar(py.NPTS_DIM, py.t_real, yerr=py.err_t_real, color="y", marker="o", markersize=10, capsize=5, label=r"python + numpy + numba")
ax0.errorbar(pp.NPTS_DIM, pp.t_real, yerr=pp.err_t_real, color="c", marker="o", markersize=10, capsize=5, label=r"python + numpy")
ax0.set(xlim=[10, 5000], xlabel=r"Number of points per dimension", ylabel=r"Time (s)", title=r"Real Numbers", xscale="log", yscale="log")
ax0.legend()

fig1, ax1 = plt.subplots()
ax1.errorbar(cc.NPTS_DIM, cc.t_complex, yerr=cc.err_t_complex, color="k", marker="o", markersize=10, capsize=5, label=r"gcc")
ax1.errorbar(co.NPTS_DIM, co.t_complex, yerr=co.err_t_complex, color="b", marker="o", markersize=10, capsize=5, label=r"gcc -O3")
ax1.errorbar(jl.NPTS_DIM, jl.t_complex, yerr=jl.err_t_complex, color="r", marker="o", markersize=10, capsize=5, label=r"julia")
ax1.errorbar(py.NPTS_DIM, py.t_complex, yerr=py.err_t_complex, color="y", marker="o", markersize=10, capsize=5, label=r"python + numpy + numba")
ax1.errorbar(pp.NPTS_DIM, pp.t_complex, yerr=pp.err_t_complex, color="c", marker="o", markersize=10, capsize=5, label=r"python + numpy")
ax1.set(xlim=[10, 5000], xlabel=r"Number of points per dimension", ylabel=r"Time (s)", title=r"Complex Numbers", xscale="log", yscale="log")
ax1.legend()

fig0.tight_layout()
fig1.tight_layout()
# fig0.savefig("real.pgf", backend="pgf")
# fig1.savefig("complex.pgf", backend="pgf")
fig0.savefig("mandelbrot_real.png")
fig1.savefig("mandelbrot_complex.png")
