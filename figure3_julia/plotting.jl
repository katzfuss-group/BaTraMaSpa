using JLD
using StatsBase
using Plots
using LaTeXStrings
Plots.scalefontsizes()
#scale fonts
Plots.scalefontsizes(1.75)
# the parameters for the unconditional problem
q = 6
n = 2 ^ q
N = n ^ 2
μ = 0.05
ν = 1.0
T = 25

# loading the unconditional chain
uncond_chain = load("./out/unconditional_chain_kmax_10000000_nu_0.04_mu_0.002_T_25_q_6_delta_0.25.jld")["chain"]

# loading the conditional chains
cond_chain_1 = load("./out/conditional_chain_p_1_kmax_10000000_nu_0.04_mu_0.002_T_25_q_6_delta_0.25.jld")["chain"]
cond_chain_2 = load("./out/conditional_chain_p_2_kmax_10000000_nu_0.04_mu_0.002_T_25_q_6_delta_0.25.jld")["chain"]
cond_chain_3 = load("./out/conditional_chain_p_3_kmax_10000000_nu_0.04_mu_0.002_T_25_q_6_delta_0.25.jld")["chain"]
cond_chain_4 = load("./out/conditional_chain_p_4_kmax_10000000_nu_0.04_mu_0.002_T_25_q_6_delta_0.25.jld")["chain"]
cond_chain_5 = load("./out/conditional_chain_p_5_kmax_10000000_nu_0.04_mu_0.002_T_25_q_6_delta_0.25.jld")["chain"]

plot_list = []
histogram_list = []

k = 606
x = -4 : 0.01 : 4
std_gaussian(x) = 1 / sqrt(2π) * exp(-x^2 / 2)

pl = plot(x, std_gaussian.(x))

heat_plot = heatmap(reshape(uncond_chain[end], n, n), aspect_ratio=:equal, size=(500,500), xticks=false, yticks=false, xaxis=false, yaxis=false, colorbar=false)
savefig(heat_plot, "./figs/cahn_hilliard_heatmap.pdf")

# plotting unconditional histogram
mn, st = mean_and_std([f[k] for f in uncond_chain])
uncond_hist = histogram([(f[k] - mn) / st for f in uncond_chain], normalize=true)
uncond_plot = plot(uncond_hist, legend=false)
plot!(uncond_plot, x, std_gaussian.(x), label=L"\mathcal{N}(0, 1)", linewidth=3, size=(500,500), xlims=(-3.0,3.0), ylims=(0, 0.45))
savefig(uncond_plot, "figs/cahn_hilliard_uncond_hist.pdf")

# plotting conditional histogram
mn, st = mean_and_std([f[k] for f in cond_chain_1])
cond_hist_1 = histogram([(f[k] - mn) / st for f in cond_chain_1], normalize=true)
cond_plot_1 = plot(cond_hist_1, legend=false)
plot!(cond_plot_1, x, std_gaussian.(x), label=L"\mathcal{N}(0, 1)", linewidth=3, size=(500,500), xlims=(-3.0,3.0), ylims=(0, 0.45))
savefig(cond_plot_1, "figs/cahn_hilliard_cond_hist_1.pdf")

# plotting conditional histogram
mn, st = mean_and_std([f[k] for f in cond_chain_2])
cond_hist_2 = histogram([(f[k] - mn) / st for f in cond_chain_2], normalize=true)
cond_plot_2 = plot(cond_hist_2, legend=false)
plot!(cond_plot_2, x, std_gaussian.(x), label=L"\mathcal{N}(0, 1)", linewidth=3, size=(500,500), xlims=(-3.0,3.0), ylims=(0, 0.45), linealpha=0.75)
savefig(cond_plot_2, "figs/cahn_hilliard_cond_hist_2.pdf")

# plotting conditional histogram
mn, st = mean_and_std([f[k] for f in cond_chain_3])
cond_hist_3 = histogram([(f[k] - mn) / st for f in cond_chain_3], normalize=true)
cond_plot_3 = plot(cond_hist_3, legend=false)
plot!(cond_plot_3, x, std_gaussian.(x), label=L"\mathcal{N}(0, 1)", linewidth=3, size=(500,500), xlims=(-3.0,3.0), ylims=(0, 0.45))
savefig(cond_plot_3, "figs/cahn_hilliard_cond_hist_3.pdf")

# plotting conditional histogram
mn, st = mean_and_std([f[k] for f in cond_chain_4])
cond_hist_4 = histogram([(f[k] - mn) / st for f in cond_chain_4], normalize=true)
cond_plot_4 = plot(cond_hist_4, legend=false)
plot!(cond_plot_4, x, std_gaussian.(x), label=L"\mathcal{N}(0, 1)", linewidth=3, size=(500,500), xlims=(-3.0,3.0), ylims=(0, 0.45))
savefig(cond_plot_4, "figs/cahn_hilliard_cond_hist_4.pdf")

# plotting conditional histogram
mn, st = mean_and_std([f[k] for f in cond_chain_5])
cond_hist_5 = histogram([(f[k] - mn) / st for f in cond_chain_5], normalize=true)
cond_plot_5 = plot(cond_hist_5, legend=false)
plot!(cond_plot_5, x, std_gaussian.(x), label=L"\mathcal{N}(0, 1)", linewidth=3, size=(500,500), xlims=(-3.0,3.0), ylims=(0, 0.45))
savefig(cond_plot_5, "figs/cahn_hilliard_cond_hist_5.pdf")