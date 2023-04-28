# +
using CSV
using DataFrames
using ReactionMechanismSimulator.PyPlot
using ReactionMechanismSimulator.YAML
using ReactionMechanismSimulator.PyCall
using ReactionMechanismSimulator.SparseArrays

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
font0 = Dict(
    "font.size" => 12,
)
merge!(rcParams, font0)


println("Loading RMS model...")
mech = YAML.load_file("/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/liquid_mechanism/chem.rms")


# +

majorspeciesnames = ["N-BUTANE", "2-BUTENE", "1,3-BUTADIENE", "CYCLOPENTADIENE", "BENZENE", "1,3-CYCLOHEXADIENE", "TOLUENE", "STYRENE"]
radicalspcnames = [spc["name"] for spc in mech["Phases"][1]["Species"] if spc["radicalelectrons"] != 0.0]
minorspcnames = [spc["name"] for spc in mech["Phases"][1]["Species"] if !(spc["name"] in majorspeciesnames)]
spcnames = [spc["name"] for spc in mech["Phases"][1]["Species"]]

rxnstrings = []
for rxn in mech["Reactions"]
    rxnstring = ""
    for spc in rxn["reactants"]
        rxnstring *= spc
        rxnstring *= "+"
    end
    rxnstring = rxnstring[1:end-1]
    rxnstring *= "<=>"
    for spc in rxn["products"]
        rxnstring *= spc
        rxnstring *= "+"
    end
    rxnstring = rxnstring[1:end-1]
    push!(rxnstrings, rxnstring)
end

source_string_index_dict = Dict(rxnstring => i for (i, rxnstring) in enumerate(rxnstrings))
source_string_index_dict["inlet"] = length(rxnstrings) + 1
source_string_index_dict["outlet"] = length(rxnstrings) + 2
source_string_index_dict["cond"] = length(rxnstrings) + 3
source_string_index_dict["evap"] = length(rxnstrings) + 4

t0 = 0.0
delta_t = 64.0
tf = 3600.0
ts = t0:delta_t:tf+delta_t

d = 2.5
h = 0.3
A = (d / 2)^2 * pi
Vliq = A * h
trays = 1:40

cmap = get_cmap(:plasma)
cs = cmap(trays / length(trays))

perturbed_target_list = ["1,3-BUTADIENE", "CYCLOPENTADIENE"]
perturbed_factor_list = [0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 1.9]
# -

# # Load simulation

liquid_mols = Dict()
for perturbed_target in perturbed_target_list
    for perturbed_factor in perturbed_factor_list
        path = "../simulation_results/1,3-BUTADIENE_$(perturbed_factor)_3600.0_$(delta_t)/simulation_vapor_liquid_yliqn_$(ts[end]).csv"
        liquid_mols[perturbed_target, perturbed_factor] = DataFrame(CSV.File(path))
    end
end

# +
liquid_rops = Dict()
perturbed_target = "1,3-BUTADIENE"
perturbed_factor = 1.0
spc_name_index_dict = Dict(spcname => i for (i, spcname) in enumerate(spcnames))
for tray in trays
    path = "../simulation_results/1,3-BUTADIENE_$(perturbed_factor)_3600.0_$(delta_t)/simulation_vapor_liquid_liqrop_$(tray).csv"
    df = DataFrame(CSV.File(path))
    spc_inds = [spc_name_index_dict[spcname] for spcname in df[!, "rop_spcname"]]
    source_inds = [source_string_index_dict[sourcestring] for sourcestring in df[!, "rop_sourcestring"]]
    liquid_rops[tray] = sparse(spc_inds, source_inds, df[!, "rop"])
end

liquid_rates = Dict()
perturbed_target = "1,3-BUTADIENE"
perturbed_factor = 1.0
path = "../simulation_results/1,3-BUTADIENE_$(perturbed_factor)_3600.0_$(delta_t)/alpha_rates.yml"
liquid_rates = YAML.load_file(path)
alpha1, alpha2, alphas, alphas_DA = liquid_rates[1]
rates = liquid_rates[2]
radical_production_rates = liquid_rates[3]
# -

rates

# # Radical concentration

ENV["COLUMNS"] = "500"

# +
perturbed_target = "1,3-BUTADIENE"
perturbed_factor = 1.0

radicalmoltray = sum(eachcol(liquid_mols[perturbed_target, perturbed_factor][!, radicalspcnames]))
radicaloutlettray = rates["R._outlet"] + rates["R._evap"]
radicalreactiontray = rates["R._Add"] + rates["R._Habs"] + rates["R._Recomb"] + rates["R._Disprop"]
radicalconsumptiontray = radicaloutlettray + radicalreactiontray
radicalproductiontray = sum(values(radical_production_rates))

fig = figure(figsize=(9, 9))
gs = fig.add_gridspec(3, 2)

element(i, j) = get(gs, (i, j))
slice(i, j) = pycall(pybuiltin("slice"), PyObject, i, j)

ax = fig.add_subplot(element(0, 0))
ax.plot(trays, radicalmoltray / Vliq, "-x",)
ax.set_ylabel("R.(liq) (mol/m^3)", fontsize=12)
ax.set_xlabel("Tray", fontsize=12)
ax.set_yscale(:log)
ax.set_title("(a)", loc="left")

ax = fig.add_subplot(element(0, 1))
radical_chemical_lifetime = radicalmoltray ./ radicalreactiontray
radical_residence_time = radicalmoltray ./ radicaloutlettray
ax.plot(trays, radical_chemical_lifetime ./ radical_residence_time, "-o",)
ax.set_yscale(:log)
ax.set_ylabel("R.(liq) " * L"\tau_\mathrm{chem}" * "/" * L"\tau_\mathrm{residence}")
ax.set_xlabel("Tray")
ax.set_title("(b)", loc="left")

ax = fig.add_subplot(element(1, slice(0, 2)))
global bottom
bottom = zeros(length(trays))
for source in ["R._outlet", "R._Add", "R._Habs", "R._Recomb", "R._Disprop", "R._evap",]
    global bottom
    contributions = rates[source] ./ radicalconsumptiontray * 100
    ax.bar(trays, bottom=bottom, contributions, label=source)
    bottom += contributions
end
ax.set_ylabel("R.(liq) consumption (%)")
ax.set_ylim([0, 100])
# ax.set_yscale(:log)
# ax.set_ylim([1e-16, 1e-2])
ax.set_title("(c)", loc="left")
ax.set_xticks([])
ax.set_xticklabels([])
ax.legend(bbox_to_anchor=(1, 1))

ax = fig.add_subplot(element(2, slice(0, 2)))
bottom = zeros(length(trays))
for source in ["R._inlet", "R._cond", "R._RevDisprop"]
    global bottom
    contributions = radical_production_rates[source] ./ radicalproductiontray * 100
    ax.bar(trays, bottom=bottom, contributions, label=source)
    bottom += contributions
end
ax.set_ylabel("R.(liq) source (%)")
ax.set_xlabel("Tray")
ax.set_ylim([0, 100])
# ax.set_yscale(:log)
# ax.set_ylim([1e-16, 1e-2])
ax.set_title("(d)", loc="left")
ax.legend(bbox_to_anchor=(1, 1))

fig.tight_layout()
fig.savefig("basecase_debutanizer_model_liquid_radical.pdf", bbox_inches="tight")

# +
fig = figure(figsize=(10, 5))
gs = fig.add_gridspec(2, 3)

ax = fig.add_subplot(element(0, 0))
ax.scatter(trays, alphas_DA, c=cs, zorder=2)
ax.plot(trays, alphas_DA, zorder=1, "grey")
ax.set_ylabel(L"\alpha_\mathrm{DA}", fontsize=13)
ax.set_xlabel("Tray")
ax.set_title("(a)", loc="left")
ax.set_yscale(:log)

ax = fig.add_subplot(element(1, 0))
ax.scatter(trays, alphas, c=cs, zorder=2)
ax.plot(trays, alphas, zorder=1, "grey")
ax.set_ylim([0, 1])
ax.set_ylabel(L"\alpha_\mathrm{R.}", fontsize=13)
ax.set_xlabel("Tray")
ax.set_title("(b)", loc="left")

ax = fig.add_subplot(element(slice(0, 2), slice(1, 3)), projection="3d")
ns = 1:15
function Wn(n, alpha)
    return n * (1 - alpha)^2 * alpha^(n - 1)
end

for n in ns
    cs = cmap(trays ./ length(trays))
    ax.bar(trays, Wn.(n, alphas), zs=n, zdir="x", color=cs)
end
ax.set_xlabel(L"k")
ax.set_ylabel("Tray", labelpad=10)
ax.set_zlabel("W(k)", labelpad=10)
ax.set_title("(c)", loc="left")

fig.tight_layout()
fig.savefig("basecase_debutanizer_model_ASF_distribution.pdf", bbox_inches="tight")
# -




