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

rms_path = ARGS[1]
model_name = ARGS[2]

println("Loading RMS model...")
mech = YAML.load_file(rms_path)

# +

majorspeciesnames = ["N-BUTANE", "2-BUTENE", "1,3-BUTADIENE", "CYCLOPENTADIENE", "BENZENE", "1,3-CYCLOHEXADIENE", "TOLUENE", "STYRENE"]
constantspecies = majorspeciesnames

if model_name == "trace_oxygen_perturbed_debutanizer_model"
    push!(constantspecies, "OXYGEN")
end

radicalspcnames = [spc["name"] for spc in mech["Phases"][1]["Species"] if spc["radicalelectrons"] == 1]
if model_name == "trace_oxygen_perturbed_debutanizer_model"
    carboncenterradicalspcnames = [spc["name"] for spc in mech["Phases"][1]["Species"] if spc["radicalelectrons"] == 1 && (occursin("[C", spc["name"]) || occursin("[c", spc["name"]))]
    peroxylradicalspcnames = [spc["name"] for spc in mech["Phases"][1]["Species"] if spc["radicalelectrons"] == 1 && (occursin("O[O]", spc["name"]) || occursin("[O]O", spc["name"]))]
    alkoxylradicalspcnames = [spc["name"] for spc in mech["Phases"][1]["Species"] if spc["radicalelectrons"] == 1 && (occursin("C[O]", spc["name"]) || occursin("[O]C", spc["name"]))]
end
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

if model_name == "basecase_debutanizer_model"
    perturbed_target_list = ["1,3-BUTADIENE", "CYCLOPENTADIENE"]
    perturbed_factor_string_list = ["0.5", "0.7", "0.9", "1.0", "1.1", "1.3", "1.5", "1.9"]
elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
    perturbed_target_list = ["OXYGEN"]
    perturbed_factor_string_list = ["0.0", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2", "1e-1", "1e0"]
end
# -

# # Load simulation
if model_name == "basecase_debutanizer_model"
    perturbed_target = "1,3-BUTADIENE"
    perturbed_factor = "1.0"
elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
    perturbed_target = "OXYGEN"
    perturbed_factor = "1e0"
end

liquid_mols = Dict()
for perturbed_target in perturbed_target_list
    for perturbed_factor in perturbed_factor_string_list
        path = "../simulation_results/$(perturbed_target)_$(perturbed_factor)_3600.0_$(delta_t)/simulation_vapor_liquid_yliqn_$(ts[end]).csv"
        liquid_mols[perturbed_target, perturbed_factor] = DataFrame(CSV.File(path))
    end
end

# +
liquid_rops = Dict()
spc_name_index_dict = Dict(spcname => i for (i, spcname) in enumerate(spcnames))
for tray in trays
    path = "../simulation_results/$(perturbed_target)_$(perturbed_factor)_3600.0_$(delta_t)/simulation_vapor_liquid_liqrop_$(tray).csv"
    df = DataFrame(CSV.File(path))
    spc_inds = [spc_name_index_dict[spcname] for spcname in df[!, "rop_spcname"]]
    source_inds = [source_string_index_dict[sourcestring] for sourcestring in df[!, "rop_sourcestring"]]
    liquid_rops[tray] = sparse(spc_inds, source_inds, df[!, "rop"])
end

liquid_rates = Dict()
path = "../simulation_results/$(perturbed_target)_$(perturbed_factor)_3600.0_$(delta_t)/alpha_rates.yml"
liquid_rates = YAML.load_file(path)
alpha1, alpha2, alphas, alphas_DA = liquid_rates[1]
rates = liquid_rates[2]
radical_production_rates = liquid_rates[3]
# -

rates

# # Radical concentration

ENV["COLUMNS"] = "500"

# +

carboncenterradicalmoltray = sum(eachcol(liquid_mols[perturbed_target, perturbed_factor][!, carboncenterradicalspcnames]))
carboncenterradicaloutlettray = rates["R._outlet"] + rates["R._evap"]
carboncenterradicalreactiontray = rates["R._Add"] + rates["R._Habs"] + rates["R._Recomb"] + rates["R._Disprop"]
if model_name == "trace_oxygen_perturbed_debutanizer_model"
    carboncenterradicalreactiontray += rates["R.+O2"] + rates["R._CycEther"]

    peroxylradicalmoltray = zeros(40)
    for (i, mols) in enumerate(eachrow(liquid_mols[perturbed_target, perturbed_factor][!, peroxylradicalspcnames]))
        peroxylradicalmoltray[i] = sum(mols[mols .> 0.0])
    end

    peroxylradicaloutlettray = rates["ROO._outlet"] + rates["ROO._evap"]
    peroxylradicalreactiontray = rates["ROO._Add"] + rates["ROO._Habs"] + rates["ROO._Recomb"] + rates["ROO._Disprop"] + rates["ROO._eli"]

    alkoxylradicalmoltray = zeros(40)
    for (i, mols) in enumerate(eachrow(liquid_mols[perturbed_target, perturbed_factor][!, alkoxylradicalspcnames]))
        alkoxylradicalmoltray[i] = sum(mols[mols .> 0.0])
    end
    alkoxylradicaloutlettray = rates["RO._outlet"] + rates["RO._evap"]
    alkoxylradicalreactiontray = rates["RO._Add"] + rates["RO._Habs"] + rates["RO._Recomb"] + rates["RO._Disprop"] + rates["RO._CycEther"]
end

fig = figure(figsize=(9, 9))
gs = fig.add_gridspec(3, 2)

element(i, j) = get(gs, (i, j))
slice(i, j) = pycall(pybuiltin("slice"), PyObject, i, j)

ax = fig.add_subplot(element(0, 0))
if model_name == "basecase_debutanizer_model"
    ax.plot(trays, carboncenterradicalmoltray / Vliq, "-x",)
elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
    ax.plot(trays, carboncenterradicalmoltray / Vliq, "-o", label="RC.")
    concs = peroxylradicalmoltray / Vliq
    ax.plot(trays, concs, "-x", label="ROO.")
    concs = alkoxylradicalmoltray / Vliq
    ax.plot(trays, concs, "-s", label="RO.")
end
ax.set_ylabel("R.(liq) (mol/m^3)", fontsize=12)
ax.set_xlabel("Tray", fontsize=12)
ax.set_yscale(:log)
ax.set_ylim([1e-18, 1e-5])
ax.set_title("(a)", loc="left")

ax = fig.add_subplot(element(0, 1))
carboncenterradical_chemical_lifetime = carboncenterradicalmoltray ./ carboncenterradicalreactiontray
carboncenterradical_residence_time = carboncenterradicalmoltray ./ carboncenterradicaloutlettray
if model_name == "basecase_debutanizer_model"
    ax.plot(trays, carboncenterradical_chemical_lifetime ./ carboncenterradical_residence_time, "-o",)
elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
    ax.plot(trays, carboncenterradical_chemical_lifetime ./ carboncenterradical_residence_time, "-o", label="RC.")
    peroxylradical_chemical_lifetime = peroxylradicalmoltray ./ peroxylradicalreactiontray
    peroxylradical_residence_time = peroxylradicalmoltray ./ peroxylradicaloutlettray
    ratios = peroxylradical_chemical_lifetime ./ peroxylradical_residence_time
    ax.plot(trays[1:23], ratios[1:23], "-x", label="ROO.")
    alkoxylradical_chemical_lifetime = alkoxylradicalmoltray ./ alkoxylradicalreactiontray
    alkoxylradical_residence_time = alkoxylradicalmoltray ./ alkoxylradicaloutlettray
    ratios = alkoxylradical_chemical_lifetime ./ alkoxylradical_residence_time
    ax.plot(trays[1:23], ratios[1:23], "-s", label="RO.")
    ax.legend(bbox_to_anchor=(1, 1))
end
ax.set_yscale(:log)
ax.set_ylabel("R.(liq) " * L"\tau_\mathrm{chem}" * "/" * L"\tau_\mathrm{residence}")
ax.set_xlabel("Tray")
ax.set_title("(b)", loc="left")

ax = fig.add_subplot(element(1, slice(0, 2)))
global bottom
bottom = zeros(length(trays))
sources = ["R._outlet", "R._evap", "R._Add", "R._Habs", "R._Recomb", "R._Disprop"]
if model_name == "trace_oxygen_perturbed_debutanizer_model"
    append!(sources, ["R.+O2", "R._CycEther",
                        "ROO._outlet", "ROO._evap", "ROO._Add", "ROO._Habs", "ROO._Recomb", "ROO._Disprop", "ROO._eli",
                        "RO._outlet", "RO._evap", "RO._Add", "RO._Habs", "RO._Recomb", "RO._Disprop", "RO._CycEther"])
end
allradicalconsumptiontray = sum(rates[source] for source in sources)
for source in sources
    global bottom
    contributions = rates[source] ./ allradicalconsumptiontray * 100
    if all(contributions .< 2.0)
        continue
    end
    label = source
    if model_name == "trace_oxygen_perturbed_debutanizer_model"
        if occursin("R.", source) && !occursin("ROO.", source)
            label = "RC."*split(source, "R.")[end]
        end
    end
    ax.bar(trays, bottom=bottom, contributions, label=label)
    bottom += contributions
end
ax.bar(trays, bottom=bottom, 100 .- bottom, label="other")
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
sources = ["R._inlet", "R._cond", "R._RevDisprop"]
if model_name == "trace_oxygen_perturbed_debutanizer_model"
   append!(sources, ["ROO._inlet", "ROO._cond", "RO._inlet", "RO._cond", "RO._BondDiss"])
end
allradicalproductiontray = sum(radical_production_rates[source] for source in sources)
for source in sources
    global bottom
    contributions = radical_production_rates[source] ./ allradicalproductiontray * 100
    if all(contributions .< 2.0)
        continue
    end
    label = source
    if model_name == "trace_oxygen_perturbed_debutanizer_model"
        if occursin("R.", source) && !occursin("ROO.", source)
            label = "RC."*split(source, "R.")[end]
        end
    end
    ax.bar(trays, bottom=bottom, contributions, label=label)
    bottom += contributions
end
ax.bar(trays, bottom=bottom, 100 .- bottom, label="other")
ax.set_ylabel("R.(liq) source (%)")
ax.set_xlabel("Tray")
ax.set_ylim([0, 100])
# ax.set_yscale(:log)
# ax.set_ylim([1e-16, 1e-2])
ax.set_title("(d)", loc="left")
ax.legend(bbox_to_anchor=(1, 1))

fig.tight_layout()
fig.savefig("$(model_name)_liquid_radical.pdf", bbox_inches="tight")

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
ns = 1:40
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
fig.savefig("$(model_name)_ASF_distribution.pdf", bbox_inches="tight")
# -




