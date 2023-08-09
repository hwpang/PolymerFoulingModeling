using CSV
using DataFrames
using ReactionMechanismSimulator.PyPlot
using ReactionMechanismSimulator.YAML

majorspeciesnames = ["N-BUTANE", "2-BUTENE", "1,3-BUTADIENE", "CYCLOPENTADIENE", "BENZENE", "1,3-CYCLOHEXADIENE", "TOLUENE", "STYRENE"]
mech = YAML.load_file("/home/gridsan/hwpang/Software/PolymerFoulingModeling/debutanizer_models/basecase/liquid_mechanism/chem.rms")
radicalspcnames = [spc["name"] for spc in mech["Phases"][1]["Species"] if spc["radicalelectrons"] != 0.0]
minorspcnames = [spc["name"] for spc in mech["Phases"][1]["Species"] if !(spc["name"] in majorspeciesnames)]

# +
t0 = 0.0
# delta_ts = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0]
delta_ts = [1024.0, 512.0, 256.0, 128.0, 64.0, 32.0]
tf = 3600.0



d = 2.5
h = 0.3
A = (d / 2)^2 * pi
Vliq = A * h

yliqn_minorspc_dict = Dict()
yliqn_radicalspc_dict = Dict()

for delta_t in delta_ts
    ts = t0:delta_t:tf+delta_t
    for t in ts
        path = "../time_step_studies/simulation_results/1,3-BUTADIENE_1.0_3600.0_$(delta_t)/simulation_vapor_liquid_yliqn_$(t).csv"
        if isfile(path)
            df = DataFrame(CSV.File(path))
            yliqn_minorspc_dict[delta_t, t] = sum(eachcol(df[!, minorspcnames]))
            yliqn_radicalspc_dict[delta_t, t] = sum(eachcol(df[!, radicalspcnames]))
        end
    end
end


# +
fig, axs = subplots(figsize=(6, 6), ncols=2, nrows=3)
trays = [1, 10, 20, 30, 40]
for delta_t in delta_ts
    ts = t0:delta_t:tf+delta_t
    for (i, tray) in enumerate(trays)
        ts_plot = []
        minorspcconcs = []
        for t in ts
            if (delta_t, t) in keys(yliqn_minorspc_dict)
                push!(minorspcconcs, yliqn_minorspc_dict[delta_t, t][tray] / Vliq)
                push!(ts_plot, t)
            end
        end
        axs[i].plot(ts_plot, minorspcconcs, label=L"\Delta" * "t=$(delta_t)")
        axs[i].set_title("Tray $(tray)")
        axs[i].set_xlim([t0, tf])
    end
end
axs[end].axis(:off)
handles, labels = axs[1].get_legend_handles_labels()
axs[end].legend(handles, labels)

axs[2].set_ylabel("Minor species (mol/m^3)", fontsize=12)
axs[3].set_xlabel("Time (sec)", fontsize=12)
fig.tight_layout()
fig.savefig("basecase_debutanizer_model_minor_species_steady_state.pdf", bbox_to_anchor="tight")

# +

fig, axs = subplots(figsize=(6, 6), ncols=2, nrows=3)
trays = [1, 10, 20, 30, 40]
for delta_t in delta_ts
    ts = t0:delta_t:tf+delta_t
    for (i, tray) in enumerate(trays)
        ts_plot = []
        radicalspcconcs = []
        for t in ts
            if (delta_t, t) in keys(yliqn_radicalspc_dict)
                push!(radicalspcconcs, yliqn_radicalspc_dict[delta_t, t][tray] / Vliq)
                push!(ts_plot, t)
            end
        end
        axs[i].plot(ts_plot, radicalspcconcs, label=L"\Delta" * "t=$(delta_t)")
        axs[i].set_title("Tray $(tray)")
        axs[i].set_xlim([t0, tf])
        axs[i].set_xlim([t0, tf])
    end
end
axs[end].axis(:off)
handles, labels = axs[1].get_legend_handles_labels()
axs[end].legend(handles, labels)

axs[2].set_ylabel("Radical (mol/m^3)", fontsize=12)
axs[3].set_xlabel("Time (sec)", fontsize=12)
fig.tight_layout()
fig.savefig("basecase_debutanizer_model_radical_steady_state.pdf", bbox_to_anchor="tight")
# -




