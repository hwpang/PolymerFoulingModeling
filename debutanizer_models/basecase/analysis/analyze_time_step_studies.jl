using CSV
using DataFrames
using ReactionMechanismSimulator.PyPlot
using ReactionMechanismSimulator.YAML

rms_path = ARGS[1]
model_name = ARGS[2]

majorspeciesnames = ["N-BUTANE", "2-BUTENE", "1,3-BUTADIENE", "CYCLOPENTADIENE", "BENZENE", "1,3-CYCLOHEXADIENE", "TOLUENE", "STYRENE"]
constantspecies = majorspeciesnames
if model_name == "trace_oxygen_perturbed_debutanizer_model"
    push!(constantspecies, "OXYGEN")
end
mech = YAML.load_file(rms_path)
carboncenterradicalspcnames = [spc["name"] for spc in mech["Phases"][1]["Species"] if occursin("[C", spc["smiles"]) || occursin("[c", spc["smiles"])]
if model_name == "trace_oxygen_perturbed_debutanizer_model"
    peroxylradicalspcnames = [spc["name"] for spc in mech["Phases"][1]["Species"] if occursin("O[O]", spc["smiles"]) || occursin("[O]O", spc["smiles"])]
    alkoxylradicalspcnames = [spc["name"] for spc in mech["Phases"][1]["Species"] if occursin("C[O]", spc["smiles"]) || occursin("[O]C", spc["smiles"])]
end
minorspcnames = [spc["name"] for spc in mech["Phases"][1]["Species"] if !(spc["name"] in constantspecies)]

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
yliqn_carboncenterradicalspc_dict = Dict()
if model_name == "trace_oxygen_perturbed_debutanizer_model"
    yliqn_peroxylradicalspc_dict = Dict()
    yliqn_alkoxylradicalspc_dict = Dict()
end

if model_name == "basecase_debutanizer_model"
    simulation_folder = "1,3-BUTADIENE_1.0"
elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
    simulation_folder = "OXYGEN_1.0"
end

for delta_t in delta_ts
    ts = t0:delta_t:tf+delta_t
    for t in ts
        path = "../time_step_studies/simulation_results/$(simulation_folder)_3600.0_$(delta_t)/simulation_vapor_liquid_yliqn_$(t).csv"
        if isfile(path)
            df = DataFrame(CSV.File(path))
            yliqn_minorspc_dict[delta_t, t] = sum(eachcol(df[!, minorspcnames]))
            yliqn_carboncenterradicalspc_dict[delta_t, t] = sum(eachcol(df[!, carboncenterradicalspcnames]))
            if model_name == "trace_oxygen_perturbed_debutanizer_model"
                yliqn_peroxylradicalspc_dict[delta_t, t] = sum(eachcol(df[!, peroxylradicalspcnames]))
                yliqn_alkoxylradicalspc_dict[delta_t, t] = sum(eachcol(df[!, alkoxylradicalspcnames]))
            end
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
fig.savefig("$(model_name)_minor_species_steady_state.pdf", bbox_to_anchor="tight")

# +

fig, axs = subplots(figsize=(6, 6), ncols=2, nrows=3)
trays = [1, 10, 20, 30, 40]
for delta_t in delta_ts
    ts = t0:delta_t:tf+delta_t
    for (i, tray) in enumerate(trays)
        ts_plot = []
        carboncenterradicalspcconcs = []
        for t in ts
            if (delta_t, t) in keys(yliqn_carboncenterradicalspc_dict)
                push!(carboncenterradicalspcconcs, yliqn_carboncenterradicalspc_dict[delta_t, t][tray] / Vliq)
                push!(ts_plot, t)
            end
        end
        axs[i].plot(ts_plot, carboncenterradicalspcconcs, label=L"\Delta" * "t=$(delta_t)")
        axs[i].set_title("Tray $(tray)")
        axs[i].set_xlim([t0, tf])
        axs[i].set_xlim([t0, tf])
    end
end
axs[end].axis(:off)
handles, labels = axs[1].get_legend_handles_labels()
axs[end].legend(handles, labels)

if model_name == "basecase_debutanizer_model"
    axs[2].set_ylabel("Radical (mol/m^3)", fontsize=12)
elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
    axs[2].set_ylabel("RC. (mol/m^3)", fontsize=12)
end
axs[3].set_xlabel("Time (sec)", fontsize=12)
fig.tight_layout()
fig.savefig("$(model_name)_carboncenterradical_steady_state.pdf", bbox_to_anchor="tight")
# -

if model_name == "trace_oxygen_perturbed_debutanizer_model"
    fig, axs = subplots(figsize=(6, 6), ncols=2, nrows=3)
    trays = [1, 10, 20, 30, 40]
    for delta_t in delta_ts
        ts = t0:delta_t:tf+delta_t
        for (i, tray) in enumerate(trays)
            ts_plot = []
            peroxylradicalspcconcs = []
            for t in ts
                if (delta_t, t) in keys(yliqn_peroxylradicalspc_dict)
                    push!(peroxylradicalspcconcs, yliqn_peroxylradicalspc_dict[delta_t, t][tray] / Vliq)
                    push!(ts_plot, t)
                end
            end
            axs[i].plot(ts_plot, peroxylradicalspcconcs, label=L"\Delta" * "t=$(delta_t)")
            axs[i].set_title("Tray $(tray)")
            axs[i].set_xlim([t0, tf])
            axs[i].set_xlim([t0, tf])
        end
    end
    axs[end].axis(:off)
    handles, labels = axs[1].get_legend_handles_labels()
    axs[end].legend(handles, labels)

    axs[2].set_ylabel("ROO. (mol/m^3)", fontsize=12)
    axs[3].set_xlabel("Time (sec)", fontsize=12)
    fig.tight_layout()
    fig.savefig("$(model_name)_peroxylradical_steady_state.pdf", bbox_to_anchor="tight")

    fig, axs = subplots(figsize=(6, 6), ncols=2, nrows=3)
    trays = [1, 10, 20, 30, 40]
    for delta_t in delta_ts
        ts = t0:delta_t:tf+delta_t
        for (i, tray) in enumerate(trays)
            ts_plot = []
            alkoxylradicalspcconcs = []
            for t in ts
                if (delta_t, t) in keys(yliqn_alkoxylradicalspc_dict)
                    push!(alkoxylradicalspcconcs, yliqn_alkoxylradicalspc_dict[delta_t, t][tray] / Vliq)
                    push!(ts_plot, t)
                end
            end
            axs[i].plot(ts_plot, alkoxylradicalspcconcs, label=L"\Delta" * "t=$(delta_t)")
            axs[i].set_title("Tray $(tray)")
            axs[i].set_xlim([t0, tf])
            axs[i].set_xlim([t0, tf])
        end
    end
    axs[end].axis(:off)
    handles, labels = axs[1].get_legend_handles_labels()
    axs[end].legend(handles, labels)

    axs[2].set_ylabel("RO. (mol/m^3)", fontsize=12)
    axs[3].set_xlabel("Time (sec)", fontsize=12)
    fig.tight_layout()
    fig.savefig("$(model_name)_alkoxylradical_steady_state.pdf", bbox_to_anchor="tight") 
end



