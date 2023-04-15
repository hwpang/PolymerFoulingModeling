# %%
using XLSX
using DataFrames
using ReactionMechanismSimulator.PyPlot
using ReactionMechanismSimulator.YAML
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
font0 = Dict(
        "font.size" => 20,
)
merge!(rcParams, font0)

# %%
ENV["COLUMNS"] = 500
ENV["ROWS"] = 500

aspen_results_path = ARGS[1]
model_name = ARGS[2]

# %%
TPFQ = DataFrame(XLSX.readtable(aspen_results_path,"TPFQ"))

# %%
liquid_molfraction = DataFrame(XLSX.readtable(aspen_results_path,"LiquidCompositionMoleBasis"))

# %%
vapor_molfraction = DataFrame(XLSX.readtable(aspen_results_path,"VaporCompositionMoleBasis"))

# %%
stream = DataFrame(XLSX.readtable(aspen_results_path,"StreamResults"))

# %%
spcnamedict = Dict(
    "N-BUT-01"=>"N-BUTANE",
    "2-BUT-01"=>"2-BUTENE",
    "1:3-B-01"=>"1,3-BUTADIENE",
    "CYCLO-01"=>"CYCLOPENTADIENE",
    "BENZE-01"=>"BENZENE",
    "1:3-C-01"=>"1,3-CYCLOHEXADIENE",
    "TOLUE-01"=>"TOLUENE",
    "STYRE-01"=>"STYRENE")
spcnames = [
    "N-BUTANE",
    "2-BUTENE",
    "1,3-BUTADIENE",
    "CYCLOPENTADIENE",
    "BENZENE",
    "1,3-CYCLOHEXADIENE",
    "TOLUENE",
    "STYRENE"
]

if model_name == "trace_oxygen_perturbed_debutanizer_model"
    spcnamedict["OXYGEN"] = "OXYGEN"
    push!(spcnames,"OXYGEN")
end

# %%
vapor_molfraction

# %%
initial_conditions = Dict()

liquid_molar_density = stream[14,"DC4-L"]*1000 #mol/m^3
liquid_concentration = liquid_molfraction[:,2:end] .* liquid_molar_density
liquid_outlet_molar_flowrate = TPFQ[2:end,"Liquid from (Mole)"]*1000 #mol/s
liquid_outlet_volumetric_flowrate = liquid_outlet_molar_flowrate ./ liquid_molar_density
initial_conditions["liquid_outlet_volumetric_flowrate"] = liquid_outlet_volumetric_flowrate

vapor_molar_density = stream[14,"DC4-V"]*1000 #mol/m^3
vapor_concentration = vapor_molfraction[:,2:end] .* vapor_molar_density
vapor_outlet_molar_flowrate = TPFQ[2:end,"Vapor from (Mole)"]*1000 #mol/s
vapor_outlet_volumetric_flowrate = vapor_outlet_molar_flowrate ./ vapor_molar_density
initial_conditions["vapor_outlet_volumetric_flowrate"] = vapor_outlet_volumetric_flowrate

T = TPFQ[2:end,"Temperature"]
initial_conditions["T"] = T
P = TPFQ[2:end,"Pressure"]
initial_conditions["P"] = P

initial_conditions["liquid_concentration"] = Dict()
initial_conditions["vapor_concentration"] = Dict()
for (aspenname, spcname) in spcnamedict
    initial_conditions["liquid_concentration"][spcname] = liquid_concentration[:,aspenname]
    initial_conditions["vapor_concentration"][spcname] = vapor_concentration[:,aspenname]
end

if model_name == "basecase_debutanizer_model"
    save_name = "aspen_conditions"
elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
    save_name = "aspen_conditions_oxygen"
end
YAML.write_file("$save_name.yml", initial_conditions)

# %%
d = 2.5
h = 0.3
A = (d/2)^2*pi
Vliq = A*h
spacing = 0.6
Vvap = A*(spacing-h)

lines = ["-v", "-^", "-<", "->", "-1", "-2", "-3", "-4", "-8", "-s", "-p", "-P", "-*", "-h", "-H", "-+", "-x", "-X", "-D", "-d",]

if model_name == "basecase_debutanizer_model"
    fig,axs = subplots(nrows=4,ncols=2,figsize=(10,14))
elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
    fig,axs = subplots(nrows=5,ncols=2,figsize=(10,17))
end

fig.subplots_adjust(wspace=0, hspace=0)

axs[1,1].plot(1:40,initial_conditions["T"],color="k","-o")
axs[1,1].set_ylabel("Temperature (K)",size=20)
axs[1,1].set_title("(a)",loc="left",size=20)
axs[1,1].set_xticklabels([])

for (i,spc) in enumerate(spcnames)
    label = spc
    axs[2,1].plot(1:40,initial_conditions["liquid_concentration"][spc],lines[i],label=label)
end

axs[2,1].set_ylabel("Conc. (mol/m^3)",size=20)
axs[2,1].set_title("(b) Liquid phase",loc="left",size=20)
axs[2,1].set_xticklabels([])
axs[2,1].plot([0,0],[-400,800],"k--")
axs[2,1].plot([0,41],[-400,-400],"k--")
axs[2,1].plot([41,41],[-400,800],"k--")
axs[2,1].plot([0,41],[800,800],"k--")

for (i,spc) in enumerate(spcnames)
    label = spc
    axs[3,1].plot(1:40,initial_conditions["liquid_concentration"][spc],lines[i],label=label)
end
axs[3,1].set_ylabel("Conc. (mol/m^3)",size=20)
axs[3,1].set_ylim([-20,600])
axs[3,1].set_title("Zoomed in on (b)",loc="left",size=20)
axs[3,1].set_xticklabels([])


for (i,spc) in enumerate(spcnames)
    label = spc
    axs[2,2].plot(1:40,initial_conditions["vapor_concentration"][spc],lines[i],label=label)
end
axs[2,2].set_title("(c) Vapor phase",loc="left",size=20)
axs[2,2].set_xticklabels([])
axs[2,2].plot([0,0],[-10,20],"k--")
axs[2,2].plot([0,41],[-10,-10],"k--")
axs[2,2].plot([41,41],[-10,20],"k--")
axs[2,2].plot([0,41],[20,20],"k--")
handles, labels = axs[2,2].get_legend_handles_labels()

axs[1,2].legend(handles=handles, loc="center",fontsize=14)
axs[1,2].axis("off")

for (i,spc) in enumerate(spcnames)
    label = spc
    axs[3,2].plot(1:40,initial_conditions["vapor_concentration"][spc],lines[i],label=label)
end
axs[3,2].set_ylim([-1,20])
axs[3,2].set_title("Zoomed in on (c)",loc="left",size=20)
axs[3,2].set_xticklabels([])

flowrate = initial_conditions["liquid_outlet_volumetric_flowrate"][:]
flowrate[end] = flowrate[end-1]
axs[4,1].plot(1:40,Vliq ./ flowrate,color="k","-o")
axs[4,1].set_ylabel("Residence time (s)",size=20)
axs[4,1].set_title("(d) Liquid phase",loc="left",size=20)
if model_name == "basecase_debutanizer_model"
    axs[4,1].set_xlabel("Tray",size=20)
elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
    axs[4,1].set_xticklabels([])
end

flowrate = initial_conditions["vapor_outlet_volumetric_flowrate"][:]
flowrate[1] = flowrate[2]
axs[4,2].plot(1:40,Vvap ./ flowrate,color="k","-o")
axs[4,2].set_title("(e) Vapor phase",loc="left",size=20)
if model_name == "basecase_debutanizer_model"
    axs[4,2].set_xlabel("Tray",size=20)
elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
    axs[4,2].set_xticklabels([])
end

if model_name == "trace_oxygen_perturbed_debutanizer_model"
    spc = "OXYGEN"
    i = findfirst(x->x==spc,spcnames)
    label = spc
    axs[5,1].plot(1:40,initial_conditions["liquid_concentration"][spc],lines[i],label=spc,color="C8")
    axs[5,1].set_ylabel("Conc. (mol/m^3)",size=20)
    axs[5,1].set_title("(f) Liquid phase",loc="left",size=20)
    axs[5,1].set_xlabel("Tray",size=20)
    axs[5,1].set_yscale("log")

    axs[5,2].plot(1:40,initial_conditions["vapor_concentration"][spc],lines[i],label=spc,color="C8")
    axs[5,2].set_title("(g) Vapor phase",loc="left",size=20)
    axs[5,2].set_xlabel("Tray",size=20)
    axs[5,2].set_yscale("log")
end

fig.tight_layout()

plt.tight_layout()
if model_name == "basecase_debutanizer_model"
    save_name = "monomer_conc_trays"
elseif model_name == "trace_oxygen_perturbed_debutanizer_model"
    save_name = "monomer_conc_trays_oxygen"
end 
savefig("$save_name.pdf",bbox_inches="tight")

# %%


# %%



