using CSV
using DataFrames
using Base.Threads
using LinearAlgebra
using ReactionMechanismSimulator
using ReactionMechanismSimulator.SciMLBase
using ReactionMechanismSimulator.Sundials
using ReactionMechanismSimulator.YAML
using ReactionMechanismSimulator.SparseArrays

########## Parse command line arguments ########
# path to liquid mechanism in .rms format
liquid_mechanism_path = ARGS[1]
println("Liquid mechanism path = $(liquid_mechanism_path)")
# perturbation factor for O2
perturb_factor_string = ARGS[2]
println("perturb_factor = $(perturb_factor_string)")
perturb_factor = parse(Float64, perturb_factor_string)
################################################

# calculating initial conditions based on experimental conditions
rho_dict = Dict()
rho_dict["Dowtherm A"] = 1056 #kg/m^3
rho_dict["methylcyclohexadiene"] = 770 #kg/m^3

Mw_dict = Dict()
Mw_dict["Dowtherm A"] = 324.463 / 1000.0 #kg/mol
Mw_dict["methylcyclohexadiene"] = 94.15 / 1000.0 #kg/mol

rho_mixture = 0.7 * rho_dict["Dowtherm A"] + 0.3 * rho_dict["methylcyclohexadiene"] #kg/m^3

initial_conc_dict = Dict()
conc_Dowtherm_A = 0.7 * rho_dict["Dowtherm A"] / Mw_dict["Dowtherm A"] #mol/m^3
initial_conc_dict["c1ccccc1"] = conc_Dowtherm_A

conc_methylcyclohexadiene = 0.3 * rho_dict["methylcyclohexadiene"] / Mw_dict["methylcyclohexadiene"] #mol/m^3
initial_conc_dict["5-methylcyclohexadiene"] = conc_methylcyclohexadiene * 3 / 100
initial_conc_dict["1-methylcyclohexadiene"] = conc_methylcyclohexadiene * 47 / 100
initial_conc_dict["2-methylcyclohexadiene"] = conc_methylcyclohexadiene * 48 / 100
initial_conc_dict["methylenecyclohexene"] = conc_methylcyclohexadiene * 2 / 100

pure_O2_sat_molfrac = 8.1 * 1e-4 #mol fraction measured experimentally from Battino, R., Rettich, T. R., & Tominaga, T. (1983). The Solubility of Oxygen and Ozone in Liquids. Journal of Physical and Chemical Reference Data, Vol. 12, pp. 163–178. https://doi.org/10.1063/1.555680
benzene_density = 876 #kg/m^3
benzene_mw = 78.11 / 1000 #kg/mol
benzene_sat_conc = benzene_density / benzene_mw #mol/m^3
pure_O2_sat_conc = pure_O2_sat_molfrac * benzene_sat_conc #mol/m^3
air_O2_sat_conc = 0.21 * pure_O2_sat_conc #mol/m^3
initial_conc_dict["[O][O]"] = air_O2_sat_conc
initial_conc_dict["[O][O]"] *= perturb_factor

for (key, value) in initial_conc_dict
    println("$(key) (mol/m^3) = $(value)")
end

Vliqdot = 150 #uL/min
Vliqdot = Vliqdot * 1e-6 * 1e-3 / 60 #m^3/s
T = 90 + 273.15
tf = 3600 * 2
Vliq = 40 * 1e-9 #m^3
abstol = 1e-20
reltol = 1e-6

phaseDict = readinput(liquid_mechanism_path)
liqspcs = phaseDict["phase"]["Species"]
liqrxns = phaseDict["phase"]["Reactions"]
solvent = phaseDict["Solvents"][1]

liq = IdealDiluteSolution(liqspcs, liqrxns, solvent; name="phase", diffusionlimited=true)

spcnames = getfield.(liq.species, :name)
liqrxnstrings = getrxnstr.(liq.reactions)

simulation_result_folder = "simulation_results"
if !isdir(simulation_result_folder)
    mkdir(simulation_result_folder)
end

save_directory = joinpath(simulation_result_folder, "O2_$(perturb_factor_string)")
if !isdir(save_directory)
    mkdir(save_directory)
end

initialconds = Dict{String,Float64}(name => conc * Vliq for (name, conc) in initial_conc_dict)
initialconds["T"] = T
initialconds["V"] = Vliq
domain, y0, p = ConstantTVDomain(phase=liq, initialconds=initialconds)


conddict = Dict{String,Float64}(name => conc / sum(values(initial_conc_dict)) for (name, conc) in initial_conc_dict)
conddict["T"] = T
conddict["P"] = 1e8
F = x -> sum(values(initial_conc_dict)) * Vliqdot
Vout = x -> Vliqdot

interfaces = [Inlet(domain, conddict, F),
    VolumetricFlowRateOutlet(domain, Vout),
]
react = Reactor(domain, y0, (0.0, tf), interfaces; p=p) #Create the reactor object
@time sol = solve(react.ode, react.recommendedsolver, abstol=abstol, reltol=reltol) #solve the ode associated with the reactor
df = DataFrame(sol)
rename!(df, names(df)[domain.indexes[1]+1:domain.indexes[2]+1] .=> spcnames)
CSV.write("$save_directory/simulation_liquid_1.csv", df)

bsol = Simulation(sol, domain)
ropmat = rops(bsol, bsol.sol.t[end])
rxnind, spcind, rop = findnz(ropmat)
rop_sourcestrings = liqrxnstrings[rxnind]
rop_spcnames = spcnames[spcind]
inlet = zeros(length(spcnames))
outlet = -interfaces[2].Vout(bsol.sol.t[end]) .* bsol.sol[end] ./ Vliq
evap = zeros(length(spcnames))
cond = zeros(length(spcnames))
rop_sourcestrings = [rop_sourcestrings; ["inlet" for _ in 1:length(spcnames)]; ["outlet" for _ in 1:length(spcnames)]; ["evap" for _ in 1:length(spcnames)]; ["cond" for _ in 1:length(spcnames)]]
rop_spcnames = [rop_spcnames; spcnames; spcnames; spcnames; spcnames]
rop = [rop; inlet; outlet; evap; cond]
df = DataFrame([:rop_sourcestring => rop_sourcestrings, :rop_spcname => rop_spcnames, :rop => rop])
CSV.write("$save_directory/simulation_liquid_liqrop_1.csv", df)
