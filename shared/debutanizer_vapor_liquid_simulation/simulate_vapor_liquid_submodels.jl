"""
This script simulates the vapor-liquid submodels in the debutanizer fouling model with Speth et al. rebalanced operator splitting.
"""

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
# path to Aspen Plus initial conditions
aspen_condition_path = ARGS[1]
println("Aspen condition path: $(aspen_condition_path)")
# path to liquid mechanism in .rms format
liquid_mechanism_rms_path = ARGS[2]
println("Liquid mechanism path: $(liquid_mechanism_rms_path)")
# species to perturb
perturb_species = ARGS[3]
println("Perturbation species: $(perturb_species)")
# perturbation factor
perturb_factor_string = ARGS[4]
println("Perturbation factor: $(perturb_factor_string)")
# time step
delta_t = parse(Float64, ARGS[5])
println("Time step: $(delta_t)")
# model name
model_name = ARGS[6]
println("Model name: $(model_name)")
################################################

println("Number of threads: $(nthreads())")

d = 2.5
h = 0.3
A = (d / 2)^2 * pi
Vliq = A * h
spacing = 0.6
Vvap = A * (spacing - h)
trays = 1:40
println("Number of trays: $(length(trays))")

t0 = 0.0
tf = 3600.0

println("Loading initial conditions")
aspen_initial_conditions = YAML.load_file(aspen_condition_path)

println("Loading mechanism")
phaseDict = readinput(liquid_mechanism_rms_path)
liqspcs = phaseDict["phase"]["Species"];
liqrxns = phaseDict["phase"]["Reactions"];
solvent = phaseDict["Solvents"][1]
vapspcs = liqspcs;

vap = IdealGas(vapspcs, []; name="vap");
liq = IdealDiluteSolution(liqspcs, liqrxns, solvent; name="liq", diffusionlimited=true);

spcnames = getfield.(liq.species, :name)
liqrxnstrings = getrxnstr.(liq.reactions)

aspenspeciesnames = ["N-BUTANE", "2-BUTENE", "1,3-BUTADIENE", "CYCLOPENTADIENE", "BENZENE", "1,3-CYCLOHEXADIENE", "TOLUENE", "STYRENE"]

if model_name == "trace_oxygen_perturbed_debutanizer_model"
    push!(aspenspeciesnames, "OXYGEN")
end

perturb_factor = parse(Float64, perturb_factor_string)

simulation_result_folder = "simulation_results"
if !isdir(simulation_result_folder)
    mkdir(simulation_result_folder)
end
save_directiry = joinpath(simulation_result_folder, "$(perturb_species)_$(perturb_factor_string)_$(tf)_$(delta_t)")
if !isdir(save_directiry)
    mkdir(save_directiry)

end
println("Saving to $(save_directiry)")

println("Initializing initial conditions")

yvap1 = zeros(length(trays), length(spcnames))
yliq1 = zeros(length(trays), length(spcnames))
yvap2 = zeros(length(trays), length(spcnames))
yliq2 = zeros(length(trays), length(spcnames))
yvap3 = zeros(length(trays), length(spcnames))
yliq3 = zeros(length(trays), length(spcnames))
yvapn = zeros(length(trays), length(spcnames))
yliqn = zeros(length(trays), length(spcnames))

for (spcind, spcname) in enumerate(spcnames)
    if spcname in aspenspeciesnames
        for tray in trays
            yvapn[tray, spcind] = aspen_initial_conditions["vapor_concentration"][spcname][tray] * Vvap
            yliqn[tray, spcind] = aspen_initial_conditions["liquid_concentration"][spcname][tray] * Vliq
            if spcname == perturb_species
                yvapn[tray, spcind] *= perturb_factor
                yliqn[tray, spcind] *= perturb_factor
            end
        end
    end
end

kLAs = zeros(length(trays), length(spcnames))
kHs = zeros(length(trays), length(spcnames))
for tray in trays
    T = aspen_initial_conditions["T"][tray]
    for (spcind, spcname) in enumerate(spcnames)
        kLAs[tray, spcind] = liq.species[spcind].liquidvolumetricmasstransfercoefficient(T=T)
        kHs[tray, spcind] = liq.species[spcind].henrylawconstant(T=T)
    end
end

shiftvectorvapn = zeros(length(trays), length(spcnames))
shiftvectorliqn = zeros(length(trays), length(spcnames))

function dydt_T_spc!(dy, y, shiftvectorvapnspc::AbstractArray, shiftvectorliqnspc::AbstractArray, Vvap::Float64, Vliq::Float64, Ts::Array{Float64,1}, Vvapdotout::Array{Float64,1}, Vliqdotout::Array{Float64,1}, kLAsspc::AbstractArray, kHsspc::AbstractArray, trays::UnitRange{Int64})
    dy .= 0.0
    @views dyvap = dy[trays]
    @views dyliq = dy[trays.+length(trays)]
    @views yvap = y[trays]
    @views yliq = y[trays.+length(trays)]
    for tray in trays
        net_evap = kLAsspc[tray] * Vliq * (yliq[tray] / Vliq - yvap[tray] / Vvap * R * Ts[tray] / kHsspc[tray])
        if tray != length(trays)
            dyvap[tray] += Vvapdotout[tray+1] * yvap[tray+1] / Vvap
        end
        dyvap[tray] += -Vvapdotout[tray] * yvap[tray] / Vvap + net_evap + shiftvectorvapnspc[tray]
        if tray != 1
            dyliq[tray] += Vliqdotout[tray-1] * yliq[tray-1] / Vliq
        end
        dyliq[tray] += -Vliqdotout[tray] * yliq[tray] / Vliq - net_evap + shiftvectorliqnspc[tray]
    end
end

function T_reactor!(yvap0::Array{Float64,2}, yliq0::Array{Float64,2}, yvapf::Array{Float64,2}, yliqf::Array{Float64,2},
    shiftvectorvapn::Array{Float64,2}, shiftvectorliqn::Array{Float64,2}, kLAs::Array{Float64,2}, kHs::Array{Float64,2},
    t0::Float64, tf::Float64, aspen_initial_conditions::Dict{Any,Any}, Vvap::Float64, Vliq::Float64,
    trays::UnitRange{Int64}, spcnames::Array{String,1}, aspenspeciesnames::Array{String,1})
    # solve for vap
    @threads for spcname in spcnames
        spcind = findfirst(isequal(spcname), spcnames)
        if !(spcname in aspenspeciesnames)
            y0 = vcat(yvap0[:, spcind], yliq0[:, spcind])
            @views kLAsspc = kLAs[:, spcind]
            @views kHsspc = kHs[:, spcind]
            @views shiftvectorvapnspc = shiftvectorvapn[:, spcind]
            @views shiftvectorliqnspc = shiftvectorliqn[:, spcind]
            dydt!(dy, y, p, t) = dydt_T_spc!(dy, y, shiftvectorvapnspc, shiftvectorliqnspc, Vvap, Vliq, aspen_initial_conditions["T"],
                aspen_initial_conditions["vapor_outlet_volumetric_flowrate"],
                aspen_initial_conditions["liquid_outlet_volumetric_flowrate"], kLAsspc,
                kHsspc, trays)
            odefcn = ODEFunction(dydt!)
            ode = ODEProblem(odefcn, y0, (t0, tf))
            sol = solve(ode, CVODE_BDF(), abstol=1e-18, reltol=1e-6)
            @views yvapf[:, spcind] .= sol[end][trays]
            @views yliqf[:, spcind] .= sol[end][trays.+length(trays)]
        else
            @views yvapf[:, spcind] .= yvap0[:, spcind]
            @views yliqf[:, spcind] .= yliq0[:, spcind]
        end
    end
end

function solve_R_reactor_liq_tray(yliq0tray::AbstractArray, shiftvectorliqntray::AbstractArray, t0::Float64, tf::Float64, T::Float64, Vliq::Float64, trays::UnitRange{Int64}, spcnames::Array{String,1}, aspenspeciesnames::Array{String,1})

    # liquid
    initialconds = Dict(["T" => T, "V" => Vliq])
    for (spcind, spcname) in enumerate(spcnames)
        initialconds[spcname] = yliq0tray[spcind]
    end

    domain, y0, p = ConstantTVDomain(phase=liq, initialconds=initialconds, constantspecies=aspenspeciesnames)

    react = Reactor(domain, y0, (t0, tf), []; p=p, shiftvector=shiftvectorliqntray)

    sol = solve(react.ode, react.recommendedsolver, abstol=1e-18, reltol=1e-6)

    return sol, domain, p
end

function dydt_R_vap_tray!(dy, y, shiftvectorvapntray::AbstractArray)
    dy .= 0.0
    dy .-= shiftvectorvapntray
end

function R_reactor!(yvap0::Array{Float64,2}, yliq0::Array{Float64,2}, yvapf::Array{Float64,2}, yliqf::Array{Float64,2},
    shiftvectorvapn::Array{Float64,2}, shiftvectorliqn::Array{Float64,2},
    t0::Float64, tf::Float64, aspen_initial_conditions::Dict{Any,Any}, Vvap::Float64, Vliq::Float64,
    trays::UnitRange{Int64}, spcnames::Array{String,1}, aspenspeciesnames::Array{String,1})
    @threads for tray in trays
        @views yliq0tray = yliq0[tray, :]
        @views shiftvectorliqntray = shiftvectorliqn[tray, :]
        T = aspen_initial_conditions["T"][tray]
        sol, domain, p = solve_R_reactor_liq_tray(yliq0tray, shiftvectorliqntray, t0, tf, T, Vliq, trays, spcnames, aspenspeciesnames)
        @views yliqf[tray, :] .= sol[end][domain.indexes[1]:domain.indexes[2]]

        @views shiftvectorvapntray = shiftvectorvapn[tray, :]
        dydt!(dy, y, p, t) = dydt_R_vap_tray!(dy, y, shiftvectorvapntray)
        odefcn = ODEFunction(dydt!)
        ode = ODEProblem(odefcn, yvap0[tray, :], (t0, tf))
        sol = solve(ode, CVODE_BDF(), abstol=1e-18, reltol=1e-6)
        @views yvapf[tray, :] .= sol[end]
    end
end

function save_rop(yliqn::Array{Float64,2}, shiftvectorliqn::Array{Float64,2}, t0::Float64, tf::Float64, aspen_initial_conditions::Dict{Any,Any}, Vvap::Float64, Vliq::Float64, kLAs::Array{Float64,2}, kHs::Array{Float64,2}, trays::UnitRange{Int64}, spcnames::Array{String,1}, aspenspeciesnames::Array{String,1}, save_directiry::String)

    for tray in trays
        @views yliq0tray = yliqn[tray, :]
        @views shiftvectorliqntray = shiftvectorliqn[tray, :]
        T = aspen_initial_conditions["T"][tray]
        sol, domain, p = solve_R_reactor_liq_tray(yliq0tray, shiftvectorliqntray, t0, tf, T, Vliq, trays, spcnames, aspenspeciesnames)
        sim = Simulation(sol, domain, p)
        ropmat = rops(sim, sol.t[1])
        rxnind, spcind, rop = findnz(ropmat)

        if tray != trays[1]
            inlet = aspen_initial_conditions["liquid_outlet_volumetric_flowrate"][tray-1] .* yliqn[tray-1, :] ./ Vliq
        else
            inlet = zeros(length(spcnames))
        end
        @views outlet = -aspen_initial_conditions["liquid_outlet_volumetric_flowrate"][tray] .* yliqn[tray, :] ./ Vliq
        T = aspen_initial_conditions["T"][tray]
        @views evap = -kLAs[tray, :] * Vliq .* yliqn[tray, :] / Vliq
        @views cond = kLAs[tray, :] * Vliq .* yvapn[tray, :] / Vvap * R * T ./ kHs[tray, :]
        sourcestrings = [liqrxnstrings[rxnind]; ["inlet" for _ in 1:length(spcnames)]; ["outlet" for _ in 1:length(spcnames)]; ["evap" for _ in 1:length(spcnames)]; ["cond" for _ in 1:length(spcnames)]]
        sourcespcnames = [spcnames[spcind]; spcnames; spcnames; spcnames; spcnames]
        sources = [rop; inlet; outlet; evap; cond]
        df = DataFrame([:rop_sourcestring => sourcestrings, :rop_spcname => sourcespcnames, :rop => sources])
        CSV.write("$(save_directiry)/simulation_vapor_liquid_liqrop_$(tray).csv", df)
    end
end

df = DataFrame(yvapn, spcnames)
CSV.write("$(save_directiry)/simulation_vapor_liquid_yvapn_0.0.csv", df)
df = DataFrame(yliqn, spcnames)
CSV.write("$(save_directiry)/simulation_vapor_liquid_yliqn_0.0.csv", df)

println("Starting simulation...")
for t in t0:delta_t:tf

    @time T_reactor!(yvapn, yliqn, yvap1, yliq1, shiftvectorvapn, shiftvectorliqn, kLAs, kHs, t, t + delta_t / 2, aspen_initial_conditions, Vvap, Vliq, trays, spcnames, aspenspeciesnames)
    vap_norm_1 = norm(abs.(yvap1 .- yvapn))
    liq_norm_1 = norm(abs.(yliq1 .- yliqn))
    println("vap_norm_1 = $(vap_norm_1)")
    println("liq_norm_1 = $(liq_norm_1)")
    # df = DataFrame(yvap1, spcnames)
    # CSV.write("$(save_directiry)/simulation_vapor_liquid_yvap1_$(t).csv",df)
    # df = DataFrame(yliq1, spcnames)
    # CSV.write("$(save_directiry)/simulation_vapor_liquid_yliq1_$(t).csv",df)

    @time R_reactor!(yvap1, yliq1, yvap2, yliq2, shiftvectorvapn, shiftvectorliqn, t, t + delta_t, aspen_initial_conditions, Vvap, Vliq, trays, spcnames, aspenspeciesnames)
    vap_norm_2 = norm(abs.(yvap2 .- yvap1))
    liq_norm_2 = norm(abs.(yliq2 .- yliq1))
    println("vap_norm_2 = $(vap_norm_2)")
    println("liq_norm_2 = $(liq_norm_2)")
    majornorm = norm((yliq2.-yliq1)[1:end, 5])
    println("major species norm = $(majornorm)")
    # df = DataFrame(yvap2, spcnames)
    # CSV.write("$(save_directiry)/simulation_vapor_liquid_yvap2_$(t).csv",df)
    # df = DataFrame(yliq2, spcnames)
    # CSV.write("$(save_directiry)/simulation_vapor_liquid_yliq2_$(t).csv",df)

    @time T_reactor!(yvap2, yliq2, yvap3, yliq3, shiftvectorvapn, shiftvectorliqn, kLAs, kHs, t + delta_t / 2, t + delta_t, aspen_initial_conditions, Vvap, Vliq, trays, spcnames, aspenspeciesnames)
    vap_norm_3 = norm(abs.(yvap3 .- yvap2))
    liq_norm_3 = norm(abs.(yliq3 .- yliq2))
    println("vap_norm_3 = $(vap_norm_3)")
    println("liq_norm_3 = $(liq_norm_3)")
    # df = DataFrame(yvap3, spcnames)
    # CSV.write("$(save_directiry)/simulation_vapor_liquid_yvap3_$(t+delta_t).csv",df)
    # df = DataFrame(yliq3, spcnames)
    # CSV.write("$(save_directiry)/simulation_vapor_liquid_yliq3_$(t+delta_t).csv",df)

    # check steady-state
    # calculate norm of the difference between the last two time steps
    vap_norm = norm(abs.(yvap3 .- yvapn))
    liq_norm = norm(abs.(yliq3 .- yliqn))

    # udpate
    shiftvectorliqn .+= 1 / (2 * delta_t) * (-yliq3 .+ 2 * yliq2 .- 2 * yliq1 .+ yliqn)
    shiftvectorvapn .+= 1 / (2 * delta_t) * (-yvap3 .+ 2 * yvap2 .- 2 * yvap1 .+ yvapn)
    yvapn .= yvap3
    yliqn .= yliq3

    df = DataFrame(yvapn, spcnames)
    CSV.write("$(save_directiry)/simulation_vapor_liquid_yvapn_$(t+delta_t).csv", df)
    df = DataFrame(yliqn, spcnames)
    CSV.write("$(save_directiry)/simulation_vapor_liquid_yliqn_$(t+delta_t).csv", df)

    println("t = $(t+delta_t)")
    println("vap_norm = $(vap_norm)")
    println("liq_norm = $(liq_norm)")
end

@time save_rop(yliqn, shiftvectorliqn, tf, tf + delta_t, aspen_initial_conditions, Vvap, Vliq, kLAs, kHs, trays, spcnames, aspenspeciesnames, save_directiry)
