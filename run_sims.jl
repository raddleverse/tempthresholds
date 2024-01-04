using Mimi
using MimiCIAM
using Query
using CSV
using StatsBase
using NetCDF
using DataFrames

include("montecarlo/ciamMonteCarlo.jl")
include("montecarlo/defmcs.jl")
include("montecarlo/run_ciam_mcs.jl")
include("montecarlo/brickLSL.jl")
include("montecarlo/processResults.jl")

#brickfile = "https://zenodo.org/record/6626335/files/sneasybrick_projections_csv.zip"
brickfile = "sneasybrick_projections_csv.zip"
outputdir = joinpath(@__DIR__, "output", "MonteCarlo")
isdir(outputdir) || mkpath(outputdir)

ssp_files = Dict(1 => "IIASAGDP_SSP1_v9_130219",
                 2 => "IIASAGDP_SSP2_v9_130219",
                 3 => "IIASAGDP_SSP3_v9_130219",
                 4 => "IIASAGDP_SSP4_v9_130219",
                 5 => "IIASAGDP_SSP5_v9_130219")
popinput = 0                          # population density input data (only 0 is supported currently)
ssp_rcp_scenarios = [(1,26), (2,45), (4,60), (5,85)]  # what combinations of SSP (first) and RCP (second)?
nensemble = 10                      # how many ensemble members for the Monte Carlo?
surgeoption = 0  # which surge data sets to use (0 = original CIAM/DINAS-COAST; 1 = GTSR-corrected D-C; 2 = GTSR nearest data points)

#for (ssp, rcp) in ssp_rcp_scenarios

    println("Running SSP",ssp, "-RCP",rcp,"...")

    # write the init file
    init_settings = Dict(
        :init_filename   => "MCdriver_init.csv",
        :subset          => false,
        :ssp             => ssp_files[ssp],
        :ssp_simplified  => ssp,
        :rcp             => rcp,
        :b => "lsl_rcp85_p50.csv" # won't matter, overwritten when doing the SLR sampling anyway
    )
    init_file = joinpath(outputdir,init_settings[:init_filename])
    textheader="run_name,lslr,subset,ssp,ssp_simplified\n"
    textstr = "base,$(init_settings[:b]),$(init_settings[:subset]),$(init_settings[:ssp]),$(init_settings[:ssp_simplified])"
    txtfile=open(init_file,"w") do io
        write(io,textheader)
        write(io,textstr)
    end

    # Define input trial parameters: BRICK model, number of trials (n), min and max BRICK percentile,
    # start and end year, and timestep. Resetting each time through the loop for consistency.
    trial_params = Dict(
        :brickfile  => brickfile,
        :n      => nensemble,
        :high   => 100,
        :low    => 0,
        :ystart => 2010,
        :yend   => 2150,
        :tstep  => 10,
    )

    # Define other eneded parameters and settings for the model
    adaptRegime1 = Dict(
        :fixed          => true,
        :t              => 15,
        :noRetreat      => false,
        :allowMaintain  => false,
        :popval         => popinput,
        :GAMSmatch      => false,
        :surgeoption    => surgeoption,
        :subset         => false
    )

    # vary SLR and CIAM parameters
    trial_params[:low] = 0
    trial_params[:high] = 100
    runname = string("SSP",init_settings[:ssp_simplified],"_BRICK",init_settings[:rcp],"_global")
    runTrials(init_settings[:rcp], init_settings[:ssp_simplified], trial_params, adaptRegime1, outputdir, init_file, vary_slr=true, vary_ciam=true, runname=runname)