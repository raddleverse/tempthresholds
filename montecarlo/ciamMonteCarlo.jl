##==============================================================================
## ciamMonteCarlo.jl
##
## Original code: Catherine Ledna (4 Feb 2021)
## Modified code: Tony Wong (16 May 2021)
##==============================================================================

# The file structure created by this process will look as follows:
#
# - a (unique) top directory created with the code:
#   joinpath(outputdir, runname, "CIAM $(Dates.format(now(), "yyyy-mm-dd HH-MM-SS")) MC$(trial_params[:n])")
#   which will hold a RawResults folder with results directly from the monte carlo
#   runs in run_ciam_mcs.jl and a PostProcessing folder written with the code at
#   the bottom of ciamMonteCarlo.jl

function runTrials(rcp, ssp, trial_params, adaptRegime, outputdir, init_filepath; vary_slr = true, vary_ciam=true, runname="default_run")


    # make an output directory for this SSP-RCP scenario
    scenario = "SSP$(ssp)_BRICK$(rcp)"
    scenario_outputdir = joinpath(outputdir, scenario)
    isdir(scenario_outputdir) || mkpath(scenario_outputdir)

    # within the scenario output directory, make one for this particular set of simulations
    #outputdir = joinpath(scenario_outputdir, runname, "CIAM $(Dates.format(now(), "yyyy-mm-dd HH-MM-SS")) MC$(trial_params[:n])")
    outputdir = joinpath(scenario_outputdir, runname, "CIAM MC$(trial_params[:n])")
    isdir(outputdir) || mkpath(outputdir)

    # Output Files: Trials, NPV, Global Time Series, Regional Spotlight

    # Load CIAM parameters from file
    if adaptRegime[:subset] == false
        segIDs = false
    else
        subs = MimiCIAM.load_subset(adaptRegime[:subset])
        sort!(subs)
        segIDs = MimiCIAM.segStr_to_segID(subs)
    end

    # Get non-climatic LSLR
    # Years are 2010, 2020, ... , 2200
    data_dir = joinpath(@__DIR__, "..", "data", "lslr")
    lsl_nonclim = CSV.read(joinpath(data_dir, "lsl_rcp0_p50.csv"), DataFrame) |> DataFrame
    years_nonclim = range(2010, stop=2200, length=20) |> collect

    # Filter according to subset segments
    # Will still have possibly too many time steps
    if adaptRegime[:subset] != false
        lsl_nonclim = lsl_nonclim[!, subs]
    end

    # Need to make sure that lsl_nonclim has same segment order as lslr (climatic, below) 
    segnames = get_segnames(segIDs)
    col_names = [i for i in names(lsl_nonclim) if string(i) in segnames]
    col_names = sort(col_names)
    lsl_nonclim = lsl_nonclim[!, col_names]

    # Load BRICK data
    # - segIDs will subset (if needed), and lsl_nonclim already subsetted
    # - years_nonclim and lsl_nonclim go into brick_lsl, to be trimmed and added
    # - lsl that comes out (lsl[1]) will include lsl_nonclim added
    lsl = brick_lsl(rcp, segIDs, trial_params[:brickfile], trial_params[:n], lsl_nonclim, 
                    years_nonclim, trial_params[:low], trial_params[:high],
                    trial_params[:ystart], trial_params[:yend],
                    trial_params[:tstep], false)
    if vary_slr
        lslr    =lsl[1]
        gmsl    =lsl[2]
        ensInds =lsl[3] # Indices of original BRICK array
        temps   =lsl[4] # Full temperature time series
        years   =lsl[5]
        # reference temperatures to 1850-1900 mean
        ref_temps = mean(eachrow(temps[findall(t->t in (1850:1900), years), :]))
        temps_norm = temps - repeat(transpose(ref_temps), size(temps)[1],1)
        temps_norm_2100 = temps_norm[findall(t->t==2100, years),:]

    elseif trial_params[:low] == trial_params[:high]
        lslr = repeat(lsl[1],outer=(trial_params[:n], 1, 1))
        gmsl = repeat(lsl[2], outer = (1,trial_params[:n]))
        ensInds = fill(lsl[3],trial_params[:n]) # only 1 element coming back so `fill` instead of `repeat`

    else
        error("Not varying SLR in Monte Carlo sampling, but the low and high quantiles requested are not equal.")
    end

    num_ens = trial_params[:n]

    m = MimiCIAM.get_model(t = adaptRegime[:t], initfile = init_filepath,
                            fixed = adaptRegime[:fixed], noRetreat = adaptRegime[:noRetreat],
                            allowMaintain = adaptRegime[:allowMaintain], popinput = adaptRegime[:popval],
                            surgeoption = adaptRegime[:surgeoption])
 
    # TODO do we want these? They are the central values as used in run_ciam_mcs.jl
    update_param!(m, :slrcost, :movefactor, 1) # From Diaz (2016); incorporates communication with Mendelsohn and Anthoff and Tol (2014)
    update_param!(m, :slrcost, :dvbm, 5.376)    # Updated to 2010USD from FUND
    update_param!(m, :slrcost, :vslel, 0.47)     # From Viscusi and Aldy (2003)
    update_param!(m, :slrcost, :vslmult, 200)      # From FUND (originally Cline (1992))
    update_param!(m, :slrcost, :wvel, 1.16)      # From Brander et al (2006)
    update_param!(m, :slrcost, :wvpdl, 0.47)     # From Brander et al (2006)

    # get the segments and their corresponding World Bank regions
    dfSR = CSV.read("../data/segments_regions_WB.csv", DataFrame)
    dfSR[!,"ids"] = [parse(Int64,replace(i, r"[^0-9]"=> "")) for i in dfSR[!,"segments"]]

    # unique World Bank regions
    wbrgns = unique(dfSR[!,"global region"])

    # array of distributions, could also pre-compute these if you want
    distribs = (
        # (movefactor, Truncated(Normal(1,1),0.5,3) # From Diaz (2016); incorporates communication with Mendelsohn and Anthoff and Tol (2014)
        (param = :dvbm, distrib = Truncated(Normal(5.376,2.688),0.0,Inf)), # Updated to 2010USD from FUND
        (param = :vslel, distrib =Truncated(Normal(0.47,0.15),0.0,Inf)), # From Viscusi and Aldy (2003)
        (param = :vslmult, distrib =Truncated(Normal(200,100),0.0,Inf)), # From FUND (originally Cline (1992))
        (param = :wvel, distrib =Truncated(Normal(1.16, 0.46),0.0,Inf)), # From Brander et al (2006)
        (param = :wvpdl, distrib =Truncated(Normal(0.47,0.12),0.0,1.0)) # From Brander et al (2006)
    )

    # preallocate arrays
    if vary_ciam
        outtrials = Array{Float64}(undef, num_ens, length(distribs)) # array of random CIAM parameter draws
    end 
    outts = DataFrame() # Dataframe of outputs, TODO performance would improve if this was an Array but fine for now
    globalNPV = Array{Float64}(undef, num_ens)
    regionNPV = Array{Float64}(undef, num_ens, length(wbrgns))
    regionNPV1 = Array{Float64}(undef, num_ens, length(wbrgns)) # for first time step optimal costs, to subtract off later

    run(m) # Run the model once so we instantiate a model instance to modify

    p = Progress(num_ens; showspeed=true)

    for i = 1:num_ens

        # Progress Meter
        next!(p; showvalues = [(:iter,i)])

        # Work with the model instance when updating parameters instead of the 
        # model def to prevent rebuilding, this would be handled with the Mimi
        # MCS framework using an empirical distribution for lslr but this keeps
        # your general structure

        update_param!(m.mi,:slrcost, :lslr,view(lslr, i,:,:)) # use a view to prevent allocation from slices

        if vary_ciam
            for (j, distrib) in enumerate(distribs)
                draw = rand(distrib.distrib)
                outtrials[i,j] = draw
                update_param!(m.mi, :slrcost, distrib.param, draw)
            end
        end

        # Run the model again, noting it will not rebuild because we updated 
        # m.mi under the hood
        run(m)

        # Peformance TODO This section is slow because the getTimeSeries
        # function is complex, and less importantly we are appending to a DataFrame
        # instead of a preallocated Array, but leavign it alone for now.
        ts = MimiCIAM.getTimeSeries(m,i,rgns=false,sumsegs="global")
        append!(outts, ts)

        globalNPV[i] = m[:slrcost,:NPVOptimalTotal]

        # Performance TODO: does this need to happen within the loop? Not a priority it is fast
        # get the segments for each region and aggregate
        for rgn in wbrgns
            segIDs_rgn = filter(:"global region" => ==(rgn), dfSR)[!,"ids"]
            idx_rgn = findall(x -> x in segIDs_rgn, dfSR[!,"ids"])
            col_rgn = findfirst(x->x==rgn, wbrgns)
            regionNPV[i,col_rgn] = sum(m[:slrcost,:NPVOptimal][idx_rgn])
            regionNPV1[i,col_rgn] = sum(m[:slrcost,:OptimalCost][1,idx_rgn]) # regional (World Bank regions) first time step optimal costs
        end
    end

    # get regional NPV as DataFrame for output
    outregionNPV = DataFrame(regionNPV, :auto)
    rename!(outregionNPV,wbrgns)

    # regional (World Bank regions) first time step optimal costs
    outregionNPV1 = DataFrame(regionNPV1, :auto)
    rename!(outregionNPV1,wbrgns)

    # Write Trials, Global NPV and Time Series
    postprocessing_outputdir = joinpath(outputdir, "PostProcessing")
    isdir(postprocessing_outputdir) || mkpath(postprocessing_outputdir)

    outtrialsname = joinpath(postprocessing_outputdir, "trials_$(runname).csv")
    outnpvname= joinpath(postprocessing_outputdir, "globalnpv_$(runname).csv")
    outrgnname = joinpath(postprocessing_outputdir, "regionnpv_$(runname).csv")
    outtsname= joinpath(postprocessing_outputdir, "globalts_$(rcp)_$(runname).csv")

    # Turn outtrials into a DataFrame and save if vary_ciam was true
    if vary_ciam 
        outtrials_df = DataFrame(outtrials,[distrib.param for distrib in distribs])
        insertcols!(outtrials_df, 1, :ens => i)
        CSV.write(outtrialsname, outtrials)
    end

    procGlobalOutput(globalNPV,gmsl,temps_norm_2100,ensInds,trial_params[:brickfile],rcp,adaptRegime[:noRetreat],outnpvname)
    CSV.write(outtsname, outts)
    CSV.write(outrgnname, outregionNPV)

    # also write out the regional (World Bank regions) first time step optimal costs
    outrgnname1= joinpath(postprocessing_outputdir, "regionts_$(runname).csv")
    CSV.write(outrgnname1, outregionNPV1)

end

##==============================================================================
## End
##==============================================================================
