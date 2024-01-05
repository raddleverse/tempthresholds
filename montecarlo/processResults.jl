##==============================================================================
## processResults.jl
##
## Original code: Catherine Ledna (4 Feb 2021)
## Modified code: Tony Wong (16 May 2021)
##==============================================================================

# Functions to process CIAM data and make plots

function procGlobalOutput(glob,gmsl,temp2100,inds,brickfile,rcp,noRetreat,outfile=false,tstart=2010,tend=2100)

    # Get Global NPV
    npv = glob

    #gmsl_fd = gmsl[10,fd_inds]
    #gmsl_nofd=gmsl[10,nofd_inds]
    gmsl2040 = gmsl[4,:]
    gmsl2060 = gmsl[6,:]
    gmsl2080 = gmsl[8,:]
    gmsl2100 = gmsl[10,:]
    # CAUTION! this is hardwired for 2010-2100

    if noRetreat==false
        lab="With Retreat"
    else
        lab="No Retreat"
    end

    df1 = DataFrame([npv gmsl2100 transpose(temp2100)], :auto)
    #df1[:brick]="Fast Dynamics"
    df1[!,:retreat] = fill(lab, size(df1)[1])

    #df2 = DataFrame([npv_nofd gmsl_nofd], :auto)
    #df2[:brick]="No Fast Dynamics"
    #df2[:retreat]=lab

    #outdf = [df1;df2]
    outdf = df1

    rename!(outdf,[:npv,:gmsl2100,:temp2100,:retreat])
    outdf[!,:brickEnsInd] = inds

    if outfile==false
        CSV.write("output/ciammcs/globalNPV_rcp$(rcp)_noRetreat$(noRetreat)_mcs_newPop_maintainFalse.csv",outdf)
    else
        CSV.write(outfile,outdf)
    end

end

function procSegResults(m,seg,lsl,inds,brickfile,rcp,noRetreat,tstart=2010,tend=2100)
    segmap = load_segmap()
    segIDs = m[:slrcost,:segID]
    segNames =segID_to_seg(segIDs,segmap)

    npv=seg[1]
    opt2050=seg[2]
    lev2050=seg[3]
    opt2100=seg[4]
    lev2100=seg[5]
  #  lsl2050=lsl45[:,5,:]
  #  lsl2100=lsl45[:,10,:]

    for k in string.(keys(dict1))
        val = dict1[k]
        md = [mode(val[:,i]) for i in 1:size(val)[2]]
        minval = minimum(val,dims=1)'
        maxval = maximum(val,dims=1)'
        df = DataFrame([segNames md minval maxval], :auto)
        names!(df,[:segment,:modalOption,:minOpt,:maxOpt])
        CSV.write("output/segAdapt_stats_rcp$(rcp)_noRetreat$(noRetreat)_$(k).csv",df)
    end
end

function getbricktime(brickfile,tstart,tend)
    years=ncread(brickfile,"time_proj")
    start_ind =findall(x->x==tstart,years)[1]
    end_ind=findall(x-> x==tend,years)[1]
    return start_ind,end_ind
end

# Function: costs as percentage of GDP
function plotMap(tabstr)
    ciamLonLat = CSV.read("data/diva_segment_latlon.csv")
    df = CSV.read(tabstr)
end

# Function: plot costs on map (5-95%, outliers as insets)

# Function: plot dist
# function plotDists(tabstr,globOrSeg="glob")
#     tab = CSV.read(tabstr)
#     if globOrSeg=="glob"
#          ## GMSL Distribution Plot
#          gmslPlot = @df tab density(:gmsl,group=(:brickOutput),label=[:brickOutput])

#     else
#     end

# end
