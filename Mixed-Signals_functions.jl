# Install and load required packages
using Pkg
neededPackages = [:StatsBase, :CSV, :DataFrames, :Plots, :JLD2, :KernelDensity, :KernelDensitySJ, :Statistics, :ProgressMeter, :ColorSchemes] 
using Pkg;
for neededpackage in neededPackages
    (String(neededpackage) in keys(Pkg.project().dependencies)) || Pkg.add(String(neededpackage))
    @eval using $neededpackage
end
using Plots.PlotMeasures

# Functions for initializing and updating soil and age arrays
initialise_soil = function(soildepth, nlayers, layer_thickness, bleaching_depth, grains_per_layer, bd = bd) # Fill the soil and age arrays based on the parameters
    soil = zeros(Float64, nlayers, 4)   # Array of the soil profile. Index 1: thickness [m]. Index 2: mass [kg]. Index 3: midpoint depth [m]. Index 4: cumulative depth [m].
    ages = fill(Int[],nlayers)          # Array of arrays of OSL age particles
    remainingsoil = soildepth
    for l in 1:nlayers
        if l==1 
            soil[l,1] = bleaching_depth
        else
            if remainingsoil > layer_thickness
                soil[l,1] = layer_thickness
            elseif remainingsoil <= 0
                soil[l,1] = 0
            else 
                soil[l,1] = remainingsoil
            end
        end
        remainingsoil -= soil[l,1]
    
        if l == nlayers && remainingsoil > 0
            soil[l,1] += remainingsoil
        end

        soil[l,2] = bd*soil[l,1] 
        ages[l] = fill(0, ceil(Int32, grains_per_layer))
    end
    depth = 0
    for l in 1:nlayers
        depth += soil[l,1]/2
        soil[l,3] = depth
        depth += soil[l,1]/2
    end
    soil[:,4]=cumsum(soil[:,1])
    return soil, ages
end

update_soil_and_luminescence = function(soil, ages, bd = bd)
    soil, ages = update_layers(soil, ages, bd)
    soil, ages = update_OSL(soil, ages)
    return(soil, ages)
end

update_layers = function(soil, ages, bd = bd) # Split or combine layers based on their thickness and the thickness tolerance (55%)
    for l in 1:(nlayers-1)

        if soil[l,2]/bd > (layer_thickness*(1+0.55))
            # Split
            if l < nlayers-2
                # merge bottom two layers and move layers down to create space
                soil[nlayers,2] += soil[nlayers-1,2]
                soil[nlayers-1,2] = 0
                soil, ages = transfer_OSL(soil, ages, nlayers-1, nlayers, 1,0)

                for ll in (nlayers-1):-1:(l+2)
                    soil[ll,2] = soil[ll-1,2]
                    soil[ll-1,2] = 0
                    soil, ages = transfer_OSL(soil, ages, ll-1, ll, 1,0)
                end
            end

            # Even splitting of the layers, except when resulting layer thickness exceeds initial layer thickness
            # In that case, the splitted layer becomes layer_thickness and the remainder is moved down for further splitting
            split_ratio = 0.5
            if (soil[l,2]/bd)/2 > layer_thickness
                split_ratio =  layer_thickness / (soil[l,2]/bd)
            end
            soil[l+1,2] = soil[l,2] * (1-split_ratio)
            soil[l,2] *= split_ratio
            soil, ages = transfer_OSL(soil, ages, l, l+1, (1-split_ratio), 0)
        end
        

        # Combine
        if l != 1 && soil[l,2]/bd < (layer_thickness*(1-0.55)) 
            # Determine thinnest neighbouring layer to merge with
            l2 = l+1
            if (l-1)>1 && soil[l-1,2] < soil[l+1,2]
                l2 = l-1
            end

            soil[l2,2] += soil[l,2]
            soil[l,2] = 0
            soil, ages = transfer_OSL(soil, ages, l, l2, 1,0)

            # Move layers up
            for ll in l:(nlayers-1)
                soil[ll,2] = soil[ll+1,2]
                soil[ll+1,2] = 0
                soil, ages = transfer_OSL(soil, ages, ll+1, ll, 1,0)
            end
        end
        soil = recalculate_layer_thicknesses(soil, bd)
    end
    depth = 0
    for l in 1:nlayers
        depth += soil[l,1]/2
        soil[l,3] = depth
        depth += soil[l,1]/2
    end
    soil[:,4]=cumsum(soil[:,1])
    return soil, ages
end

update_OSL = function(soil, ages) # Update OSL properties: thickness of bleaching layer and OSL ages
    # update thickness top layer to bleaching depth
    if soil[1,2]/bd > bleaching_depth # Layer too thick, give to layer below
        sep_fraction = (soil[1,2]/bd-bleaching_depth)/(soil[1,2]/bd)
        transport = sep_fraction * soil[1,2]
        soil[1,2] -= transport
        soil[2,2] += transport
        soil, ages = transfer_OSL(soil, ages, 1,2,sep_fraction,0)
    end

    if soil[1,2]/bd < bleaching_depth # Layer too thin, take from layer below
        sep_fraction = (bleaching_depth-soil[1,2]/bd)/(soil[2,2]/bd)
        transport = sep_fraction * soil[2,2]
        soil[2,2] -= transport
        soil[1,2] += transport
        soil, ages = transfer_OSL(soil, ages, 2,1,sep_fraction,0)
    end
    soil = recalculate_layer_thicknesses(soil)

    ages[1] = ages[1] .*0 # Bleach particles in the top layer
    for l in 2:nlayers # Add a year to particles in all other layers
        ages[l] = ages[l] .+ 1 
    end
    return soil, ages
end

transfer_OSL = function(soil, ages, layer, otherlayer, P_layer, P_otherlayer) # Transfer OSL particles in between layers, based on transfer probability
    ages_from = ages[layer]
    ages_to = ages[otherlayer]
    if isnan(P_layer) || isinf(P_layer)
        P_Layer = 0
    end
    if isnan(P_otherlayer) || isinf(P_otherlayer)
        P_otherlayer = 0
    end

    ind_l = sample(0:1,ProbabilityWeights([1-P_layer,P_layer]),length(ages_from)).==1
    ind_ol = sample(0:1,ProbabilityWeights([1-P_otherlayer,P_otherlayer]),length(ages_to)).==1
    
    ages[layer] = vcat(ages_from[map(!,ind_l)], ages_to[ind_ol])
    ages[otherlayer] = vcat(ages_from[ind_l], ages_to[map(!,ind_ol)]) 
    return soil, ages
end

recalculate_layer_thicknesses = function(soil, bd = bd) # Update layer thicknesses after transport of material
    for ll in 1:nlayers
        soil[ll,1] = soil[ll,2] / bd
    end
    return soil
end

# Bioturbation functions
BT_mixing = function(soil, ages, _BT_pot, _depth_function, _dd, _dd_exch = dd_exch) # Bioturbation by subsurface mixing
    if _BT_pot > 0
        depth = 0

        for l in 1:nlayers
            # Calculate BT occurring in this layer
            BT_layer = _BT_pot * BT_layer_index(soil, depth, depth + soil[l,1], _depth_function, _dd)
            if BT_layer > soil[l,2]
                BT_layer = soil[l,2]
            end

            depth += soil[l,1]/2

            if BT_layer > 0
                # Calculate exchange with other layers

                #exchange index
                tot_exch_index = -1/_dd_exch * (exp(-_dd_exch*(depth))-1)+ -1/_dd_exch*(exp(-_dd_exch*(soildepth-depth))-1)
                otherdepth = 0

                for ol in 1:nlayers
                    z_upp = otherdepth
                    z_low = otherdepth + soil[ol,1]
                    if(ol<l) # above
                        ol_exch_index = -1/_dd_exch*(exp(-_dd_exch*(depth-z_upp)) - exp(-_dd_exch * (depth - z_low)));
                    end

                    if(ol>l) # below
                        ol_exch_index = -1/_dd_exch*(exp(-_dd_exch*(z_low - depth)) - exp(-_dd_exch * (z_upp - depth)));
                    end

                    otherdepth += soil[ol,1]

                    if(l!=ol) # calculate exchange
                        interlayer_exchange = BT_layer * ol_exch_index / tot_exch_index
                        dmass_l = min(soil[l,2],interlayer_exchange / 2)
                        dmass_ol = min(soil[ol,2],interlayer_exchange / 2)

                        P_transfer_l = dmass_l / soil[l,2]
                        P_transfer_ol = dmass_ol / soil[ol,2]
                        transfer_OSL(soil, ages, l, ol, P_transfer_l, P_transfer_ol)
                    end
                end
            end
            depth += soil[l,1]/2
        end
    end
    return soil, ages
end

BT_mounding = function(soil, ages, _BT_pot, _depth_function, _dd) # Bioturbation by mounding
    if _BT_pot > 0
        depth = 0
        for l in 2:nlayers
            # Calculate BT occurring in this layer
            BT_layer = _BT_pot * BT_layer_index(soil, depth, depth + soil[l,1], _depth_function, _dd)
            if BT_layer > soil[l,2]
                BT_layer = soil[l,2]
            end

            if BT_layer > 0
                # Transport to upper layer
                P_transfer_l = BT_layer / soil[l,2]
                soil[l,2] -= BT_layer
                soil[1,2] += BT_layer
                soil, ages = transfer_OSL(soil, ages, l, 1, P_transfer_l, 0)
            end
            depth += soil[l,1]
        end
    end
    return soil, ages
end

BT_upheaval = function(md, freq, t) # Bioturbation by upheaval
    if t % freq==0 & t != ntime # determine is there is an upheaval event in this timestep

        # Pepare for mixing
        particles_contributed = zeros(Int64, nlayers)
        mixed_particles = Array{Int64}(undef, 0)
        tot_particles = 0
        depth_mixed = 0

        # Take up particles from source layers
        for l in 1:nlayers
            if depth_mixed < md
                if depth_mixed + soil[l,1] < md
                    particles_contributed[l] = length(ages[l])
                    mixed_particles = vcat(mixed_particles, ages[l])
                    ages[l] = []
                else
                    particles_contributed[l] = round(Int32,length(ages[1]) * (md - depth_mixed) / soil[l,1])
                    ages[l] = shuffle(ages[l])
                    mixed_particles = vcat(mixed_particles, ages[l][1:particles_contributed[l]])
                    ages[l] = ages[l][(particles_contributed[l]+1):length(ages[l])]
                end
                depth_mixed += soil[l,1]
            end
            if depth_mixed >= md
                break
            end
        end
        
        # Mix OSL particles
        mixed_particles = shuffle(mixed_particles)

        # And give back to layers
        tot_particles = sum(particles_contributed)

        count = 1
        for l in 1:nlayers
            ages[l] = vcat(ages[l],mixed_particles[count:(count+particles_contributed[l]-1)])
            count += particles_contributed[l]
            if count >= tot_particles
                break
            end
        end
    end
    return(soil, ages)
end

BT_layer_index = function(soil, depth, depth_low, _depth_function, _dd) # Determine fraction of bioturbation occuring in a layer, depending on the depth function
    tot_thickness = sum(soil[1:nlayers,1])

    if _depth_function == "grd" # Gradational / linear decline in rate
        if depth < (1/_dd)
            if depth_low > (1/_dd)
                depth_low = 1/_dd
            end
            layer_index = -_dd/2 * (depth_low^2 - depth^2) + (depth_low - depth)
            total_index = -_dd/2 * (1/_dd)^2 + 1/_dd
        else
            layer_index = 0
            total_index = 1
        end
    end
    
    if _depth_function == "exp" # Exponential decline in rate
        layer_index = exp(-_dd*depth) - exp(-_dd * (depth_low))
        total_index = 1 - exp(-_dd*tot_thickness)
    end

    if _depth_function == "abr" # Constant rate, or abrupt decline in rate
        if depth < _dd
            if depth_low > _dd
                depth_low = _dd
            end
            layer_index = depth_low - depth
            total_index = _dd
        else
            layer_index = 0
            total_index = 1
        end
    end

    return layer_index / total_index
end

# Calibration functions
create_calibration_parameters = function(cal_depthfunction, cal_BT_pot, cal_active_mixing_depth, cal_rel_process)
    calib_pars = zeros(Float64, length(cal_depthfunction)*length(cal_BT_pot)*length(cal_active_mixing_depth)*length(cal_rel_process), 8) 
    ct = 0
    for i1 in eachindex(cal_depthfunction)
        for i2 in eachindex(cal_BT_pot)
            for i3 in eachindex(cal_active_mixing_depth)
                for i4 in eachindex(cal_rel_process)
                    ct += 1
                    calib_pars[ct,2] = cal_BT_pot[i2]
                    if cal_depthfunction[i1]=="grd"
                        calib_pars[ct,1] = 1
                        calib_pars[ct,3] = 1 / cal_active_mixing_depth[i3]
                    end
                    if cal_depthfunction[i1]=="exp"
                        calib_pars[ct,1] = 2
                        calib_pars[ct,3] = -log(0.0025) / (cal_active_mixing_depth[i3])
                    end
                    
                    if cal_depthfunction[i1]=="abr"
                        calib_pars[ct,1] = 3
                        calib_pars[ct,3] = cal_active_mixing_depth[i3]
                    end
                    calib_pars[ct,4] = cal_rel_process[i4]
                end
            end
        end
    end
    return(calib_pars)
end  

BT_calibration_run = function(run_parameters, calib_results, i) # Function for bioturbation runs, which reads its parameters from the bioturbation scenarios
    ## Objects and parameters 
    soil, ages = initialise_soil(run_parameters[2], Int64(run_parameters[3]), run_parameters[4], run_parameters[5], Int64(run_parameters[6]), run_parameters[7])
    ntime = run_parameters[1]
    depthfunctions = ["grd", "exp", "abr"]
    # Simulations
    for t in (1:ntime)
        #BT_mounding(BT_pot * cal_BT_pot[i2], cal_depthfunction[i1], dd_ref * cal_dd[i3])
        # soil, ages = BT_mounding(soil, ages, calib_results[i,2], cal_depthfunction[Int32(calib_results[i,1])],calib_results[i,3])
        soil, ages = BT_mounding(soil, ages, calib_results[i,2] * calib_results[i,4], depthfunctions[Int32(calib_results[i,1])],calib_results[i,3])
        soil, ages = BT_mixing(soil, ages, calib_results[i,2] * (1-calib_results[i,4]), depthfunctions[Int32(calib_results[i,1])],calib_results[i,3], run_parameters[8])

        update_soil_and_luminescence(soil, ages, run_parameters[7])
    end
    return soil, ages
end

BT_calibration_errors = function(soil, ages, calibration_data, cal_mode=true, cal_iqr=true, cal_fbio=true, _ntime = ntime)
    # Calibration
    n_scores = 0
    n_cal_pars = 0
    error_mode = 0
    error_fbio = 0
    error_iqr = 0
    error_total = 0
    depths = unique(calibration_data[:,"depth"])
    for d in 1:length(depths)
        # println(d)
        it_layer = 1
        depth_ref = 0
        depth_dist = 9999
        
        # determine corresponding model layer
        for l in 1:nlayers
            depth_ref += soil[l,1]/2
            if abs(depth_ref-depths[d])<depth_dist
                depth_dist = abs(depth_ref-depths[d])
                it_layer +=1
            else
                break
            end
            depth_ref += soil[l,1]/2
        end
        
        if it_layer > nlayers
            it_layer = nlayers
        end
        # extract field ages
        calibration_ages = calibration_data[calibration_data[:,"depth"].==depths[d],"age"]
        model_ages = copy(ages[it_layer])

        # determine fbio
        fbio_exp = calibration_data[calibration_data[:,"depth"].==depths[d],"fbio"][1]
        fbio_mod = 1 - length(findall(model_ages .>= _ntime)) / length(model_ages)

        # remove non-bioturbated grains and normalize
        filter!(e-> e < _ntime,calibration_ages)
        filter!(e-> e < _ntime,model_ages)

        calibration_ages ./= ntime
        
        model_ages = model_ages ./ _ntime

        if length(model_ages) > 2 && length(calibration_ages) > 2 # if there is sufficient data available in simulated layer that is bioturbated
            # determine IQR and normalize
            if cal_iqr
                error = 1
                if std(model_ages) > 0.00001
                    iqr_exp = (Statistics.quantile(calibration_ages, 0.75)  - Statistics.quantile(calibration_ages, 0.25))
                    iqr_mod = (Statistics.quantile(model_ages, 0.75)  - Statistics.quantile(model_ages, 0.25))
                    error = (iqr_exp - iqr_mod)^2 
                end
                error_iqr += error
            end

            if cal_mode
                # Sheather and Jones bandwidth estimation
                error = 1
            
                if std(model_ages)>0.0000001 && std(calibration_ages) > 0.0000001
                    bw1=bwsj(calibration_ages)
                    bw2=bwsj(model_ages)
                    
                    if !(isnan(bw1) || isnan(bw2))
                        kde1=kde(calibration_ages,bandwidth=bw1,boundary=(0,1),npoints=2001)
                        kde2=kde(model_ages,bandwidth=bw2,boundary=(0,1),npoints=2001)
                
                        mode_field = kde1.x[kde1.density.==maximum(kde1.density)][1]
                        mode_model = kde2.x[kde2.density.==maximum(kde2.density)][1]
                        error = (mode_field - mode_model)^2
                    end
                end
                error_mode += error 
            end

            if cal_fbio
                error_fbio += (fbio_exp - fbio_mod)^2 
            end
        else
            error_mode += 1
            error_iqr += 1
            error_fbio += 1
        end
        n_scores += 1
    end

    if cal_mode
        error_total += error_mode
        n_cal_pars += 1
    end
    if cal_iqr
        error_total += error_iqr
        n_cal_pars += 1
    end
    if cal_fbio
        error_total += error_fbio
        n_cal_pars += 1
    end

    return(error_mode / n_scores, error_fbio / n_scores, error_iqr / n_scores, error_total / (n_cal_pars * n_scores))
end

# Visualization functions
plot_age_depth_profile = function(soil_all, ages_all, aggr_layers, plot_IQR = true, every_nth_layer = 1, xlim = [0,ntime], ylim = [-0.05,soildepth*1.1])
    # if input are two-dimensional arrays, add another dimension to facilitate the plotting, which is also able to plot multi-dimensional input
    if length(size(soil_all)) == 2
        soil_all = reshape(soil_all, (1, size(soil_all)...))
        ages_all = reshape(ages_all, (1, size(ages_all)...))
    end

    ndims = size(soil_all)[1]
    if ndims == 1
        colscheme = [colorant"black"]
    else
        colscheme = palette(:viridis, ndims)
    end

    plt = plot(0,  xlims = xlim, ylims = ylim, xlabel="Age [a]", ylabel="Depth [m]",plot_title="Bioturbated fraction [-]", plot_titlefontsize = 12, grid=false); yflip!(true)
    
    for i in 1:ndims
        soil = soil_all[i,:,:]
        ages = ages_all[i,:,:]

        modes = Array{Float64}(undef, 0)
        depths = Array{Float64}(undef, 0)
        fbios = Array{Float64}(undef, 0)
        IQRs = Array{Float64}(undef, 0)
        for l in 1:round(Int32,nlayers/aggr_layers)
            lay_sel = range((l-1)*aggr_layers+1, (l)*aggr_layers)
            lay_ind = lay_sel .<= nlayers
            lay_sel = lay_sel[lay_ind]
            data = vcat(copy(ages[lay_sel])...)
            l1=length(data)
            filter!(e->e<ntime,data)
            fbio = (length(data)/l1)*ntime
            depth = mean(soil[lay_sel,3])

            IQR = NaN
            mode = NaN
            if std(data) > 1
                IQR = iqr(data)

                bw = bwsj(data)
                mode = NaN
                if(!isnan(bw))
                    kde_ages = KernelDensity.kde(data,bandwidth=bw,boundary=(0,ntime*1.1),npoints=round(Int64,ntime*1.1/5+1))
                    mode = kde_ages.x[kde_ages.density.==maximum(kde_ages.density)]
                    kde_ages.density = kde_ages.density./maximum(kde_ages.density)
                    kde_ages.density[kde_ages.density.<0.01] .= NaN # remove entries with P<0.05. Can also be set to other values
                    kde_ages.density = kde_ages.density.*-0.1 .+ depth

                    if !plot_IQR
                        if (l-1) %every_nth_layer == 0
                            plot!(plt, kde_ages.x,kde_ages.density, c=colscheme[i], label = false, linestyle =:dash)    
                        end
                    end
                end
            end
            modes=vcat(modes,mode)
            depths = vcat(depths,depth)
            fbios = vcat(fbios,fbio)
            IQRs = vcat(IQRs, IQR)
        end
        plot!(plt, modes, depths,label = false,lw=2,c=colscheme[i])
        if plot_IQR
            plot!(plt, IQRs, depths,label = false,lw=1,c=colscheme[i], linestyle =:dash)
        end
        plot!(plt, fbios, depths, label = false, c=colscheme[i])
    end

    # Create the legend by plotting invisible lines and descriptions
    plot!(plt, [0], [1], label = "Mode", c=:black, lw=2)
    if plot_IQR
        plot!(plt, [0], [1], label = "Interquartile range", c=:black, lw=1, linestyle=:dash)
    else
        plot!(plt, [0], [1], label = "Age distributions", c=:black, lw=1, linestyle=:dash)
    end
    plot!(plt, [0], [1], [-10,-10], label = "Bioturbated fraction", c=:black)
    
    if ndims > 1
        for i in 1:ndims
            plot!(plt, [0], [1], [-10,-10], label = string("Scenario ",i), c=colscheme[i])
        end
    end
  
    # Plot the secondary axis for the bioturbated fraction and manually place the legend in the plot margins. Automatic placement doesn't work with secondary axes, as it will mess up the plot lay-out.
    axis2 = twiny()
    plot!(axis2,[0],[1], xlims=[0,1], ylims=ylim, label=false, grid=false);yflip!(true)
    plot!(right_margin=50mm)
    plot!(legend=(1.15,1))

    return(plt)
end

add_calibration_PDFs = function(plt, calibration_data, plot_IQR = true, ylim = [-0.05,soildepth*1.1])
    depths = unique(calibration_data[:,"depth"])
    modes = Array{Float64}(undef, 0)
    fbios = Array{Float64}(undef, 0)
    IQRs = Array{Float64}(undef, 0)
    for d in depths
        calibration_ages = calibration_data[calibration_data[:,"depth"].==d, "age"]
        l1=length(calibration_ages)
        filter!(e->e < ntime,calibration_ages)
        fbio = calibration_data[calibration_data[:,"depth"].==d, "fbio"][1]
        IQR = iqr(calibration_ages)
        # Sheather and Jones bandwidth estimation
        bw=bwsj(calibration_ages)

        kde_ages = KernelDensity.kde(calibration_ages,bandwidth=bw,boundary=(0,ntime*1.1),npoints=round(Int64,ntime*1.1/5+1))
        mode = kde_ages.x[kde_ages.density.==maximum(kde_ages.density)]
        kde_ages.density = kde_ages.density./maximum(kde_ages.density)
        kde_ages.density[kde_ages.density.<0.01] .= NaN # remove entries with P<0.05. Can also be set to other values
        kde_ages.density = kde_ages.density.*-0.1 .+ d
        
        if !plot_IQR
            plot!(plt, kde_ages.x, kde_ages.density, label = false, legend = false, c=:red, linestyle=:dash)
        end
        
        modes=vcat(modes,mode)
        fbios = vcat(fbios, fbio)
        IQRs = vcat(IQRs, IQR)
    end
    plot!(plt, modes, depths, label = false, legend = false, c=:red,lw=2)
    if plot_IQR
        plot!(plt, IQRs, depths, label = false, c=:red, legend = false, linestyle=:dash)
    end    
    plot!(plt, fbios * ntime, depths, label = false, c=:red, legend=:bottomleft)
    plot!(plt, [0], [1], label = "Data sources", linecolor = invisible())
    plot!(plt, [0], [1], label = "Simulated", c=:black)
    plot!(plt, [0], [1], label = "Experimental", c=:red)

    # Plot the secondary axis for the bioturbated fraction and manually place the legend in the plot margins. Automatic placement doesn't work with secondary axes, as it will mess up the plot lay-out.
    axis2 = twiny()
    plot!(axis2,[0],[1], xlims=[0,1], ylims=ylim, label=false, grid=false);yflip!(true)
    plot!(right_margin=50mm)
    plot!(legend=(1.15,1))
    
    return(plt)
end

plot_calibration_curves = function(cal_data)
    error_range = [0, maximum(vec(sum.(eachrow(cal_data[:, 5:7])) ./ 3))]
    plt = plot([0],[1], xlims = [0,1], ylims = error_range, label = false, xlab = "Fraction mounding", ylab = "Calibration error")

    depthfunctions = Int64.(unique(cal_data[:,1]))
    depthfunctions_labels = ["Gradational", "Exponential", "Abrupt"]
    BT_rates = unique(cal_data[:,2])

    colscheme = palette(:viridis, length(BT_rates))
    linestyles = [:solid, :dash, :dot]

    for i in depthfunctions # depth functions
        for j in eachindex(BT_rates) # potential bioturbation
            cal_temp = cal_data[cal_data[:,1] .== i .&& cal_data[:,2] .== BT_rates[j],:]
            error = vec(sum.(eachrow(cal_temp[:, 5:7])) ./ 3)
            plot!(plt, cal_temp[:,4], error, c=colscheme[j], linestyle=linestyles[i], label=false)# label = string("df ",depthfunctions[i],". BT ", BT_rates[j]))
        end
    end
    plot!(plt, [0], [1], linecolor = invisible(), label = "Depth functions")
    for i in depthfunctions
        plot!(plt, [0], [1], linestyle=linestyles[i], c=:black, label = depthfunctions_labels[i])
    end
    plot!(plt, [0], [1], linecolor = invisible(), label = "Bioturbation rate\n[kg m⁻² a⁻¹]")
    for i in eachindex(BT_rates)
        plot!(plt, [0], [1], c=colscheme[i], label = BT_rates[i])
    end
    plot!(legend=:outerright)
    return(plt)
end

# Functions for reading and writing
write_soil_ages_CSV = function(soil, ages, scenario_name, folder)
    # if input are two-dimensional arrays, add another dimension to facilitate the export, which is designed for input from different scenarios
    if length(size(soil)) == 2
        soil = reshape(soil, (1, size(soil)...))
        ages = reshape(ages, (1, size(ages)...))
    end

    # Convert ages array-of-arrays to a two-dimensional array, with each row representing a luminescence particles
    ages_matrix = Array{Int64}(undef, 0, 3)
    soil_matrix = Array{Int64}(undef, 0, 6)
    n_scenarios = size(soil)[1]
    n_layers = size(soil)[2]
   
    for i in 1:n_scenarios
        soil_matrix = vcat(soil_matrix, hcat(fill(i, n_layers),1:n_layers, soil[i,:,:]))
        for j in 1:n_layers
            ages_temp = ages[i,j]
            ages_matrix = vcat(ages_matrix, hcat(fill(i, length(ages_temp)), fill(j, length(ages_temp)), ages_temp))
        end
    end

    soil = DataFrame(soil_matrix,:auto)
    rename!(soil, ["scenario", "layer", "thickness_m", "mass_kg", "midpoint_depth_m", "cumulative_depth_m"])
    
    ages = DataFrame(ages_matrix,:auto)
    rename!(ages, ["scenario", "layer", "age_a"])

    # Create folder, if it not already exists, and write the output
    if(length(folder)) > 0
        isdir(folder) || mkdir(folder) 
    end
    CSV.write(joinpath(folder, string("soil_",scenario_name,".csv")), soil)
    CSV.write(joinpath(folder, string("ages_",scenario_name,".csv")), ages)
end

read_soil_ages_CSV = function(soil_directory, ages_directory)
    soil_read = CSV.read(soil_directory, DataFrame)
    ages_read = CSV.read(ages_directory, DataFrame)

    n_scenarios = Int(maximum(soil_read[:,"scenario"]))
    n_layers = Int(maximum(soil_read[:,"layer"]))
    
    # Create empty soil and age arrays
    soil = zeros(Float64, n_scenarios, n_layers, 4)
    ages = fill(Int[],n_scenarios, n_layers)

    for i in 1:n_scenarios
        soil[i,:,:] = Matrix(soil_read[soil_read[:,1].==i,3:6])
        for j in 1:n_layers
            ages[i,j] = ages_read[ages_read[:,1].==i .&& ages_read[:,2].==j,3]
        end
    end
    return soil, ages
end

write_soil_ages_JLD2 = function(soil, ages, scenario_name, folder)
       # if input are two-dimensional arrays, add another dimension to facilitate the plotting, which is also able to plot multi-dimensional input
       if length(size(soil)) == 2
        soil = reshape(soil, (1, size(soil)...))
        ages = reshape(ages, (1, size(ages)...))
    end

    # Create folder, if it not already exists, and write the output
    if(length(folder)) > 0
        isdir(folder) || mkdir(folder) 
    end
    save_object(joinpath(folder, string("soil_",scenario_name,".jld2")), soil)
    save_object(joinpath(folder, string("ages_",scenario_name,".jld2")), ages)
end

read_soil_ages_JLD2 = function(soil_directory, ages_directory)
    soil = load_object(soil_directory)
    ages = load_object(ages_directory)
    return soil, ages
end

