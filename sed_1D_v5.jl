#################################################################################
#
# sed_1d_v5.jl - a Julia script for modeling 1-D sediment transport
#
# concept is to exchange sediments in water column with sediment bed layers
# based on differences in river sediment-carrying capacity
#
#################################################################################


### constants

const rho_w = 1000.		# density of suspending medium (~water)
const rho_s = 2650.     # density of sediment particles
const rho_b = 1600.     # bulk sediment density
const g = 9.807         # gravitational acceleration
const nu = 0.9e-6       # kinematic viscosity of water (m^2/sec)
const s = rho_s/rho_w  	# sediment particle specific gravity


### sediment physical characteristics

type Sediment
	name::AbstractString
    d50::Float64
    d90::Float64
    floc::Float64 					# multiplier for increasing d50 to account for floccing
    upgradient::Float64 			# fractional abundance of this sediment class included in the upgradient boundary condition
	dnd::Float64 					# dimensionless particle size
    vs::Float64 					# settling velocity
end

function AssignSeds()
	# read sediment characteristics file and assign individual sediments to respective type
	sediment = Sediment[]
	data = readdlm("sediments.txt", '\t', header=true)
	for i = 1:size(data[1], 1)
		name = data[1][i, 1]
		d50 = Float64(data[1][i, 2])
		d90 = Float64(data[1][i, 3])
		floc = Float64(data[1][i, 4])
		upgradient = Float64(data[1][i, 5])
		dnd = Dnd(d50, 1.0)
		vs = SettlingVel(d50, floc)
		push!(sediment, Sediment(name, d50, d90, floc, upgradient, dnd, vs))
	end
	return sediment
end

function Dnd(d50::Float64, floc::Float64)
	# non-dimensional particle diameter
    return ((s-1.)*g/nu^2)^(1./3.) * (d50 * floc)
end

function SettlingVel(d50::Float64, floc::Float64)
	# particle settling velocity (Cheng, 1997 approximation)
    return nu/(d50*floc) * ((25. + 1.2*Dnd(d50,floc)^2)^0.5 - 5.)^1.5
end


### dredging history

type Dredge
	t_ex::Array{Float64, 1}			# time of dredging event
	dz::Array{Float64, 2} 			# dredging depth; rows = time events, columns = cells
end

function AssignDredge()
	# read dredge event file and assign all to Dredge type
	data = readdlm("dredge.txt", '\t', header=true)
	t_ex = data[1][:,1]
	dz = data[1][:,2:end]
	dredge = Dredge(t_ex, dz)
	return dredge
end


### model parameters and supporting functions

type Model
    num_cells::Int64 			# number of river cells included in the model
	w_up::Float64 				# wdith of river section upgradient of first model cell
    hw::Float64 				# height (thickness) of water column
    elev_w::Float64 			# elevation of top of water column above datum
    dt::Float64 				# time step size (fixed)
    t_max::Float64 				# simulation end time
    t_profile::Float64 			# time step size for writing profile output files
    t_series::Float64 			# time step size for monitoring output
    num_slices::Int64 			# (initial) number of slices for sediment bed
    mixing_depth::Float64 		# scour/mixing depth (fixed)
    mode_Q::Int64 				# 0 = fixed discharge rate; 1 = variable (probabilistic) discharge rate
    log_min_Q::Float64 			# minimum discharge rate (log scale)
    log_max_Q::Float64 			# maximum discharge rate (log scale)
    log_avg_Q::Float64 			# average discharge rate (log scale)
    log_stdev_Q::Float64 		# standard deviation of the log discharge rate
    dt_Q::Float64 				# time interval after which to select a new discharge rate if mode_Q == 1
end

function ModelParams()
	# read model parameters file and assign values to Model type
	data = readdlm("model_params.txt", '\t', header=false)
	num_cells = Int64(data[1, 2])
	w_up = Float64(data[2, 2])
	hw = Float64(data[3, 2])
	elev_w = Float64(data[4, 2])
	dt = Float64(data[5, 2])
	t_max = Float64(data[6, 2])
	t_profile = Float64(data[7, 2])
	t_timeseries = Float64(data[8, 2])
	num_slices = Int64(data[9, 2])
	mixing_depth = Float64(data[10, 2])
	mode_Q = Int64(data[11, 2])
	log_min_Q = Float64(data[12, 2])
	log_max_Q = Float64(data[13, 2])
	log_avg_Q = Float64(data[14, 2])
	log_stdev_Q = Float64(data[15, 2])
	dt_Q = Float64(data[16, 2])
	model = Model(num_cells, w_up, hw, elev_w, dt, t_max, t_profile, t_timeseries, num_slices, mixing_depth, mode_Q, log_min_Q, log_max_Q, log_avg_Q, log_stdev_Q, dt_Q)
	return model
end

function Velocity(Q::Float64, h::Float64, w::Float64)
    # mean flow velocity
    return Q/(w * h)
end
		
function Me(q::Float64, ucr::Float64, sed::Sediment)
    # mobility parameter for simplified load models (for all class sizes)
    if q > ucr
		me = (q-ucr)/((s-1)*g*sed.d50)^0.5
    else
		me = 0.
	end
    return me
end		

function Ucr(h::Float64, sed::Sediment)
    # critical velocity for currents (for all class sizes)
    if sed.d50 < 0.0005
		ucr = 0.19*sed.d50^0.1 * log10(12*h/(3.*sed.d90))
    else
		ucr = 8.5*sed.d50^0.6 * log10(12*h/(3.*sed.d90))
	end
    return ucr
end
		
function TransBL(q::Float64, h::Float64, sed::Sediment)
    # simplified bed load transport formula for steady-state current flow (for all class sizes)
    ucr = Ucr(h, sed)
    me = Me(q, ucr, sed)
    return 0.015 * rho_s * q * h * (sed.d50/h)^1.2 * me^1.5      # dimensions = mass length^-1 time^-1
end

function TransSL(q::Float64, h::Float64, sed::Sediment)
	# simplified suspended load transport formula for steady-state current flow
    ucr = Ucr(h, sed)
    M_e = Me(q, ucr, sed)                            
    return 0.012 * rho_s * q * sed.d50 * M_e^2.4 * sed.dnd^(-0.6)    # dimensions = mass length^-1 time^-1
end


### stratigraphy: data and methods

type Strat
    dz::Float64 									# model individual slice thickness (fixed at the start of the simulation)
    active_dz::Float64 								# thickness of the top (active) slice, the only slice which can have a thickness != dz
    lith_frac::Array{Array{Float64,1},1}			# array of volume proportions of different sediment classes represented in each slice
end

function AssignStrat(model::Model)

	# read stratigraphy file (per model cell), which specifies sediment layer characteristics
	data = readdlm("stratigraphy.txt", '\t', header=true)
	strat = Strat[]
	num_seds = size(data[1], 2)-3 					# note: sediment ordering is assumed to be preserved in input file; no check is made ...	
	start_line = zeros(Int64, model.num_cells) 		# starting point for cell[i] in the table
	num_layers = zeros(Int64, model.num_cells) 		# number of layers corresponding to cell[i] in the table 
	cell_indices = data[1][:,1]
	for i = 1:model.num_cells
		start_line[i] = minimum(findin(cell_indices,i))
		num_layers[i] = length(findin(cell_indices,i))		
	end
	
	for icell = 1:model.num_cells 					# for each cell
		bottom = Float64[]
		top = Float64[]
		lith_frac_matrix = Array[] 					# this is just a container for sets of Float64 arrays
		for i = start_line[icell]:(start_line[icell] + num_layers[icell]) - 1		# for each defined layer
			push!(bottom, Float64(data[1][i, 2]))
			push!(top, Float64(data[1][i, 3]))
			lith_frac_layer = Float64[]		
			for j = 1:num_seds
				push!(lith_frac_layer, Float64(data[1][i, 3+j]))
			end
			push!(lith_frac_matrix, copy(lith_frac_layer))
		end
		
		# distribute layer characteristics across appropriate model slices
		dz = (maximum(top) - minimum(bottom))/model.num_slices     			# model slice thickness will always remain fixed, even for new added layers
		active_dz = dz                            							# initial active slice thickness
		lith_frac = Array[] 												# just a container for sets of Float64 arrays
		for i = 1:model.num_slices
			midpoint = (i+0.5) * dz
			for j = 1:num_layers[icell]
				if midpoint >= bottom[j] && midpoint < top[j]
					push!(lith_frac, copy(lith_frac_matrix[j]))  			# this 2D array is all that ever needs to be tracked subsequently, since dz is already specified
					break
				end
			end
		end
		push!(strat, Strat(dz, active_dz, lith_frac))
	
	end
	
	return strat
	
end

function Dig(strat::Strat, water_depth::Float64, dredge_depth::Float64)
    # remove sediment to match dredging event authorized depth
	dig = max(dredge_depth - water_depth, 0.0)
    num_dig_layers = Int64(floor(dig/strat.dz))
	for i = 1:num_dig_layers
		deleteat!(strat.lith_frac, 1)
	end
    strat.active_dz = strat.dz
	return dig
end

function Blend(strat::Strat, d::Float64)
	# determine size class composition of mixing zone at top of sediment bed
	blend_frac = zeros(strat.lith_frac[1]) 			# composite sediment class composition of mixed zone (initialize array)
	blend_w = Float64[] 							# weighting fractions associated with active slices
	d_1 = max(d - strat.active_dz,0.)        		# extent blending zone beneath active layer
	d_2 = d_1 % strat.dz 							# division remainder
	num_mix_layer = Int64(floor(d_1/strat.dz)) 		# division quotient
	# weight the influences of blending layers on average sediment composition
	push!(blend_w,min(strat.active_dz/d, 1.))       # active layer
	for i = 1:num_mix_layer
		push!(blend_w, strat.dz/d) 					# complete layers underneath active layer
	end
	if d_2 > 0.
		push!(blend_w, d_2/d)                  		# partial layer on bottom
	end	
	for j = 1:length(blend_frac), i = 1:length(blend_w)
		blend_frac[j] += blend_w[i] * strat.lith_frac[i][j]
	end
	return blend_frac
end

function Deposit(strat::Strat, d::Float64, blended_lith::Array{Float64,1}, dh::Array{Float64,1})         

    # mix in deposited material (positive) or remove (negative) from blended zone
	for j = 1:length(dh)
		dh[j] = max(dh[j],-d*blended_lith[j])          # net change in sediment thickness
	end	
    delta_d = sum(dh)

    # (1) determine updated average composition of sediment mixture down to the blending depth
	for j = 1:length(dh)
		blended_lith[j] = ((d-delta_d)*blended_lith[j] + dh[j])/d
	end
	
	# (2) determine any layering changes
	if delta_d + strat.active_dz > strat.dz
		# need to add new layer(s)
		
		global num_added_layers = Int64(floor((delta_d + strat.active_dz)/strat.dz))		# division quotient	
		
		strat.active_dz = (delta_d + strat.active_dz) % strat.dz 					# division remainder
		for i = 1:num_added_layers
			insert!(strat.lith_frac, 2, copy(blended_lith)) 						# insert just beneath active layer
		end
	elseif (delta_d < 0.) && (delta_d + strat.active_dz <= 0.)
		# need to remove layer(s)
		del_active_dz = abs(delta_d + strat.active_dz) % strat.dz 						# division remainder
		num_del_layers = Int64(floor(abs(delta_d + strat.active_dz)/strat.dz)) + 1		# division quotient (plus one layer)
		strat.active_dz = strat.dz - del_active_dz
		for i = 1:num_del_layers
			deleteat!(strat.lith_frac, 1)
		end
	else
		# active layer thickness changes, but no change in layer numbering
		strat.active_dz += delta_d
	end
	
	# (3) apply compositional changes to top layers from surface down to depth not compensated by layering changes
    strat.lith_frac[1] = copy(blended_lith) 										# active layer
    comp_depth = max(d - strat.active_dz, 0.)         								# depth into sediment pile, beyond active layer
	bottom_blend_d = comp_depth % strat.dz 											# division remainder
	num_comp_layers = Int64(floor(comp_depth/strat.dz))								# division quotient
    for i = 1:num_comp_layers
		strat.lith_frac[i+1] = copy(blended_lith) 									# layers beneath active layer
	end
	# bottom boundary of mixing zone
	for j = 1:length(dh)
		strat.lith_frac[num_comp_layers+2][j] = (bottom_blend_d/strat.dz)*blended_lith[j] + (1.-bottom_blend_d/strat.dz)*strat.lith_frac[num_comp_layers+2][j]
	end
	
end


### external sources

type Source
	t_src::Array{Array{Float64,1},1} 					# list of source time arrays (one array per active-source cell-sediment combo)
	J::Array{Array{Float64,1},1} 						# matching list of flux arrays
	index::Array{Tuple{Int64, Int64}, 1} 				# list of tuples (cell, material index)
    load::Array{Float64, 1} 							# cumulative loading, per sediment size class (i.e., 1-D array)
	sed_indices::Array{Int64,1} 						# set of sediment list indices employed as sources
	sed_names::Array{AbstractString,1} 					# list of sediment names accompanying sed_indices
end

function AssignSources(sediment::Array{Sediment,1})

	# assign t_src and J arrays to source arrays; create corresponding list of tuples (cell, material index) to track source arrays

	t_src  = Array[] 									# initial multi-source time and flux array container
	J = Array[]
	index_list = Tuple{Int64, Int64}[] 					# initialize list of tuples to hold cell and sediment indices for source(s)
	load = zeros(Float64,length(sediment)) 				# initialize cumulative loading array (size = number of sediment particle classes)
	sed_name = AbstractString[] 						# names of all sediment classes in sediment
	sed_names = AbstractString[] 						# subset of names of those classes with external sources
	sed_indices = Int64[] 								# sediment class indices associated with external sources
	sed_names = AbstractString[] 						# subset of names of those classes with external sources
	time_array = Float64[] 								# temporary holders for housing time and flux arrays
	flux_array = Float64[]	
	
	for i = 1:length(sediment)
		push!(sed_name, sediment[i].name)
	end	
	data = readdlm("sources.txt", '\t', header=true)
	
    for i = 1:1:size(data[1], 1)

		cell_index = data[1][i, 4]	
		sed_index = findin(sed_name,[data[1][i, 1]])[1]
		push!(sed_indices, sed_index)
		
		if i == 1
			# this is a new cell and/or sediment type (in source definition)		
			push!(index_list, (cell_index, sed_index)) 
			push!(time_array, data[1][i, 2])
			push!(flux_array, data[1][i, 3])			
	
		elseif (i > 1) && ((cell_index, sed_index) != index_list[end])
			# this is a new cell and/or sediment type (in source definition); data exists for a prior source term 
			push!(t_src, time_array)
			push!(J, flux_array)
			push!(index_list, (cell_index, sed_index)) 			
			time_array = Float64[] 
			flux_array = Float64[]	
			push!(time_array, data[1][i, 2])
			push!(flux_array, data[1][i, 3])
			
		else
			# this is more of the same source term; add to time and mass flux arrays
			push!(time_array, data[1][i, 2])
			push!(flux_array, data[1][i, 3])			
			
		end
		
	end
	
	# capture final read
	push!(t_src, time_array) #####
	push!(J, flux_array)	
	
	# match names on sed_name list to sed_indices
	sed_indices = sort(collect(Set(sed_indices)))	
	for i in sed_indices
		push!(sed_names, sed_name[i])
	end
	
	# populate source type
	source = Source(t_src, J, index_list, load, sed_indices, sed_names)	


	#println(source.t_src[1])


	
	return source
	
end

function SourceTerm(icell::Int64, ised::Int64, t::Float64, source::Source)
	# return a matching source flux for ised in icell at t
	loc = (icell, ised)
	find_flag = false
	source_index = 0
	for (i, value) in enumerate(source.index)
		if loc == value
			find_flag = true 							# note that i is the index to use in the t and J lists
			source_index = i
			break
		end
	end
	if find_flag == true
		# compute flux from arrays
		J_index = 1 + round(Int, (t - source.t_src[source_index][1]) / (source.t_src[source_index][end] - source.t_src[source_index][1]) * length(source.t_src[source_index]))
		if (J_index < 1) || (J_index > length(source.t_src[source_index]))
			flux = 0.
		else
			flux = source.J[source_index][J_index]
		end
	else
		flux = 0.
	end
	return flux
end


### cell geometry

type Cell
    w::Float64 					# river cell width
    area::Float64 				# areal footprint of the river cell
end

function CellGeom(num_cells::Int64)
	# read cell geometry specs file
	cell = Cell[]
	data = readdlm("cells.txt", '\t', header=true)
	for i = 1:num_cells
		w = Float64(data[1][i, 2])
		L = Float64(data[1][i, 3])
		area = w * L
		push!(cell, Cell(w, area))
	end
	return cell
end


### utility functions

function SettlingDep(J0_SL::Array{Float64}, J0_BL::Array{Float64}, J_SL::Array{Float64,2}, J_BL::Array{Float64,2}, sediment::Array{Sediment,1}, source::Source, model::Model, Q::Float64, t::Float64)
	# compute net deposition per time step using the settling velocity method (for comparison to model algorithm)
	dS_2 = zeros(Float64, model.num_cells) 							# net deposition per time step, settling velocity method
	num_seds = length(sediment)
	dep = zeros(Float64, num_seds) 							# sediment height change, per sediment type	
	for icell = 1:model.num_cells
		if icell == 1
			for i = 1:num_seds
				C = ((0.5*J0_SL[i] + 0.5*J_SL[icell, i]) + (0.5*J0_BL[i] + 0.5*J_BL[icell, i])  +  SourceTerm(icell, i, t, source)) / Q
				flux = C * sediment[i].vs
				dep[i] = model.dt * flux / rho_b
			end
		else
			for i = 1:num_seds
				C = ((0.5*J_SL[icell-1, i] + 0.5*J_SL[icell, i]) + (0.5*J_BL[icell-1, i] + 0.5*J_BL[icell, i])  +  SourceTerm(icell, i, t, source)) / Q
				flux = C * sediment[i].vs
				dep[i] = model.dt * flux / rho_b
			end			
		end
		dS_2[icell] = sum(dep)
	end
	return dS_2
end

function RandNorm(avg::Float64, stdev::Float64, lower_bound::Float64, upper_bound::Float64)
	# clipped,normally-distributed random number generation
	r = 0.
	accept = 0
	while accept == 0
		r = avg + stdev*randn(1)[1]
		if (r >= lower_bound) && (r <= upper_bound)
			accept = 1
		end
	end
	return r
end

function WaterDepth(model::Model, strat::Array{Strat, 1})
	# return updated array of water depths (i.e., water-column thicknesses) across all model cells
	w_depth = zeros(Float64, model.num_cells)
	for icell = 1:model.num_cells
		w_depth[icell] = model.elev_w - strat[icell].dz*(length(strat[icell].lith_frac)-1) + strat[icell].active_dz	
	end
	return w_depth
end


### output functions

function WriteProfile(fname::AbstractString, strat::Strat, sediment::Array{Sediment, 1})
	# write sediment class volume proportions, by layer, to output file
	csvfile = open(fname, "w")
	line_out = "z"
	for i = 1:length(sediment)
		line_out = line_out * "," * sediment[i].name
	end
	println(csvfile,line_out)	
	num_layers = length(strat.lith_frac)
	for i = 1:num_layers
		z = (num_layers - (i+0.5))*strat.dz
		line_out = string(z)
		for j = 1:length(sediment)
			line_out = line_out * "," * string(strat.lith_frac[i][j])
		end
		println(csvfile,line_out)
	end
	close(csvfile)
end

function WriteThickSeries(t::Float64, Q::Float64, strat::Array{Strat, 1}, init_flag::Bool)
    # write select variables to time series file
	fname = "sed_thick.csv"
	num_cells = length(strat)
	if init_flag == true
		# create time series output file and header row
		csvfile = open(fname,"w")
		line_out = "t" * "," * "Q"
		for i = 1:num_cells
			line_out = line_out * "," * "cell_" * string(i)
		end
		println(csvfile,line_out)	
	else
		csvfile = open(fname,"a")
	end
	# update ...
	line_out = string(t) * "," * string(Q)
	for i = 1:num_cells
		thick = (length(strat[i].lith_frac)-1)*strat[i].dz + strat[i].active_dz
		line_out = line_out * "," * string(thick)	
	end
	println(csvfile,line_out)
	close(csvfile)		
end

function WriteDepSeries(t::Float64, Q::Float64, dS_1::Array{Float64}, dS_2::Array{Float64}, init_flag::Bool)		
	# write deposition thicknesses, per time step, to output file
	fname = "dep_compare.csv"
	num_cells = length(dS_1)
	if init_flag == true
		# create header row
		csvfile = open(fname,"w")
		line_out = "t" * "," * "Q"
		for j = 1:2
			for i = 1:num_cells
				line_out = line_out * "," * "dS_" * string(j) * "_cell_" * string(i)
			end
		end
		println(csvfile,line_out)	
	else
		csvfile = open(fname,"a")
	end
	# update ...
	line_out = string(t) * "," * string(Q)
	# depositional calc method #1 (employed): differences in sediment carrying capacity
	for i = 1:num_cells
		line_out = line_out * "," * string(dS_1[i])	
	end
	# depositional calc method #2: settling velocity X concentration
	for i = 1:num_cells
		line_out = line_out * "," * string(dS_2[i])	
	end
	println(csvfile,line_out)
	close(csvfile)		
end

function WriteSourceSeries(t::Float64, strat::Array{Strat, 1}, source::Source, cell::Array{Cell,1}, init_flag::Bool)
    # write select variables to time series file
	fname = "external_residual.csv"
	num_cells = length(strat)
	if init_flag == true
		# create time series output file and header row
		csvfile = open(fname,"w")
		line_out = "t"
		for sed in source.sed_names
			line_out = line_out * "," * sed * "_input"		
			for icell = 1:num_cells
				line_out = line_out * "," * sed * "_cell_" * string(icell)
			end
		end
		println(csvfile,line_out)	
	else
		csvfile = open(fname,"a")
	end
	
	# update ...
	line_out = string(t)
	for ised in source.sed_indices
		line_out = line_out * "," * string(source.load[ised])
		for icell = 1:num_cells
			sum_mass = strat[icell].lith_frac[1][ised] * strat[icell].active_dz
			for islice = 2:length(strat[icell].lith_frac)
				sum_mass += strat[icell].lith_frac[islice][ised] * strat[icell].dz
			end
			sum_mass *= cell[icell].area * rho_b
			line_out = line_out * "," * string(sum_mass)
		end
	end
	
	println(csvfile,line_out)
	close(csvfile)		
end

function WriteSpoils(t::Float64, tot_spoils::Float64, dredge_count::Int64)		
	# write dredging spoils uncorrected volume estimates, per dredge evenet, to output file
	fname = "spoils.csv"
	if dredge_count == 1
		# create header row
		csvfile = open(fname,"w")
		line_out = "t" * "," * "spoils_vol"
		println(csvfile, line_out)	
	else
		csvfile = open(fname,"a")
	end
	# update ...
	line_out = string(t) * "," * "," * string(tot_spoils)	
	println(csvfile,line_out)
	close(csvfile)		
end

function WriteSedFluxes(t::Float64, J_influx::Array{Float64, 1}, J_outflux::Array{Float64, 1}, sediment::Array{Sediment, 1}, init_flag::Bool)
	# update cumulative sediment influxes and outflows, by sediment class
	fname = "sed_fluxes.csv"
	
	if init_flag == true
		# create time series output file and header row
		csvfile = open(fname,"w")
		line_out = "t"
		for i = 1:length(sediment)
			line_out = line_out * "," * sediment[i].name * "_in"
		end
		for i = 1:length(sediment)
			line_out = line_out * "," * sediment[i].name* "_out"
		end		
		println(csvfile,line_out)	
	else
		csvfile = open(fname,"a")
	end
	
	# update ...
	line_out = string(t)
	for i = 1:length(sediment)
		line_out = line_out * "," * string(J_influx[i])
	end	
	for i = 1:length(sediment)
		line_out = line_out * "," * string(J_outflux[i])
	end	
	println(csvfile,line_out)
	close(csvfile)	
	
end


### main script

function SedTran()

	println("Initializing sediment-1D model ...")

	sediment = AssignSeds()
	println("... read sediment characteristics")

	model = ModelParams()
	println("... read model parameters")

	cell = CellGeom(model.num_cells)
	println("... read cell geometry")

	strat = AssignStrat(model)
	println("... read initial stratigraphy")

	source = AssignSources(sediment)
	println("... read external source terms.")

	dredge = AssignDredge()
	dredge_count = 0
	println("... read dredge history")

	# initialize flux containers
	t = 0.
	num_seds = length(sediment)
	J0_SL = zeros(Float64, num_seds) 								# incoming (upstream) suspended load fluxes
	J0_BL = zeros(Float64, num_seds)								# incoming (upstream) bed load fluxes
	dh = zeros(Float64,num_seds) 									# change in sediment thickness, per size class
	J_SL = zeros(Float64,(model.num_cells, num_seds)) 				# outgoing (current-cell) suspended load fluxes; 1st = cell index, 2nd = sediment size class
	J_BL = zeros(Float64,(model.num_cells, num_seds)) 				# outgoing (current-cell) bed load fluxes; 1st = cell index, 2nd = sediment size class
	dS_1 = zeros(Float64, model.num_cells) 							# net deposition per time step, sediment-carrying capacity method
	J_influx = zeros(Float64, num_seds)								# cumulative sediment fluxes, by size class, across upstream boundary
	J_outflux = zeros(Float64, num_seds)							# cumulative sediment fluxes, by size class, out of last cell
	
	# initialize discharge
	if model.mode_Q == 1
		Q = 10.^RandNorm(model.log_avg_Q, model.log_stdev_Q, model.log_min_Q, model.log_max_Q)
	else
		Q = 10.^model.log_avg_Q
	end

	# time-step loop ...

	while t < model.t_max

		t += model.dt
		
		# remove sediment thickness for a simplified dredging event
		dredge_match = findin(dredge.t_ex,[t]) 						# check if current time corresponds to a dredging event
		if length(dredge_match) == 1
			dredge_count += 1	
			spoils = zeros(Float64, model.num_cells) 						# depth excavated (per cell)
			w_depth = WaterDepth(model, strat)   							# compute/update water depth array  
			println("Modeling dredge event # ", dredge_count, ", time = ", t)
			for i = 1:model.num_cells				
				dig = Dig(strat[i], w_depth[i], dredge.dz[dredge_count, i])
				spoils[i] = dig * cell[i].area
			end
			WriteSpoils(t, sum(spoils), dredge_count) 						# append to dredging spoils volume file
		end	
		
		# boundary condition for first cell 
		q = Velocity(Q, model.hw, model.w_up)
		for i = 1:num_seds
			J0_SL[i] = sediment[i].upgradient * TransSL(q, model.hw, sediment[i]) * model.w_up
			J0_BL[i] = sediment[i].upgradient * TransBL(q, model.hw, sediment[i]) * model.w_up
		end
		
		# sediment fluxes per model cell
		blended_lith = Array[]
		w_depth = WaterDepth(model, strat)   							# compute/update water depth array       		
		for icell = 1:model.num_cells

			# average current sediment composition across blending depth for icell
			push!(blended_lith, Blend(strat[icell], model.mixing_depth))	
			
			# outgoing water-column fluxes
			if icell == 1
				q = 0.5*Velocity(Q, model.hw, model.w_up) + 0.5*Velocity(Q, w_depth[icell], cell[icell].w)
			else
				q = 0.5*Velocity(Q, w_depth[icell-1], cell[icell-1].w) + 0.5*Velocity(Q, w_depth[icell], cell[icell].w)
			end
			
			for i = 1:num_seds
				J_SL[icell,i] = blended_lith[icell][i] * TransSL(q, w_depth[icell], sediment[i]) * cell[icell].w
				J_BL[icell,i] = blended_lith[icell][i] * TransBL(q, w_depth[icell], sediment[i]) * cell[icell].w
			end	
			
		end
		
		# track cumulative sediment fluxes across boundaries
		for i = 1:num_seds
			J_influx[i] += model.dt * (J0_SL[i] + J0_BL[i])	
			J_outflux[i] += model.dt * (J_SL[end, i] + J_BL[end, i])			
		end
		
		# calculate implied net deposition from differences in sediment carrying capacity + source term(s)
		for icell = 1:model.num_cells
			if icell == 1 				# first cell in river sees incoming sediment fluxed from upgradient boundary condition
				for i = 1:num_seds
					dh[i] = model.dt * ((J0_SL[i]-J_SL[icell, i]) + (J0_BL[i]-J_BL[icell, i])  +  SourceTerm(icell, i, t, source)) / (rho_b * cell[icell].area)
				end
			else 						# cell is downstream from boundary condition
				for i = 1:num_seds
					dh[i] = model.dt * ((J_SL[icell-1, i]-J_SL[icell, i]) + (J_BL[icell-1, i]-J_BL[icell, i])  +  SourceTerm(icell, i, t, source)) / (rho_b * cell[icell].area)
				end	
			end
			dS_1[icell] = sum(dh)
			if dS_1[icell] > model.mixing_depth
				println("Time step is too large; dS_1 = ", dS_1[icell])
				assert(1==0) 		# throw an exception
			end
			
			# modify active layer (by mixing), or change active layer by strat layer addition or subtraction
			Deposit(strat[icell], model.mixing_depth, blended_lith[icell], dh)	
		end

		# summarize source loading
		for loc in source.index
			source.load[loc[2]] += SourceTerm(loc[1], loc[2], t, source) * model.dt
		end
		
		# write profiles to output files, as warranted
		if (t % model.t_profile == 0) || (dredge_match == 1)
			println("Writing output at t = " * string(t) * ", Q = " * string(Q))
			for icell = 1:model.num_cells
				fname = string(icell) * "_" * string(t) * "_out.csv"
				WriteProfile(fname, strat[icell], sediment)
			end
		end

		# write/append to time series to output files, as warranted
		if t % model.t_series == 0
			init_flag = (t == model.t_series)
			dS_2 = SettlingDep(J0_SL, J0_BL, J_SL, J_BL, sediment, source, model, Q, t)
			WriteDepSeries(t, Q, dS_1, dS_2, init_flag)		
			WriteThickSeries(t, Q, strat, init_flag)
			WriteSourceSeries(t, strat, source, cell, init_flag)
			WriteSedFluxes(t, J_influx, J_outflux, sediment, init_flag)
		end
		
		# select a new value for Q, if warranted by t (and mode_Q selection)
		if (model.mode_Q == 1) && (t % model.dt_Q == 0)
			Q = 10.^RandNorm(model.log_avg_Q, model.log_stdev_Q, model.log_min_Q, model.log_max_Q)
		end

	end

	println("Done.")
	
end


### run script

SedTran()
