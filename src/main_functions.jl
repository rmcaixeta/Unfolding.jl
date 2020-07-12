using ImageMorphology:thinning
using NearestNeighbors
using LightGraphs:dijkstra_shortest_paths,connected_components
using SimpleWeightedGraphs
using StatsBase
using MultivariateStats
using Optim
using LinearAlgebra
using Random

### MAIN FUNCTIONS ###

# Isomap
function _isomap(ref_coords::AbstractArray{<:Number,2};neigh_type="knn",neigh_val=15,anchor=1500)

	@assert neigh_type in ["knn","inrange"] "Invalid neighborhood type"

    idxs, dists = _isomap_neighbors(ref_coords, neigh_type, neigh_val)

    sources,destinations,weights = [Int64[],Int64[],Float64[]]

	if neigh_type=="knn"
    	for i in 1:length(idxs)
	        for j in 2:neigh_val
	            push!(sources, i)
	            push!(destinations, idxs[i][j])
	            push!(weights, dists[i][j])
	        end
        end
	elseif neigh_type=="inrange"
		for i in 1:length(idxs)
			for j in idxs[i][idxs[i].!=i]
				push!(sources, i)
				push!(destinations, j)
				push!(weights, euclidean(ref_coords[:,i],ref_coords[:,j]))
			end
		end
	end

    g = SimpleWeightedGraph(convert(Array{Int64,1}, sources), convert(Array{Int64,1}, destinations), weights)
	comps = length(connected_components(g))
	@assert comps==1 string("There are ",comps," subgroups of isolated vertices. Need to increase number of neighbors")

	use_anchors = length(idxs)>anchor
	wgt = neigh_type=="inrange" ? Weights([1/length(x) for x in idx]) : Weights([mean(x) for x in dists])
    anchor_ids = use_anchors ? sort!(sample(1:length(idxs), wgt, anchor, replace=false)) : collect(1:length(idxs))
	anchor = minimum([anchor,length(idxs)])

	ADM = zeros(Float64,(anchor,anchor))
	for i in 1:anchor
	    d = dijkstra_shortest_paths(g,anchor_ids[i]).dists
	    for j in 1:anchor
	        ADM[i,j] = d[anchor_ids[j]]
	    end
	end

    G = dmat2gram(ADM)
    F = eigen(G)
    EM = (F.vectors[:,sortperm(F.values,rev=true)])[:,[1,2]]
    sq_eigenvals = sort(F.values,rev=true)[1:2].^0.5
    AM = Diagonal(sq_eigenvals)
    anchor_coords = EM*AM
    G = nothing

	transf_ref_coords = zeros(Float64,size(ref_coords))
	transf_ref_coords[1:2,anchor_ids] .= permutedims(anchor_coords)

	if use_anchors
		# allocate another points
		other_ids = setdiff(Array(1:length(idxs)),anchor_ids)
		other_coords = zeros(Float64,(length(other_ids),2))
		M1 = permutedims(EM)
		M1[1,:] ./= sq_eigenvals[1]
		M1[2,:] ./= sq_eigenvals[2]
		M3 = zeros(Float64,anchor,1)
		mean!(M3,ADM.^2)

		for i in 1:length(other_ids)
			M2 = zeros(Float64,(anchor,1))
			d = dijkstra_shortest_paths(g,other_ids[i]).dists
			for x in 1:anchor
				M2[x,1] = d[anchor_ids[x]]^2
			end
			out = -0.5*M1*(M2-M3)
			other_coords[i,1] = out[1,1]
			other_coords[i,2] = out[2,1]
		end

		transf_ref_coords[1:2,other_ids] .= permutedims(other_coords)
	end

    return transf_ref_coords
end

# Optimization coordinates
function _opt(points_known_true::AbstractArray{<:Number,2},
	points_known_transf::AbstractArray{<:Number,2},
	points_to_transf::AbstractArray{<:Number,2};
	xyzguess=[0],opt_neigh=8)

    tree = BallTree(points_known_true)
    idxs, dists = knn(tree, points_to_transf, opt_neigh, true)

    out_coords = zeros(Float64,3,length(idxs))
    remake_known = []
    remake_transf = []

    for i in 1:length(idxs)
        locs = points_known_transf[:,idxs[i]]
        d = dists[i]
		initial_guess = points_known_transf[:,idxs[i][1]]
		if length(xyzguess)>1
			initial_guess = xyzguess[:,i]
		end

        opt = optimize(x->_mse_coords(x, locs, d), initial_guess)#, LBFGS())

        res = Optim.minimizer(opt)
        out_coords[:,i]=res

    end

   return out_coords

end

# Skeletonization
function ref_surface_from_blocks(blks::AbstractArray{<:Number,2};
	 axis=["X","Y","Z"])

	if typeof(axis)==String
		axis = [axis]
	end
	coord_ids = Dict("X"=>1, "Y"=>2, "Z"=>3)
	ref = [Float64[],Float64[],Float64[]]
	cells = zeros(Float64,3)
	for c in 1:3
		coords =  sort!(unique(blks[c,:]))
		cell_ = sort!(unique([coords[x]-coords[x-1] for x in 2:length(coords)]))
		#@assert length(cell_)==1 "Block model must have regular cell sizes"
		cells[c] = cell_[1]
	end

	#@assert axis in ["X","Y","Z"] "Invalid axis"
	for c in axis
		ax = coord_ids[c]
	    axis_coords = sort!(unique(blks[ax,:]))
		sec = setdiff([1,2,3],ax)

	    for s in axis_coords
	        section = blks[:,blks[ax,:] .== s]

	        min_i = minimum(section[sec[1],:])
	        min_j = minimum(section[sec[2],:])
	        i = section[sec[1],:] .- min_i
	        j = section[sec[2],:] .- min_j

	        i = convert(Array{Int,1}, i/cells[sec[1]]) .+ 1
	        j = convert(Array{Int,1}, j/cells[sec[2]]) .+ 1

	        img = zeros(Bool,maximum(i)+1,maximum(j)+1)

	        for x in 1:length(i)
	            img[i[x],j[x]]=true
	        end

	        img = thinning(img)

	        for x in 1:length(i)
	            if img[i[x],j[x]]==true
	                io = (i[x]-1)*cells[sec[1]]+min_i
	                jo = (j[x]-1)*cells[sec[2]]+min_j
	                push!(ref[ax], s)
	                push!(ref[sec[1]], io)
	                push!(ref[sec[2]], jo)
	            end
	        end
	    end
	end
	out_surf = permutedims(hcat(ref[1],ref[2],ref[3]))

	if length(axis)==1
		return out_surf
	else
		return _remove_duplicates(out_surf)
	end
end

function unfold(refsurf_true_coords::AbstractArray{<:Number,2},
	input_blocks::AbstractArray{<:Number,2}, input_samps=nothing;
	neigh_type="knn",neigh_val=15,seed=1234567890)

	# random seed
	Random.seed!(seed)

	# Getting normals
	normals_neigh = neigh_type=="knn" ? neigh_val : 25
	ref_normals = _normals(refsurf_true_coords,normals_neigh)

	# Doing Isomap to get reference surface points
	ref_surf_transf = _isomap(refsurf_true_coords,neigh_type=neigh_type,neigh_val=neigh_val)
	good, bad = unfold_error(refsurf_true_coords, ref_surf_transf, 5, 5, false)

	# Get points to allocate
	xyz_finals = _xyzguess(input_blocks, refsurf_true_coords[:,good], ref_surf_transf[:,good], ref_normals[:,good])

	# Allocating blocks in chunks
	nb_chunks = 3
	shuffled_ids = shuffle(1:size(input_blocks)[2])
	ids_to_loop = collect(Iterators.partition(shuffled_ids, Int(floor(length(shuffled_ids)/nb_chunks))))

	for (i,ids) in enumerate(ids_to_loop)

		known_coords = refsurf_true_coords[:,good]
		known_tcoords = ref_surf_transf[:,good]
		ids_to_opt = ids

		if i>1

			ref_ids = [ids_to_loop[x][y] for x in 1:(i-1) for y in 1:length(ids_to_loop[x])]
			good2, bad2 = unfold_error(input_blocks[:,ref_ids], xyz_finals[:,ref_ids], 5, 5, false)
			println("good and bad: ",length(good2)," ",length(bad2))

			known_coords = hcat(known_coords,input_blocks[:,ref_ids][:,good2])
			known_tcoords = hcat(known_tcoords,xyz_finals[:,ref_ids][:,good2])
			ids_to_opt = unique(append!(Array{Int}(ids),Array{Int}(view(ref_ids,bad2))))
		end

		println("Pass ",i)

		out_transf_coords = _opt(known_coords, known_tcoords, input_blocks[:,ids_to_opt], xyzguess=xyz_finals[:,ids_to_opt])
		xyz_finals[:,ids_to_opt] .= out_transf_coords
	end

	if input_samps==nothing
		return xyz_finals
	else

		xyz_guess = _xyzguess(input_samps, input_blocks, xyz_finals)
		xyz_dh_finals = _opt(input_blocks, xyz_finals, input_samps, xyzguess=xyz_guess)

		return xyz_finals, xyz_dh_finals
	end
end
