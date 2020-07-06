using ImageMorphology:thinning
using NearestNeighbors
using LightGraphs:dijkstra_shortest_paths
using SimpleWeightedGraphs
using StatsBase
using MultivariateStats
using Optim
using LinearAlgebra
using Random
#using Mmap

### MAIN FUNCTIONS ###

# Isomap
function p1_isomap(ref_coords;neighs=16,anchor=800)

    tree = BallTree(ref_coords)
    idxs, dists = knn(tree, ref_coords, neighs, true)

    sources,destinations,weights = [Int64[],Int64[],Float64[]]

    for i in 1:length(idxs)
        for j in 2:neighs
            push!(sources, i)
            push!(destinations, idxs[i][j])
            push!(weights, dists[i][j])
        end
    end

    g = SimpleWeightedGraph(convert(Array{Int64,1}, sources), convert(Array{Int64,1}, destinations), weights)

	use_anchors = length(idxs)>anchor
    anchor_ids = use_anchors ? sort!(sample(1:length(idxs), anchor, replace=false)) : collect(1:length(idxs))
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
	transf_ref_coords[1:2,anchor_ids] .= transpose(anchor_coords)

	if use_anchors
		# allocate another points
		other_ids = setdiff(Array(1:length(idxs)),anchor_ids)
		other_coords = zeros(Float64,(length(other_ids),2))
		M1 = transpose(EM)
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

		transf_ref_coords[1:2,other_ids] .= transpose(other_coords)
	end

    return transf_ref_coords
end

# Optimization coordinates
function p2_opt(points_known_true,points_known_transf,points_to_transf;xyzguess=[0],opt_neigh=8)

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

        opt = optimize(x->a2_mse_coords(x, locs, d), initial_guess)#, LBFGS())

        res = Optim.minimizer(opt)
        out_coords[:,i]=res

    end

   return out_coords

end

# Skeletonization
function ref_surface_from_blocks(blks; axis="X")

	coord_ids = Dict("X"=>1, "Y"=>2, "Z"=>3)
	@assert axis in ["X","Y","Z"] "Invalid axis"
	ax = coord_ids[axis]
    axis_coords = sort!(unique(blks[ax,:]))
	sec = setdiff([1,2,3],ax)
	cells = zeros(Float64,3)

	cellsize = sort!(unique([axis_coords[x]-axis_coords[x-1] for x in 2:length(axis_coords)]))
	#@assert length(cellsize)==1 "Block model must have regular cell sizes"
	cells[ax] = cellsize[1]

	for c in sec
		coords =  sort!(unique(blks[c,:]))
		cell_ = sort!(unique([coords[x]-coords[x-1] for x in 2:length(coords)]))
		#@assert length(cell_)==1 "Block model must have regular cell sizes"
		cells[c] = cell_[1]
	end

    ref = [Float64[],Float64[],Float64[]]

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

	out_surf = transpose(hcat(ref[1],ref[2],ref[3]))

    return out_surf
end

function unfold(refsurf_true_coords,input_blocks,input_samps=nothing)

	# Getting normals
	ref_normals = a1_normals(refsurf_true_coords)

	# Doing Isomap to get reference surface points
	ref_surf_transf = p1_isomap(refsurf_true_coords)
	good, bad = unfold_error(refsurf_true_coords, ref_surf_transf, 5, 5, false)

	# Get points to allocate
	xyz_finals = a3_xyzguess(input_blocks, refsurf_true_coords[:,good], ref_surf_transf[:,good], ref_normals[:,good])

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

		out_transf_coords = p2_opt(known_coords, known_tcoords, input_blocks[:,ids_to_opt], xyzguess=xyz_finals[:,ids_to_opt])
		xyz_finals[:,ids_to_opt] .= out_transf_coords
	end

	if input_samps==nothing
		return xyz_finals
	else

		xyz_guess = a3_xyzguess(input_samps, input_blocks, xyz_finals)
		xyz_dh_finals = p2_opt(input_blocks, xyz_finals, input_samps, xyzguess=xyz_guess)

		return xyz_finals, xyz_dh_finals
	end
end
