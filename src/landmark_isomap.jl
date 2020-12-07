
"""
	landmark_isomap(input_coords; isomap_search="knn", isomap_neigh=16, anchors=1500)

This function runs Landmark Isomap at 3-D points and returns a coordinate matrix
of shape (3,:) with the unfolded points (Z=0 for all points).

## Parameters:

* `input_coords`  - coordinate matrix of shape (3,:) of the reference points.
* `isomap_search` - search type to build neighbors graph for Isomap ("knn" for
  k-nearest neighbor or "inrange" for radius search)
* `isomap_neigh`  - number of neighbors (for `isomap_search`="knn") or radius
  distance (for `isomap_search`="inrange") to build neighbors graph for Isomap.
* `anchors`       - number of anchors/landmark points for the dimensionality reduction.
"""
function landmark_isomap(input_coords::AbstractArray{<:Number,2};isomap_search="knn",isomap_neigh=16,anchors=1500)

	nb_points = size(input_coords)[2]
	g, anchor_ids = _make_graph_and_set_anchors(input_coords,isomap_search,isomap_neigh,anchors)
	anchors = minimum([anchors,nb_points])
	use_anchors = nb_points>anchors
	anchor_coords, other_ids, M1, M3 = _anchors_mds(nb_points,anchors,g,anchor_ids,use_anchors)

	if use_anchors

		#other_coords = @distributed (hcat) for i in 1:length(other_ids)
		# 	_triang(i,g,anchors,anchor_ids,other_ids,M1,M3)
		#end

		other_coords = zeros(Float64,(2,length(other_ids)))
		Threads.@threads for i in 1:length(other_ids)
			out = _triang(i,g,anchors,anchor_ids,other_ids,M1,M3)
			other_coords[1,i] = out[1,1]
			other_coords[2,i] = out[2,1]
		end

		transf_ref_coords = zeros(Float64,size(input_coords))
		transf_ref_coords[1:2,anchor_ids] .= permutedims(anchor_coords)
		transf_ref_coords[1:2,other_ids] .= other_coords
		return transf_ref_coords
	else
		transf_ref_coords = zeros(Float64,size(input_coords))
		transf_ref_coords[1:2,anchor_ids] .= permutedims(anchor_coords)
		return transf_ref_coords
	end
end

# Isomap neighborhood
function _get_neighbors(ref_coords::AbstractArray{<:Number,2}, neigh_type::String, neigh_val)

	ref_coords = typeof(ref_coords)<:AbstractArray{Float64} ? ref_coords : convert(Array{Float64}, ref_coords)
	tree = BallTree(ref_coords)
	if neigh_type=="knn"
		idxs, dists = knn(tree, ref_coords, neigh_val, true)
		return idxs, dists
	else
		idxs = inrange(tree, ref_coords, neigh_val, true)
		return idxs, nothing
	end
end

function _make_graph_and_set_anchors(ref_coords::AbstractArray{<:Number,2},neigh_type,neigh_val,anchor)

	nb_points = size(ref_coords)[2]
	@assert neigh_type in ["knn","inrange"] "Invalid neighborhood type"
	idxs, dists = _get_neighbors(ref_coords, neigh_type, neigh_val)
	sources,destinations,weights = [Int64[],Int64[],Float64[]]

	if neigh_type=="knn"
		for i in 1:nb_points
			for j in 2:neigh_val
				push!(sources, i)
				push!(destinations, idxs[i][j])
				push!(weights, dists[i][j])
			end
		end
	elseif neigh_type=="inrange"
		for i in 1:nb_points
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

	use_anchors = nb_points>anchor
	wgt = neigh_type=="inrange" ? Weights([1/length(x) for x in idxs]) : Weights([mean(x) for x in dists])
	anchor_ids = use_anchors ? sort!(sample(1:nb_points, wgt, anchor, replace=false)) : collect(1:nb_points)
	return g,anchor_ids
end

function _anchors_mds(nb_points,anchor,g,anchor_ids,use_anchors)
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

	if use_anchors
		# allocate another points
		other_ids = setdiff(Array(1:nb_points),anchor_ids)
		M1 = permutedims(EM)
		M1[1,:] ./= sq_eigenvals[1]
		M1[2,:] ./= sq_eigenvals[2]
		M3 = zeros(Float64,anchor,1)
		mean!(M3,ADM.^2)
		return anchor_coords,other_ids,M1,M3
	else
		return anchor_coords,nothing,nothing,nothing
	end
end

function _triang(i,g,anchor,anchor_ids,other_ids,M1,M3)
	M2 = zeros(Float64,(anchor,1))
	d = dijkstra_shortest_paths(g,other_ids[i]).dists
	for x in 1:anchor
		M2[x,1] = d[anchor_ids[x]]^2
	end
	out = -0.5*M1*(M2-M3)
	return out
end
