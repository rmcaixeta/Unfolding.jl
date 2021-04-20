
"""
	landmark_isomap(coords; search="knn", nneigh=16, anchors=1500)

This function runs Landmark Isomap at 3-D points and returns a coordinate matrix
of shape (3,:) with the unfolded points (Z=0 for all points).

## Parameters:

* `coords`  - coordinate matrix of the reference points.
* `search` - search type to build neighbors graph for Isomap ("knn" for
  k-nearest neighbor or "inrange" for radius search)
* `nneigh`  - number of neighbors (for `search`="knn") or radius
  distance (for `search`="inrange") to build neighbors graph for Isomap.
* `anchors`       - number of anchors/landmark points for the dimensionality reduction.
"""
function landmark_isomap(coords::AbstractMatrix; search="knn", nneigh=16,
	                     anchors=1500, dim=2)

	n = size(coords, 2)
	n < anchors && (anchors = n)

	g, ianchors = graph_and_anchors(coords, search, nneigh, anchors)
	ADM = Array{Float64}(undef, (anchors, anchors))
	dissmatrix!(ADM, g, ianchors)

	if n > anchors
		iothers = setdiff(1:n, ianchors)
		atcoords, M1, M3 = anchors_mds(ADM, maxoutdim)
		dim = size(atcoords, 1)
		otcoords = Array{Float64}(undef, (dim, length(iothers)))
		Threads.@threads for i in 1:length(iothers)
			otcoords[:,i] .= triangulation(D,metric,iothers[i],ianchors,M1,M3)
		end
		tcoords = Array{Float64}(undef, (dim,n))
		tcoords[:,ianchors] .= atcoords
		tcoords[:,iothers] .= otcoords
	else
		M = fit(MDS, ADM, maxoutdim=maxoutdim, distances=true)
		tcoords = transform(M)
		#println("Explained variance: $(sum(sortλ[1:maxdim])/sum(sortλ[λ .> 0]))")
	end
 	tcoords
end

# nb_points = size(coords,2)
# g, anchor_ids = _make_graph_and_set_anchors(coords,search,nneigh,anchors)
# anchors = minimum([anchors,nb_points])
# use_anchors = nb_points>anchors
# anchor_coords, other_ids, M1, M3 = _anchors_mds(nb_points,anchors,g,anchor_ids,use_anchors)
#
# if use_anchors
#
# 	other_coords = zeros(Float64,(2,length(other_ids)))
# 	Threads.@threads for i in 1:length(other_ids)
# 		out = _triang(i,g,anchors,anchor_ids,other_ids,M1,M3)
# 		other_coords[1,i] = out[1,1]
# 		other_coords[2,i] = out[2,1]
# 	end
#
# 	transf_ref_coords = zeros(Float64,size(coords))
# 	transf_ref_coords[1:2,anchor_ids] .= permutedims(anchor_coords)
# 	transf_ref_coords[1:2,other_ids] .= other_coords
# 	return transf_ref_coords
# else
# 	transf_ref_coords = zeros(Float64,size(coords))
# 	transf_ref_coords[1:2,anchor_ids] .= permutedims(anchor_coords)
# 	return transf_ref_coords
# end

# Isomap neighborhood
# function _get_neighbors(ref_coords::AbstractMatrix, neigh_type::String, neigh_val)
#
# 	ref_coords = typeof(ref_coords)<:AbstractArray{Float64} ? ref_coords : convert(Array{Float64}, ref_coords)
# 	tree = BallTree(ref_coords)
# 	if neigh_type=="knn"
# 		idxs, dists = knn(tree, ref_coords, neigh_val, true)
# 		return idxs, dists
# 	else
# 		idxs = inrange(tree, ref_coords, neigh_val, true)
# 		return idxs, nothing
# 	end
# end
#
function graph_and_anchors(ref_coords::AbstractMatrix,neigh_type,neigh_val,anchor)

	nb_points = size(ref_coords,2)
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
				push!(weights, euclidean(view(ref_coords,:,i),view(ref_coords,:,j)))
			end
		end
	end

	g = SimpleWeightedGraph(convert(Array{Int64,1}, sources), convert(Array{Int64,1}, destinations), weights)
	comps = length(connected_components(g))
	@assert comps==1 string("There are ",comps," subgroups of isolated vertices. Need to increase number of neighbors")

	use_anchors = nb_points>anchor
	wgt = neigh_type=="inrange" ? Weights([1/length(x) for x in idxs]) : Weights([mean(x) for x in dists])
	anchor_ids = use_anchors ? sort!(sample(1:nb_points, wgt, anchor, replace=false)) : collect(1:nb_points)
	g,anchor_ids
end
#
# function _anchors_mds(nb_points,anchor,g,anchor_ids,use_anchors)
# 	ADM = zeros(Float64,(anchor,anchor))
# 	for i in 1:anchor
# 		d = dijkstra_shortest_paths(g,anchor_ids[i]).dists
# 		for j in 1:anchor
# 			ADM[i,j] = d[anchor_ids[j]]
# 		end
# 	end
#
# 	G = dmat2gram(ADM)
#     F = eigen(Symmetric(G))
#     EM = (F.vectors[:,sortperm(F.values,rev=true)])[:,[1,2]]
#     sq_eigenvals = sort(F.values,rev=true)[1:2].^0.5
#     AM = Diagonal(sq_eigenvals)
#     anchor_coords = EM*AM
#
# 	if use_anchors
# 		# allocate another points
# 		other_ids = setdiff(Array(1:nb_points),anchor_ids)
# 		M1 = permutedims(EM)
# 		M1 ./= reshape(sq_eigenvals,2,1)
# 		M3 = zeros(Float64,anchor,1)
# 		mean!(M3,ADM.^2)
# 		return anchor_coords,other_ids,M1,M3
# 	else
# 		return anchor_coords,nothing,nothing,nothing
# 	end
# end
#
# function _triang(i,g,anchor,anchor_ids,other_ids,M1,M3)
# 	M2 = zeros(Float64,(anchor,1))
# 	d = dijkstra_shortest_paths(g,other_ids[i]).dists
# 	for x in 1:anchor
# 		M2[x,1] = d[anchor_ids[x]]^2
# 	end
# 	out = -0.5*M1*(M2-M3)
# 	out
# end



function dissmatrix!(ADM, g, iax::Vector{Int})
	n = size(ADM,1)
	for (i,ia) in enumerate(iax)
		dcols = dijkstra_shortest_paths(g, ia).dists[iax[i:n]]
		for (d,j) in zip(dcols,i:n)
			ADM[i,j] = ADM[j,i] = d
		end
	end
end

function anchors_mds(ADM, maxoutdim)
	nx = size(ADM,1)
	G = dmat2gram(ADM)
    F = eigen(Symmetric(G))
	λ = F.values
	maxdim = minimum([maxoutdim, nx-1, sum(λ .> 0)])
	sorti = sortperm(F.values,rev=true)
	sortλ = λ[sorti]
    EM = (F.vectors[:,sorti])[:,1:maxdim]
	#println("Explained variance: $(sum(sortλ[1:maxdim])/sum(sortλ[sortλ .> 0]))")
    sq_eigenvals = sortλ[1:maxdim].^0.5
    AM = Diagonal(sq_eigenvals)
    atcoords = permutedims(EM*AM)

	# to use later for another points allocations
	M1 = permutedims(EM)
	M1 ./= reshape(sq_eigenvals,maxdim,1)
	M3 = Array{Float64}(undef,nx) #Array{Float64}(undef,(nx,1))
	mean!(M3,ADM.^2)

	atcoords, M1, M3
end

function triangulation(D::LocalGeoData, metric::Type{<:LocalMetric},i,j,M1,M3)
	M2 = colwise(D, metric, i, j) .^ 2
	-0.5*M1*(M2-M3)
end
