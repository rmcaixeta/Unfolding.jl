
"""
	landmark_isomap(coords; search="knn", neighval=16, anchors=1500)

This function runs Landmark Isomap at 3-D points and returns a coordinate matrix
with the unfolded points.

## Parameters:

* `coords`  - coordinate matrix of the reference points.
* `search`  - search type to build neighbors graph for Isomap ("knn" for
  k-nearest neighbor or "inrange" for radius search)
* `neighval`  - number of neighbors (for `search`="knn") or radius
  distance (for `search`="inrange") to build neighbors graph for Isomap.
* `anchors` - number of anchors/landmark points for the dimensionality reduction.
"""
function landmark_isomap(coords::AbstractMatrix; search="knn", neighval=16,
	                     anchors=1500, dim=2)
	n = size(coords, 2)
	n < anchors && (anchors = n)

	g, ianchors = graph_and_anchors(coords, search, neighval, anchors)
	ADM = Array{Float64}(undef, (anchors, anchors))
	dissmatrix!(ADM, g, ianchors)

	if n > anchors
		iothers = setdiff(1:n, ianchors)
		atcoords, M1, M3 = anchors_mds(ADM, dim)

		otcoords = Array{Float64}(undef, (3, length(iothers)))
		Threads.@threads for i in 1:length(iothers)
			otcoords[:,i] .= triangulation(g, iothers[i], ianchors, M1, M3)
		end
		tcoords = Array{Float64}(undef, (3,n))
		tcoords[:,ianchors] .= atcoords
		tcoords[:,iothers] .= otcoords
	else
		M = fit(MDS, ADM, maxoutdim=dim, distances=true)
		tcoords = dim == 2 ? vcat(transform(M),zeros(1, n)) : transform(M)
	end
 	tcoords
end

function graph_and_anchors(ref_coords::AbstractMatrix, nhood, neigh_val, nanchors)
	n = size(ref_coords,2)
	idxs, dists = get_neighbors(ref_coords, nhood, neigh_val)
	src, dest, dwgt = [Int64[], Int64[], Float64[]]

	if nhood=="knn"
		for i in 1:n
			for j in 2:neigh_val
				push!(src, i)
				push!(dest, idxs[i][j])
				push!(dwgt, dists[i][j])
			end
		end
	elseif nhood=="inrange"
		for i in 1:n
			for j in idxs[i][idxs[i].!=i]
				push!(src, i)
				push!(dest, j)
				push!(dwgt, euclidean(view(ref_coords,:,i),view(ref_coords,:,j)))
			end
		end
	end

	g = SimpleWeightedGraph(Int.(src), Int.(dest), dwgt)
	comps = length(connected_components(g))
	@assert comps==1 string("$comps subgroups of isolated vertices, need to increase number of neighbors")

	wgt = nhood=="inrange" ? [1/length(x) for x in idxs] : [mean(x) for x in dists]
	ianchors = n > nanchors ? sort!(sample(1:n, Weights(wgt), nanchors, replace=false)) : collect(1:n)
	g, ianchors
end

function dissmatrix!(ADM, g, iax::Vector{Int})
	n = size(ADM,1)
	for (i,ia) in enumerate(iax)
		dcols = dijkstra_shortest_paths(g, ia).dists[iax[i:n]]
		for (d,j) in zip(dcols,i:n)
			ADM[i,j] = ADM[j,i] = d
		end
	end
end

function anchors_mds(ADM, dims)
	nx = size(ADM,1)
	G = dmat2gram(ADM)
    F = eigen(Symmetric(G))
	λ = F.values
	sorti = sortperm(F.values,rev=true)
	sortλ = λ[sorti]
    EM = (F.vectors[:,sorti])[:,1:dims]
    sq_eigenvals = sortλ[1:dims].^0.5
    AM = Diagonal(sq_eigenvals)
    atcoords = dims == 2 ? vcat(permutedims(EM*AM),zeros(1,nx)) : permutedims(EM*AM)

	# to use later for another points allocations
	M1 = permutedims(EM)
	M1 ./= reshape(sq_eigenvals,dims,1)
	M3 = Array{Float64}(undef,nx)
	mean!(M3,ADM.^2)

	atcoords, M1, M3
end

function triangulation(g::SimpleWeightedGraph, i, j, M1, M3)
	M2 = dijkstra_shortest_paths(g, i).dists[j] .^ 2
	size(M1,1) == 2 ? vcat(-0.5*M1*(M2-M3),0) : -0.5*M1*(M2-M3)
end
