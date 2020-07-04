using CSV
using ImageMorphology:thinning
using NearestNeighbors
using LightGraphs:dijkstra_shortest_paths
using SimpleWeightedGraphs
using StatsBase
using MultivariateStats
using Optim
using Distances
using DelimitedFiles
using LinearAlgebra
using Random


### AUXILIAR FUNCTIONS ###

# Get surface normals
function a1_normals(ref_surf)
	tree = BallTree(ref_surf)
	idxs, dists = knn(tree, ref_surf, 15, true)

	ref_normals = zeros(Float64,size(ref_surf))
	bad_id = Int64[]

	for i in 1:length(idxs)

		M = fit(PCA, ref_surf[:,idxs[i]], maxoutdim=3, pratio=1)
		ev = projection(M)

		if size(ev)[2]>=3
			ref_normals[:,i] .= ev[:,3]
		elseif size(ev)[2]==2
			ref_normals[:,i] .= cross(ev[:,1],ev[:,2])
		#else
		#	ref_normals[:,i] .= [-99,-99,-99]
		#	append!(bad_id,i)
		end
	end

	# Iterate through all normals and check for consistency
	ok = false
	pre_loop = Array(1:length(idxs))
	to_loop = Int64[]
	visited = Int64[]

	xi = pre_loop[floor(Int,length(pre_loop)/2)]
	append!(visited,xi)
	to_loop = union(to_loop, idxs[xi])
	extra_loop = setdiff(to_loop,visited)
	while ok==false
		xi = extra_loop[1]
		append!(visited,xi)
		to_loop = union(to_loop, idxs[xi])
		extra_loop = setdiff(to_loop,visited)
		if length(extra_loop)==0
			ok=true
		end
	end

	missing_pts = setdiff(pre_loop,visited)
	if length(missing_pts)>0
		println("bad",length(missing_pts))
		#deal with not visited points
	end

	for i in unique(to_loop)
		X = ref_normals[:,idxs[i]]
		Y = fill(ref_normals[1,i],size(X))
		Y[2,:] .= ref_normals[2,i]
		Y[3,:] .= ref_normals[3,i]
		CD = colwise(CosineDist(), X, Y)
		ref_normals[:,idxs[i][CD .> 1]] .*= -1.0
	end

	return ref_normals
end

# Error function
function a2_mse_coords(x, locations, distances)
    mse = 0.0
    for k in 1:length(distances)
        loc, d = [locations[:,k],distances[k]]
        distance_calculated = euclidean([x[1],x[2],x[3]],loc)
        mse += (distance_calculated - d)^2
    end
    return mse / length(distances)
end

# Guess initial XYZ
function a3_xyzguess(ref_surf_coords, ref_surf_tcoords, ref_surf_normals, coords_to_allocate)

    tree = BallTree(ref_surf_coords)
	idxs, dists = knn(tree, coords_to_allocate, 1, true)

	xyzguess = zeros(Float64,size(coords_to_allocate))

	for x in 1:length(idxs)
		filt = idxs[x][1]
		refs_coords = ref_surf_coords[:,filt]
		refn_coords = ref_surf_normals[:,filt]
		ref_vect = coords_to_allocate[:,x]-refs_coords
		dotn_val = dot(refn_coords,ref_vect)
		if dotn_val<0
			xyzguess[3,x] = -1*dists[x][1]
		else
			xyzguess[3,x] = dists[x][1]
		end

		xyzguess[1,x] = ref_surf_tcoords[1,filt[1]]
		xyzguess[2,x] = ref_surf_tcoords[2,filt[1]]
	end

	return xyzguess
end

# Check error
function unfold_error(true_coords, transf_coords, nneigh, max_error;main_return=true)

    tree = BallTree(true_coords)
	idxs, dists = knn(tree, true_coords, nneigh, true)

	tdists = zeros(Float64,size(dists[1]))
	bad_ids = Int[]
	error = Float64[]

	for x in 1:length(idxs)
		for y in length(idxs[x])
			tdists[y] = euclidean(transf_coords[:,x],transf_coords[:,idxs[x][y]])
		end

		tdists .-= dists[x]
		tdists .= abs.(tdists)
		check = findall(tdists .> max_error)

		if main_return==false
			append!(error,tdists[2:end])
		elseif length(check)>0
			push!(bad_ids,x)
			append!(bad_ids,idxs[x][check])
		end
	end

	if main_return==false
		return error
	else
		bad = unique(bad_ids)
		good = setdiff(Array(1:size(true_coords)[2]),bad)
		return good,bad
	end
end
