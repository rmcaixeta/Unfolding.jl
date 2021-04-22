

# Optimization coordinates
function opt(known_pts::AbstractMatrix, known_unf::AbstractMatrix,
	          to_unf::AbstractMatrix; guess=[0], nneigh=16)

	idxs, dists = get_neighbors(known_pts, to_unf, "knn", nneigh)
	out_coords  = zeros(Float64, size(to_unf))

	Threads.@threads for i in 1:length(idxs)
		locs = view(known_unf, :, idxs[i])
		initguess = length(guess) > 1 ? guess[:,i] : known_unf[:,idxs[i][1]]
        opt = optimize(x->_mse_coords(x, locs, dists[i]), initguess)
        res = Optim.minimizer(opt)
        out_coords[:,i] .= res
	end
    out_coords
end



# Get surface normals
function getnormals(ref::AbstractMatrix, search::String, nneigh::Int, good)
	ref_surf = view(ref, :, good)
	idxs, dists = get_neighbors(ref_surf, search, nneigh)
	n = length(idxs)

	ref_normals = zeros(Float64, size(ref_surf))
	ignore = Int[]

	for i in 1:n
		M  = fit(PCA, view(ref_surf, :, idxs[i]), maxoutdim=3, pratio=1)
		ev = projection(M)
		ndims = size(ev,2)

		if ndims <= 1
			push!(ignore,i)
		elseif ndims == 2
			ref_normals[:,i] .= cross(view(ev,:,1),view(ev,:,2))
		else
			ref_normals[:,i] .= view(ev,:,3)
		end
	end

	# Get a sequential path to visit all nodes
	pre_loop = setdiff(1:n, ignore)
	ok = false
	n  = length(pre_loop)
	to_loop, visited, refx = Int64[], Int64[], Int64[]

	xi = pre_loop[floor(Int,n/2)]
	append!(visited,xi)
	append!(to_loop, idxs[xi])
	append!(refx, [xi for x in idxs[xi]])
	extra_loop = setdiff(to_loop,visited)
	while ok==false
		xj = extra_loop[1]
		append!(visited,xj)
		append!(to_loop, idxs[xj])
		append!(refx, [xj for x in idxs[xj]])
		extra_loop = setdiff(to_loop,visited)
		if length(extra_loop)==0
			ok=true
		end
	end

	ids_unique = unique(i -> to_loop[i], 1:length(to_loop))
	path = [p=>q for (p,q) in zip(refx[ids_unique],to_loop[ids_unique])]

	# Iterate through all normals and check for consistency
	for x in 1:n
		i, j = path[x]
		d = evaluate(CosineDist(), ref_normals[:,i], ref_normals[:,j])
		d >= 1.5 && (ref_normals[:,j] .*= -1)
	end

	ref_normals[:,pre_loop], good[pre_loop]
end

# Error function
function _mse_coords(x::AbstractArray{Float64,1}, locations::AbstractMatrix, distances::AbstractArray{Float64,1})
    mse = 0.0
    for k in 1:length(distances)
        loc, d = [view(locations,:,k),distances[k]]
        distance_calculated = euclidean([x[1],x[2],x[3]],loc)
        mse += (distance_calculated - d)^2
    end
    mse / length(distances)
end

# Guess initial XYZ
function firstguess(to_unf::AbstractMatrix, ref_pts::AbstractMatrix,
    unf_ref::AbstractMatrix, ref_normals=nothing)
	# get the closest point as initial guess
	idxs, dists = get_neighbors(ref_pts, to_unf, "knn", 1)
	n   = length(idxs)
	ids = [idxs[x][1] for x in 1:n]
	guess = unf_ref[:,ids]

	# if reference surface normals are informed, define z as the signed distance
	if !isnothing(ref_normals)
		for (x, id) in enumerate(ids)
			refcoords = view(ref_pts, :, id)
			refnormal = view(ref_normals, :, id)
			ref_vect = view(to_unf, :, x) - refcoords
			dotn_val = dot(refnormal, ref_vect)
			guess[3,x] = dotn_val < 0 ? -1*dists[x][1] : dists[x][1]
		end
	end

	guess
end


"""
	error_ids(true_coords, transf_coords; nneigh=16, max_error=5)

Unfolding distorts the original distances between neighbor points. This
function give the IDs of the points above and below the `max_error` threshold.
Returns a tuple of two arrays. The first array with the ID of points that passed
the tests. The second array with the ID of points that failed during the tests.

## Parameters:

* `true_coords`   - coordinate matrix of the points before unfolding.
* `transf_coords` - coordinate matrix of the points after unfolding.
* `nneigh`        - number of nearest neighbors used to make the validations.
* `max_error`     - the maximum accepted absolute difference of the distances
  for the closest neighbors after unfolding.
"""
function error_ids(true_coords::AbstractMatrix,
	 transf_coords::AbstractMatrix;
	 nneigh=16, max_error=5)

	!(true_coords[1] isa Float64) && (true_coords = Float64.(true_coords))
    tree = KDTree(true_coords)
	idxs, dists = knn(tree, true_coords, nneigh, true)

	bad_ids = Int[]

	for x in 1:length(idxs)
		tdists = zeros(Float64,size(idxs[x]))
		for y in length(idxs[x])
			tdists[y] = euclidean(view(transf_coords,:,x),view(transf_coords,:,idxs[x][y]))
		end
		tdists .-= dists[x]
		tdists .= abs.(tdists)
		ids = findall(tdists .> max_error)

		if length(ids)>0
			push!(bad_ids,x)
			append!(bad_ids,idxs[x][ids])
		end
	end

	bad = unique(bad_ids)
	good = setdiff(Array(1:size(true_coords,2)),bad)
	good,bad
end

"""
	error_dists(true_coords, transf_coords; nneigh=16)

Unfolding distorts the original distances between neighbor points. This
function output the difference of the expected distance for each pair analyzed.
Can be used as input to boxplot to verify distortions.

## Parameters:

* `true_coords`   - coordinate matrix of the points before unfolding.
* `transf_coords` - coordinate matrix of the points after unfolding.
* `nneigh`        - number of nearest neighbors used to make the validations.
"""
function error_dists(true_coords::AbstractMatrix, transf_coords::AbstractMatrix; nneigh=16)
	!(true_coords[1] isa Float64) && (true_coords = Float64.(true_coords))
    tree = KDTree(true_coords)
	idxs, dists = knn(tree, true_coords, nneigh, true)

	error = Float64[]

	for x in 1:length(idxs)
		tdists = zeros(Float64,size(idxs[x]))
		for y in length(idxs[x])
			tdists[y] = euclidean(view(transf_coords,:,x),view(transf_coords,:,idxs[x][y]))
		end

		tdists .-= dists[x]
		tdists .= abs.(tdists)

		append!(error,tdists[2:end])
	end

	error
end
