

# Optimization coordinates
function _opt(points_known_true::AbstractMatrix,
	points_known_transf::AbstractMatrix,
	points_to_transf::AbstractMatrix;
	xyzguess=[0],opt_neigh=16)

	points_known_true = typeof(points_known_true)<:AbstractArray{Float64} ? points_known_true : convert(Array{Float64}, points_known_true)
	points_to_transf = typeof(points_to_transf)<:AbstractArray{Float64} ? points_to_transf : convert(Array{Float64}, points_to_transf)

    tree = KDTree(points_known_true)
    idxs, dists = knn(tree, points_to_transf, opt_neigh, true)

	out_coords = zeros(Float64,size(points_to_transf))
	Threads.@threads for i in 1:length(idxs)
		locs = view(points_known_transf,:,idxs[i])
        d = dists[i]
		initial_guess = points_known_transf[:,idxs[i][1]]
		if length(xyzguess)>1
			initial_guess = xyzguess[:,i]
		end

        opt = optimize(x->_mse_coords(x, locs, d), initial_guess)#, LBFGS())
        res = Optim.minimizer(opt)
        out_coords[:,i] .= res
	end
    out_coords
end



# Get surface normals
function _normals(ref_surf::AbstractMatrix,nneigh::Number)

	ref_surf = typeof(ref_surf)<:AbstractArray{Float64} ? ref_surf : convert(Array{Float64}, ref_surf)
	tree = KDTree(ref_surf)
	idxs, dists = knn(tree, ref_surf, nneigh, true)

	ref_normals = zeros(Float64,size(ref_surf))

	last = nothing

	for i in 1:length(idxs)
		#pcapts = standardize(ZScoreTransform, ref_surf[:,idxs[i]], dims=2, scale=false)
		M = fit(PCA, view(ref_surf, :, idxs[i]), maxoutdim=3, pratio=1)
		ev = projection(M)

		if size(ev,2)>=3
			last = view(ev,:,3)
			ref_normals[:,i] .= last
		elseif size(ev,2)==2
			last = cross(view(ev,:,1),view(ev,:,2))
			ref_normals[:,i] .= last
		else # not best way to deal with 1 PC
			if last==nothing
				ref_normals[:,i] .= [0.,0.,1.]
			else
				ref_normals[:,i] .= last
			end
		end
	end

	# Get a sequential path to visit all nodes
	n = length(idxs)
	ok = false
	pre_loop = Array(1:n)
	to_loop = Int64[]
	visited = Int64[]
	refx = Int64[]

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

	ref_normals
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
function _xyzguess(coords_to_allocate::AbstractMatrix, ref_coords::AbstractMatrix,
	 ref_transf_coords::AbstractMatrix, ref_surf_normals=nothing)

	ref_coords = typeof(ref_coords)<:AbstractArray{Float64} ? ref_coords : convert(Array{Float64}, ref_coords)
	coords_to_allocate = typeof(coords_to_allocate)<:AbstractArray{Float64} ? coords_to_allocate : convert(Array{Float64}, coords_to_allocate)
	tree = KDTree(ref_coords)
	idxs, dists = knn(tree, coords_to_allocate, 1, true)

	if ref_surf_normals==nothing
		filt = [idxs[x][1] for x in 1:length(idxs)]
		return ref_transf_coords[:,filt]
	else
		xyzguess = zeros(Float64,size(coords_to_allocate))

		for x in 1:length(idxs)
			filt = idxs[x][1]
			refs_coords = view(ref_coords,:,filt)
			refn_coords = view(ref_surf_normals,:,filt)
			ref_vect = view(coords_to_allocate,:,x)-refs_coords
			dotn_val = dot(refn_coords,ref_vect)
			if dotn_val<0
				xyzguess[3,x] = -1*dists[x][1]
			else
				xyzguess[3,x] = dists[x][1]
			end

			xyzguess[1,x] = ref_transf_coords[1,filt[1]]
			xyzguess[2,x] = ref_transf_coords[2,filt[1]]
		end

		return xyzguess
	end
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

	true_coords = typeof(true_coords)<:AbstractArray{Float64} ? true_coords : convert(Array{Float64}, true_coords)
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
function error_dists(true_coords::AbstractMatrix,
	 transf_coords::AbstractMatrix; nneigh=16)

	true_coords = typeof(true_coords)<:AbstractArray{Float64} ? true_coords : convert(Array{Float64}, true_coords)
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
