

# Optimization coordinates
function opt(known_pts, known_unf, to_unf, search, neigh, guess=nothing)
	idxs, dists = get_neighbors(known_pts, to_unf, search, neigh, true)
	out_coords  = zeros(Float64, size(to_unf))

	Threads.@threads for i in 1:length(idxs)
		locs = view(known_unf, :, idxs[i])
		initguess = isnothing(guess) ? known_unf[:,idxs[i][1]] : guess[:,i]
        opt = optimize(x->mse_coords(x, locs, dists[i]), initguess)
        res = Optim.minimizer(opt)
        out_coords[:,i] .= res
	end
    out_coords
end



# Get surface normals
function getnormals(ref_surf::AbstractMatrix, search, neigh)
	idxs, dists = get_neighbors(ref_surf, search, neigh)
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

	# get a sequential path to visit all nodes
	# might be possible to use unf indices to create a path in a simpler way in the future
	good = setdiff(1:n, ignore)
	ok = false
	n  = length(good)
	to_loop, visited, refx = Int64[], Int64[], Int64[]

	xi = good[floor(Int,n/2)]
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
	prev = Int[]
	for x in 1:n
		i, j = path[x]

		prev = union(prev,i)
		chk = intersect(prev,idxs[j])
		nx = length(chk)

		d = [evaluate(CosineDist(), ref_normals[:,k], ref_normals[:,j]) for k in chk]
		sum(d .>= 1.0) >= ceil(nx/2) && (ref_normals[:,j] .*= -1)
		prev = union(prev,j)
	end

	ref_normals, good
end

# Error function
function mse_coords(x, locations, distances)
    mse = 0.0
    for k in 1:length(distances)
        loc, d = view(locations,:,k), distances[k]
        distance_calculated = euclidean([x[1], x[2], x[3]], loc)
        mse += (distance_calculated - d)^2
    end
    mse / length(distances)
end

# Guess initial XYZ
function firstguess(to_unf, ref_pts, unf_ref, ref_normals=nothing)
	# get the closest point as initial guess
	idxs, dists = get_neighbors(ref_pts, to_unf, :knn, 1)
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
	errors(mode, coords, unf_coords; search=:knn, neigh=16, maxerr=5)

Unfolding distorts the original distances between neighbor points. This
function returns two types of results depending on `mode`.

If `mode = :ids`, this gives the IDs of the points above and below the `maxerr`
threshold as a tuple such that (ids_passed_the_tests, ids_not_passed_the_tests)

If `mode = :dists`, this function output the difference of the expected distance
for each pair analyzed in the neighborhood. Can be used as input in a boxplot to
verify distortions.

## Parameters:

* `mode`       - type of error, :ids or :dists.
* `coords`     - coordinate matrix of the points before unfolding.
* `unf_coords` - coordinate matrix of the points after unfolding.
* `search`     - type of neighborhood, :knn or :radius.
* `neigh`      - number of neighbors or radius used to make the validations.
* `maxerr`     - the maximum accepted absolute difference of the distances
  if mode = :ids.
"""
function errors(mode::Symbol, coords::AbstractMatrix, unf_coords::AbstractMatrix;
	            search=:knn, neigh=16, maxerr=5)
	idxs, dists = get_neighbors(coords, search, neigh, true)
	if mode == :ids
		error_ids(coords, unf_coords, search, neigh, maxerr)
	elseif mode == :dists
		error_dists(coords, unf_coords, search, neigh)
	else
		throw(ArgumentError("mode should be :ids or :dists"))
	end
end

function error_ids(coords, unf_coords, search=:knn, neigh=16, maxerr=5)
	idxs, dists = get_neighbors(coords, search, neigh, true)
	bad_ids = Int[]
	n = length(idxs)

	for x in 1:n
		p0 = view(unf_coords,:,x)
		tdists = [euclidean(p0, view(unf_coords,:,y)) for y in idxs[x]]
		tdists .-= dists[x]
		ids = findall(abs.(tdists) .> maxerr)

		if length(ids) > 0
			push!(bad_ids,x)
			append!(bad_ids,idxs[x][ids])
		end
	end

	bad  = unique(bad_ids)
	good = setdiff(collect(1:n), bad)
	good, bad
end

function error_dists(coords, unf_coords, search=:knn, neigh=16)
	!(coords[1] isa Float64) && (coords = Float64.(coords))
	idxs, dists = get_neighbors(coords, search, neigh, true)
	error = Float64[]

	for x in 1:length(idxs)
		p0 = view(unf_coords,:,x)
		tdists = [euclidean(p0, view(unf_coords,:,y)) for y in idxs[x]]
		tdists .-= dists[x]
		tdists .= abs.(tdists)
		mask = idxs[x] .!= x

		append!(error,tdists[mask])
	end

	error
end
