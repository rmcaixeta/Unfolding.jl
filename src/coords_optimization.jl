

# Optimization coordinates
function _opt(points_known_true::AbstractArray{<:Number,2},
	points_known_transf::AbstractArray{<:Number,2},
	points_to_transf::AbstractArray{<:Number,2};
	xyzguess=[0],opt_neigh=16)

	points_known_true = typeof(points_known_true)<:AbstractArray{Float64} ? points_known_true : convert(Array{Float64}, points_known_true)
	points_to_transf = typeof(points_to_transf)<:AbstractArray{Float64} ? points_to_transf : convert(Array{Float64}, points_to_transf)

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



# Get surface normals
function _normals(ref_surf::AbstractArray{<:Number,2},nneigh::Number)

	ref_surf = typeof(ref_surf)<:AbstractArray{Float64} ? ref_surf : convert(Array{Float64}, ref_surf)
	tree = BallTree(ref_surf)
	idxs, dists = knn(tree, ref_surf, nneigh, true)

	ref_normals = zeros(Float64,size(ref_surf))

	last = nothing

	for i in 1:length(idxs)
		pcapts = standardize(ZScoreTransform, ref_surf[:,idxs[i]], dims=2, scale=false)
		M = fit(PCA, pcapts, maxoutdim=3, pratio=1)
		ev = projection(M)

		if size(ev)[2]>=3
			last = ev[:,3]
			ref_normals[:,i] .= ev[:,3]
		elseif size(ev)[2]==2
			last = cross(ev[:,1],ev[:,2])
			ref_normals[:,i] .= cross(ev[:,1],ev[:,2])
		else # not best way to deal with 1 PC
			if last==nothing
				ref_normals[:,i] .= [0.,0.,1.]
			else
				ref_normals[:,i] .= last
			end
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
function _mse_coords(x::Array{Float64,1}, locations::AbstractArray{<:Number,2}, distances::Array{Float64,1})
    mse = 0.0
    for k in 1:length(distances)
        loc, d = [locations[:,k],distances[k]]
        distance_calculated = euclidean([x[1],x[2],x[3]],loc)
        mse += (distance_calculated - d)^2
    end
    return mse / length(distances)
end

# Guess initial XYZ
function _xyzguess(coords_to_allocate::AbstractArray{<:Number,2}, ref_coords::AbstractArray{<:Number,2},
	 ref_transf_coords::AbstractArray{<:Number,2}, ref_surf_normals=nothing)

	ref_coords = typeof(ref_coords)<:AbstractArray{Float64} ? ref_coords : convert(Array{Float64}, ref_coords)
	coords_to_allocate = typeof(coords_to_allocate)<:AbstractArray{Float64} ? coords_to_allocate : convert(Array{Float64}, coords_to_allocate)
	tree = BallTree(ref_coords)
	idxs, dists = knn(tree, coords_to_allocate, 1, true)

	if ref_surf_normals==nothing
		filt = [idxs[x][1] for x in 1:length(idxs)]
		return ref_transf_coords[:,filt]
	else
		xyzguess = zeros(Float64,size(coords_to_allocate))

		for x in 1:length(idxs)
			filt = idxs[x][1]
			refs_coords = ref_coords[:,filt]
			refn_coords = ref_surf_normals[:,filt]
			ref_vect = coords_to_allocate[:,x]-refs_coords
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


# Check error
function unfold_error_ids(true_coords::AbstractArray{<:Number,2},
	 transf_coords::AbstractArray{<:Number,2};
	 nneigh=16, max_error=5)

	true_coords = typeof(true_coords)<:AbstractArray{Float64} ? true_coords : convert(Array{Float64}, true_coords)
    tree = BallTree(true_coords)
	idxs, dists = knn(tree, true_coords, nneigh, true)

	bad_ids = Int[]

	for x in 1:length(idxs)
		tdists = zeros(Float64,size(idxs[x]))
		for y in length(idxs[x])
			tdists[y] = euclidean(transf_coords[:,x],transf_coords[:,idxs[x][y]])
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
	good = setdiff(Array(1:size(true_coords)[2]),bad)
	return good,bad

end


function unfold_error_dists(true_coords::AbstractArray{<:Number,2},
	 transf_coords::AbstractArray{<:Number,2};
	 nneigh=16, plotname=nothing)

	true_coords = typeof(true_coords)<:AbstractArray{Float64} ? true_coords : convert(Array{Float64}, true_coords)
    tree = BallTree(true_coords)
	idxs, dists = knn(tree, true_coords, nneigh, true)

	error = Float64[]

	for x in 1:length(idxs)
		tdists = zeros(Float64,size(idxs[x]))
		for y in length(idxs[x])
			tdists[y] = euclidean(transf_coords[:,x],transf_coords[:,idxs[x][y]])
		end

		tdists .-= dists[x]
		tdists .= abs.(tdists)

		append!(error,tdists[2:end])
	end

	if plotname!=nothing
		_make_boxplot(error,plotname)
	end

	return error
end
