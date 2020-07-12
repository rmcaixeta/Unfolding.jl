using NearestNeighbors
using StatsBase
using MultivariateStats
using Distances
using DelimitedFiles
using LinearAlgebra
using WriteVTK
using StatsPlots
using Clustering
using Suppressor

### AUXILIAR FUNCTIONS ###

# Adjust matrix format
function julia_matrix(x;columns=nothing)

	x = columns==nothing ? x : x[:,[Symbol(v) for v in columns]]
	shape = size(x)
	@assert (shape[1]==3 || shape[2]==3) "Three coordinate matrix should be informed"

	x = typeof(x)<:Matrix{<:Number} ? x : Matrix{Float64}(x)
	x = (shape[2]==3 && shape[1]!=3) ? permutedims(x) : x
	return x
end

# Average duplicate points
function _remove_duplicates(coords::AbstractArray{<:Number,2},tol=0.5)
	groups = [union(x.core_indices,x.boundary_indices) for x in dbscan(coords, tol)]
	outcoords = zeros(Float64,3,length(groups))
	for g in 1:length(groups)
		idx = groups[g]
		if length(idx)==1
			outcoords[:,g] .= dropdims(coords[:,idx],dims=2)
		else
			outcoords[:,g] .= mean!(outcoords[:,g],coords[:,idx])
		end
	end
	return outcoords
end

# Isomap neighborhood
function _isomap_neighbors(ref_coords::AbstractArray{<:Number,2}, neigh_type::String, neigh_val)
	tree = BallTree(ref_coords)
	if neigh_type=="knn"
		idxs, dists = knn(tree, ref_coords, neigh_val, true)
		return idxs, dists
	else
		idxs = inrange(tree, ref_coords, neigh_val, true)
		return idxs, nothing
	end
end

# Get surface normals
function _normals(ref_surf::AbstractArray{<:Number,2},nneigh::Number)
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
function unfold_error(true_coords::AbstractArray{<:Number,2},
	 transf_coords::AbstractArray{<:Number,2},
	 nneigh::Int, max_error=5, outfunc=true; plotname="error")

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

		if outfunc==true
			append!(error,tdists[2:end])
		elseif length(check)>0
			push!(bad_ids,x)
			append!(bad_ids,idxs[x][check])
		end
	end

	if outfunc==true
		@suppress fig = boxplot(["Unfolding error"], error)
		savefig(plotname)
		return error
	else
		bad = unique(bad_ids)
		good = setdiff(Array(1:size(true_coords)[2]),bad)
		return good,bad
	end
end

function data_to_csv(coords::AbstractArray{<:Number,2},outname::String,colnames)

	out = length(size(coords))==1 ? coords : coords'
	open(string(outname,".csv"); write=true) do f
		write(f, string(join(colnames,","),"\n"))
		writedlm(f, out, ',')
	end
end

function data_to_vtk(coords::AbstractArray{<:Number,2},outname::String,extra_props=nothing)
	verts = [MeshCell( VTKCellTypes.VTK_VERTEX, [i]) for i in 1:size(coords)[2] ]
	outfiles = vtk_grid(outname,coords,(verts)) do vtk
		if extra_props!=nothing
			for (name,prop) in extra_props
				vtk[name] = prop
			end
		end
	end

end
