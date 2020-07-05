using NearestNeighbors
using StatsBase
using MultivariateStats
using Distances
using DelimitedFiles
using LinearAlgebra
using WriteVTK
using StatsPlots


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
function a3_xyzguess(coords_to_allocate, ref_coords, ref_transf_coords, ref_surf_normals=nothing)

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
function unfold_error(true_coords, transf_coords, nneigh, max_error=5, outfunc=true; plotname="error")

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
		fig = boxplot(["Unfolding error"], error)
		savefig(plotname)
		return error
	else
		bad = unique(bad_ids)
		good = setdiff(Array(1:size(true_coords)[2]),bad)
		return good,bad
	end
end

function data_to_csv(coords,outname,colnames)

	out = length(size(coords))==1 ? coords : coords'
	open(string(outname,".csv"); write=true) do f
		write(f, string(join(colnames,","),"\n"))
		writedlm(f, out, ',')
	end
end

function data_to_vtk(coords,outname,extra_props=nothing)
	verts = [MeshCell( VTKCellTypes.VTK_VERTEX, [i]) for i in 1:size(coords)[2] ]
	outfiles = vtk_grid(outname,coords,(verts)) do vtk
		if extra_props!=nothing
			for (name,prop) in extra_props
				vtk[name] = prop
			end
		end
	end

end
