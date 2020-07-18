
function unfold(refsurf_true_coords::AbstractArray{<:Number,2},
	input_blocks::AbstractArray{<:Number,2}, input_samps=nothing;
	isomap_search="knn",isomap_neigh=15,seed=1234567890,
	max_error=5, neighs_to_valid=16, nb_chunks=4)

	# random seed
	Random.seed!(seed)

	# Doing Isomap to get reference surface points
	ref_surf_transf = landmark_isomap(refsurf_true_coords,neigh_type=isomap_search,neigh_val=isomap_neigh)
	good, bad = unfold_error_ids(refsurf_true_coords, ref_surf_transf, nneigh=8, max_error=max_error)

	# Get initial guess
	normals_neigh = isomap_search=="knn" ? isomap_neigh : 25
	ref_normals = _normals(refsurf_true_coords,normals_neigh)
	xyz_finals = _xyzguess(input_blocks, refsurf_true_coords[:,good], ref_surf_transf[:,good], ref_normals[:,good])

	# Allocating points in random chunks
	shuffled_ids = shuffle(1:size(input_blocks)[2])
	ids_to_loop = collect(Iterators.partition(shuffled_ids, Int(floor(length(shuffled_ids)/nb_chunks))))

	for (i,ids) in enumerate(ids_to_loop)

		known_coords = refsurf_true_coords[:,good]
		known_tcoords = ref_surf_transf[:,good]
		ids_to_opt = ids

		if i>1

			ref_ids = [ids_to_loop[x][y] for x in 1:(i-1) for y in 1:length(ids_to_loop[x])]
			good2, bad2 = unfold_error_ids(input_blocks[:,ref_ids], xyz_finals[:,ref_ids], nneigh=neighs_to_valid, max_error=max_error)

			known_coords = hcat(known_coords,input_blocks[:,ref_ids][:,good2])
			known_tcoords = hcat(known_tcoords,xyz_finals[:,ref_ids][:,good2])
			ids_to_opt = unique(append!(Array{Int}(ids),Array{Int}(view(ref_ids,bad2))))
		end

		out_transf_coords = _opt(known_coords, known_tcoords, input_blocks[:,ids_to_opt], xyzguess=xyz_finals[:,ids_to_opt], opt_neigh=neighs_to_valid)
		xyz_finals[:,ids_to_opt] .= out_transf_coords
	end

	if input_samps==nothing
		return xyz_finals
	else

		xyz_guess = _xyzguess(input_samps, input_blocks, xyz_finals)
		xyz_dh_finals = _opt(input_blocks, xyz_finals, input_samps, xyzguess=xyz_guess)

		return xyz_finals, xyz_dh_finals
	end
end
