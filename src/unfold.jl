"""
	unfold(ref_pts, input_domain, input_samps=nothing;
	isomap_search="knn",isomap_neigh=16, seed=1234567890,
	max_error=5, neighs_to_valid=16, nb_chunks=4, reftol=0.01)

Unfold the input points based on the reference points informed. Returns
a coordinate matrix with the unfolded domain points. Or a tuple of two
matrices (unfolded domain and unfolded samples points).

## Parameters:

* `ref_pts`         - coordinate matrix with the reference points
  for unfolding
* `input_domain`    - coordinate matrix with domain points for
  unfolding (blocks or mesh points)
* `input_samps`     - coordinate matrix of the sample points (optional)
* `isomap_search`   - search type to build neighbors graph for Isomap ("knn" for
  k-nearest neighbor or "inrange" for radius search).
* `isomap_neigh`    - number of neighbors (for `isomap_search`="knn") or radius distance
  (for `isomap_search`="inrange") to build neighbors graph for Isomap.
* `seed`            - seed for random values used during the process.
* `neighs_to_valid` - number of nearest neighbors to use for validations during
  the process.
* `max_error`       - the maximum accepted absolute difference of the distances
  for the closest neighbors after deformation.
* `nb_chunks`       - number of rounds of optimization.
"""
function unfold(ref_pts::AbstractMatrix, input_domain::AbstractMatrix,
	input_samps=nothing; isomap_search="knn", isomap_neigh=16, seed=1234567890,
	max_error=5, neighs_to_valid=16, nb_chunks=4, reftol=0.01)

	# random seed
	Random.seed!(seed)

	# pre-process reference points
	ref_pts = remove_duplicates(ref_pts, tol=reftol)
	resol   = get_resolution(ref_pts)

	# do landmark isomap at reference points
	ref_surf_transf = landmark_isomap(ref_pts,isomap_search=isomap_search,isomap_neigh=isomap_neigh)
	good, bad = error_ids(ref_pts, ref_surf_transf, nneigh=8, max_error=resol)
	if length(bad)>length(good)
		good, bad = error_ids(ref_pts, ref_surf_transf, nneigh=8, max_error=2*resol)
	end

	# Get initial guess
	normals_neigh = isomap_search=="knn" ? isomap_neigh : 25
	ref_normals = _normals(ref_pts,normals_neigh)
	xyz_finals = _xyzguess(input_domain, view(ref_pts,:,good), view(ref_surf_transf,:,good), view(ref_normals,:,good))

	# Allocating points in random chunks
	shuffled_ids = shuffle(1:size(input_domain,2))
	ids_to_loop = collect(Iterators.partition(shuffled_ids, Int(floor(length(shuffled_ids)/nb_chunks))))

	for (i,ids) in enumerate(ids_to_loop)

		known_coords = ref_pts[:,good]
		known_tcoords = ref_surf_transf[:,good]
		ids_to_opt = ids

		if i>1

			ref_ids = [ids_to_loop[x][y] for x in 1:(i-1) for y in 1:length(ids_to_loop[x])]
			good2, bad2 = error_ids(view(input_domain,:,ref_ids), view(xyz_finals,:,ref_ids), nneigh=neighs_to_valid, max_error=max_error)

			known_coords = hcat(known_coords,view(view(input_domain,:,ref_ids),:,good2))
			known_tcoords = hcat(known_tcoords,view(view(xyz_finals,:,ref_ids),:,good2))
			ids_to_opt = unique(append!(Array{Int}(ids),Array{Int}(view(ref_ids,bad2))))
		end

		out_transf_coords = _opt(known_coords, known_tcoords, view(input_domain,:,ids_to_opt), xyzguess=view(xyz_finals,:,ids_to_opt), opt_neigh=neighs_to_valid)
		xyz_finals[:,ids_to_opt] .= out_transf_coords
	end

	if input_samps==nothing
		xyz_finals
	else
		xyz_guess = _xyzguess(input_samps, input_domain, xyz_finals)
		xyz_dh_finals = _opt(input_domain, xyz_finals, input_samps, xyzguess=xyz_guess)
		xyz_finals, xyz_dh_finals
	end
end
