"""
	unfold(ref_pts, domain, samps=nothing; search="knn", neighval=16,
	seed=1234567890, max_error=5, neighs_to_valid=16, nb_chunks=4, reftol=0.01)

Unfold the input points based on the reference points informed. Returns
a coordinate matrix with the unfolded domain points. Or a tuple of two
matrices (unfolded domain and unfolded samples points).

## Parameters:

* `ref_pts`         - coordinate matrix with the reference points
  for unfolding
* `domain`    - coordinate matrix with domain points for
  unfolding (blocks or mesh points)
* `samps`     - coordinate matrix of the sample points (optional)
* `search`   - search type to build neighbors graph for Isomap ("knn" for
  k-nearest neighbor or "inrange" for radius search).
* `neighval`    - number of neighbors (for `search`="knn") or radius distance
  (for `search`="inrange") to build neighbors graph for Isomap.
* `seed`            - seed for random values used during the process.
* `neighs_to_valid` - number of nearest neighbors to use for validations during
  the process.
* `max_error`       - the maximum accepted absolute difference of the distances
  for the closest neighbors after deformation.
* `nb_chunks`       - number of rounds of optimization.
"""
function unfold(ref_pts::AbstractMatrix, domain::AbstractMatrix,
	samps=nothing; search="knn", neighval=16, seed=1234567890,
	max_error=5, neighs_to_valid=16, nb_chunks=4, reftol=0.01)

	@assert search in ["knn","inrange"] "invalid neighborhood type"

	# random seed
	Random.seed!(seed)

	# pre-process reference points
	ref_pts = remove_duplicates(ref_pts, tol=reftol)
	resol   = get_resolution(ref_pts)

	# do landmark isomap at reference points
	unf_ref = landmark_isomap(ref_pts, search=search, neighval=neighval)
	good, bad = error_ids(ref_pts, unf_ref, nneigh=8, max_error=resol)
	if length(bad)>length(good)
		good, bad = error_ids(ref_pts, unf_ref, nneigh=8, max_error=2*resol)
	end

	# Get initial guess
	search == "inrange" && (neighval = 25)
	normals = _normals(ref_pts, neighval)
	unf_dom = _xyzguess(domain, view(ref_pts,:,good), view(unf_ref,:,good), view(normals,:,good))

	# Allocating points in random chunks
	shuffled_ids = shuffle(1:size(domain,2))
	ids_to_loop = collect(Iterators.partition(shuffled_ids, Int(floor(length(shuffled_ids)/nb_chunks))))

	for (i,ids) in enumerate(ids_to_loop)

		known_orig = ref_pts[:,good]
		known_unf = unf_ref[:,good]
		ids_to_opt = ids

		if i>1
			ref_ids = [ids_to_loop[x][y] for x in 1:(i-1) for y in 1:length(ids_to_loop[x])]
			good2, bad2 = error_ids(view(domain,:,ref_ids), view(unf_dom,:,ref_ids), nneigh=neighs_to_valid, max_error=max_error)

			known_orig = hcat(known_orig,view(view(domain,:,ref_ids),:,good2))
			known_unf = hcat(known_unf,view(view(unf_dom,:,ref_ids),:,good2))
			ids_to_opt = unique(append!(Array{Int}(ids),Array{Int}(view(ref_ids,bad2))))
		end

		unf_dom[:,ids_to_opt] .= _opt(known_orig, known_unf, view(domain,:,ids_to_opt),
		    xyzguess=view(unf_dom,:,ids_to_opt), opt_neigh=neighs_to_valid)
	end

	if samps==nothing
		unf_dom
	else
		xyz_guess = _xyzguess(samps, domain, unf_dom)
		unf_samps = _opt(domain, unf_dom, samps, xyzguess=xyz_guess)
		unf_dom, unf_samps
	end
end
