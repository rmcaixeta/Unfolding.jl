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

	# conversions if necessary
	!(ref_pts[1] isa Float64) && (ref_pts = Float64.(ref_pts))
	!(domain[1] isa Float64) && (domain = Float64.(domain))
	!isnothing(samps) && !(samps[1] isa Float64) && (samps = Float64.(samps))

	# random seed
	Random.seed!(seed)

	# pre-process reference points
	ref_pts = remove_duplicates(ref_pts, tol=reftol)
	resol   = get_resolution(ref_pts)

	# do landmark isomap at reference points
	unf_ref = landmark_isomap(ref_pts, search=search, neighval=neighval)
	good, bad = error_ids(ref_pts, unf_ref, nneigh=8, max_error=resol)
	if length(bad) > length(good)
		good, bad = error_ids(ref_pts, unf_ref, nneigh=8, max_error=2*resol)
	end

	# Get initial guess
	normals, good = getnormals(ref_pts, search, neighval, good)
	unf_dom = firstguess(domain, view(ref_pts,:,good), view(unf_ref,:,good), normals)

	# Allocating points in random chunks
	ndom = size(domain,2)
	shuffled_ids = shuffle(1:ndom)
	ids_to_loop = collect(Iterators.partition(shuffled_ids, ceil(Int, ndom/nb_chunks)))

	for (i, ids) in enumerate(ids_to_loop)

		known_orig = ref_pts[:,good]
		known_unf = unf_ref[:,good]
		ids_to_opt = ids

		if i>1
			ref_ids = vcat(ids_to_loop[1:(i-1)]...)
			good2, bad2 = error_ids(view(domain,:,ref_ids), view(unf_dom,:,ref_ids), nneigh=neighs_to_valid, max_error=max_error)

			known_orig = hcat(known_orig,view(view(domain,:,ref_ids),:,good2))
			known_unf  = hcat(known_unf,view(view(unf_dom,:,ref_ids),:,good2))
			ids_to_opt = Int.(union(ids, view(ref_ids, bad2)))
		end

		unf_dom[:,ids_to_opt] .= opt(known_orig, known_unf, view(domain,:,ids_to_opt),
		    guess=view(unf_dom,:,ids_to_opt), nneigh=neighs_to_valid)
	end

	# unfold samples if informed; otherwise, return just the domain unfolded
	if !isnothing(samps)
		initguess = firstguess(samps, domain, unf_dom)
		unf_samps = opt(domain, unf_dom, samps, guess=initguess)
		unf_dom, unf_samps
	else
		unf_dom
	end
end
