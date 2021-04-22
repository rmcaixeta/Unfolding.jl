"""
	unfold(ref_pts, domain, samps=nothing; search=:knn, neigh=16,
	seed=1234567890, maxerr=5, neighs_to_valid=16, nchunks=4, reftol=0.01)

Unfold the input points based on the reference points informed. Returns
a coordinate matrix with the unfolded domain points. Or a tuple of two
matrices (unfolded domain and unfolded samples points).

## Parameters:

* `ref_pts`         - coordinate matrix with the reference points
  for unfolding
* `domain`    - coordinate matrix with domain points for
  unfolding (blocks or mesh points)
* `samps`     - coordinate matrix of the sample points (optional)
* `search`   - search type to build neighbors graph for Isomap (:knn for
  k-nearest neighbor or :radius for radius search).
* `neigh`    - number of neighbors (for `search`=:knn) or radius distance
  (for `search`=:radius) to build neighbors graph for Isomap.
* `seed`            - seed for random values used during the process.
* `neighs_to_valid` - number of nearest neighbors to use for validations during
  the process.
* `maxerr`       - the maximum accepted absolute difference of the distances
  for the closest neighbors after deformation.
* `nchunks`       - number of rounds of optimization.
"""
function unfold(ref_pts::AbstractMatrix, domain::AbstractMatrix, samps=nothing;
	isomap=:default, optim=:default)

	# read isomap and optim parameters or assign the defaults below
	ipars = (search=:knn, neigh=16, anchors=1500, seed=1234567890, reftol=0.01)
	opars = (search=:knn, neigh=16, maxerr=5, nchunks=4)

	ipars = updatepars(ipars, isomap)
	opars = updatepars(opars, optim)

	# conversions if necessary
	!(ref_pts[1] isa Float64) && (ref_pts = Float64.(ref_pts))
	!(domain[1] isa Float64) && (domain = Float64.(domain))
	!isnothing(samps) && !(samps[1] isa Float64) && (samps = Float64.(samps))

	# random seed
	Random.seed!(ipars.seed)

	# pre-process reference points
	ref_pts = remove_duplicates(ref_pts, tol=ipars.reftol)
	resol   = get_resolution(ref_pts)

	# do landmark isomap at reference points
	unf_ref = landmark_isomap(ref_pts, ipars.search, ipars.neigh, ipars.anchors)
	good, bad = error_ids(ref_pts, unf_ref, nneigh=8, maxerr=resol)
	if length(bad) > length(good)
		good, bad = error_ids(ref_pts, unf_ref, nneigh=8, maxerr=2*resol)
	end

	# get initial guess
	normals, good = getnormals(ref_pts, ipars.search, ipars.neigh, good)
	unf_dom = firstguess(domain, view(ref_pts,:,good), view(unf_ref,:,good), normals)

	# unfold points in random chunks
	ndom = size(domain,2)
	shuffled_ids = shuffle(1:ndom)
	nchunks = ceil(Int, ndom/opars.nchunks)
	chunks  = collect(Iterators.partition(shuffled_ids, nchunks))

	for (i, ids) in enumerate(chunks)
		# only reference surface as conditioner
		known_org = ref_pts[:,good]
		known_unf = unf_ref[:,good]
		ids_to_unf = ids

		# add unfolded from previous chunks
		if i > 1
			extra_ids = vcat(chunks[1:(i-1)]...)
			extra_org = view(domain,:,extra_ids)
			extra_unf = view(unf_dom,:,extra_ids)
			good_, bad_ = error_ids(extra_org, extra_unf, nneigh=opars.neigh,
			                      maxerr=opars.maxerr)

			known_org = hcat(known_org, view(extra_org, :, good_))
			known_unf = hcat(known_unf, view(extra_unf, :, good_))
			ids_to_unf = Int.(union(ids, view(extra_ids, bad_)))
		end

		# unfold points
		initguess = view(unf_dom, :, ids_to_unf)
		to_unf    = view(domain, :, ids_to_unf)
		unf_dom[:,ids_to_unf] .= opt(known_org, known_unf, to_unf, opars.search,
		                             opars.neigh, initguess)
	end

	# unfold samples if informed; otherwise, return just the domain unfolded
	if !isnothing(samps)
		initguess = firstguess(samps, domain, unf_dom)
		unf_samps = opt(domain, unf_dom, samps, opars.search, opars.neigh, initguess)
		unf_dom, unf_samps
	else
		unf_dom
	end
end


function updatepars(default, newpars)
	if newpars != :default
		for (k,v) in zip(keys(newpars), values(newpars))
			default = Setfield.setindex(default, v, k)
		end
	end
	default
end
