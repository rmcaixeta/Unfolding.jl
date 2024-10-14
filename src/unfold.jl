"""
	unfold(ref, to_unf, samps=nothing; isomap=:default, optim=:default)

Unfold the input points `to_unf` based on the reference points informed in `ref`.
Returns a coordinate matrix with the unfolded points. Or a tuple of matrices if
`to_unf` is a list of points.

## Parameters:

* `to_unf`  - coordinate matrix with the points that will be unfolded
* `ref`     - coordinate matrix with the reference points for unfolding
* `isomap`  - additional parameters for isomap step; either :default or a
  NamedTuple with the keys shown below
* `optim`   - additional parameters for optimization step; either :default or a
  NamedTuple with the keys shown below

### `isomap` parameters:

*Default:* isomap = (search=:knn, neigh=16, anchors=1500, reftol=0.01)
* `search`  - search type to build neighbors graph for Isomap (:knn for
  k-nearest neighbor or :radius for radius search).
* `neigh`   - number of neighbors (for search=:knn) or radius distance
  (for search=:radius) to build neighbors graph for Isomap.
* `anchors` - number of anchors for landmark isomap
* `reftol`  - distance to which reference points are considered duplicates

### `optim` parameters:

*Default:* optim = (search=:knn, neigh=50, maxerr=5, nchunks=4)
* `search`  - search type to optimize locally (:knn for k-nearest neighbor or
  :radius for radius search).
* `neigh`   - number of neighbors (for search=:knn) or radius distance
  (for search=:radius) to local optimization.
* `maxerr`  - the maximum accepted distortion for unfolding optimization.
* `nchunks` - number of rounds of optimization.
"""
function unfold(
    to_unf::AbstractMatrix,
    ref::AbstractMatrix;
    isomap = :default,
    optim = :default,
)

    # read isomap and optim parameters or assign the defaults below
    ipars = (search = :knn, neigh = 16, anchors = 1500, reftol = 0.01)
    opars = (search = :knn, neigh = 50, maxerr = 5, nchunks = 4)

    ipars = updatepars(ipars, isomap)
    opars = updatepars(opars, optim)

    # conversions if necessary
    !(ref[1] isa Float64) && (ref = Float64.(ref))
    !(to_unf[1] isa Float64) && (to_unf = Float64.(to_unf))

    # random seed
    Random.seed!(1234567890)

    # pre-process reference points
    ref = remove_duplicates(ref, tol = ipars.reftol)
    resol = get_resolution(ref)

    # do landmark isomap at reference points
    unf_ref = landmark_isomap(ref, ipars.search, ipars.neigh, ipars.anchors)

    # get initial guess
    normals, good = getnormals(ref, ipars.search, ipars.neigh)
    unf = firstguess(
        to_unf,
        view(ref, :, good),
        view(unf_ref, :, good),
        view(normals, :, good),
    )

    # unfold points in random chunks
    ndom = size(to_unf, 2)
    shuffled_ids = shuffle(1:ndom)
    nchunks = ceil(Int, ndom / opars.nchunks)
    chunks = collect(Iterators.partition(shuffled_ids, nchunks))

    for (i, ids) in enumerate(chunks)
        # only reference surface as conditioner
        known_org = ref[:, good]
        known_unf = unf_ref[:, good]
        ids_to_unf = ids

        # add unfolded from previous chunks
        if i > 1
            extra_ids = vcat(chunks[1:(i-1)]...)
            extra_org = view(to_unf, :, extra_ids)
            extra_unf = view(unf, :, extra_ids)
            good_, bad_ =
                error_ids(extra_org, extra_unf, opars.search, opars.neigh, opars.maxerr)

            known_org = hcat(known_org, view(extra_org, :, good_))
            known_unf = hcat(known_unf, view(extra_unf, :, good_))
            ids_to_unf = Int.(union(ids, view(extra_ids, bad_)))
        end

        # unfold points
        initguess = view(unf, :, ids_to_unf)
        chunk_to_unf = view(to_unf, :, ids_to_unf)
        unf[:, ids_to_unf] .=
            opt(known_org, known_unf, chunk_to_unf, opars.search, opars.neigh, initguess)
    end

    unf
end

function unfold(
    to_unf::AbstractVector,
    ref::AbstractMatrix;
    isomap = :default,
    optim = :default,
)

    merged_to_unf = reduce(hcat, to_unf)
    unf = unfold(merged_to_unf, ref, isomap = isomap, optim = optim)
    j_ = cumsum(size.(to_unf, 2))
    i_ = [i == 1 ? 1 : j_[i-1] + 1 for i = 1:length(j_)]
    [unf[:, i:j] for (i, j) in zip(i_, j_)]
end

"""
	unfold(to_unf, ref, unf_ref; optim=:default)

Unfold the input points `to_unf` based on points already unfolded `unf_ref` and
their correspondent in the original space `ref`. Useful to unfold samples after
the domain is already unfolded, for example. Returns a coordinate matrix with
the unfolded points. Or a tuple of matrices if `to_unf` is a list of points.

## Parameters:

* `to_unf`  - coordinate matrix with the points that will be unfolded
* `ref`     - coordinate matrix with the reference points for unfolding
* `unf_ref` - coordinate matrix with the unfolded reference points
* `optim`   - additional parameters for optimization step; either :default or a
  NamedTuple with the keys shown below

### `optim` parameters:

*Default:* optim = (search=:knn, neigh=50, maxerr=5, nchunks=1)
* `search`  - search type to optimize locally (:knn for k-nearest neighbor or
  :radius for radius search).
* `neigh`   - number of neighbors (for search=:knn) or radius distance
  (for search=:radius) to local optimization.
* `maxerr`  - the maximum accepted distortion for unfolding optimization.
* `nchunks` - number of rounds of optimization.
"""
function unfold(
    to_unf::AbstractMatrix,
    ref::AbstractMatrix,
    unf_ref::AbstractMatrix;
    optim = :default,
)

    # read optim parameters or assign the defaults below
    opars = (search = :knn, neigh = 50, maxerr = 5, nchunks = 2)
    opars = updatepars(opars, optim)

    # conversions if necessary
    !(to_unf[1] isa Float64) && (to_unf = Float64.(to_unf))

    # unfold points in random chunks
    ndom = size(to_unf, 2)
    shuffled_ids = shuffle(1:ndom)
    nchunks = ceil(Int, ndom / opars.nchunks)
    chunks = collect(Iterators.partition(shuffled_ids, nchunks))
    unf = firstguess(to_unf, ref, unf_ref)

    for (i, ids) in enumerate(chunks)
        # only reference surface as conditioner
        known_org = ref
        known_unf = unf_ref
        ids_to_unf = ids

        # add unfolded from previous chunks
        if i > 1
            extra_ids = vcat(chunks[1:(i-1)]...)
            extra_org = view(to_unf, :, extra_ids)
            extra_unf = view(unf, :, extra_ids)
            good_, bad_ =
                error_ids(extra_org, extra_unf, opars.search, opars.neigh, opars.maxerr)

            known_org = hcat(known_org, view(extra_org, :, good_))
            known_unf = hcat(known_unf, view(extra_unf, :, good_))
            ids_to_unf = Int.(union(ids, view(extra_ids, bad_)))
        end

        # unfold points
        initguess = view(unf, :, ids_to_unf)
        chunk_to_unf = view(to_unf, :, ids_to_unf)
        unf[:, ids_to_unf] .=
            opt(known_org, known_unf, chunk_to_unf, opars.search, opars.neigh, initguess)
    end
    unf
end

function unfold(
    to_unf::AbstractVector,
    ref::AbstractMatrix,
    unf_ref::AbstractMatrix;
    optim = :default,
)

    merged_to_unf = reduce(hcat, to_unf)
    unf = unfold(merged_to_unf, ref, unf_ref, optim = optim)
    j_ = cumsum(size.(to_unf, 2))
    i_ = [i == 1 ? 1 : j_[i-1] + 1 for i = 1:length(j_)]
    [unf[:, i:j] for (i, j) in zip(i_, j_)]
end

function updatepars(default, newpars)
    if newpars != :default
        for (k, v) in zip(keys(newpars), values(newpars))
            default = Setfield.setindex(default, v, k)
        end
    end
    default
end
