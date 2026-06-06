"""
	unfold(to_unf, ref; isomap=:default, optim=:default, normals=nothing)

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
* `normals` - normals matrix of the reference points; auto-calculated if not informed

### `isomap` parameters:

*Default:* isomap = (search=:knn, neigh=16, anchors=1500, reftol=0.01, align_axis=:x)
* `search`  - search type to build neighbors graph for Isomap (:knn for
  k-nearest neighbor or :radius for radius search).
* `neigh`   - number of neighbors (for search=:knn) or radius distance
  (for search=:radius) to build neighbors graph for Isomap.
* `anchors` - number of anchors for landmark isomap
* `reftol`  - distance to which reference points are considered duplicates
* `align_axis` - make :x or :y unfolded axis correlate well with original one

### `optim` parameters:

*Default:* optim = (search=:knn, neigh=50, maxerr=5, nchunks=4)
* `search`  - search type to optimize locally (:knn for k-nearest neighbor or
  :radius for radius search).
* `neigh`   - number of neighbors (for search=:knn) or radius distance
  (for search=:radius) to local optimization.
* `maxerr`  - the maximum accepted distortion for unfolding optimization.
* `nchunks` - number of rounds of optimization.
"""
function unfold(to_unf::AbstractMatrix, ref::AbstractMatrix; isomap=:default, optim=:default, normals=nothing)

  # read isomap and optim parameters or assign the defaults below
  ipars = (search=:knn, neigh=16, anchors=1500, reftol=0.01, align_axis=:x)
  opars = (search=:knn, neigh=50, maxerr=5, nchunks=4)

  ipars = updatepars(ipars, isomap)
  opars = updatepars(opars, optim)

  # conversions if necessary
  !(ref[1] isa Float64) && (ref = Float64.(ref))
  !(to_unf[1] isa Float64) && (to_unf = Float64.(to_unf))

  # random seed
  Random.seed!(1234567890)

  # pre-process reference points
  ref = remove_duplicates(ref, tol=ipars.reftol)

  # do landmark isomap at reference points
  unf_ref = landmark_isomap(ref, ipars.search, ipars.neigh, ipars.anchors, ipars.align_axis)

  # unfold other pts
  displ = unf_ref .- ref
	ttab = georef((I=displ[1,:], J=displ[2,:], K=displ[3,:],
	               x=ref[1,:], y=ref[2,:], z=ref[3,:]), (:x,:y,:z))
  
  dom = georef((x=to_unf[1,:], y=to_unf[2,:], z=to_unf[3,:]), (:x,:y,:z))
  dom = domain(dom)
	
	# displacement
	uvw = mapreduce(hcat, [:I,:J,:K]) do ev
		ttab |> Select(ev) |> InterpolateNeighbors(dom, model=IDW(1), maxneighbors=40)
	end |> values
	
  unf = to_unf .+ [uvw.I uvw.J uvw.K]'
  unf
end

function unfold(to_unf::AbstractVector, ref::AbstractMatrix; isomap=:default, optim=:default, normals=nothing)
  merged_to_unf = reduce(hcat, to_unf)
  unf = unfold(merged_to_unf, ref; isomap, optim, normals)
  j_ = cumsum(size.(to_unf, 2))
  i_ = [i == 1 ? 1 : j_[i - 1] + 1 for i in 1:length(j_)]
  [unf[:, i:j] for (i, j) in zip(i_, j_)]
end

"""
	unfold(to_unf, ref, unf_ref; optim=:default, normals=nothing)

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
* `normals` - normals matrix of the reference points; auto-calculated if not informed

### `optim` parameters:

*Default:* optim = (search=:knn, neigh=50, maxerr=5, nchunks=2)
* `search`  - search type to optimize locally (:knn for k-nearest neighbor or
  :radius for radius search).
* `neigh`   - number of neighbors (for search=:knn) or radius distance
  (for search=:radius) to local optimization.
* `maxerr`  - the maximum accepted distortion for unfolding optimization.
* `nchunks` - number of rounds of optimization.
"""
function unfold(to_unf::AbstractMatrix, ref::AbstractMatrix, unf_ref::AbstractMatrix; optim=:default, normals=nothing)

  # read optim parameters or assign the defaults below
  opars = (search=:knn, neigh=50, maxerr=5, nchunks=4)
  opars = updatepars(opars, optim)

  # conversions if necessary
  !(to_unf[1] isa Float64) && (to_unf = Float64.(to_unf))
  
  # unfold other pts
  displ = unf_ref .- ref
	ttab = georef((I=displ[1,:], J=displ[2,:], K=displ[3,:],
	               x=ref[1,:], y=ref[2,:], z=ref[3,:]), (:x,:y,:z))
  
  dom = georef((x=to_unf[1,:], y=to_unf[2,:], z=to_unf[3,:]), (:x,:y,:z))
  dom = domain(dom)
	
	# displacement
	uvw = mapreduce(hcat, [:I,:J,:K]) do ev
		ttab |> Select(ev) |> InterpolateNeighbors(dom, model=IDW(1), maxneighbors=40)
	end |> values
	
  unf = to_unf .+ [uvw.I uvw.J uvw.K]'
  unf
end

function unfold(to_unf::AbstractVector, ref::AbstractMatrix, unf_ref::AbstractMatrix; optim=:default, normals=nothing)
  merged_to_unf = reduce(hcat, to_unf)
  unf = unfold(merged_to_unf, ref, unf_ref; optim, normals)
  j_ = cumsum(size.(to_unf, 2))
  i_ = [i == 1 ? 1 : j_[i - 1] + 1 for i in 1:length(j_)]
  [unf[:, i:j] for (i, j) in zip(i_, j_)]
end

function unfold(ref::AbstractMatrix; isomap=:default)

  # read isomap and optim parameters or assign the defaults below
  ipars = (search=:knn, neigh=16, anchors=1500, reftol=0.01, align_axis=:x)

  ipars = updatepars(ipars, isomap)

  # conversions if necessary
  !(ref[1] isa Float64) && (ref = Float64.(ref))

  # random seed
  Random.seed!(1234567890)

  # pre-process reference points
  ref = remove_duplicates(ref, tol=ipars.reftol)

  # do landmark isomap at reference points
  unf_ref = landmark_isomap(ref, ipars.search, ipars.neigh, ipars.anchors, ipars.align_axis)
  ref, unf_ref
end

function updatepars(default, newpars)
  if newpars != :default
    for (k, v) in zip(keys(newpars), values(newpars))
      default = Setfield.setindex(default, v, k)
    end
  end
  default
end
