
"""
	coords(X; columns=nothing)

Reshape and filter dataframe or matrix and return the matrix in right format
for unfolding.

## Parameters:

* `X`       - the input Dataframe or matrix with the coordinates.
* `columns` - if the input is a Dataframe with more than 3 columns, the column
  names of the three coordinates must be informed here (e.g. ["X","Y","Z"]).
"""
function coords(x; columns=nothing)
  !isnothing(columns) && (x = x[:, Symbol.(columns)])
  shape = size(x)
  @assert (3 in shape) "3-D coordinate matrix should be informed"
  !(x isa Matrix) && (x = Matrix{Float64}(x))
  (shape[2] == 3 && shape[1] != 3) && (x = permutedims(x))
  x
end

"""
	to_vtk(outname, coordinates, props)

Export the `coordinates` as VTK points named `outname`.vtu. Optionally, extra
properties are passed via `props` as `NamedTuple` or `Dict`, where the first item
is the property name and the second is an array with the property values.
"""
function to_vtk(outname::String, coords::AbstractMatrix, props=nothing)
  verts = [MeshCell(VTKCellTypes.VTK_VERTEX, [i]) for i in 1:size(coords, 2)]
  outfiles = vtk_grid(outname, coords, (verts)) do vtk
    if !isnothing(props)
      for (name, prop) in zip(keys(props), values(props))
        vtk[string(name)] = prop
      end
    end
  end
end

# get neighbors
function get_neighbors(origin::AbstractMatrix, target::AbstractMatrix, nhood::Symbol, neigh_val, calcdists=false)
  @assert nhood in [:knn, :radius] "invalid neighborhood type"
  tree = KDTree(origin)
  if nhood == :knn
    nneigh = minimum([neigh_val, size(origin, 2)])
    idxs, dists = knn(tree, target, nneigh, true)
    idxs, dists
  elseif nhood == :radius
    idxs = inrange(tree, target, neigh_val)
    if calcdists
      dists = Vector{Vector{Float64}}(undef, 0)
      for i in 1:length(idxs)
        p1 = view(target, :, i)
        d = [euclidean(p1, view(origin, :, j)) for j in idxs[i]]
        push!(dists, d)
      end
    else
      dists = nothing
    end
    idxs, dists
  end
end

get_neighbors(ref_coords::AbstractMatrix, nhood, neigh_val, calcdists=false) =
  get_neighbors(ref_coords, ref_coords, nhood, neigh_val, calcdists)
