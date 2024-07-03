
"""
    getreference(blocks; axis=["X","Y","Z"])

Extract the mid-surface points that split the informed block model in two.
Returns a coordinate matrix with the reference points for unfolding.

## Parameters:

* `blocks` - coordinate matrix of the regular blocks centroids.
* `axis`   - axis or list of axis that will be scanned to extract reference
  points. Defaults to ["X","Y","Z"].
"""
function getreference(blocks::AbstractMatrix; axis=["X","Y","Z"])
	typeof(axis) in (String, Symbol) && (axis = [axis])
	axis[1] isa Symbol && (axis = String.(axis))
	@assert axis ⊆ ["X","Y","Z"] "invalid axis: $axis"
	coord_ids = Dict("X"=>1, "Y"=>2, "Z"=>3)

	# get cell sizes
	cells = zeros(Float64,3)
	for c in 1:3
		coords   = sort!(unique(blocks[c,:]))
		cells[c] = minimum([coords[x]-coords[x-1] for x in 2:length(coords)])
	end

	# scan sections and save skeleton points
	ref = Float64[]
	for c in axis
		ax = coord_ids[c]
	    axis_coords = sort!(unique(blocks[ax,:]))
		sec = setdiff([1,2,3],ax)

	    for s in axis_coords
	        section = view(blocks,:,view(blocks,ax,:) .== s)

	        min_i = minimum(section[sec[1],:])
	        min_j = minimum(section[sec[2],:])
	        i = section[sec[1],:] .- min_i
	        j = section[sec[2],:] .- min_j

			i = Int.(i/cells[sec[1]]) .+ 1
			j = Int.(j/cells[sec[2]]) .+ 1

	        img = zeros(Bool,maximum(i)+1,maximum(j)+1)

	        for x in 1:length(i)
	            img[i[x],j[x]]=true
	        end

	        img = thinning(img)

	        for x in 1:length(i)
	            if img[i[x],j[x]]==true
	                io = (i[x]-1)*cells[sec[1]]+min_i
	                jo = (j[x]-1)*cells[sec[2]]+min_j
					order = sortperm([ax,sec[1],sec[2]])
					append!(ref, [s, io, jo][order])
	            end
	        end
	    end
	end

	reshape(ref, 3, Int(length(ref)/3))
end

# Average duplicate points
function remove_duplicates(coords::AbstractMatrix; tol=0.01)
	groups = [union(x.core_indices,x.boundary_indices) for x in dbscan(coords, tol).clusters]
	outcoords = zeros(Float64,3,length(groups))
	for g in 1:length(groups)
		idx = groups[g]
		if length(idx)==1
			outcoords[:,g] .= dropdims(view(coords,:,idx),dims=2)
		else
			outcoords[:,g] .= mean!(view(outcoords,:,g),view(coords,:,idx))
		end
	end
	outcoords
end

function get_resolution(ref_pts::AbstractMatrix)
	idxs, dists = get_neighbors(ref_pts, :knn, 2)
	closest_pt = [sum(x) for x in dists]
	quantile(closest_pt,0.75)
end
