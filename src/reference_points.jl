
# Skeletonization
function ref_surface_from_blocks(blks::AbstractArray{<:Number,2};
	 axis=["X","Y","Z"])

	if typeof(axis)==String
		axis = [axis]
	end
	coord_ids = Dict("X"=>1, "Y"=>2, "Z"=>3)
	ref = [Float64[],Float64[],Float64[]]
	cells = zeros(Float64,3)
	for c in 1:3
		coords =  sort!(unique(blks[c,:]))
		cell_ = sort!(unique([coords[x]-coords[x-1] for x in 2:length(coords)]))
		#@assert length(cell_)==1 "Block model must have regular cell sizes"
		cells[c] = cell_[1]
	end

	#@assert axis in ["X","Y","Z"] "Invalid axis"
	for c in axis
		ax = coord_ids[c]
	    axis_coords = sort!(unique(blks[ax,:]))
		sec = setdiff([1,2,3],ax)

	    for s in axis_coords
	        section = blks[:,blks[ax,:] .== s]

	        min_i = minimum(section[sec[1],:])
	        min_j = minimum(section[sec[2],:])
	        i = section[sec[1],:] .- min_i
	        j = section[sec[2],:] .- min_j

	        i = convert(Array{Int,1}, i/cells[sec[1]]) .+ 1
	        j = convert(Array{Int,1}, j/cells[sec[2]]) .+ 1

	        img = zeros(Bool,maximum(i)+1,maximum(j)+1)

	        for x in 1:length(i)
	            img[i[x],j[x]]=true
	        end

	        img = thinning(img)

	        for x in 1:length(i)
	            if img[i[x],j[x]]==true
	                io = (i[x]-1)*cells[sec[1]]+min_i
	                jo = (j[x]-1)*cells[sec[2]]+min_j
	                push!(ref[ax], s)
	                push!(ref[sec[1]], io)
	                push!(ref[sec[2]], jo)
	            end
	        end
	    end
	end
	out_surf = permutedims(hcat(ref[1],ref[2],ref[3]))

	if length(axis)==1
		return out_surf
	else
		return _remove_duplicates(out_surf)
	end
end


# Average duplicate points
function _remove_duplicates(coords::AbstractArray{<:Number,2},tol=0.5)
	groups = [union(x.core_indices,x.boundary_indices) for x in dbscan(coords, tol)]
	outcoords = zeros(Float64,3,length(groups))
	for g in 1:length(groups)
		idx = groups[g]
		if length(idx)==1
			outcoords[:,g] .= dropdims(coords[:,idx],dims=2)
		else
			outcoords[:,g] .= mean!(outcoords[:,g],coords[:,idx])
		end
	end
	return outcoords
end
