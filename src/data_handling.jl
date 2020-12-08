

"""
	coordinate_matrix(X; columns=nothing)

Reshape and filter dataframe or matrix and return the matrix in right format
for unfolding.

## Parameters:

* `X`       - the input Dataframe or matrix with the coordinates.
* `columns` - if the input is a Dataframe with more than 3 columns, the column
  names of the three coordinates must be informed here (e.g. ["X","Y","Z"]).
"""
function coordinate_matrix(x; columns=nothing)

	x = columns==nothing ? x : x[:,[Symbol(v) for v in columns]]
	shape = size(x)
	@assert (shape[1]==3 || shape[2]==3) "Three coordinate matrix should be informed"

	x = typeof(x)<:Matrix{<:Number} ? x : Matrix{Float64}(x)
	x = (shape[2]==3 && shape[1]!=3) ? permutedims(x) : x
	x
end

"""
	data_to_csv(input_matrix, outname, colnames)

Export the `input_matrix` as CSV table file named `outname`.csv. The column
names are passed in `colnames` array.
"""
function data_to_csv(input_matrix::AbstractArray,outname::String,colnames)

	out = length(size(input_matrix))==1 ? input_matrix : input_matrix'
	open(string(outname,".csv"); write=true) do f
		write(f, string(join(colnames,","),"\n"))
		writedlm(f, out, ',')
	end
end

"""
	data_to_vtk(coords, outname, extra_props)

Export the `input_matrix` as VTK points named `outname`.vtu. Optionally, extra
properties are passed via an array of tuples `extra_props`, where the first item
of the tuple is the column name and the second is an array with the property values.
"""
function data_to_vtk(coords::AbstractArray{<:Number,2},outname::String,extra_props=nothing)
	verts = [MeshCell( VTKCellTypes.VTK_VERTEX, [i]) for i in 1:size(coords,2) ]
	outfiles = vtk_grid(outname,coords,(verts)) do vtk
		if extra_props!=nothing
			for (name,prop) in extra_props
				vtk[name] = prop
			end
		end
	end

end
