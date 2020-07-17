

# Adjust matrix format
function coordinate_matrix(x;columns=nothing)

	x = columns==nothing ? x : x[:,[Symbol(v) for v in columns]]
	shape = size(x)
	@assert (shape[1]==3 || shape[2]==3) "Three coordinate matrix should be informed"

	x = typeof(x)<:Matrix{<:Number} ? x : Matrix{Float64}(x)
	x = (shape[2]==3 && shape[1]!=3) ? permutedims(x) : x
	return x
end

function data_to_csv(coords::AbstractArray,outname::String,colnames)

	out = length(size(coords))==1 ? coords : coords'
	open(string(outname,".csv"); write=true) do f
		write(f, string(join(colnames,","),"\n"))
		writedlm(f, out, ',')
	end
end

function data_to_vtk(coords::AbstractArray{<:Number,2},outname::String,extra_props=nothing)
	verts = [MeshCell( VTKCellTypes.VTK_VERTEX, [i]) for i in 1:size(coords)[2] ]
	outfiles = vtk_grid(outname,coords,(verts)) do vtk
		if extra_props!=nothing
			for (name,prop) in extra_props
				vtk[name] = prop
			end
		end
	end

end

function _make_boxplot(error,plotname)
	fig = boxplot(["Unfolding error"], error)
	savefig(plotname)
end
