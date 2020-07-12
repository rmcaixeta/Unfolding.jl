using Unfolding
using CSV
using Test

@testset "Unfolding.jl" begin
    # Input info
    block_file = "./block_model.csv"
    block_coords = ["XC","YC","ZC"]

    samp_file = "./samples.csv"
    samp_coords = ["X","Y","Z"]

    # Reading file
    df_samp = CSV.read( samp_file, select=samp_coords )
    df_block = CSV.read( block_file, select=block_coords )

    # Get coordinate points as matrix of shape (3 x nb_points)
    input_block = Matrix{Float64}(df_block)'
    input_samp = Matrix{Float64}(df_samp)'

    # Get reference surface points for unfolding. Otherwise, read it from out
    ref_surface = ref_surface_from_blocks(input_block, axis="Y")
    unf_block, unf_samp = unfold(ref_surface,input_block,input_samp)

    # Write new XT, YT and ZT columns with the transformed coordinates
    for (i,c) in enumerate([:XT,:YT,:ZT])
        df_samp[c] = unf_samp[i,:]
        df_block[c] = unf_block[i,:]
    end

    good, bad = unfold_error(refsurf_true_coords, ref_surf_transf, 5, 5, false)
    @test length(good)>2*length(bad)
end
