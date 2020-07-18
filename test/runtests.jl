using Unfolding
using CSV
using Test

@testset "Unfolding.jl" begin
    # Reading data
    df_samp = CSV.read("samples.csv")
    df_block = CSV.read("block_model.csv")

    # Get coordinate points as matrix
    input_block = coordinate_matrix( df_block, columns=["XC","YC","ZC"] )
    input_samp = coordinate_matrix( df_samp, columns=["X","Y","Z"] )

    # Get reference surface points for unfolding
    ref_surface = ref_surface_from_blocks(input_block, axis="Y")
    # Get transformed coordinates of blocks and samples after unfolding
    unf_block, unf_samp = unfold(ref_surface,input_block,input_samp)

    # Write new XT, YT and ZT columns with the transformed coordinates
    for (i,c) in enumerate([:XT,:YT,:ZT])
        df_samp[:,c] = unf_samp[i,:]
        df_block[:,c] = unf_block[i,:]
    end

    good, bad = unfold_error_ids(hcat(input_samp,input_block), hcat(unf_samp,unf_block), nneigh=5, max_error=5)
    @test length(good)>90*length(bad)
end
