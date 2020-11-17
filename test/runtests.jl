using Unfolding
using CSV
using DataFrames
using Test

@testset "Unfolding.jl" begin
    # Reading data
    df_samp = CSV.read("samples.csv", DataFrame)
    df_block = CSV.read("block_model.csv", DataFrame)
    println("- Samples loaded")

    # Get coordinate points as matrix
    input_block = coordinate_matrix( df_block, columns=["XC","YC","ZC"] )
    input_samp = coordinate_matrix( df_samp, columns=["X","Y","Z"] )

    # Get reference surface points for unfolding
    ref_surface = ref_surface_from_blocks(input_block, axis=["X","Y"])
    println("- Reference surface extracted")

    # Get transformed coordinates of blocks and samples after unfolding
    unf_block, unf_samp = unfold(ref_surface,input_block,input_samp)
    println("- Unfolding finished")

    # Write new XT, YT and ZT columns with the transformed coordinates
    for (i,c) in enumerate([:XT,:YT,:ZT])
        df_samp[:,c] = unf_samp[i,:]
        df_block[:,c] = unf_block[i,:]
    end

    data_to_vtk(unf_samp,"test_out_vtk",[(string("test",x),[rand() for i in 1:size(unf_samp)[2]]) for x in 1:4])
    data_to_csv(unf_samp,"test_out_csv",["X","Y","Z"])
    #error = unfold_error_dists(hcat(input_samp,input_block), hcat(unf_samp,unf_block), nneigh=8, plotname="test_boxplot")

    good, bad = unfold_error_ids(hcat(input_samp,input_block), hcat(unf_samp,unf_block))
    @test length(good)>50*length(bad)
end
