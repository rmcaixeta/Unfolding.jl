using Unfolding
using CSV
using DataFrames
using Test

@testset "Unfolding.jl" begin
    # Reading data
    df_samps = CSV.read("samples.csv", DataFrame)
    df_block = CSV.read("block_model.csv", DataFrame)
    println("- Samples loaded")

    # Get coordinate points as matrix
    input_block = coords(df_block, columns=["XC","YC","ZC"])
    input_samps = coords(df_samps, columns=["X","Y","Z"])

    # Get reference surface points for unfolding
    ref_surface = getreference(input_block, axis=["X","Y"])
    println("- Reference surface extracted")

    # Get transformed coordinates of blocks and samples after unfolding
    unf_block, unf_samps = unfold(ref_surface, input_block, input_samps)
    println("- Unfolding finished")

    # Write new XT, YT and ZT columns with the transformed coordinates
    for (i,c) in enumerate([:XT,:YT,:ZT])
        df_samps[:,c] = unf_samps[i,:]
        df_block[:,c] = unf_block[i,:]
    end

    to_vtk(unf_samps,"test_out_vtk",[(string("test",x),[rand() for i in 1:size(unf_samps,2)]) for x in 1:4])
    to_csv(unf_samps,"test_out_csv",["X","Y","Z"])
    error = error_dists(hcat(input_samps,input_block), hcat(unf_samps,unf_block), nneigh=8)

    good, bad = error_ids(hcat(input_samps,input_block), hcat(unf_samps,unf_block))
    @test length(good)>50*length(bad)
end
