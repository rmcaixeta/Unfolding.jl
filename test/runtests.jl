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
    optim = (neigh=15,)
    unf_block, unf_samps = unfold(ref_surface, input_block, input_samps, optim=optim)
    println("- Unfolding finished")

    # Write new XT, YT and ZT columns with the transformed coordinates
    for (i,c) in enumerate([:XT,:YT,:ZT])
        df_samps[:,c] = unf_samps[i,:]
        df_block[:,c] = unf_block[i,:]
    end

    to_vtk("test_out_vtk", unf_samps, (test=rand(size(unf_samps,2)),))

    all_input = hcat(input_samps,input_block)
    all_unf   = hcat(unf_samps,unf_block)
    error     = errors(:dists, all_input, all_unf)
    good, bad = errors(:ids, all_input, all_unf)

    @test length(good) > (30 * length(bad))
end
