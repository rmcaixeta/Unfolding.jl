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
  input_block = coords(df_block, columns=["XC", "YC", "ZC"])
  input_samps = coords(df_samps, columns=["X", "Y", "Z"])

  # Get reference surface points for unfolding
  ref_surface = getreference(input_block, axis=["X", "Y"])
  println("- Reference surface extracted")

  # Sample blocks to make tests feasible for CI
  input_block = input_block[:, df_block.ZC .< -90]
  input_samps = input_samps[:, df_samps.Z .< -90]
  ref_surface = ref_surface[:, ref_surface[3, :] .< -90]

  # Get transformed coordinates of blocks and samples after unfolding
  optim = (neigh=20,)
  isomap = (anchors=100,)
  unf_block = unfold(input_block, ref_surface; optim, isomap)
  unf_samps = unfold(input_samps, input_block, unf_block; optim)
  println("- Unfolding finished")

  to_vtk("test_out_vtk", unf_samps, (test=rand(size(unf_samps, 2)),))

  all_input = hcat(input_samps, input_block)
  all_unf = hcat(unf_samps, unf_block)
  error = errors(:dists, all_input, all_unf)
  good, bad = errors(:ids, all_input, all_unf)

  @test length(good) > (2 * length(bad))

  println("- Quick test of unfold with multiple inputs")
  s1, s2 = unfold([input_samps[:, 1:5], input_samps[:, 5:end]], ref_surface; optim, isomap)
  s1, s2 = unfold([input_samps[:, 1:5], input_samps[:, 5:end]], input_block, unf_block; optim)
end
