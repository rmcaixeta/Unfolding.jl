module Unfolding

using Clustering: dbscan
using DelimitedFiles
using Distances
using ImageMorphology: thinning
using LightGraphs: connected_components, dijkstra_shortest_paths
using LinearAlgebra: Diagonal, Symmetric, cross, dot, eigen
using MultivariateStats: PCA, dmat2gram, fit, projection
using NearestNeighbors
using Optim
using Random
using Setfield
using SimpleWeightedGraphs
using StatsBase: Weights, mean, mean!, quantile, sample
using WriteVTK

include("coords_optimization.jl")
include("data_handling.jl")
include("landmark_isomap.jl")
include("reference_points.jl")
include("unfold.jl")

export
    coords,
    errors,
    getreference,
    landmark_isomap,
    to_vtk,
    unfold
end
