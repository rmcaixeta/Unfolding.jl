module Unfolding

using Clustering:dbscan
using DelimitedFiles
using ImageMorphology:thinning
using LinearAlgebra:Diagonal,eigen,cross,dot
using MultivariateStats:dmat2gram,PCA,fit,projection
using NearestNeighbors
using Random
using SimpleWeightedGraphs
using StatsBase:mean,mean!,Weights,sample,standardize,ZScoreTransform
using StatsPlots:boxplot,savefig
using WriteVTK

using Distributed
using Distances
using LightGraphs:dijkstra_shortest_paths,connected_components
using Optim

include("coords_optimization.jl")
include("data_handling.jl")
include("landmark_isomap.jl")
include("reference_points.jl")
include("unfold.jl")

export
    unfold,
    ref_surface_from_blocks,
    unfold_error_ids,
    unfold_error_dists,
    coordinate_matrix,
    data_to_vtk,
    data_to_csv,
    landmark_isomap

end
