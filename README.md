# Unfolding.jl

[![Build Status](https://travis-ci.com/rmcaixeta/Unfolding.jl.svg?branch=master)](https://travis-ci.com/rmcaixeta/Unfolding.jl)
[![Coverage](https://codecov.io/gh/rmcaixeta/Unfolding.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rmcaixeta/Unfolding.jl)

`Unfolding.jl` is a package written in Julia to perform unfolding of 3-D geometries. It was developed for geostatistical cases where complex 3-D domains are modeled and need to be unfolded for appropriate variography, estimations and simulations.

Julia was used due to its high performance and easy coding. This package was successfully tested with some big mining datasets but is still under development, so please enter in contact if you have some issue or feel free to contribute to the code (it is open source!).

## Installation

It is necessary to install Julia to run this code. Installation instructions for Windows, Linux and macOS are available [here](https://julialang.org/downloads/platform/).

## Usage

There a two possible workflows:

* Unfolding based on a block model
* Unfolding based on a points of a reference surface

The example below is based on a block model within a folded domain. It is only necessary the X, Y and Z coordinates of the block centroids. Optionally, samples coordinates can also be informed.

`FIGURE HERE
(blocks 2x2x2m; blocks outside the domain are not informed)`

### Julia example

The usage in Julia is detailed in the code below.

`Julia code here`

### Python example

Julia is not a widespread language yet. For those more familiar with Python, the Julia code can be called inside Python script and the data can be informed as Numpy arrays. Julia still need to be installed before calling it in Python. Additionally, an extra Python library shoulde be installed in this case:

`pip install julia`

The data for the unfolding functions should be informed as Numpy arrays.

`Python code here`

### Results

The output information are the transformed coordinates of the original data. It can be saved in CSV or VTK format for further uses. The example below loads the output VTK data using Python [`pyvista`](https://www.pyvista.org/) library.

`Output figure`

## Unfolding documentation

Three main functions are available for use:

### ref_surface_from_blocks()
### unfold()
### data_to_vtk()
