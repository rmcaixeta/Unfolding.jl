# Unfolding.jl

[![Build Status][build-img]][build-url] [![Coverage][codecov-img]][codecov-url]

`Unfolding.jl` is a package written in Julia to perform unfolding of 3-D geometries. It was developed for geostatistical cases where complex 3-D domains are modeled and need to be unfolded for appropriate spatial analysis, estimations and simulations.

Julia was used due to its high performance and easy coding. This package was successfully tested with some big mining datasets but is still under development, so please enter in contact if you have some issue or feel free to contribute to the code.

## Installation

It is necessary to install Julia to run this code. Installation instructions for Windows, Linux and macOS are available [here](https://julialang.org/downloads/platform/).

After that, it is necessary to install the Unfolding package. Open a terminal, type `julia` to open the REPL and then install the package with the following command. Additionally, the CSV and DataFrames packages are also installed to run the examples in the sequence.

```julia
using Pkg; Pkg.add(["Unfolding","CSV","DataFrames"])
```

## Usage

The algorithm is based on [Caixeta & Costa (2021)](https://doi.org/10.1016/j.cageo.2021.104844). There are three possible workflows:

* Extract reference points from a given block model. Run unfolding for any points using them.
* Reference points are already available from other sources. Run unfolding for any points using them.
* Extract reference points from a given mesh. Run unfolding for any points using them _(not available yet)_.

The example below is based on a block model within a folded domain. It is only necessary the X, Y and Z coordinates of the block centroids. Optionally, samples coordinates can also be informed.

<p align="center">
  <img src="docs/images/unfolding_usage.png">
</p>

### Julia example

The usage in Julia is detailed in the code below.

```julia
# Julia code
using CSV
using DataFrames
using Unfolding

# Reading data
df_samp = CSV.read("samples.csv", DataFrame)
df_block = CSV.read("block_model.csv", DataFrame)

# Get coordinate points as matrix
input_block = coords( df_block, columns=["XC","YC","ZC"] )
input_samp = coords( df_samp, columns=["X","Y","Z"] )

# Get reference surface points for unfolding
ref_surface = getreference(input_block)
# Get transformed coordinates of blocks and samples after unfolding
unf_block, unf_samp = unfold([input_block,input_samp], ref_surface)

# Write new XT, YT and ZT columns with the transformed coordinates
for (i,c) in enumerate([:XT,:YT,:ZT])
    df_samp[:,c] = unf_samp[i,:]
    df_block[:,c] = unf_block[i,:]
end

# Write output to CSV
CSV.write("out_dh.csv", df_samp )
CSV.write("out_blks.csv", df_block )

# Write output to VTK format
to_vtk("out_blks", unf_block)
to_vtk("out_dh", unf_samp)
```

The code can be saved in a textfile with `.jl` extension and be called in a terminal: `julia file.jl` or `julia -t 4 file.jl` to run faster using 4 threads (or any number of threads you want). Or you can organize it in notebooks.

### Python example (experimental)

Julia is not a widespread language yet. For those more familiar with Python, the Julia code can be called inside Python scripts.

Julia still need to be installed before calling it in Python. Additionally, an extra Python library must be installed in this case:

```python
pip install julia
import julia
julia.install()
```

The data for the unfolding functions should be informed as Numpy arrays. Note that the input arrays are transposed because multidimensional arrays in Julia are stored in column-major order.

```python
# Python example
from julia import Julia
Julia(compiled_modules=False) # excluding these first lines makes it run faster; but may crash in some systems

import pandas as pd
from julia import Unfolding as unf

# Reading data
df_block = pd.read_csv("block_model.csv",usecols=["XC","YC","ZC"])
df_samp = pd.read_csv("samples.csv",usecols=["X","Y","Z"])

# Get coordinate points as matrix. Data points must be informed in columns instead rows
input_block = df_block.to_numpy().T
input_samp = df_samp.to_numpy().T

# Get reference surface points for unfolding
ref_surface = unf.getreference(input_block)
# Get transformed coordinates of blocks and samples after unfolding
unf_block, unf_samp = unf.unfold([input_block,input_samp], ref_surface)

# Write new XT, YT and ZT columns with the transformed coordinates
for i,c in enumerate(["XT","YT","ZT"]):
    df_block[c] = unf_block[i,:]
    df_samp[c] = unf_samp[i,:]

# Write output to CSV
df_block.to_csv("out_blks.csv",index=False)
df_samp.to_csv("out_samp.csv",index=False)

# Write output to VTK format
unf.to_vtk("out_blks", unf_block)
unf.to_vtk("out_dh", unf_samp)
```

### Results

The output information are the transformed coordinates of the original data. It can be saved in CSV or VTK format for further uses. The example below loads the output VTK data using Python [`pyvista`](https://www.pyvista.org/) library.

```python
# Python code
import pyvista as pv

p = pv.Plotter(notebook=False)
p.add_mesh( pv.read('out_blks.vtu'), point_size=5, opacity=0.15, color='white' )
p.add_mesh( pv.read('out_dh.vtu'), point_size=8, color='red' )
p.show()
```

<p align="center">
  <img src="docs/images/result.png">
</p>


## Unfolding documentation

The documentation of the main functions are available as [docstrings](https://juliahub.com/docs/Unfolding)

[build-img]: https://img.shields.io/github/actions/workflow/status/rmcaixeta/Unfolding.jl/CI.yml?branch=master
[build-url]: https://github.com/rmcaixeta/Unfolding.jl/actions

[codecov-img]: https://codecov.io/gh/rmcaixeta/Unfolding.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/rmcaixeta/Unfolding.jl
