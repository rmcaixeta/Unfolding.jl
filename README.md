# Unfolding.jl

[![Build Status](https://travis-ci.com/rmcaixeta/Unfolding.jl.svg?branch=master)](https://travis-ci.com/rmcaixeta/Unfolding.jl)
[![Coverage](https://codecov.io/gh/rmcaixeta/Unfolding.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rmcaixeta/Unfolding.jl)

`Unfolding.jl` is a package written in Julia to perform unfolding of 3-D geometries. It was developed for geostatistical cases where complex 3-D domains are modeled and need to be unfolded for appropriate variography, estimations and simulations.

Julia was used due to its high performance and easy coding. This package was successfully tested with some big mining datasets but is still under development, so please enter in contact if you have some issue or feel free to contribute to the code (it is open source!).

## Installation

It is necessary to install Julia to run this code. Installation instructions for Windows, Linux and macOS are available [here](https://julialang.org/downloads/platform/).

After that, it is necessary to install the Unfolding package. Open a terminal, type `julia` to open the REPL and then install the package with the following command:

```julia
using Pkg; Pkg.add("Unfolding")
```

## Usage

There a two possible workflows:

* Unfolding based on a block model
* Unfolding based on a points of a reference surface

The example below is based on a block model within a folded domain. It is only necessary the X, Y and Z coordinates of the block centroids. Optionally, samples coordinates can also be informed.

`FIGURE HERE
(blocks 2x2x2m; blocks outside the domain are not informed)`
![](docs/images/result.png)


### Julia example

The usage in Julia is detailed in the code below.

```julia
# Julia code
using CSV
using Unfolding

# Reading data
df_samp = CSV.read("samples.csv")
df_block = CSV.read("block_model.csv")

# Get coordinate points as matrix
input_block = coordinate_matrix( df_block, columns=["XC","YC","ZC"] )
input_samp = coordinate_matrix( df_samp, columns=["X","Y","Z"] )

# Get reference surface points for unfolding
ref_surface = ref_surface_from_blocks(input_block)
# Get transformed coordinates of blocks and samples after unfolding
unf_block, unf_samp = unfold(ref_surface,input_block,input_samp)

# Write new XT, YT and ZT columns with the transformed coordinates
for (i,c) in enumerate([:XT,:YT,:ZT])
    df_samp[:,c] = unf_samp[i,:]
    df_block[:,c] = unf_block[i,:]
end

# Write output to CSV
CSV.write( "out_dh", df_samp )
CSV.write( "out_blks", df_block )

# Write output to VTK format
data_to_vtk(unf_block,"out_blks")
data_to_vtk(unf_samp,"out_dh")
```

The code can be saved in a textfile with `.jl` extension and be called in a terminal: `julia file.jl`. Or you can organize it in Jupyter notebooks (see instructions [here](https://github.com/JuliaLang/IJulia.jl)).

### Python example

Julia is not a widespread language yet. For those more familiar with Python, the Julia code can be called inside Python scripts.

Julia still need to be installed before calling it in Python. Additionally, an extra Python library must be installed in this case:

```python
pip install julia
```

The data for the unfolding functions should be informed as Numpy arrays.

```python
# This part is only necessary for the first run; after that, this part can be deleted
import julia
julia.install()

# Python example
from julia import Julia
Julia(compiled_modules=False) # excluding this line makes it run faster; but may crash in some systems

import pandas as pd
from julia import Unfolding as unf

# Reading data
df_block = pd.read_csv("block_model.csv",usecols=["XC","YC","ZC"])
df_samp = pd.read_csv("samples.csv",usecols=["X","Y","Z"])

# Get coordinate points as matrix
input_block = df_block.to_numpy().T
input_samp = df_samp.to_numpy().T

# Get reference surface points for unfolding
ref_surface = unf.ref_surface_from_blocks(input_block)
# Get transformed coordinates of blocks and samples after unfolding
unf_block, unf_samp = unf.unfold(ref_surface,input_block,input_samp)

# Write new XT, YT and ZT columns with the transformed coordinates
for i,c in enumerate(["XT","YT","ZT"]):
    df_block[:,c] = unf_block[i,:]
    df_samp[:,c] = unf_samp[i,:]

# Write output to CSV
df_block.to_csv("out_blks.csv",index=False)
df_samp.to_csv("out_samp.csv",index=False)

# Write output to VTK format
unf.data_to_vtk(unf_block,"out_blks")
unf.data_to_vtk(unf_samp,"out_dh")
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

Three main functions are available for use:

### Main functions
#### ref_surface_from_blocks()
a
#### unfold()
a

### Data Handling
#### coordinate_matrix()
a
#### data_to_vtk()
a
#### data_to_csv()
a
#### unfold_error_ids()
a
#### unfold_error_dists()
a

### Extra
#### landmark_isomap()
a
