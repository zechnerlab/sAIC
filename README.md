# Shared Antithetic Integral Controller

Julia code for shared antithetic integral control in dynamic cell populations

## Usage
### Dependencies
This repository already contains `Project.toml` and `Manifest.toml` files specifying
is dependencies in the form of a `Pkg` environment.  
To use it:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()
```
or, in the Julia REPL:
```julia
] activate .
] instantiate
] precompile
```
This will download, install and precompile all the dependencies.

### General usage and project structure
The code is structured as a main module `SAIC` which directly exposes
the main data structures and SSA features and additional submodules `Models` and `Figures`.

#### Sequential execution
It is possible to run the trajectories sequentially by running (in the repository directory):
```
julia --project
```
```julia
import SAIC
```

#### Parallel execution
Running SSA trajectories in parallel is supported through Julia's `Distributed` 
features.  
Please note that it requires specifying the number of worker processes when launching Julia:
```
julia -p auto
```
Running simulations in parallel then requires the following lines of code to be executed
before running any simulation:
```julia
using Distributed
@everywhere import Pkg
@everywhere Pkg.activate(".")
import SAIC
```

### Generating the paper figures
The functions generating the figures for the paper are available in the 
`SAIC.Figures` submodule.  
They can be generated with the following:
```
julia -p auto
```
and
```julia
using Distributed
@everywhere import Pkg
@everywhere Pkg.activate(".")
import SAIC

SAIC.Figures.MassControl()
SAIC.Figures.NumberControl()
```

The simulation results are also serialized and saved into `*.jser` files, so that 
tweaking of the plots does not require running the SSAs all over again.  
To load the results from these saves, use the `loadFromDump=true` flag:
```julia
SAIC.Figures.MassControl(; loadFromDump=true)
SAIC.Figures.NumberControl(; loadFromDump=true)
```
The figure generation is a sequential operation, so it is fine to launch julia
as a single process (see Sequential execution).
