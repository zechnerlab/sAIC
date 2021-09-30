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
using SAIC
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
@everywhere using SAIC	# Repeat this if it fails on the first attempt
```

Note that depending on the Julia version and if Revise is in use, the 
`@everywhere using SAIC` command may fail throwing errors related to compilecache.
It is usually sufficient to run `@everywhere using SAIC` a second time to succeed
in loading the package on all workers.

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
@everywhere using SAIC	# Repeat this if it fails on the first attempt

Figures.MassControl()
Figures.NumberControl()
```

The simulation results are also serialized and saved into `*.jser` files, so that 
tweaking of the plots does not require running the SSAs all over again.  
To load the results from these saves, use the `loadFromDump=true` flag:
```julia
Figures.MassControl(; loadFromDump=true)
Figures.NumberControl(; loadFromDump=true)
```
The figure generation is a sequential operation, so it is fine to launch julia
as a single process (see Sequential execution).
