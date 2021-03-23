module SAIC

using Distributions
using StatsBase
using Random
using InlineExports			# Provides the @export macro
using Distributed

export Models, Figures

include("SAIC_structs.jl")
include("SAIC_SSA.jl")

module Models
	using ..SAIC	# Use exported function from the parent module
	using InlineExports			# Provides the @export macro
	using Random, Distributions
	include("SAIC_models.jl")
end # submodule

module Figures
	using ..SAIC	# Use exported function from the parent module
	using InlineExports		# Provides the @export macro
	using LaTeXStrings		# For using LaTeX in plot labels (e.g. L"X^2")
	using Measures			# Required to use `mm` or `pt` as units for plot margins
	import KernelDensity 	# Required for KernelDensityEstimation
	using Serialization 	# Required for saving/loading simulation results
	using Parameters 		# Provides the @unpack macro
	using Plots; pyplot()	# This first pyplot() call is required as it needs to
							# precompile some generated functions that will be
							# required when compiling __init__()

	"""
		__init__()

	Force the use of TrueType fonts in plots (ref: http://phyletica.org/matplotlib-fonts/ )
	This is to avoid errors because of Type 3 fonts.
	"""
	function __init__()
		pyplot()	# This second pyplot() call is required to actually load the
					# pyplot backend at runtime
		Plots.pyrcparams["pdf.fonttype"] = 42
		Plots.pyrcparams["ps.fonttype"] = 42
		Plots.pyrcparams["axes.labelpad"] = 10.0 # 4.0 is the default
	end

	include("SAIC_figures.jl")
end # submodule

end  # module
