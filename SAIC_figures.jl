#=
Definition of figures used in the paper, so that they are reproducible.

Uses Plots with pyplot backend.
=#

"""
    MassControl()

Executing this function without arguments runs the simulations and then uses their
results to assemble the MassControl figure of the paper, with the same parameters.

Translation table between the parameter names used in the code vs. in the paper:
| Code      | Paper |
| --------- | ----- |
| omega     | ω     |
| theta     | θ     |
| eta       | η     |
| k_a       | kₐ    |
| k_1       | k₁    |
| k_2       | k₂    |
| k_3       | k₃    |
| k_4       | k₄    |
| k_div     | k_div |
| k_0       | k₀    |
"""
@export function MassControl(; omega=20.0, theta=0.02, eta=0.005,
                                             k_1=0.1, k_a=1.0, k_2=1.0,
                                             k_3=2.0, k_4=1.0,
                                             k_div=0.01, k_0=0.1,
                                             durations=[200.0;200;200], changes=[1.0,4.0,0.5],
                                             N_SSA::Int64=100, seed::Union{Nothing,Int64}=nothing,
                                             N_Cells::Int=20,
                                             pdfBandwidth::Union{AbstractFloat,Symbol}=:auto,   # This should be right, but can be adjusted
                                             fillalpha::AbstractFloat=0.4,
                                             outdir::String="Figures",
                                             fname::Union{Nothing,String}="MassControl_x$(N_SSA)", # The filename should NOT contain the extension
                                             dumpfname::Union{Nothing,String}=fname, # The filename should NOT contain the extension
                                             readFromDump::Bool=false,
                                             # figSize=(1280,920), # (W, H) size for the figure,
                                             figSize=(600,1150), # (W, H) size for the figure,
                                             fontsizeScale=1.6, # to be used with scalefontsizes(...)
                                          )
    mkpath(outdir)
    fname = outdir*"/"*fname
    dumpfname = outdir*"/"*dumpfname
    res = nothing
    # Compute or load the results
    if readFromDump
        # Read the simulation output from a serialized dump
        @assert !isnothing(dumpfname) "FATAL: Please provide a dump filename (dumpfname=\"foobar.jser\")"
        res = deserialize(dumpfname*".jser")
        @unpack t_ssa, mom_ssa, ts, vs, time_all, Moms, Vars, t_setpoint, val_setpoint, n = res
    else
        # Run a new set of SSA simulations and serialize its results
        seed!=nothing ? Random.seed!(seed) : nothing
        n0 = zeros(Int64, 5, N_Cells+1) # Rows=species, Cols=cells+environment
        n0[1,2:end] .= 1 # cells
        S = Models.MassControl(omega=omega, theta=theta, eta=eta,
                                 k_1=k_1, k_a=k_a, k_2=k_2,
                                 k_3=k_3, k_4=k_4,
                                 k_div=k_div, k_0=k_0)
        t_ssa, mom_ssa, ts, vs, _ = SSA_perturbations(deepcopy(S), n0, durations, changes)
        @time time_all, Moms, Vars, t_setpoint, val_setpoint, n = SSA_perturbations(deepcopy(S),
                                                                                n0,
                                                                                durations,
                                                                                changes,
                                                                                N_SSA;
                                                                                exportRawOutput=true)
        res = FigureResults(t_ssa, mom_ssa, ts, vs, time_all, Moms, Vars, t_setpoint, val_setpoint, n)
        # Save the results
        if !isnothing(dumpfname)
            serialize(dumpfname*".jser", res)
        end
    end

    closeall() # Close any existing figure window

    scalefontsizes(fontsizeScale) # Scale all the plot fonts by the given factor

    pTot = nothing
    try
        # Panel 1: One SSA trajectory
        p1 = plot(; title=L"(a)", titlepos=:left)
        plot!(t_ssa, mom_ssa[2,:], color="magenta", label=L"Z_1", lw=1., legend=:topleft)
        plot!(t_ssa, mom_ssa[3,:], color="orange", label=L"Z_2", lw=1., legend=:topleft)
        plot!(t_ssa, mom_ssa[6,:], color="green", label=L"Q_2", lw=1., legend=:topleft)
        plot!(t_setpoint, val_setpoint, color="black", linestyle=:dash, label=L"\omega^*", legend=:topleft)
        xlabel!("Time [a.u.]")
        ylabel!("Abundance")

        # Panel 2: Averages
        p2 = plot(; title=L"(b)", titlepos=:left)
        plot!(time_all, Moms[2,:], ribbon=( sqrt.(Vars[2,:]), sqrt.(Vars[2,:]) ),
                                fillalpha=fillalpha,
                                color="magenta", legend=:topleft, label=L"\langle Z_1 \rangle")
        plot!(time_all, Moms[3,:], ribbon=( sqrt.(Vars[3,:]), sqrt.(Vars[3,:]) ),
                                fillalpha=fillalpha,
                                color="orange",
                                legend=:topleft, label=L"\langle Z_2 \rangle")
        plot!(time_all, Moms[6,:], ribbon=( sqrt.(Vars[6,:]), sqrt.(Vars[6,:]) ),
                                     fillalpha=fillalpha,
                                     color="green",
                                     legend=:topleft, label=L"\langle Q_2 \rangle")
        plot!(t_setpoint, val_setpoint, color="black", linestyle=:dash, legend=:topleft, label=L"\omega^*")
        xlabel!("Time [a.u.]")
        ylabel!("Abundance")

        # Panel 3: Number of cells
        p3 = plot(; title=L"(c)", titlepos=:left)
        local N = Moms[4,:]
        local ylims = minimum(N)==maximum(N) ? (0,2*maximum(N)) : (0,Inf)
        plot!(time_all, N, ribbon=( sqrt.(Vars[4,:]), sqrt.(Vars[4,:]) ),
                                fillalpha=fillalpha,
                                color="black",
                                label=L"\langle N \rangle",
                                legend=:topleft,
                                ylim=ylims,
                                )
        xlabel!("Time [a.u.]")
        ylabel!("Number of cells")

        # Panel 4: pooled PDF of X2
        K = computePooledPDF(n, 5; bandwidth=pdfBandwidth)
        numCellsMean = Moms[4, :]
        expectedMean = val_setpoint[end] / numCellsMean[end]
        p4 = plot(; title=L"(d)", titlepos=:left)
        plot!(K.x, K.density; color="green",
                                label=L"pdf(Q_2)",
                                xlim=(max(0,minimum(K.x)),Inf),
                                )
        vline!([expectedMean]; color="black",
                                linestyle=:dash,
                                label=L"\frac{ ω^* }{ \langle N \rangle }")
        xlabel!("Protein mass per cell")
        ylabel!("Probability density")

        # # Panel 4b: summed PDF of X2
        # K = computeSummedPDF(n, 5; bandwidth=pdfBandwidth)
        # expectedMean = val_setpoint[end]
        # p4 = plot(; title=L"(d)", titlepos=:left)
        # plot!(K.x, K.density; color="green",
        #                           label=L"pdf(Q_2)",
        #                           xlim=(max(0,minimum(K.x)),Inf),
        #                           )
        # vline!([expectedMean]; color="black",
        #                           linestyle=:dash,
        #                           label=L"ω^*")
        # xlabel!("Protein mass")
        # ylabel!("Probability density")

        ### Compose the subplots together
        # pTot = plot(p1, p3, p2, p4, layout=(2,2), size=figSize, margin=5mm)
        pTot = plot(p1, p2, p3, layout=(3,1), size=figSize)

        if !isnothing(fname)
          savefig(fname*".pdf")
          savefig(fname*".png")
        end

    finally
        scalefontsizes(1/fontsizeScale) # Set it back to the default one
    end

    return pTot
end

"""
    MassControl_ConstantCellNumber()

Executing this function without arguments runs the simulations and then uses their
results to assemble a supplementary figure, not used in the paper.
It is the same as the MassControl figure, but forcing the number of cells to 
stay constant throughout the simulation.
"""
@export MassControl_ConstantCellNumber(; 
            N_SSA::Int64=100,
            fname::Union{Nothing,String}="MassControl_ConstCellNum_x$(N_SSA)",
            kwargs...) = MassControl(; k_div=0.0, k_0=0.0,
                                      N_SSA=N_SSA,
                                      N_Cells=100,
                                      fname=fname,
                                      kwargs...)

"""
    NumberControl()

Executing this function without arguments runs the simulations and then uses their
results to assemble the NumberControl figure of the paper, with the same parameters.

Translation table between the parameter names used in the code vs. in the paper:
| Code      | Paper |
| --------- | ----- |
| omega     | ω     |
| theta     | θ     |
| eta       | η     |
| k_a       | kₐ    |
| k_1       | k₁    |
| k_2       | k₂    |
| k_3       | k₃    |
| k_4       | k₄    |
| k_div     | k_div |
| k_0       | k₀    |
"""
@export function NumberControl(; omega=200.0, theta=4.0, eta=0.01,
                                 k_1=0.0, k_a=1.0, k_2=1.0,
                                 k_3=0.25, k_4=1.0,
                                 k_div=0.01, k_0=0.1,
                                 durations=[200.0;200;200], changes=[1.,6,2/3],
                                 N_SSA::Int64=100, seed::Union{Nothing,Int64}=nothing,
                                 N_Cells::Int=50,
                                 pdfBandwidth::Union{AbstractFloat,Symbol}=:auto,   # This should be right, but can be adjusted
                                 fillalpha::AbstractFloat=0.4,
                                 outdir::String="Figures",
                                 fname::Union{Nothing,String}="NumberControl_x$(N_SSA)", # The filename should NOT contain the extension
                                 dumpfname::Union{Nothing,String}=fname, # The filename should NOT contain the extension
                                 readFromDump::Bool=false,
                                 figSize=(600,850), # (W, H) size for the figure
                                 fontsizeScale=1.6, # to be used with scalefontsizes(...)
                              )
    mkpath(outdir)
    fname = outdir*"/"*fname
    dumpfname = outdir*"/"*dumpfname
    res = nothing
    # Compute or load the results
    if readFromDump
        # Read the simulation output from a serialized dump
        @assert !isnothing(dumpfname) "FATAL: Please provide a dump filename (dumpfname=\"foobar.jser\")"
        res = deserialize(dumpfname*".jser")
        @unpack t_ssa, mom_ssa, ts, vs, time_all, Moms, Vars, t_setpoint, val_setpoint, n = res
    else
        # Run a new set of SSA simulations and serialize its results
        seed!=nothing ? Random.seed!(seed) : nothing
        n0 = zeros(Int64, 5, N_Cells+1) # Rows=species, Cols=cells+environment
        n0[1,2:end] .= 1 # cells
        S = Models.NumberControl(omega=omega, theta=theta, eta=eta,
                                 k_1=k_1, k_a=k_a, k_2=k_2,
                                 k_3=k_3, k_4=k_4,
                                 k_div=k_div, k_0=k_0)
        t_ssa, mom_ssa, ts, vs, _ = SSA_perturbations(deepcopy(S), n0, durations, changes)
        @time time_all, Moms, Vars, t_setpoint, val_setpoint, n = SSA_perturbations(deepcopy(S),
                                                                                n0,
                                                                                durations,
                                                                                changes,
                                                                                N_SSA;
                                                                                exportRawOutput=true)
        res = FigureResults(t_ssa, mom_ssa, ts, vs, time_all, Moms, Vars, t_setpoint, val_setpoint, n)
        # Save the results
        if !isnothing(dumpfname)
            serialize(dumpfname*".jser", res)
        end
    end

    closeall() # Close any existing figure window

    scalefontsizes(fontsizeScale) # Scale all the plot fonts by the given factor

    pTot = nothing
    try
        # Panel 1: One SSA trajectory
        p1 = plot(; title=L"(a)", titlepos=:left)
        plot!(t_ssa, mom_ssa[4,:], color="#585858", label=L"N", lw=1.,
                                          # legend=:topright,
                                          legend=:topleft,
                                          ylim=(0.0,Inf),
                                          )
        plot!(t_setpoint, val_setpoint, color="black", linestyle=:dash, label=L"N^*",
                                          # legend=:topright,
                                          legend=:topleft,
                                          )
        ylabel!("Number of cells")
        xlabel!("Time [a.u.]")
        y2 = twinx(p1)
        plot!(y2, t_ssa, mom_ssa[2,:];
                      color="magenta",
                      guidefontcolor="magenta",
                      # ytickfontcolor="magenta",
                      ylim=(0.0,Inf),
                      label=L"Z_1",
                      lw=1.,
                      legend=:bottomright,
                      formatter=:plain,
                      # yformatter=:scientific,
                      )
        ## Don't show Z2 since it's not visible (always close to 0)
        # plot!(y2, t_ssa, mom_ssa[3,:];
        #                 color="orange",
        #                 ylim=(0.0,Inf),
        #                 label=L"Z_2",
        #                 lw=1.,
        #                 legend=:bottomright,
        #                 formatter=:plain,
        #                 )
        ylabel!(y2, "Abundance")

        # Panel 2: Averages
        p2 = plot(; title=L"(b)", titlepos=:left)
        plot!(time_all, Moms[4,:], ribbon=( sqrt.(Vars[4,:]), sqrt.(Vars[4,:]) ),
                                 fillalpha=fillalpha,
                                 color="black",
                                 label=L"\langle N \rangle",
                                 # legend=:topright,
                                 legend=:topleft,
                                 ylim=(0.0,Inf),
                                 )
        plot!(t_setpoint, val_setpoint, color="black",
                                  linestyle=:dash,
                                  label=L"N^*",
                                  # legend=:topright,
                                  legend=:topleft,
                                  ylim=(0.0,Inf),
                                  )
        xlabel!("Time [a.u.]")
        ylabel!("Number of cells";
                  # right_margin=10mm,
                  )
        y2 = twinx(p2)
        plot!(y2, time_all, Moms[2,:], ribbon=( sqrt.(Vars[2,:]), sqrt.(Vars[2,:]) ),
                                 fillalpha=fillalpha,
                                 color="magenta",
                                 guidefontcolor="magenta",
                                 # ytickfontcolor="magenta",
                                 lw=1.0,
                                 label=L"\langle Z_1 \rangle",
                                 legend=:bottomright,
                                 ylim=(0.0,Inf),
                                 formatter=:plain,
                                 # yformatter=:scientific,
                                 )
        ## Don't show Z2 since it's not visible (always close to 0)
        # plot!(y2, time_all, Moms[3,:], ribbon=( sqrt.(Vars[3,:]), sqrt.(Vars[3,:]) ),
        #                            fillalpha=fillalpha,
        #                            color="orange",
        #                            lw=1.0,
        #                            label=L"\langle Z_2 \rangle",
        #                            legend=:bottomright,
        #                            ylim=(0.0,Inf),
        #                            formatter=:plain,
        #                            )

        ylabel!(y2, "Abundance")

        # Panel 4: PDF of numCells across trajectories
        K = computeSummedPDF(n, 1; bandwidth=pdfBandwidth)
        expectedMean = val_setpoint[end]
        p3 = plot(; title=L"(c)", titlepos=:left)
        plot!(K.x, K.density; color="green",
                                label=L"pdf(N)",
                                xlim=(max(0,minimum(K.x)),Inf),
                                )
        vline!([expectedMean]; color="black",
                                linestyle=:dash,
                                label=L"N^*")
        xlabel!("Number of cells")
        ylabel!("Probability density")

        ### Compose the subplots together
        pTot = plot(p1, p2, layout=(2,1), size=figSize, right_margin=25mm)

        if !isnothing(fname)
          savefig(fname*".pdf")
          savefig(fname*".png")
        end

    finally
        scalefontsizes(1/fontsizeScale) # Set it back to the default one
    end

    return pTot
end

### Auxiliary structs and functions

"""
    computePooledPDF(n, species; bandwidth=:auto)

Compute the PDF of the distribution of values of `species` pooling all cells
from all trajectories in `n`.

Return: K with fields K.x and K.density ready to be plotted.
"""
function computePooledPDF(n::Array{Array{Int64,2},1}, species::Int;
                                bandwidth=:auto,
                                avoidEnvironment::Bool=true,
                                trajectoryTransform=identity)
    n_slice = []
    if avoidEnvironment
        n_slice = [ trajectoryTransform(x[species, 2:end]) for x in n ]
    else
        n_slice = [ trajectoryTransform(x[species, :]) for x in n ]
    end
    pooledValues = reshape(vcat(n_slice...), :)
    # @show pooledValues #debug
    K = nothing
    if bandwidth==:auto
        K = KernelDensity.kde(pooledValues,
                                # boundary=(minimum(pooledValues), maximum(pooledValues)),
                                )
    else
        K = KernelDensity.kde(pooledValues,
                                # boundary=(minimum(pooledValues), maximum(pooledValues)),
                                bandwidth=bandwidth)
    end
    return K
    # return x->KernelDensity.pdf(K, x), extrema(pooledValues)
end

"""
    computeSummedPDF(args...; kwargs...) = computePooledPDF(args...; trajectoryTransform=sum, kwargs...)

Compute the PDF of the distribution of the total mass of `species` across all
trajectories in `n`.
"""
computeSummedPDF(args...; kwargs...) = computePooledPDF(args...; trajectoryTransform=sum, kwargs...)

struct FigureResults
    t_ssa
    mom_ssa
    ts
    vs
    time_all
    Moms
    Vars
    t_setpoint
    val_setpoint
    n
end

#eof
