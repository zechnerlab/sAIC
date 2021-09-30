#=
Definition of figures used in the paper, so that they are reproducible.

Uses Plots with pyplot backend.
=#

"""
    FigureResults(...)

Wrapper around various simulation result data to be serialized/deserialized 
and used for figures.
"""
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
| k_5       | k₅    |
| k_div     | k_div |
| k_0       | k₀    |
"""
@export function MassControl(; omega=10.0, theta=0.01, eta=0.005,
                                k_1=0.1, k_a=1.0/100, k_2=1.0,
                                k_3=2.0, k_4=1.0,
                                k_5=4.0, k_div=π/4000, k_0=3e-4,
                                durations=[[200.0;200;200], [300.0;300.0]], 
                                changes=[[1.0,4.0,0.5], [1.0,10.0]],
                                N_SSA::Int64=100, 
                                seed::Union{Nothing,Int64}=nothing,
                                N_Cells::Int=50,
                                fillalpha::AbstractFloat=0.4,
                                outdir::String="Figures",
                                fname::Union{Nothing,String}="MassControl_x$(N_SSA)", # The filename should NOT contain the extension
                                dumpfname::Union{Nothing,String}=fname, # The filename should NOT contain the extension
                                readFromDump::Bool=false,
                                figSize=(650*1.16,950*0.83), # (W, H) size for the figure,
                                fontsizeScale=1.6, # to be used with scalefontsizes(...)
                                plotTrajectory::Bool=true,
                                )
    mkpath(outdir) # Make sure the output directory exists
    fname = outdir*"/"*fname
    dumpfname = outdir*"/"*dumpfname
    res = nothing

    # Compute or load the results
    if readFromDump
        # Read the simulation output from a serialized dump
        @assert !isnothing(dumpfname) "FATAL: Please provide a dump filename (dumpfname=\"foobar.jser\")"
        res = deserialize(dumpfname*"_1.jser")
        res2 = deserialize(dumpfname*"_2.jser")
        @unpack t_ssa, mom_ssa, ts, vs, time_all, Moms, Vars, t_setpoint, val_setpoint, n = res2
        t_ssa2, mom_ssa2, ts2, vs2, time_all2, Moms2, Vars2, t_setpoint2, val_setpoint2, n2 = t_ssa, mom_ssa, ts, vs, time_all, Moms, Vars, t_setpoint, val_setpoint, n
        @unpack t_ssa, mom_ssa, ts, vs, time_all, Moms, Vars, t_setpoint, val_setpoint, n = res
    else
        # Run a new set of SSA simulations and serialize its results
        seed!=nothing ? Random.seed!(seed) : nothing
        n0 = zeros(Int64, 6, N_Cells+1) # Rows=species, Cols=cells+environment
        n0[1,2:end] .= 1 # cells
        S = Models.MassControl(omega=omega, theta=theta, eta=eta,
                                 k_1=k_1, k_a=k_a, k_2=k_2,
                                 k_3=k_3, k_4=k_4,
                                 k_5=k_5,
                                 k_div=k_div, k_0=k_0)
        println("> MassControl::Panel(a)::SingleTrajectory...")
        @time t_ssa, mom_ssa, ts, vs, _ = SSA_perturbations(deepcopy(S), n0, 
                                                durations[1], changes[1],
                                                )
        println("> MassControl::Panel(a)::Averages...")
        @time time_all, Moms, Vars, t_setpoint, val_setpoint, n = SSA_perturbations(deepcopy(S),
                                                                                n0,
                                                                                durations[1],
                                                                                changes[1],
                                                                                N_SSA;
                                                                                exportRawOutput=true)
        S2 = Models.MassControl(omega=2*omega, theta=theta, eta=eta,
                                 k_1=k_1, k_a=k_a, k_2=k_2,
                                 k_3=k_3, k_4=k_4,
                                 k_5=k_5,
                                 k_div=k_div, k_0=k_0)
        println("> MassControl::Panel(b)::SingleTrajectory...")
        @time t_ssa2, mom_ssa2, ts2, vs2, _ = SSA_perturbations(deepcopy(S2), n0, 
                                                durations[2], changes[2],
                                                transitionClassIndex=9, # Cell division
                                                )
        println("> MassControl::Panel(b)::Averages...")
        @time time_all2, Moms2, Vars2, t_setpoint2, val_setpoint2, n2 = SSA_perturbations(deepcopy(S2),
                                                                n0,
                                                                durations[2],
                                                                changes[2],
                                                                N_SSA;
                                                                transitionClassIndex=9, # Cell division
                                                                exportRawOutput=true)

        res = FigureResults(t_ssa, mom_ssa, ts, vs, time_all, Moms, Vars, t_setpoint, val_setpoint, n)
        res2 = FigureResults(t_ssa2, mom_ssa2, ts2, vs2, time_all2, Moms2, Vars2, t_setpoint2, val_setpoint2, n2)
        # Save the results
        if !isnothing(dumpfname)
            println("> MassControl::Serializing...")
            serialize(dumpfname*"_1.jser", res)
            serialize(dumpfname*"_2.jser", res2)
        end
    end

    closeall() # Close any existing figure window

    scalefontsizes(fontsizeScale) # Scale all the plot fonts by the given factor

    legendPos = :topleft
    legendPos2 = :topright
    chemTicks = 1000:1000:4000
    cellTicks(X) = (M = round(1.2*maximum(X), sigdigits=1); s = round(M/3, sigdigits=1); s:s:M)

    pTot = nothing
    try
        local lmargin, tmargin, rmargin = 5mm, 5mm, 18mm
        local tx, ty = -27mm, 4.5mm

        # Panel 1: Averages
        p1 = plot(; left_margin=lmargin, top_margin=tmargin, right_margin=rmargin)
        ## Left y axis
        ### <Z1>
        plot!(time_all, Moms[2,:], ribbon=( sqrt.(Vars[2,:]), sqrt.(Vars[2,:]) ),
                                fillalpha=fillalpha,
                                color="magenta", legend=legendPos, label=L"\langle Z_1 \rangle")
        ### ssa(Z1)
        plotTrajectory && plot!(t_ssa, mom_ssa[2,:], linestyle=:dot, linealpha=0.9*0.7,
                                color="magenta", legend=legendPos, label=false)
        ### <Q2>
        plot!(time_all, Moms[6,:], ribbon=( sqrt.(Vars[6,:]), sqrt.(Vars[6,:]) ),
                                     fillalpha=fillalpha,
                                     color="green",
                                     ylim=(-Inf,Inf),
                                     xlim=(0.0, Inf),
                                     legend=legendPos, label=L"\langle Q_2 \rangle")
        ### ssa(Q2)
        plotTrajectory && plot!(t_ssa, mom_ssa[6,:],
                                     linestyle=:dot, linealpha=0.7,
                                     color="green",
                                     ylim=(-Inf,Inf),
                                     xlim=(0.0, Inf),
                                     legend=legendPos, label=false)
        ### setpoint
        plot!(t_setpoint, val_setpoint, color="black", linestyle=:dash, legend=legendPos, 
                label=L"\omega^*")
        
        yticks!(chemTicks)
        ylabel!("Abundance")
        
        ## Right y axis
        p1y2 = twinx(p1)
        ### <N>
        plot!(p1y2, time_all, Moms[4,:], ribbon=( sqrt.(Vars[4,:]), sqrt.(Vars[4,:]) ),
                                 fillalpha=fillalpha,
                                 color="black",
                                 label=L"\langle N \rangle",
                                 legend=legendPos2,
                                 ylim=(0.0,5.5*maximum(Moms[4,:])),
                                 xlim=(0.0, Inf),
                                 xticks=false, # Avoid drawing xticks twice
                                 )
        ### ssa(N)
        plotTrajectory && plot!(p1y2, t_ssa, mom_ssa[4,:],
                                 color="#585858",
                                 label=false,
                                 linestyle=:dot, linealpha=0.7,
                                 legend=legendPos2,
                                 ylim=(0.0,5.5*maximum(mom_ssa[4,:])),
                                 xlim=(0.0, Inf),
                                 xticks=false, # Avoid drawing xticks twice
                                 )

        yticks!(p1y2, [100,300,600,900,1200])
        ylabel!(p1y2, "Number of cells")
        
        ## Panel label
        plot!(; # Float title
            title=L"(a)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=3,
            bg_inside=nothing,
            )

        # Panel 2: One SSA trajectory, perturbation of the system parameters
        p2 = plot(; left_margin=lmargin, top_margin=tmargin, right_margin=rmargin)
        ## Left y axis
        ### <Z1>
        plot!(time_all2, Moms2[2,:], ribbon=( sqrt.(Vars2[2,:]), sqrt.(Vars2[2,:]) ),
                                fillalpha=fillalpha,
                                color="magenta", legend=legendPos, label=L"\langle Z_1 \rangle")
        ### <Q2>
        plot!(time_all2, Moms2[6,:], ribbon=( sqrt.(Vars2[6,:]), sqrt.(Vars2[6,:]) ),
                                     fillalpha=fillalpha,
                                     color="green",
                                     ylim=(-Inf,1.25*maximum(Moms2[6,:])),
                                     xlim=(0.0, Inf),
                                     legend=legendPos, label=L"\langle Q_2 \rangle")
        ### ssa(Z1)
        plotTrajectory && plot!(t_ssa2, mom_ssa2[2,:], 
                                    color="magenta", 
                                    label=false, 
                                    linestyle=:dot, linealpha=0.9*0.7,
                                    lw=1., 
                                    legend=legendPos)
        ### ssa(Q2)
        plotTrajectory && plot!(t_ssa2, mom_ssa2[6,:], 
                                    color="green", 
                                    label=false, 
                                    lw=1., 
                                    linestyle=:dot, linealpha=0.7,
                                    legend=legendPos,
                                    ylim=(-Inf,1.25*maximum(mom_ssa2[6,:])),
                                    xlim=(0.0, Inf),
                                    )
        ### setpoint
        plot!(ts2, vs2, color="black", linestyle=:dash, label=L"\omega^*", legend=legendPos)
        
        xlabel!("Time")
        yticks!(chemTicks)
        ylabel!("Abundance")
        
        ## Right y axis
        p2y2 = twinx(p2)
        ### <N>
        plot!(p2y2, time_all, Moms2[4,:], ribbon=( sqrt.(Vars2[4,:]), sqrt.(Vars2[4,:]) ),
                                 fillalpha=fillalpha,
                                 color="black",
                                 label=L"\langle N \rangle",
                                 legend=legendPos2,
                                 ylim=(0.0,3.5*maximum(Moms2[4,:])),
                                 xlim=(0.0, Inf),
                                 xticks=false, # Avoid drawing xticks twice
                                 )
        ### ssa(N)
        plotTrajectory && plot!(p2y2, t_ssa2, mom_ssa2[4,:], color="#585858", 
                                        # label=L"N", 
                                        label=false, 
                                        lw=1.,
                                        linestyle=:dot, linealpha=0.7,
                                        legend=legendPos2,
                                        ylim=(0.0,3.5*maximum(mom_ssa2[4,:])),
                                        xlim=(0.0, Inf),
                                        xticks=false, # Avoid drawing xticks twice
                                        )
        ## Shaded gray background
        vspan!(p2y2, cumsum(durations[2])[1:2], color="#585858", alpha=0.1, label=false)
        
        yticks!(p2y2, 300:600:2800)
        ylabel!(p2y2, "Number of cells")
        
        ## Panel label
        plot!(; # Float title
            title=L"(b)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=3,
            bg_inside=nothing,
            )

        # Compose the two panels into the whole figure
        pTot = plot(p1, p2, layout=(2,1), size=figSize)

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
                kwargs...
        ) = MassControl(; 
                k_div=0.0, k_0=0.0,
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
@export function NumberControl(;
                                 omega=60.0, theta=0.4, eta=0.005,
                                 k_1=0.1, k_a=1.0/100, k_2=1.0,
                                 k_3=2.0, k_4=1.0,
                                 k_div=π/4000, k_0=3e-4,
                                 durations=[[500.0;400;400], [650.0,650.0]], 
                                 changes=[[1.,2.,2/3], [1.0, 1/5.0]],
                                 N_SSA::Int64=100, 
                                 seed::Union{Nothing,Int64}=nothing,
                                 N_Cells::Int=50,
                                 fillalpha::AbstractFloat=0.4,
                                 outdir::String="Figures",
                                 fname::Union{Nothing,String}="NumberControl_x$(N_SSA)", # The filename should NOT contain the extension
                                 dumpfname::Union{Nothing,String}=fname, # The filename should NOT contain the extension
                                 readFromDump::Bool=false,
                                 figSize=(650*1.16,850*0.83), # (W, H) size for the figure
                                 fontsizeScale=1.6, # to be used with scalefontsizes(...)
                                 plotTrajectory::Bool=true,
                              )
    mkpath(outdir) # Make sure the output directory exists
    fname = outdir*"/"*fname
    dumpfname = outdir*"/"*dumpfname
    res = nothing

    # Compute or load the results
    if readFromDump
        # Read the simulation output from a serialized dump
        @assert !isnothing(dumpfname) "FATAL: Please provide a dump filename (dumpfname=\"foobar.jser\")"
        res = deserialize(dumpfname*"_1.jser")
        res2 = deserialize(dumpfname*"_2.jser")
        @unpack t_ssa, mom_ssa, ts, vs, time_all, Moms, Vars, t_setpoint, val_setpoint, n = res2
        t_ssa2, mom_ssa2, ts2, vs2, time_all2, Moms2, Vars2, t_setpoint2, val_setpoint2, n2 = t_ssa, mom_ssa, ts, vs, time_all, Moms, Vars, t_setpoint, val_setpoint, n
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
        println("> NumberControl::Panel(a)::SingleTrajectory...")
        @time t_ssa, mom_ssa, ts, vs, _ = SSA_perturbations(deepcopy(S), n0, 
                                                durations[1], changes[1],
                                                )
        println("> NumberControl::Panel(a)::Averages...")
        @time time_all, Moms, Vars, t_setpoint, val_setpoint, n = SSA_perturbations(deepcopy(S),
                                                                                n0,
                                                                                durations[1],
                                                                                changes[1],
                                                                                N_SSA;
                                                                                exportRawOutput=true)
        S2 = Models.NumberControl(omega=2*omega, theta=theta, eta=eta,
                                 k_1=k_1, k_a=k_a, k_2=k_2,
                                 k_3=k_3, k_4=k_4,
                                 k_div=k_div, k_0=k_0)
        println("> NumberControl::Panel(b)::SingleTrajectory...")
        @time t_ssa2, mom_ssa2, ts2, vs2, _ = SSA_perturbations(deepcopy(S2), n0, 
                                                durations[2], changes[2],
                                                transitionClassIndex=7, # Prod X2 (k3)
                                                )
        println("> NumberControl::Panel(b)::Averages...")
        @time time_all2, Moms2, Vars2, t_setpoint2, val_setpoint2, n2 = SSA_perturbations(deepcopy(S2),
                                                                n0,
                                                                durations[2],
                                                                changes[2],
                                                                N_SSA;
                                                                transitionClassIndex=7, # Prod X2 (k3)
                                                                exportRawOutput=true)

        res = FigureResults(t_ssa, mom_ssa, ts, vs, time_all, Moms, Vars, t_setpoint, val_setpoint, n)
        res2 = FigureResults(t_ssa2, mom_ssa2, ts2, vs2, time_all2, Moms2, Vars2, t_setpoint2, val_setpoint2, n2)
        # Save the results
        if !isnothing(dumpfname)
            println("> NumberControl::Serializing...")
            serialize(dumpfname*"_1.jser", res)
            serialize(dumpfname*"_2.jser", res2)
        end
    end

    closeall() # Close any existing figure window

    scalefontsizes(fontsizeScale) # Scale all the plot fonts by the given factor

    legendPos = :topleft
    legendPos2 = :topright
    chemTicks(X) = (M = round(1.1*maximum(X), sigdigits=1); s = round(M/4, sigdigits=1); s:s:M)
    cellTicks(X) = (M = round(1.2*maximum(X), sigdigits=1); s = round(M/4, sigdigits=1); s:s:M)

    pTot = nothing
    try
        local lmargin, tmargin, rmargin = 5mm, 5mm, 0mm
        local tx, ty = -25mm, 4.5mm

        # Panel 1: Averages
        p1 = plot(; left_margin=lmargin, top_margin=tmargin, right_margin=rmargin)
        ## Left y axis
        ### <N>
        plot!(time_all, Moms[4,:], ribbon=( sqrt.(Vars[4,:]), sqrt.(Vars[4,:]) ),
                                 fillalpha=fillalpha,
                                 color="black",
                                 label=L"\langle N \rangle",
                                 # legend=:topright,
                                 legend=legendPos,
                                 ylim=(0.0, 1.5*maximum(Moms[4,:]) ),
                                 xlim=(0.0, Inf),
                                 )
        ### ssa(N)
        plotTrajectory && plot!(t_ssa, mom_ssa[4,:],
                                 linestyle=:dot, linealpha=0.7,
                                 color="#585858",
                                 label=false,
                                 # legend=:topright,
                                 legend=legendPos,
                                 ylim=(0.0, 1.5*maximum(mom_ssa[4,:]) ),
                                 xlim=(0.0, Inf),
                                 )
        ### setpoint
        plot!(t_setpoint, val_setpoint, color="black",
                                  linestyle=:dash,
                                  label=L"\omega^*",
                                  # legend=:topright,
                                  legend=legendPos,
                                  # ylim=(0.0,Inf),
                                  )
        
        yticks!(cellTicks(Moms[4,:]))
        ylabel!("Number of cells";
                  # right_margin=10mm,
                  )

        ## Right y axis
        p1y2 = twinx(p1)
        ### <Z1>
        plot!(p1y2, time_all, Moms[2,:], ribbon=( sqrt.(Vars[2,:]), sqrt.(Vars[2,:]) ),
                                 fillalpha=0.7*fillalpha,
                                 color="magenta",
                                 guidefontcolor="magenta",
                                 lw=1.0,
                                 label=L"\langle Z_1 \rangle",
                                 legend=legendPos2,
                                 xticks=false,
                                 ylim=(0.0,1.0*( maximum(Moms[2,:]) + maximum(sqrt.(Vars[2,:])) ) ),
                                 xlim=(0.0, Inf),
                                 formatter=:plain,
                                 # yformatter=:scientific,
                                 )
        ### ssa(Z1)
        plotTrajectory && plot!(p1y2, t_ssa, mom_ssa[2,:],
                                 linestyle=:dot, linealpha=0.9*0.7,
                                 color="magenta",
                                 guidefontcolor="magenta",
                                 lw=1.0,
                                 label=false,
                                 legend=legendPos2,
                                 xticks=false,
                                 ylim=(0.0,1.0*( maximum(Moms[2,:]) + maximum(sqrt.(Vars[2,:])) ) ),
                                 xlim=(0.0, Inf),
                                 formatter=:plain,
                                 )
        
        yticks!(p1y2, chemTicks(Moms[2,:]), formatter=:plain)
        ylabel!(p1y2, "Abundance")

        ## Panel label
        plot!(; # Float title
            title=L"(a)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=3,
            bg_inside=nothing,
            )

        # Panel 2: One SSA trajectory, perturbation of the system parameters
        p2 = plot(; left_margin=lmargin, top_margin=tmargin, right_margin=rmargin)
        ## Left y axis
        ### <N>
        plot!(time_all2, Moms2[4,:], ribbon=( sqrt.(Vars2[4,:]), sqrt.(Vars2[4,:]) ),
                                 fillalpha=fillalpha,
                                 color="black",
                                 label=L"\langle N \rangle",
                                 legend=legendPos,
                                 ylim=(0.0, 1.4*(maximum(Moms2[4,:]) + maximum(sqrt.(Vars2[4,:]))) ),
                                 xlim=(0.0, Inf),
                                 )
        ### ssa(N)
        plotTrajectory && plot!(t_ssa2, mom_ssa2[4,:], color="#585858", 
                                          # label=L"N", 
                                          label=false, 
                                          lw=1.,
                                          linestyle=:dot, linealpha=0.7,
                                          legend=legendPos,
                                          ylim=(0.0, 1.4*maximum(mom_ssa2[4,:]) ),
                                          xlim=(0.0, Inf),
                                          )
        ### setpoint
        plot!(ts2, vs2, color="black", linestyle=:dash, label=L"\omega^*",
                                          legend=legendPos,
                                          )

        xlabel!("Time")
        yticks!(cellTicks(mom_ssa2[4,:]))
        ylabel!("Number of cells")
        
        ## Right y axis
        p2y2 = twinx(p2)
        ### <Z1>
        plot!(p2y2, time_all2, Moms2[2,:], ribbon=( sqrt.(Vars2[2,:]), sqrt.(Vars2[2,:]) ),
                                 fillalpha=0.7*fillalpha,
                                 color="magenta",
                                 guidefontcolor="magenta",
                                 lw=1.0,
                                 label=L"\langle Z_1 \rangle",
                                 legend=legendPos2,
                                 xticks=false,
                                 ylim=(0.0,1.5*( maximum(Moms2[2,:]) + maximum(sqrt.(Vars2[2,:])) ) ),
                                 xlim=(0.0, Inf),
                                 formatter=:plain,
                                 )
        ### ssa(Z1)
        plotTrajectory && plot!(p2y2, t_ssa2, mom_ssa2[2,:];
                      color="magenta",
                      guidefontcolor="magenta",
                      ylim=(0.0, 1.5*maximum(mom_ssa2[2,:]) ),
                      label=false,
                      lw=1.,
                      linestyle=:dot, linealpha=0.9*0.7,
                      legend=legendPos2,
                      xticks=false,
                      xlim=(0.0, Inf),
                      formatter=:plain,
                      )

        ## Shaded gray background
        vspan!(p2y2, cumsum(durations[2])[1:2], color="#585858", alpha=0.1, label=false)

        yticks!(p2y2, chemTicks(mom_ssa2[2,:]), formatter=:plain)
        ylabel!(p2y2, "Abundance")

        ## Panel label
        plot!(; # Float title
            title=L"(b)",
            titlepos=:left,
            grid=false, showaxis=false, ticks=false,
            inset=(1, bbox(tx,ty, 20mm,20mm, :top, :left)),
            subplot=3,
            bg_inside=nothing,
            )

        # Compose the two panels into the whole figure
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

#eof
