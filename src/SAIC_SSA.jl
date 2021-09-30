#=
    SAIC_SSA.jl
    This module contains all the SSA machinery to simulate SAIC systems.
=#

"""
    SSA_perturbations(S, n0, durations, changes, Nsamples=100;
                             timestep=0.1, seed=nothing)

Run `Nsamples` SSA trajectories with parameter `changes` (perturbations) at the
given `durations` timepoints.
"""
@export function SSA_perturbations(S::System,
                                    n0::Matrix{Int64},
                                    durations::Vector{Float64},
                                    changes::Vector{Float64},
                                    Nsamples::Int64;
                                    transitionClassIndex::Int=1,
                                    timestep::Float64=0.1,
                                    seed::Union{Nothing,Int64}=nothing,
                                    exportRawOutput::Bool=false)
    time_all,mm,t_setpoint,val_setpoint = SSA_perturbations(deepcopy(S),
                                                                n0, durations, changes,
                                                                timestep=timestep,
                                                                transitionClassIndex=transitionClassIndex,
                                                                asserting=true)
    # Prepare futures array and spawn the processes
    F = []
    for i=1:Nsamples
        fut = @spawnat :any SSA_perturbations(deepcopy(S),
                                                        n0, durations, changes,
                                                        timestep=timestep,
                                                        transitionClassIndex=transitionClassIndex,
                                                        asserting=false)
        push!(F, fut)
    end
    
    # Now prepare to collect the results
    MM = zeros(length(keys(S.MomDict)),length(time_all))
    MM2 = copy(MM)
    
    # Take care of the raw output
    n = nothing
    if exportRawOutput
        # The n layout is (Species, Cells)(Trajectories)
        n = [ zeros(Int64, size(n0)) for i=1:Nsamples ]
    end
    
    # Collect the results
    for i=1:Nsamples
        tt, Mi, ts, vs, n_local = fetch(F[i])
        println("SSA simulation # $i")
        if exportRawOutput
            n[i] = deepcopy(n_local)
        end
        MM .+= Mi
        MM2 .+= Mi.^2
    end

    MMav = MM/Nsamples
    MMvar = MM2/(Nsamples-1)-(MM.^2)/(Nsamples*(Nsamples-1))
    return time_all, MMav, MMvar, t_setpoint, val_setpoint, n
end

"""
    SSA_perturbations(S, n0, durations, changes;
                             timestep=0.1, seed=nothing, asserting=true)

Run one SSA trajectory with parameter `changes` (perturbations) at the given
`durations` timepoints.
"""
@export function SSA_perturbations(S::System,
                                    n0::Matrix{Int64},
                                    durations::Vector{Float64},
                                    changes::Vector{Float64};
                                    transitionClassIndex::Int=1,
                                    timestep::Float64=0.1,
                                    seed::Union{Nothing,Int64}=nothing,
                                    asserting::Bool=true)
    seed!=nothing ? Random.seed!(seed) : nothing
    asserting ? assert_model(S,n0) : nothing
    @assert length(durations)==length(changes) "Invalid setpoint input"
    LL=length(changes)
    n_start = deepcopy(n0)
    n_out   = deepcopy(n0) # We want to be able to export the raw output (to compute PDFs)
    t_start=0.
    time_all=Vector{Float64}()
    moments_all=zeros(Int64,length(compute_moments(S,n0)),0)
    t_setpoint=Vector{Float64}()
    val_setpoint=Vector{Float64}()
    for i=1:LL
        S.transition_classes[transitionClassIndex].k *= changes[i]
        timepoints=collect(t_start:timestep:t_start+durations[i])
        push!(t_setpoint,t_start); push!(val_setpoint,S.transition_classes[1].k/S.transition_classes[2].k)
        push!(t_setpoint,timepoints[end]); push!(val_setpoint,S.transition_classes[1].k/S.transition_classes[2].k)
        n_out, mom_out = SSA(S, n_start, timepoints, full_story=false, asserting=false)
        time_all = [time_all;timepoints]
        moments_all = [moments_all mom_out]
        t_start += durations[i]
        n_start = deepcopy(n_out)
    end
    return time_all, moments_all, t_setpoint, val_setpoint, n_out
end


# returns the average and variance of moments of system S with initial condition n0 at times timepoints over Nsamples stochastic realizations
@export function SSA(S::System, n0::Matrix{Int64}, timepoints::Vector{T}, Nsamples::Int64) where T <: Real
    assert_model(S,n0)
    MM=zeros(length(keys(S.MomDict)),length(timepoints))
    MM2=copy(MM)
    MM, MM2 = @distributed (+) for i=1:Nsamples
        println("SSA simulation # $i")
        n,Mi=SSA(S,n0,timepoints,asserting=false)
        [Mi, Mi.^2] # Reduction will occur on these values
    end
    MMav=MM/Nsamples
    MMvar=MM2/(Nsamples-1)-(MM.^2)/(Nsamples*(Nsamples-1))
    return MMav, MMvar
end



# performs a stochastic simulation of system S with initial condition n0
## if timepoints is a vector, it returns moments at each provided time
## if timepoints is a single value, it returns moments at each event up to final time
## if keyword full_story = true, it returns population state at each timepoint or event,
## else only the final population state is returned
@export function SSA(S::System, n0::Matrix{Int64}, timepoints::Union{Vector{Float64},Float64};
                maxiters::Int64=100000, full_story::Bool=false,
                seed::Union{Nothing,Int64}=nothing,asserting::Bool=true)
    seed!=nothing ? Random.seed!(seed) : nothing
    asserting ? assert_model(S,n0) : nothing
    DD,Ncomp=size(n0)
    CL=length(S.transition_classes)
    n=copy(n0)
    rates = [S.transition_classes[i].k for i=1:CL]
    Mom=compute_moments(S,n)
    NcompM=[Ncomp]
    if length(timepoints) > 1
        path_flag=false
        TT=length(timepoints)
        MM=zeros(Int64,length(Mom),TT)
        MM[:,1]=Mom
        simtime=timepoints[1]
        full_story ? ( n_story=[zeros(Int64,DD,0) for i=1:TT]; n_story[1]=deepcopy(n0) ) : nothing
    else
        path_flag=true
        reaction_times=zeros(maxiters)
        MMt=zeros(Int64,length(Mom),maxiters)
        MMt[:,1]=Mom
        simtime=0.
        full_story ? ( n_story=[zeros(Int64,DD,0) for i=1:maxiters]; n_story[1]=deepcopy(n0) ) : nothing
    end
    r_indices=zeros(Int64,2)
    xc=[zeros(Int64,DD),zeros(Int64,DD)]
    yc=[zeros(Int64,DD),zeros(Int64,DD)]
    H_weights=zeros(Int64,CL)
    H_classes=zeros(CL)
    time_index=2
    simulation_flag=true
    class_counter=zeros(Int64,length(rates))
    while simulation_flag
        #Evaluate global propensity
        for c=1:CL
            H_weights[c] = S.transition_classes[c].H(n,Mom)
            H_classes[c] = rates[c]*float(H_weights[c])
        end
        Htot = sum(H_classes)
        #Compute time of next reaction
        Htot>0.0 ? simtime -= log(1-rand())/Htot : simtime=Inf  #if total propensity is zero, just end
        if !path_flag
            while timepoints[time_index]<simtime
                MM[:,time_index]=Mom
                full_story ? n_story[time_index]=n[:,1:Mom[1]] : nothing
                time_index+=1
                time_index>TT ? (simulation_flag=false;break) : nothing
            end
        end
        (path_flag && (simtime > timepoints)) ? simulation_flag=false : nothing
        if simulation_flag
            #Compute type and detail of next reaction
            next_class = climbtower(rand()*Htot,H_classes)
            class_counter[next_class]+=1
            if S.transition_classes[next_class].rc > 0
                draw_reactant_indices!(S.transition_classes[next_class].fast_sample_reactants!,r_indices,H_weights[next_class],xc,S,next_class,n,Mom)
            end
            if S.transition_classes[next_class].pc > 0
                draw_product_compartments!(S.transition_classes[next_class].parameters,yc,xc,S,next_class)#,n,Mom)
            end
            if Mom[1] == NcompM[1]
                if S.transition_classes[next_class].DeltaN > 0
                    n = [n zeros(Int64,size(n))]
                    NcompM[1]=size(n,2)
                end
            end
            update_all!(S,next_class,r_indices,xc,yc,n,NcompM,Mom)
            if path_flag
                reaction_times[time_index]=simtime
                MMt[:,time_index]=Mom
                full_story ? n_story[time_index]=n[:,1:Mom[1]] : nothing
                time_index+=1
                time_index>maxiters ? (println("OVERFLOW! Increase maxiters");simulation_flag=false;break) : nothing
            end
            #println("t=",simtime," Z1tot=",Mom[2]," Z2tot=",Mom[3]," N=", Mom[4]," X2tot=",Mom[6])
        end
    end
    #println(class_counter)
    full_story ? nothing : n_story=n[:,1:Mom[1]]
    if path_flag
        time_index-=1
        if full_story
            return  n_story[1:time_index], reaction_times[1:time_index], MMt[:,1:time_index]
        else
            return  n[:,1:Mom[1]], reaction_times[1:time_index], MMt[:,1:time_index]
        end
    else
        if full_story
            return  n_story, MM
        else
            return  n[:,1:Mom[1]], MM
        end
    end
end ## END SSA



#function to compute index of tower sampling from a vector
function climbtower(rr::Float64,vect::Vector{T}) where T <: Real
    i=1
    cumul=vect[1]
    while rr>cumul
        i+=1
        cumul+=vect[i]
    end
    return i
end


# if necessary, expands the state matrix in order to accomodate an increase of compartment number
function check_state_boundary!(n::Matrix{Int64},Ncomp::Int64,NcompM::Vector{Int64})
    if Ncomp == NcompM[1]
        n = [n zeros(Int64,size(n))]
        NcompM[1]=2*Ncomp
    end
end

# if the fast_sample_reactants! function is provided, samples the next reacting compartments efficiently
function draw_reactant_indices!(fast_sample!::Function,r_indices::Vector{Int64},propensity_weight::Int64,
                                xc::Vector{Vector{Int64}},
                                S::System,next_class::Int64,
                                n::Matrix{Int64},Mom::Vector{Int64})
    r_indices[2]=0  # necessary to prevent troubles with two compartments ...
    fast_sample!(r_indices,n,Mom)
    for i=1:S.transition_classes[next_class].rc
        for d=1:S.n_species xc[i][d] = n[d,r_indices[i]] end
    end
end

# if the fast_sample_reactants! function is not provided, samples the next reacting compartments of class c with tower sampling
function draw_reactant_indices!(fast_sample!::Nothing,r_indices::Vector{Int64},propensity_weight::Int64,
                                xc::Vector{Vector{Int64}},
                                S::System,next_class::Int64,
                                n::Matrix{Int64},Mom::Vector{Int64})
    #fill!(r_indices,0)
    Ncomp = Mom[1]
    DD=size(n,1)
    if S.transition_classes[next_class].rc == 1
        r_indices[1]=1
        for d=1:DD xc[1][d]=n[d,1] end
        val = S.transition_classes[next_class].g(xc)
        rv=rand()*propensity_weight
        while val < rv
            r_indices[1]+=1
            for d=1:DD xc[1][d]=n[d,r_indices[1]] end
            val += S.transition_classes[next_class].g(xc)
        end
    elseif S.transition_classes[next_class].rc == 2
        r_indices[1]=1 ; for d=1:DD xc[1][d]=n[d,1] end
        r_indices[2]=2 ; for d=1:DD xc[2][d]=n[d,2] end
        val = S.transition_classes[next_class].g(xc)
        rv=rand()*propensity_weight
        while val < rv
            if r_indices[2] != Ncomp
                r_indices[2] += 1
            else
                r_indices[1] += 1
                r_indices[2] = r_indices[1]+1
                for d=1:DD xc[1][d]=n[d,r_indices[1]] end
            end
            for d=1:DD xc[2][d]=n[d,r_indices[2]] end
            val += S.transition_classes[next_class].g(xc)
        end
    end
end

# sample product compartments
function draw_product_compartments!(param::Nothing,yc::Vector{Vector{Int64}},xc::Vector{Vector{Int64}},
                                    S::System,next_class::Int64)
    S.transition_classes[next_class].pi(yc,xc)
end

# sample product compartments, when outcome distribution pi needs some additional parameters
function draw_product_compartments!(param,yc::Vector{Vector{Int64}},xc::Vector{Vector{Int64}},
                                    S::System,next_class::Int64)
    S.transition_classes[next_class].pi(yc,xc,param)
end

# efficient implementation of x^exponent for low indices
function recursive_exponentiation(x::Int64,exponent::Int64)
    if exponent==0
        return 1
    else
        return x*recursive_exponentiation(x,exponent-1)
    end
end

# updates population state n and moments vector Mom with the drawn class and the reactant and product compartments
function update_all!(S::System, next_class::Int64,
                    r_indices::Vector{Int64},
                    xc::Vector{Vector{Int64}}, yc::Vector{Vector{Int64}},
                    n::Matrix{Int64}, NcompM::Vector{Int64}, Mom::Vector{Int64})
    Ncomp=Mom[1]
    for r=1:S.transition_classes[next_class].rc
        for l=1:length(Mom)
            val=1
            for d=1:S.n_species
                val *= recursive_exponentiation(xc[r][d],S.MomDict[l][d]) # xc[r][d]^S.MomDict[l][d]
            end
            Mom[l] -= val
        end
    end
    for p=1:S.transition_classes[next_class].pc
        for l=1:length(Mom)
            val=1
            for d=1:S.n_species
                val *= recursive_exponentiation(yc[p][d],S.MomDict[l][d])  # yc[p][d]^S.MomDict[l][d]
            end
            Mom[l] += val
        end
    end
    if S.transition_classes[next_class].DeltaN == -1
        pos_overwrite=max(r_indices[1],r_indices[2])
        for j=1:S.n_species
            n[j,pos_overwrite]=n[j,Ncomp]
        end
        if S.transition_classes[next_class].rc == 2
            for j=1:S.n_species
                n[j,r_indices[1]] = yc[1][j]
            end
        end
    elseif S.transition_classes[next_class].DeltaN == 0
        for i=1:S.transition_classes[next_class].rc
            for j=1:S.n_species
                n[j,r_indices[i]] = yc[i][j]
            end
        end
    elseif S.transition_classes[next_class].DeltaN == 1
        for j=1:S.n_species
            n[j,Ncomp+1] = yc[1][j]
        end
        if S.transition_classes[next_class].rc == 1
            for j=1:S.n_species
                n[j,r_indices[1]] = yc[2][j]
            end
        end
    else
        error()
    end
end
