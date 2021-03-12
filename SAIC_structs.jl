
# TRANSITION CLASS
# object to define the structure and the rate law of a transition class
@export mutable struct TransitionClass
    rc::Int64                  # Number of reactant compartments
    pc::Int64                  # Number of product compartments
    DeltaN::Int64              # for conveniency = pc - rc
    k::Float64                 # Rate constant
    parameters::Union{Nothing,Float64,Vector{Int64}}
    g::Union{Function,Nothing}         # Reactant compartments kernel
    pi::Union{Function,Nothing}        # Product compartments distribution
    H::Union{Function,Nothing}         # Total propensity
    fast_sample_reactants!::Union{Function,Nothing} # optional, optimized for specific class
    function TransitionClass(rc::Int64, pc::Int64, k::Float64)
            @assert rc >= 0 "Number of reactant compartments can't be negative!"
            @assert pc >= 0 "Number of product compartments can't be negative!"
            @assert rc <= 2 "The simulator doesn't support more than 2 reactant compartments!"
            @assert pc <= 2 "The simulator doesn't support more than 2 product compartments!"
            @assert k >= 0. "Rate constant can't be negative!"
        return new(rc,pc,pc-rc,k,nothing,nothing,nothing,nothing,nothing)
    end
end

# shortcut definition for one-compartment chemical reaction
## change_vector is the integer update vector of the chemical species
@export function new_chemical_reaction_class(change_vector::Vector{Int64}, rate::Float64)
    class=TransitionClass(1,1,rate)
    class.parameters = change_vector
    class.pi = pi_chemical_reaction!
    return class
end
# shortcut definition for compartment update upon chemical reaction
@export function pi_chemical_reaction!(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}}, change_vector::Vector{Int64})
    for d=1:length(yc[1])
        yc[1][d] = xc[1][d] + change_vector[d]
    end
end


# SYSTEM OBJECT
# object to collect a set of transition classes and speficy model details
# the fields moment_equations and init_moments can be optionally assigned to the system's moment equations (not required for SSA)
@export mutable struct System
    name::String
    n_species::Int64                                  # Number of chemical species
    transition_classes::Vector{TransitionClass}       # Defining model dynamics (as used in SSA simulations)
    MomDict::Dict{Int64,Vector{Int64}}                # link moment index to its gamma exponents
    moment_equations::Union{Function,Nothing}         # Moment Equations f(dM,M,S,t)
    init_moments::Union{Function,Nothing}             # Based on implementation of Moment Equations
    function System(name::String, n_species::Int64, MomDict::Dict{Int64,Vector{Int64}})
            for exponents in values(MomDict)
               @assert size(exponents,1) == n_species "Moment exponents do not match the dimensionality of the compartment content!"
            end
            @assert MomDict[1] == zeros(Int64,n_species) "The moment of index 1 must map to the zero-order moment! (i.e. $(zeros(Int64,n_species))) "
        return new(name,n_species,Vector{TransitionClass}(),MomDict,nothing,nothing)
    end
end

# adds some transition classes to a System object S
@export function add_transition_class(S::System,c...)
    @assert prod(map(class -> typeof(class) == TransitionClass, c))
    S.transition_classes = [S.transition_classes; c...]
end

# computes moments of population state n as given by the index -> exponents dictionary in the system S
@export function compute_moments(S::System, n::Matrix{Int64})
DD,Ncomp=size(n)
M=zeros(Int64,length(keys(S.MomDict)))
for k=1:length(M)
    n_vals=ones(Int64,size(n,2))
    for i=1:length(S.MomDict[k])
        n_vals .*= n[i,:].^S.MomDict[k][i]
    end
    M[k]=sum(n_vals)
end
return M
end


# check that initial condition matches the system
@export function assert_model(S::System,n0::Matrix{Int64})
    @assert size(n0,1) == S.n_species "Initial condition doesn't match the model dimensionality"
    Mom0=compute_moments(S,n0)
    for c in S.transition_classes
	c.H == nothing ? error("Incomplete transition class: specify the class propensity function 'H' ") : nothing
        try
            c.H(n0,Mom0)
        catch y
            error("The moment dictornary of the provided model is incomplete!")
        end
	c.rc > 0 && c.g == nothing && c.fast_sample_reactants! == nothing ? error("Incomplete transition class: speficy 'g' or 'fast_sample_reactants!' ") : nothing
	c.pc > 0 && c.pi == nothing ? error("Incomplete transition class: speficy the outcome distribution 'pi' ") : nothing
    end
end
