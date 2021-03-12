### Models for the paper

"""
    Paper_MassControl(;[params...])

Mass control model for the paper.
Includes cell division.

The species are:
- Z1, Z2 in the medium
- X1, X2 within each cell
"""
@export function MassControl(;
                            omega = 10.0, theta = 0.01, eta = 0.005, k_a = 1.0, # Rates for Z1, Z2
                            k_1 = 0.0, k_2 = 1.0,  # Rates for X1
                            k_3 = 2.0, k_4 = 1.0,  # Rates for X2
                            k_div = 0.005,         # Rate  for cell division
                            k_0 = 0.0025)          # Rate  for cell death
    S = System("MassControl", 1+2+2, # environment or cell , Z1 & Z2, X1 & X2
            Dict(1 =>[0,0,0,0,0], # 
                 2 =>[0,1,0,0,0], # Total Z1 mass (in medium)
                 3 =>[0,0,1,0,0], # Total Z2 mass (in medium)
                 4 =>[1,0,0,0,0], # Number of cells
                 5 =>[1,0,0,1,0], # Total X1 mass (in cell)
                 6 =>[1,0,0,0,1], # Total X2 mass (in cell)
                )
            )

    # SAIC reactions
    reference = new_chemical_reaction_class([0,1,0,0,0], omega)
    reference.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> 1
    reference.fast_sample_reactants! = fast_sample_first

    measurement = new_chemical_reaction_class([0,0,1,0,0], theta)
    measurement.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[6] # This depends on X2
    measurement.fast_sample_reactants! = fast_sample_first

    comparison = new_chemical_reaction_class([0,-1,-1,0,0], eta)
    comparison.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2]*Mom[3] # This depends on Z1+Z2
    comparison.fast_sample_reactants! = fast_sample_first

    actuation = new_chemical_reaction_class([0,0,0,1,0], k_a)
    actuation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2] # Catalyzed by Z1
    actuation.fast_sample_reactants! = fast_sample_uniform_cell

    # cells' reaction network
    prodX1 = new_chemical_reaction_class([0,0,0,1,0], k_1)
    prodX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[4]
    prodX1.fast_sample_reactants! = fast_sample_uniform_cell

    deathX1 = new_chemical_reaction_class([0,0,0,-1,0], k_2)
    deathX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[5]
    deathX1.fast_sample_reactants! = fast_sample_mass_x1

    prodX2 = new_chemical_reaction_class([0,0,0,0,1], k_3)
    prodX2.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[5] # Catalyzed by X1
    prodX2.fast_sample_reactants! = fast_sample_mass_x1

    deathX2 = new_chemical_reaction_class([0,0,0,0,-1], k_4)
    deathX2.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[6]
    deathX2.fast_sample_reactants! = fast_sample_mass_x2

    # Compartment reactions
    cell_division = TransitionClass(1, 2, k_div) # Cell division
    # cell_division.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[5] # Catalyzed by X1
    cell_division.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[6] # Catalyzed by X2
    # The below is how the content of the input compartments (xc) is split
    # into the output compartments (yc)
    cell_division.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                      yc[1][1] , yc[2][1] = 1 , 1        # Here we get 2 new compartments
                      yc[1][2] , yc[2][2] = 0 , 0        # Nothing changes for Z1
                      yc[1][3] , yc[2][3] = 0 , 0        # ... nor for Z2
                      # yc[1][4] = ceil(Int64, xc[1][4]/2) # Half goes in the first daughter
                      yc[1][4] = rand(Binomial(xc[1][4], 0.5)) # Binomial distrib for the first daugther
                      yc[2][4] = xc[1][4] - yc[1][4]     # and half goes in the second one
                      # yc[1][5] = ceil(Int64, xc[1][5]/2) # Same here ...
                      yc[1][5] = rand(Binomial(xc[1][5], 0.5)) # Same here ...
                      yc[2][5] = xc[1][5] - yc[1][5]     # ...
                end
    cell_division.fast_sample_reactants! = fast_sample_cell_div

    cell_death = TransitionClass(1, 0, k_0) # This is a cell death
    cell_death.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[4] # Depends on num cells
    cell_death.fast_sample_reactants! = fast_sample_uniform_cell

    add_transition_class(S, reference, measurement, comparison, actuation, prodX1, deathX1, prodX2, deathX2, cell_division, cell_death)

    return S
end

"""
    Paper_NumberControl(;[params...])

Cell number control model for the paper.
Includes cell division.

The species are:
- Z1, Z2 in the medium
- X1, X2 within each cell
"""
@export function NumberControl(;
                                omega = 10.0, theta = 0.01, eta = 0.005, k_a = 1.0,   # Rates for Z1, Z2
                                k_1 = 0.0, k_2 = 1.0,  # Rates for X1
                                k_3 = 2.0, k_4 = 1.0,  # Rates for X2
                                k_div = 0.005,         # Rate  for cell division
                                k_0 = 0.0025)          # Rate  for cell death
    S = System("MassControl", 1+2+2, # environment or cell , Z1 & Z2, X1 & X2
            Dict(1 =>[0,0,0,0,0], # Number of WHAT???
                 2 =>[0,1,0,0,0], # Total Z1 mass (in medium)
                 3 =>[0,0,1,0,0], # Total Z2 mass (in medium)
                 4 =>[1,0,0,0,0], # Number of cells
                 5 =>[1,0,0,1,0], # Total X1 mass (in cell)
                 6 =>[1,0,0,0,1], # Total X2 mass (in cell)
                )
            )

    # SAIC reactions
    reference = new_chemical_reaction_class([0,1,0,0,0], omega)
    reference.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> 1
    reference.fast_sample_reactants! = fast_sample_first

    measurement = new_chemical_reaction_class([0,0,1,0,0], theta)
    measurement.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[4] # This only depends on the num of cells
    measurement.fast_sample_reactants! = fast_sample_first

    comparison = new_chemical_reaction_class([0,-1,-1,0,0], eta)
    comparison.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2]*Mom[3] # This depends on Z1+Z2
    comparison.fast_sample_reactants! = fast_sample_first

    actuation = new_chemical_reaction_class([0,0,0,1,0], k_a)
    actuation.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[2] # Catalyzed by Z1
    actuation.fast_sample_reactants! = fast_sample_uniform_cell

    # cells' reaction network
    prodX1 = new_chemical_reaction_class([0,0,0,1,0], k_1)
    prodX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[4]
    prodX1.fast_sample_reactants! = fast_sample_uniform_cell

    deathX1 = new_chemical_reaction_class([0,0,0,-1,0], k_2)
    deathX1.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[5]
    deathX1.fast_sample_reactants! = fast_sample_mass_x1

    prodX2 = new_chemical_reaction_class([0,0,0,0,1], k_3)
    prodX2.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[5] # Catalyzed by X1
    prodX2.fast_sample_reactants! = fast_sample_mass_x1

    deathX2 = new_chemical_reaction_class([0,0,0,0,-1], k_4)
    deathX2.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[6]
    deathX2.fast_sample_reactants! = fast_sample_mass_x2

    # Compartment reactions
    cell_division = TransitionClass(1, 2, k_div) # Cell division
    # cell_division.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[5] # Catalyzed by X1
    cell_division.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[6] # Catalyzed by X2
    # The below is how the content of the input compartments (xc) is split
    # into the output compartments (yc)
    cell_division.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                      yc[1][1] , yc[2][1] = 1 , 1        # Here we get 2 new compartments
                      yc[1][2] , yc[2][2] = 0 , 0        # Nothing changes for Z1
                      yc[1][3] , yc[2][3] = 0 , 0        # ... nor for Z2
                      # yc[1][4] = ceil(Int64, xc[1][4]/2) # Half goes in the first daughter
                      yc[1][4] = rand(Binomial(xc[1][4], 0.5)) # Binomial distrib for the first daugther
                      yc[2][4] = xc[1][4] - yc[1][4]     # and half goes in the second one
                      # yc[1][5] = ceil(Int64, xc[1][5]/2) # Same here ...
                      yc[1][5] = rand(Binomial(xc[1][5], 0.5)) # Same here ...
                      yc[2][5] = xc[1][5] - yc[1][5]     # ...
                end
    cell_division.fast_sample_reactants! = fast_sample_cell_div

    cell_death = TransitionClass(1, 0, k_0) # This is a cell death
    cell_death.H = (n::Matrix{Int64}, Mom::Vector{Int64}) -> Mom[4] # Depends on num cells
    cell_death.fast_sample_reactants! = fast_sample_uniform_cell

    add_transition_class(S, reference, measurement, comparison, actuation, prodX1, deathX1, prodX2, deathX2, cell_division, cell_death)

    return S
end



### SUPPORTING FUNCTIONS FOR EFFICIENT SIMULATION ###


function fast_sample_first(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64})
    r_indices[1] = 1
end

function fast_sample_uniform_cell(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64})
    r_indices[1] = rand(2:Mom[1])  # first one is environment ! use Mom[1] to have actual index!
end

function fast_sample_mass_x1(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64})
    rv = rand()*Mom[5]
    r_indices[1] = 2  # first one is environment !
    val = 1.0*n[4,2]
    while val < rv
        r_indices[1]+= 1
        val += n[4, r_indices[1]]
    end
end

function fast_sample_mass_x2(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64})
    rv = rand()*Mom[6]
    r_indices[1] = 2 # first one is environment!
    val = 1.0*n[5,2]
    while val < rv
        r_indices[1]+= 1
        val += n[5, r_indices[1]]
    end
end


function fast_sample_cell_div(r_indices::Vector{Int64}, n::Matrix{Int64}, Mom::Vector{Int64})
    rv = rand()*Mom[6]
    r_indices[1] = 2 # the first one is environment!
    val = 1.0*n[5,2]
    while val < rv
        r_indices[1]+= 1
        val += n[5, r_indices[1]]
    end
end
