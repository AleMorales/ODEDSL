using DataFrames
using Iterators

#######
# Fuse models
# We overload the + operator and dispatch depending on the different possible combinations
#######
function +(x::MESource, y::MESource)
    new_x = deepcopy(x)
    merge!(new_x.Constants, y.Constants)
    merge!(new_x.Parameters, y.Parameters)
    merge!(new_x.Forcings, y.Forcings)
    merge!(new_x.States, y.States)
    merge!(new_x.Species, y.Species)
    merge!(new_x.Components, y.Components)
    merge!(new_x.Reactions, y.Reactions)
    merge!(new_x.MEReactions, y.MEReactions)
    merge!(new_x.MEEquations, y.MEEquations)
    for (key,val) in y.Equations
        if haskey(new_x.Equations, key)
            new_x.Equations[key].Expr =  parse("(" * string(new_x.Equations[key].Expr) * ")" * " + " * "(" * string(val.Expr) * ")")
        else
            new_x.Equations[key] = val
        end
    end
    return new_x
end
function +(x::MESource, y::ReactionSource)
    new_x = deepcopy(x)
    merge!(new_x.Constants, y.Constants)
    merge!(new_x.Parameters, y.Parameters)
    merge!(new_x.Forcings, y.Forcings)
    merge!(new_x.States, y.States)
    merge!(new_x.Species, y.Species)
    merge!(new_x.Reactions, y.Reactions)
    for (key,val) in y.Equations
        if haskey(new_x.Equations, key)
            new_x.Equations[key].Expr =  parse("(" * string(new_x.Equations[key].Expr) * ")" * " + " * "(" * string(val.Expr) * ")")
        else
            new_x.Equations[key] = val
        end
    end
    return new_x
end
+(x::ReactionSource, y::MESource) = y + x
function +(x::MESource, y::OdeSource)
    new_x = deepcopy(x)
    merge!(new_x.Constants, y.Constants)
    merge!(new_x.Parameters, y.Parameters)
    merge!(new_x.Forcings, y.Forcings)
    merge!(new_x.States, y.States)
    for (key,val) in y.Equations
        if haskey(new_x.Equations, key)
            new_x.Equations[key].Expr =  parse("(" * string(new_x.Equations[key].Expr) * ")" * " + " * "(" * string(val.Expr) * ")")
        else
            new_x.Equations[key] = val
        end
    end
    return new_x
end
+(x::OdeSource, y::MESource) = y + x
function +(x::OdeSource, y::OdeSource)
    new_x = deepcopy(x)
    merge!(new_x.Constants, y.Constants)
    merge!(new_x.Parameters, y.Parameters)
    merge!(new_x.Forcings, y.Forcings)
    merge!(new_x.States, y.States)
    for (key,val) in y.Equations
        if haskey(new_x.Equations, key)
            new_x.Equations[key].Expr =  parse("(" * string(new_x.Equations[key].Expr) * ")" * " + " * "(" * string(val.Expr) * ")")
        else
            new_x.Equations[key] = val
        end
    end
    return new_x
end
function +(x::ReactionSource, y::ReactionSource)
    new_x = deepcopy(x)
    merge!(new_x.Constants, y.Constants)
    merge!(new_x.Parameters, y.Parameters)
    merge!(new_x.Forcings, y.Forcings)
    merge!(new_x.States, y.States)
    merge!(new_x.Species, y.Species)
    merge!(new_x.Reactions, y.Reactions)
    for (key,val) in y.Equations
        if haskey(new_x.Equations, key)
            new_x.Equations[key].Expr =  parse("(" * string(new_x.Equations[key].Expr) * ")" * " + " * "(" * string(val.Expr) * ")")
        else
            new_x.Equations[key] = val
        end
    end
    return new_x
end
function +(x::ReactionSource, y::OdeSource)
    new_x = deepcopy(x)
    merge!(new_x.Constants, y.Constants)
    merge!(new_x.Parameters, y.Parameters)
    merge!(new_x.Forcings, y.Forcings)
    merge!(new_x.States, y.States)
    for (key,val) in y.Equations
        if haskey(new_x.Equations, key)
            new_x.Equations[key].Expr =  parse("(" * string(new_x.Equations[key].Expr) * ")" * " + " * "(" * string(val.Expr) * ")")
        else
            new_x.Equations[key] = val
        end
    end
    return new_x
end
+(x::OdeSource, y::ReactionSource) = y + x


#######
# Promote master equation to reaction
# Convert mereaction, mevariable and components into a reaction model
#######

# 1. From components, generate the Species
# 2. For each reaction, identify the species on the origin and destination of the transition and generate all the necessary reactions
# 3. Expand species in and mevariable
function convert_master_equation(model::MESource)
  new_model = deepcopy(model)
  components_by_species = classify_components(model.Components)
  combinations_per_species = component_form_combinations(components_by_species)
  reactions = expand_master_equation_rules(model, components_by_species, combinations_per_species)
  species = generate_species_from_components(model, components_by_species, combinations_per_species);
  equations = convert_master_equation_equations(model, components_by_species, combinations_per_species);
  merge!(new_model.Species, species)
  for i in keys(components_by_species)
    delete!(new_model.Species, i)
  end
  merge!(new_model.Equations, equations)
  merge!(new_model.Reactions, reactions)
  return ReactionSource(new_model.Constants, new_model.Parameters, new_model.Forcings,
                        new_model.States, new_model.Species, new_model.Compartments, new_model.Equations,
                        new_model.Reactions)
end

function classify_components(components::Dict{String, Component})
  components_by_species = Dict{String, Dict{String, Component}}()
  for (key,val) in components
    !haskey(components_by_species, val.Species) && (components_by_species[val.Species] = Dict{String, Component}())
    components_by_species[val.Species][key] = val
  end
  return components_by_species
end

function total_combinations_per_species(components_by_species::Dict{String, Dict{String, Component}})
  number_of_combinations = Dict{String, Int}()
  for (key1,val1) in components_by_species
    number_of_combinations[key1] = 1
    for (key2,val2) in val1
      number_of_combinations[key1] *= length(val2.Forms)
    end
  end
  return number_of_combinations
end

function component_form_combinations(components_by_species::Dict{String, Dict{String, Component}})
  number_of_combinations = total_combinations_per_species(components_by_species)
  combinations_per_species = Dict{String, Any}()
  for (key1,val1) in components_by_species
    # Get all form-component combinations for each form in the Species
    names_forms = {}
    for (key2,val2) in val1
      push!(names_forms, val2.Forms)
    end
    # Generate matrix with all possible combinations of form-components in the species
    combination = convert(DataFrame, Array(String, (number_of_combinations[key1], length(val1))))
    names!(combination, convert(Array{ASCIIString, 1}, collect(keys(val1))))
    small = names_forms[1]
    combination[:,1] = repeat(small; outer = [int(number_of_combinations[key1]/length(small))])
    for i in 2:length(names_forms)
        small = repeat(names_forms[i]; inner = [length(small)])
        combination[:,i] = repeat(small; outer = [int(number_of_combinations[key1]/length(small))])
    end
    combinations_per_species[key1] = combination
  end
  return combinations_per_species
end

# This return a Dict of Species that may be fused with the existing Species Dict in the model
function generate_species_from_components(model, components_by_species, combinations_per_species)
  # Create the species component of the reaction model (this actually correspons to the component-forms)
  species = Dict{String, Species}() 
  # Iterate for each species
  for (key1,val1) in components_by_species
    # Get the units and compartment from the species (this species will actually be deleted)
    units = model.Species[key1].Units # Units in Species and Components are actually redundant
    compartment =  model.Species[key1].Compartment
    names_components = collect(keys(combinations_per_species[key1].dicts[2]))
    # Iterate over the rows of the matrix containing all possible combinations
    for i in 1:size(combinations_per_species[key1])[1]
      # Construct the names of the species using "_" to append the different component-forms
        name_species = key1*"_"*names_components[1]*combinations_per_species[key1][i,1]
        location_form = findin(model.Components[names_components[1]].Forms, [combinations_per_species[key1][i,1]])
        value = model.Components[names_components[1]].Values[location_form]
      for j in 2:size(combinations_per_species[key1])[2]
        name_species *= "_"*names_components[j]*combinations_per_species[key1][i,j]
        location_form = findin(model.Components[names_components[j]].Forms, [combinations_per_species[key1][i,j]])
        value .*= model.Components[names_components[j]].Values[location_form]
      end
      value .*= value[1]*model.Species[key1].Value
      species[name_species] = Species(value[1], compartment, units)
    end
  end
  return species
end

# This will expand the transitions in a master equation model to the different reactions
function expand_master_equation_rules(model, components_by_species, combinations_per_species)
  reactions = Dict{String, Reaction}()
  # Iterate for each species
  for (key1,val1) in model.MEReactions
    # Identify the rows in the combinations_per_species that correspond to the "from" elements in the master equation reaction
    from_rows = findin(combinations_per_species[val1.Species][:,convert(ASCIIString,val1.Component[1])], [val1.From[1]])
    for i in 2:length(val1.Component)
      from_rows = intersect(from_rows, findin(combinations_per_species[val1.Species][:,convert(ASCIIString,val1.Component[i])], [val1.From[i]]))
    end
    to_rows = findin(combinations_per_species[val1.Species][:,convert(ASCIIString, val1.Component[1])], [val1.To[1]])
    for i in 2:length(val1.Component)
      to_rows = intersect(to_rows, findin(combinations_per_species[val1.Species][:,convert(ASCIIString,val1.Component[i])], [val1.To[i]]))
    end
    # All possible combinations
    combinations_from_to = collect(product(from_rows, to_rows))
    # Only those combinations where the unspecified components are equal can be matched.
    good_match = {} 
    for i = 1:length(combinations_from_to)
        select_columns = intersect(names(combinations_per_species[val1.Species]), [symbol("$i") for i in val1.Component])
        other_comps_from = combinations_per_species[val1.Species][combinations_from_to[i][1],select_columns]
        other_comps_to = combinations_per_species[val1.Species][combinations_from_to[i][2],select_columns]
        if other_comps_from == other_comps_to
            push!(good_match, combinations_from_to[i])
        end
    end      
    # For every combination, generate a Reaction where the species are create using the same notation as for the "generate_species_from_components" function
    names_components = collect(keys(combinations_per_species[val1.Species].dicts[2]))
    for i in 1:size(good_match)[1]
      name_substrate = val1.Species*"_"*names_components[1]*combinations_per_species[val1.Species][good_match[i][1],1]
      name_product = val1.Species*"_"*names_components[1]*combinations_per_species[val1.Species][good_match[i][2],1]
      for j in 2:size(combinations_per_species[val1.Species])[2]
        name_substrate *= "_"*names_components[j]*combinations_per_species[val1.Species][good_match[i][1],j]
        name_product *= "_"*names_components[j]*combinations_per_species[val1.Species][good_match[i][2],j]
      end
      # The reaction is named as v_namesubstrate_nameproduct
      # Note that I assume first order kinetics (expr cannot be edited directly)
      name_reaction = "v_"*name_substrate*"_"*name_product
      expr = parse(name_substrate*"*("*string(val1.Expr)*")")
      reactions[name_reaction] = Reaction([Reactant(name_substrate, 1)], [Reactant(name_product, 1)], expr,
                            model.Species[val1.Species].Compartment, false, val1.Dim)
    end
  end
  return reactions
end

# This will expand the MEEquations into normal equations
# The rule is that any Species[[Component[Form]]] is expanded into the sum of all corresponding forms
function convert_master_equation_equations(model, components_by_species, combinations_per_species)
  equations = Dict{String, Equation}()
  for (key,val) in model.MEEquations
    substituted_expression = expand_master_equation_expression(val.Expr, combinations_per_species)
    equations[key] = Equation(parse(substituted_expression), val.Exported, val.Dim)
  end
  return equations
end

function expand_master_equation_expression(expr, combinations_per_species)
  # Extract all species, components and forms
  species = matchall(r"[a-zA-Z0-9_\s]+(?=\[\[)", string(expr))
  species = [strip(x) for x in species]
  component_forms = matchall(r"\[\s*\[\s*[a-zA-Z0-9_\s\]\[]*\s*\]\s*\]",string(expr))
  # An ragged array of arrays. Each subarrays contains the forms associated to each species mentioned
  forms = [matchall(r"\[[a-zA-Z0-9_\s]+\]", x) for x in component_forms]
  cleaned_forms = {}
  c = 1
  for i in forms
    push!(cleaned_forms, [matchall(r"(?<=\[)[a-zA-Z0-9_\s](?=\])", x) for x in i])
    for j in 1:length(cleaned_forms[c])
      cleaned_forms[c][j][1] = strip(cleaned_forms[c][j][1])
    end
    c += 1
  end
  # An ragged array of arrays. Each subarrays contains the components associated to each species mentioned
  components = [matchall(r"[a-zA-Z0-9_]+(?=\[)", x) for x in component_forms]
  for i in 1:length(components)
    for j in 1:length(components[i])
      components[i][j] = strip(components[i][j])
    end
  end
  # for each species included in the expression, identify all rows in the combinations matrix (same as in expand_master_equation_rules)
  substituted_expression = string(expr)
  for i in 1:length(species)
    rows = findin(combinations_per_species[species[i]][:,convert(ASCIIString,components[i][1])], [cleaned_forms[i][1]])
    for j in 2:length(cleaned_forms[i])
      rows = intersect(rows, findin(combinations_per_species[species[i]][:,convert(ASCIIString,components[i][j])], [cleaned_forms[i][j]]))
    end
    names_components = collect(keys(combinations_per_species[species[i]].dicts[2]))
    expanded_name = "("
    for j in 1:size(rows)[1]
      if j == 1
        expanded_name *= species[i]
      else 
        expanded_name *= " + "*species[i]
      end
      for h in 1:size(combinations_per_species[species[i]])[2]
        expanded_name *= "_"*names_components[h]*combinations_per_species[species[i]][rows[j],h]
      end
    end
    expanded_name *= ")"
    substituted_expression = replace(substituted_expression, species[i]*component_forms[i], expanded_name)
  end
  return substituted_expression
end

#######
# Promote reaction to ode
# Convert species, reaction, and generate time derivatives
#######
function convert_reaction_model(model::ReactionSource)
  new_model = deepcopy(model)
  new_states = species_to_states(new_model.Species)
  compartment_to_input!(new_model.Compartments, new_model.Parameters, new_model.States)
  new_equations = reaction_to_derivative(new_model);
  merge!(new_model.States, new_states)
  #merge!(new_model.Equations, new_equations)
  for (key,val) in new_equations
      if haskey(new_model.Equations, key)
          new_model.Equations[key].Expr =  parse("(" * string(new_model.Equations[key].Expr) * ")" * " +
                                                    " * "(" * string(val.Expr) * ")")
      else
          new_model.Equations[key] = val
      end
  end
  return OdeSource(new_model.Constants, new_model.Parameters, new_model.Forcings,
                        new_model.States, new_model.Equations)
end

function species_to_states(species)
  # Convert species into states
  states = OrderedDict{String, Parameter}()
  for (key,val) in species
      states[key] = Parameter(val.Value, val.Units)
  end
  return states
end

function compartment_to_input!(compartments, parameters, states)
  for (key,val) in compartments
      if val.InputType == "parameter"
          parameters[key] = Parameter(val.Value, val.Units)
      elseif val.InputType == "state"
          states[key] = Parameter(val.Value, val.Units)
      end
  end
end

function reaction_to_derivative(model)
  states_as_reactants = Dict{String,Any}()
  for (key,val) in model.Reactions
    for i in val.Substrates
      !haskey(states_as_reactants, i.Name) && (states_as_reactants[i.Name] = {})
      push!(states_as_reactants[i.Name], (key, i.Stoichiometry, -1))
    end
    for i in val.Products
      !haskey(states_as_reactants, i.Name) && (states_as_reactants[i.Name] = {})
      push!(states_as_reactants[i.Name], (key, i.Stoichiometry, 1))
    end    
  end
  new_equations = Dict{String, Equation}()
  for (key,val) in model.Reactions
      new_equations[key] = Equation(val.Expr, val.Exported, val.Dim)
  end
  for (key,val) in states_as_reactants
    if !haskey(model.Species, key)
      delete!(states_as_reactants, key)
    end
  end
  new_states = collect(keys(states_as_reactants))
  for key in new_states
      time_derivative = ""
      unit = ""
      for i in states_as_reactants[key]
          if time_derivative != ""
            time_derivative *= " + " * "(" * string(i[3]) * "*" * string(i[2]) * "(" * string(new_equations[i[1]].Expr) * ")" * ")" *
                               model.Reactions[i[1]].Compartment * "/" * model.Species[key].Compartment
          else 
            time_derivative *= "(" * string(i[3]) * "*" * string(i[2]) * "(" * string(new_equations[i[1]].Expr) * ")" * ")" *
                               model.Reactions[i[1]].Compartment * "/" * model.Species[key].Compartment
          end
          # Assume that all generate the same units (otherwise there is an error!!)
          unit = new_equations[i[1]].Dim * Dimension(model.Parameters[model.Reactions[i[1]].Compartment].Units) / 
                 Dimension(model.Parameters[model.Species[key].Compartment].Units)
      end
      new_equations["d_"*key*"_dt"] = Equation(parse(time_derivative), true, unit)
  end 
  return new_equations
end

#######
# Sort equations
# Generate a datatype that is similar to the source types by including an array of sorted equations instead of equations
#######
function generate_level0(model)
  equations = Dict{String, Equation}()
  c = 1
  for (key,val) in model.Parameters
    equations[key] = Equation(parse(key * "= params[$c]"), false, val.Units.d)
    c += 1
  end
  c = 1
  for (key,val) in model.States
    equations[key] = Equation(parse(key * "= states[$c]"), false, val.Units.d)
    c += 1
  end  
  c = 1
  for (key,val) in model.Forcings
    equations[key] = Equation(parse(key * "= forcs[$c]"), false, val.Units.d)
    c += 1
  end   
  for (key,val) in model.Constants
    rhs = float(val.Value * val.Units.f)
    equations[key] = Equation(parse(key * "= $(rhs)"), false, val.Units.d)
    c += 1
  end        
  return equations
end

function sort_equations(model::OdeSource)
    checked_ode_model, rhs_to_lhs, lhs_to_rhs = check_lhs_rhs(model);
    check_derivatives(checked_ode_model);
    equations = Dict{String, Equation}[]
    push!(equations, generate_level0(model))
    unsorted_equations = collect(keys(checked_ode_model.Equations))
    while length(unsorted_equations) > 0
        nsorted = length(unsorted_equations)
        unsorted_keys = unsorted_equations
        for j in unsorted_keys
          rhs = Set(lhs_to_rhs[j])
          level = 2  
          for i in 1:length(equations)
            if length(rhs) == 0
              break
            end
            rhs = setdiff(rhs, intersect(rhs, Set(collect(keys(equations[i])))))
            level = i + 1
          end
          if length(rhs) == 0
            deleteat!(unsorted_equations, findfirst(unsorted_equations, string(j))[1])
            if level <= length(equations)
              equations[level][j] = checked_ode_model.Equations[j]
            else
              tmp = Dict{String, Equation}()
              tmp[j] =  checked_ode_model.Equations[j]
              push!(equations, tmp)
            end
          end
        end
        if nsorted == length(unsorted_equations)
          error("There seem to be some circular dependencies")
        end
    end
    return OdeSorted(checked_ode_model.Constants, checked_ode_model.Parameters, 
                  checked_ode_model.Forcings, checked_ode_model.States, equations)
end

function sort_equations(model::MESource)
  error("Before sorting the equations of the model you first need to convert into an ODE model")
end

function sort_equations(model::ReactionSource)
  error("Before sorting the equations of the model you first need to convert into an ODE model")
end

#######
# Compress model
# Reduces any model to level 2 (i.e. only input instruction and calculatio of observers and time derivatives)
#######
function compress_model(model::OdeSorted)
  equations = deepcopy(model.SortedEquations)
  if length(equations) < 3
    return model 
  end
  for i in linrange(length(equations), 3, length(equations) - 2)
    for (key,val) in equations[i]
      equations[i][key].Expr = expand_expression(val.Expr, equations)
      equations[i-1][key] = equations[i][key]
    end
  end
  new_model = deepcopy(model)
  new_model.SortedEquations = equations[1:2]
  lhs_names = collect(keys(new_model.SortedEquations[2]))
  for i in lhs_names
    if !new_model.SortedEquations[2][i].Exported
      delete!(new_model.SortedEquations[2], i)
    end
  end
  return new_model
end

# Go through the AST and expand symbols into expressions
function expand_expression(expression::Expr, equations)
  for i in 1:length(expression.args)
      expression.args[i] = expand_expression(expression.args[i], equations)
  end
  return expression
end
function expand_expression(variable::Symbol, equations)
  if in(string(variable), list_of_functions) | in(string(variable), collect(keys(equations[1])))
    return variable
  else 
    for i in 2:length(equations)
      if in(string(variable), collect(keys(equations[i])))
        return :(($(equations[i][string(variable)].Expr))) 
      end
    end
    error("I cannot find dependency for $variable")
  end
end
function expand_expression(variable::Number, equations)
  return variable
end


function get_lhs(model)
    states = collect(keys(model.States))
    parameters = collect(keys(model.Parameters))
    forcings = collect(keys(model.Forcings))
    constants = collect(keys(model.Constants)) 
    equations = collect(keys(model.Equations))    
    return states, parameters, forcings, constants, equations
end

function get_rhs(ex::Expr)
  output = Set()
  for i in ex.args
    if isa(i, Expr) || isa(i, Symbol)
      union!(output, get_rhs(i))
    end
  end
  return(output)
end

function get_rhs(ex::Symbol)
  output = Set()
  if string(ex) ∉ list_of_functions
    union!(output, {string(ex)})
  end
  return(output)
end

function get_rhs(ex::Number)
  output = Set()
end

function get_rhs(Equations)
    rhs_vars = Set()
    rhs_to_lhs = Dict{String, Array{String, 1}}()
    lhs_to_rhs = Dict{String, Array{String, 1}}()
    for (key,val) in Equations
        new_rhs_vars = get_rhs(val.Expr)
        union!(rhs_vars,new_rhs_vars)
        for j in new_rhs_vars
          if haskey(lhs_to_rhs, key)
                push!(lhs_to_rhs[key], j)
          else
                lhs_to_rhs[key] = [j]
          end            
          if haskey(rhs_to_lhs, j)
                push!(rhs_to_lhs[j], key)
          else
                rhs_to_lhs[j] = [key]
          end
        end
        if !haskey(lhs_to_rhs, key)
          lhs_to_rhs[key] = String[]
        end
   end
    rhs_vars_array = String[]
    for i in rhs_vars
        push!(rhs_vars_array, string(i))
    end
   return rhs_vars_array, rhs_to_lhs, lhs_to_rhs
end #function get_rhs_var{Dict} 

function check_lhs_rhs(model)
    lhs_states, lhs_parameters, lhs_forcings, lhs_constants, lhs_equations = get_lhs(model)
    all_rhs, rhs_to_lhs, lhs_to_rhs = get_rhs(model.Equations)
    new_model = deepcopy(model)
    for i in all_rhs
        i ∉ [lhs_states, lhs_parameters, lhs_forcings, lhs_constants, lhs_equations] ? error("The input $i is missing from the model") : nothing
    end
    for i in lhs_parameters
        if i ∉ all_rhs 
            warn("The output is not dependent on the parameter $i. It has been removed from the model.")
            delete!(new_model.Parameters, i)
        end
    end
    for i in lhs_constants
        if i ∉ all_rhs 
            warn("The output is not dependent on the constant $i. It has been removed from the model.")
            delete!(new_model.Constants, i)
        end
    end  
    for i in lhs_forcings
        if i ∉ all_rhs 
            warn("The output is not dependent on the forcing $i. It has been removed from the model.")
            delete!(new_model.Forcings, i)
        end
    end  
    for i in lhs_states
        if i ∉ all_rhs 
            warn("The output is not dependent on the state $i.")
        end
    end 
    return new_model, rhs_to_lhs, lhs_to_rhs
end


function check_derivatives(model)
    new_model = deepcopy(model)
    names_derivatives = collect(keys(model.States))
    for i in 1:length(names_derivatives)
        names_derivatives[i] = "d_"*names_derivatives[i]*"_dt"
    end
    lhs = collect(keys(model.Equations))
    for i in names_derivatives
        if i ∉ lhs
            error("The time derivative $i is not presented in the model.") 
        end
        new_model.Equations[i].Exported = true
    end
    return new_model
end


