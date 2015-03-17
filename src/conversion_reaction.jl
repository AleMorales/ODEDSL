# Convert ReactionSource into ODESource
function OdeSource(model::dt.ReactionSource)
  # We will store the modifications in new_model
  new_model = deepcopy(model)

  # Generate the list of states that can be derived from species
  # If there are components, then there is one state from species and a variable
  # With the sum of all the species.
  new_states, new_variables, new_species, all_comp_forms = species_to_states(new_model)
  merge!(new_model.States, new_states)
  merge!(new_model.Variables, new_variables)
  merge!(new_model.Species, new_species)

  # Convert reactions into variables, expanding by components if necessary
  # This will generate the intermediate variables as well as the time derivatives
  new_variables, new_derivatives = reactions_to_variables(new_model, new_states, new_variables, all_comp_forms)
  merge!(new_model.Variables, new_variables)
  merge!(new_model.Derivatives, new_derivatives)

  # Go through all derivatives and expand any pattern expression into a sum of states
  for (key,val) in new_model.Derivatives
    val.Expr = expand_components(val.Expr, all_comp_forms)
  end

  # Go through all variables and expand any pattern expression
  for (key,val) in new_model.Variables
    val.Expr = expand_components(val.Expr, all_comp_forms)
  end

  # Go through all observed and expand any pattern expression
  for (key,val) in new_model.Observed
    val.Expr = expand_components(val.Expr, all_comp_forms)
  end

  return dt.OdeSource(new_model.DynamicType, new_model.Constants,
                      new_model.Parameters, new_model.Forcings,
                      new_model.States, new_model.Derivatives,
                      new_model.Variables, new_model.Observed)
end

function expand_components(expression::ASCIIString, all_comp_forms)
    if !contains(expression, "{")
      return expression
    else
      # ID the species being matched by patterns
      species = matchall(r"[a-zA-Z0-9_]+(?=\{)", expression)
      contents = matchall(r"(?<=\{)[a-zA-Z0-9_\[\]]+(?=\})", expression)
      for k in 1:length(species)
        # Extract the component forms for a particular species
        allcompforms = all_comp_forms[species[k]]
        # Extract the elements of the patterns
        compforms = extract_pattern("$(species[k]){$(contents[k])}")
        # Iterate over all columns and keep track of matched rows
        h = [1:size(allcompforms, 1)]
        for i in compforms
          for j in 1:size(allcompforms, 2)
            if i in allcompforms[:,j]
              # Reduce h to those rows that match
              h = intersect(h, findin(allcompforms[:,j], [i]))
              break
            else
              continue
            end
          end
        end
        final_forms = allcompforms[h,:]
        expanded = "($(species[k])_" * paste("", vec(final_forms[1,:]))
        for i in 2:size(final_forms, 1)
          expanded *= " + $(species[k])_" * paste("", vec(final_forms[i,:]))
        end
        expanded *= ")"
        expression = replace(expression, "$(species[k]){$(contents[k])}", expanded)
      end
      return expression
    end
end

function species_to_states(new_model::dt.ReactionSource)
  new_states = OrderedDict{ASCIIString, dt.State}()
  new_species = OrderedDict{ASCIIString, dt.Species}()
  new_variables = OrderedDict{ASCIIString, dt.Variable}()
  all_mat_comp_forms = OrderedDict{ASCIIString, Array{ASCIIString,2}}()
  # Iterate over species
  for (key,val) in new_model.Species
    if length(val.Components) == 0
      new_states[key] = dt.State(val.Value, val.Units, val.QSS)
    # If the chemical species contains components, calculate all possible forms and generate the names
    else
      # Total number of states that will be generated
      number = prod([length(i.Forms) for i in val.Components])
      # Jagged array of component-forms
      components_forms = [ ASCIIString[i.Name*j for j in i.Forms ] for i in val.Components]
      # Jagged array of component-values
      components_values = [ i.Values for i in val.Components]
      # Matrix with all the component_forms
      # Generate all possible combinations using repeat and the correct inner/outer combinations
      mat_comp_forms = Array(ASCIIString, number, length(components_forms))
      small = components_forms[1]
      mat_comp_forms[:,1] = repeat(small, outer = [integer(number/length(small))])
      for i in 2: length(components_forms)
        try
          small = repeat(components_forms[i], inner = [length(small)])
          mat_comp_forms[:,i] = repeat(small, outer = [integer(number/length(small))])
        catch
          println("Error in generating component-forms matrix: species was $key.")
          println("Component_forms were $components_forms")
          println("New assignment was $(repeat(small, outer = [integer(number/length(small))]))")
        end
      end
      # Matrix with all the component_values
      # Generate all possible combinations using repeat and the correct inner/outer combinations
      mat_comp_values = Array(Float64, number, length(components_values))
      small = components_values[1]
      mat_comp_values[:,1] = repeat(small, outer = [integer(number/length(small))])
      for i in 2: length(components_values)
        try
          small = repeat(components_values[i], inner = [length(small)])
          mat_comp_values[:,i] = repeat(small, outer = [integer(number/length(small))])
        catch
          println("Error in generating component-values matrix: species was $key.")
          println("Component_values were $components_values")
          println("New assignment was $(repeat(small, outer = [integer(number/length(small))]))")
        end
      end
      # Reduce the forms and numbers to vectors. Concatenate without intermediate symbol
      forms = Array(ASCIIString, number)
      values = Array(Float64, number)
      for i in 1:number
        forms[i] = key * "_" * prod(mat_comp_forms[i,:])
        values[i] = prod(mat_comp_values[i,:])
        new_states[forms[i]] = dt.State(values[i], val.Units, val.QSS)
        new_species[forms[i]] = dt.Species(values[i], val.Compartment, val.Units, {}, val.QSS)
      end
      # Create a variable with the name of the species that is equal to the sum of all components
      new_variables[key] = dt.Variable(paste(" + ", forms)[1:(end-2)] , val.Units)
      all_mat_comp_forms[key] = mat_comp_forms
    end
  end
  # Return the new_states and new_variables, as well as intermediate calculations that may be useful later on
  return new_states, new_variables, new_species, all_mat_comp_forms

end

function reactions_to_variables(new_model::dt.ReactionSource,
                                new_states::OrderedDict{ASCIIString,dt.State},
                                new_variables::OrderedDict{ASCIIString,dt.Variable},
                                all_mat_comp_forms::OrderedDict{ASCIIString, Array{ASCIIString,2}})
  # Create empty Derivative objects for each element in new_states
  # Use d_state_dt syntax to generate the new derivatives
  new_derivatives = OrderedDict{ASCIIString, dt.Derivative}()
  for (key,val) in new_states
    if "d_$(key)_dt" in collect(keys(new_model.Derivatives))
      new_derivatives["d_$(key)_dt"] = new_model.Derivatives["d_$(key)_dt"]
    else
      new_derivatives["d_$(key)_dt"] = dt.Derivative("", key, val.Units)
    end
  end

  # Iterate over all reactions and populate associated new_variables and new_derivatives
  # It depends very much on whether we are using patterns in the stoichiometry or not
  for (key,val) in new_model.Reactions

    # In the case that we have patterns in the reaction stoichiometry
    if reaction_contains_components(new_model, val)
      # Id chemical species on either side of the reaction that contain components
      substrate_pattern, substrate_pattern_species::ASCIIString = who_has_components(new_model, val.Substrates)
      product_pattern, product_pattern_species::ASCIIString = who_has_components(new_model, val.Products)
      # Only allow reactions where there is one chemical species per side of reaction and they are the same chemical species
      (!isa(substrate_pattern_species, ASCIIString) || !isa(product_pattern_species, ASCIIString)) && error("Error in reaction $key: I cannot process this reaction. Too many chemical species with components.")
      substrate_pattern_species != product_pattern_species && error("Error in reaction $key: The species with components on either side of the reaction must be the same.")
      # Create matrix of substrate/product pairs from substrates and products that have patterns
      states_reactants_forms = states_from_master_equation(substrate_pattern, product_pattern, substrate_pattern_species, all_mat_comp_forms)
      # Every row will generate a new reaction
      # Any mention of a pattern in the expression is substited by the expanded state
      new_reactions = OrderedDict{ASCIIString, dt.Reaction}()
      for i in 1:size(states_reactants_forms, 1)
        expanded_form = vec(states_reactants_forms[i,:])
        # Replace references to the substrate in the mathematical expression
        expanded_expression = replace(val.Expr, substrate_pattern.Name, expanded_form[1] )
        expanded_expression = replace(expanded_expression, product_pattern.Name, expanded_form[2] )
        # Create the substrates and products for the new reaction
        substrates = dt.Reactant[dt.Reactant(expanded_form[1], substrate_pattern.Stoichiometry)]
        for i in val.Substrates
          (i == substrate_pattern) && continue
          push!(substrates, i)
        end
        products = dt.Reactant[dt.Reactant(expanded_form[2], product_pattern.Stoichiometry)]
        for i in val.Products
          (i == product_pattern) && continue
          push!(products, i)
        end
        new_reactions["v_$(expanded_form[1])_$(expanded_form[2])"] = dt.Reaction(substrates, products, expanded_expression, val.Compartment, val.Units)
      end
      # Iterate over new reactions and execute the code for no-pattern reactions
      # to generate the new reactions from extending the components
      # Need to add a Variable that captures the sum of all this individual reactions
      # so that we can still make a reference to the original reaction with components
      sum_all = ""
      for (key2,val2) in new_reactions
        process_reaction_no_pattern!(new_model, val2, key2, new_variables, new_derivatives)
        sum_all *= "($(val2.Expr))+"
      end
      sum_all = sum_all[1:(end-1)]
      new_variables[key] = dt.Variable(sum_all, val.Units)
    # In the case that we don't have patterns in the reaction stoichiometry
    else
      process_reaction_no_pattern!(new_model, val, key, new_variables, new_derivatives)
    end
  end

  return new_variables, new_derivatives
end


# Check if a reaction stoichiometry contains species with components, regardless
# of whether there are patterns or not
function reaction_contains_components(model::dt.ReactionSource, reaction::dt.Reaction)
  # Do substrates contain a chemical species with components?
  for i in reaction.Substrates
    (contains(i.Name, "{") || length(model.Species[i.Name].Components) > 0) && return true
  end
  # Do products contain a chemical species with components?
  for i in reaction.Products
    (contains(i.Name, "{") || length(model.Species[i.Name].Components) > 0) && return true
  end
  return false
end

# Return the species that contains components in a set of reactants or products,
# regardless of whether there are patterns or not
function who_has_components(new_model::dt.ReactionSource, reactants::Vector{dt.Reactant})
  for i in reactants
    if contains(i.Name, "{") || (length(model.Species[i.Name].Components) > 0)
      species = match(r"[a-zA-Z0-9_]*(?=\{)", i.Name).match
      return i, species
    end
  end
end

# Calculate a 2-column matrix where each row is a substrate/product combination
function states_from_master_equation(substrate::dt.Reactant, product::dt.Reactant, species::ASCIIString,
                                     comp_forms::OrderedDict{ASCIIString, Array{ASCIIString,2}})
  # Extract the matrices for the substrates and products containing all the combinations of component forms
  substrate_components = deepcopy(comp_forms[species])
  product_components = deepcopy(comp_forms[species])
  # Extract patterns from substrate/products -> Array with the individual component-form
  substrate_patterns = extract_pattern(substrate.Name)
  product_patterns = extract_pattern(product.Name)
  # Select rows that meet the pattern (the components targeted by the different patterns should be the same)
  target_columns_s, substrate_components = subset_with_patterns(substrate_components, substrate_patterns)
  target_columns_p, product_components = subset_with_patterns(product_components, product_patterns)
  sort(target_columns_s) != sort(target_columns_p) &&
            error("Error in states_from_master_equation. The patterns on either side don't target the same components")
  # Which columns are not targeted by patterns?
  non_target_columns = filter(x -> !in(x, target_columns_s),  1:size(substrate_components, 2))
  # Every row in substrate_components is combined with each row in product_components where
  # the non-target columns have the same elements as the corresponding non-target columns in the substrate_components
  products = [product_from_substrate(substrate_components[i,:], non_target_columns, product_components) for i in 1:size(substrate_components, 1)]
  # Generate a matrix where the first column is the string representing the state for the substrate
  # and the second column is the string of the product
  substrate_product = ASCIIString[""  ""]
  for i in 1:length(products)
    sub_name = species * "_" * paste("", vec(substrate_components[i,:]))
    prod_name = ASCIIString[species * "_" * paste("", vec(products[i][j,:])) for j in 1:size(products[i], 1)]
    tmp = reshape(repeat([sub_name],outer = [length(prod_name)]),length(prod_name), 1)
    tmp = hcat(tmp, prod_name)
    substrate_product = vcat(substrate_product, tmp)
  end
  return substrate_product[2:end,:]
end

# Extract rows of component-form matrix in the product where non-target columns
# coincide with the non-target columns of the substrate
function product_from_substrate(substrate, non_target_columns, product_components)
  size(product_components, 1) == 1 && (return product_components)
  # Extract non-target columns of both substrates and products
  non_target_substrate = substrate[:,non_target_columns]
  non_target_product   = product_components[:,non_target_columns]
  # Intersection between substrate and product non-target component-forms
  matches = Int64[]
  for i in 1:size(non_target_product, 1)
    sort(vec(non_target_product[i,:])) == sort(vec(non_target_substrate)) && (push!(matches, i))
  end
  return product_components[matches,:]

end

# Given a pattern, extract the components in it
function extract_pattern(pat::ASCIIString)
  compforms = match(r"(?<=\{)[a-zA-Z0-9_\[\]]*(?=\})", pat).match
  forms = matchall(r"(?<=\[)[a-zA-Z0-9_]+(?=\])", compforms)
  components = matchall(r"[a-zA-Z0-9_]+(?=\[)", compforms)
  ASCIIString[components[i]*forms[i] for i in 1:length(forms)]
end

# Given a matrix with all possible combinations of component forms
# and a series of component forms, extract the subset of the matrix
function subset_with_patterns(substrate_components, substrate_patterns)
  target_columns = Int64[]
  # Iterate over all patterns
  for i in substrate_patterns
    # Find the column that contains the component
    col = 0
    for j in 1:size(substrate_components, 2)
      if i in substrate_components[:,j]
        # Store the column where the pattern is found and subset the matrix of component-forms
        # by the row where the match happened
        push!(target_columns, j)
        substrate_components = substrate_components[findin(substrate_components[:,j], [i]),:]
        break
      else
        continue
      end
    end
  end
  return target_columns, substrate_components
end


function process_reaction_no_pattern!(new_model::dt.ReactionSource, reaction::dt.Reaction, key::ASCIIString, new_variables, new_derivatives)
  new_variables[key] = dt.Variable(reaction.Expr, reaction.Units)
  for i in reaction.Substrates
    # Capture current expression for the associated time derivative
    previous_expr = new_derivatives["d_$(i.Name)_dt"].Expr
    # Extract the compartments of the reaction and the species that participates in the reaction
    compartment_reaction = reaction.Compartment
    compartment_reactant = new_model.Species[i.Name].Compartment
    # If the compartments are the same, don't generate conversion. If they are not, generate the conversion
    compartment_reaction != compartment_reactant ?
                            (ratio_comp = "*$compartment_reaction/$compartment_reactant") :
                            (ratio_comp = "")
    haskey(new_model.Parameters, compartment_reaction) ?
          compartment_reaction_units = new_model.Parameters[compartment_reaction].Units :
          compartment_reaction_units = new_model.States[compartment_reaction].Units
    haskey(new_model.Parameters, compartment_reactant) ?
          compartment_reactant_units = new_model.Parameters[compartment_reactant].Units :
          compartment_reactant_units = new_model.States[compartment_reactant].Units
    ratio_comp_units = compartment_reaction_units/compartment_reactant_units
    # The new expression where the rate equation is substracted from the time derivative (important to use ()!!)
    new_expr = previous_expr * "- $(i.Stoichiometry)($key)$ratio_comp"
    new_derivatives["d_$(i.Name)_dt"] = dt.Derivative(new_expr, i.Name, reaction.Units*ratio_comp_units)
  end
  for i in reaction.Products
    # Capture current expression for the associated time derivative
    previous_expr = new_derivatives["d_$(i.Name)_dt"].Expr
    # Extract the compartments of the reaction and the species that participates in the reaction
    compartment_reaction = reaction.Compartment
    compartment_reactant = new_model.Species[i.Name].Compartment
    # If the compartments are the same, don't generate conversion. If they are not, generate the conversion
    compartment_reaction != compartment_reactant ?
                            (ratio_comp = "*$compartment_reaction/$compartment_reactant") :
                            (ratio_comp = "")
    haskey(new_model.Parameters, compartment_reaction) ?
          compartment_reaction_units = new_model.Parameters[compartment_reaction].Units :
          compartment_reaction_units = new_model.States[compartment_reaction].Units
    haskey(new_model.Parameters, compartment_reactant) ?
          compartment_reactant_units = new_model.Parameters[compartment_reactant].Units :
          compartment_reactant_units = new_model.States[compartment_reactant].Units
    ratio_comp_units = compartment_reaction_units/compartment_reactant_units
    # The new expression where the rate equation is substracted from the time derivative (important to use ()!!)
    new_expr = previous_expr * "+ $(i.Stoichiometry)($key)$ratio_comp"
    new_derivatives["d_$(i.Name)_dt"] = dt.Derivative(new_expr, i.Name, reaction.Units*ratio_comp_units)
  end

end
