####################################################################################
####################################################################################
############################       CHECK RHS      ##################################
####################################################################################
####################################################################################

function check_lhs_rhs(model)
    new_model = deepcopy(model)
    lhs_states, lhs_parameters, lhs_forcings, lhs_constants, lhs_equations = get_lhs(new_model)
    all_rhs, rhs_to_lhs, lhs_to_rhs = get_rhs(new_model.Equations)
    for i in all_rhs
        i ∉ [lhs_states, lhs_parameters, lhs_forcings, lhs_constants, lhs_equations] &&
          error("The input $i is missing from the model")
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
            warn("The output is not dependent on the state $i, but it has NOT been removed.")
        end
    end
    for i in lhs_equations
        if (i ∉ all_rhs) && (i ∉ [model.NamesObserved, model.NamesDerivatives])
            #warn("The output is not dependent on the intermediate variable $i. It has been removed from the model.")
            delete!(new_model.Equations, i)
        end
    end
    return new_model, rhs_to_lhs, lhs_to_rhs
end


####################################################################################
####################################################################################
###########################       CHECK UNITS      #################################
####################################################################################
####################################################################################

function calculate_units(exp::Expr, units::Dict{Symbol, dt.Dimension}, c::Int64)
  ex = deepcopy(exp)
  for i in 1:length(ex.args)
    if isa(ex.args[i], Expr)
      ex.args[i], c = calculate_units(ex.args[i], units, c)
    elseif isa(ex.args[i], Symbol) && haskey(units, symbol(ex.args[i]))
      c += 1
      name = parse("__unit__$(c)")
      @eval global $name = ($(units[ex.args[i]]))
      ex.args[i] = name
    end
  end
  ex, c
end

function check_units(sorted_model::dt.OdeSorted)
  given_units = Dict{Symbol, dt.Dimension}()
  for (key,val) in sorted_model.Constants  given_units[symbol(key)] = val.Units.d end
  for (key,val) in sorted_model.Parameters  given_units[symbol(key)] = val.Units.d end
  for (key,val) in sorted_model.States  given_units[symbol(key)] = val.Units.d end
  for (key,val) in sorted_model.Forcings  given_units[symbol(key)] = val.Units.d end
  for level in sorted_model.SortedEquations
    for (key,val) in level  given_units[symbol(key)] = val.Units.d end
  end
  # Add units of time (as input that is always available)
  given_units[symbol("time")] = dt.s.d
  for level in sorted_model.SortedEquations[2:end]
    for (key,val) in level
      expected = given_units[symbol(key)]
      if isa(val.Expr, Symbol)
            infered = eval(given_units[val.Expr])
      else
        infered_expr, c = calculate_units(val.Expr, given_units, 0)
        infered = try eval(infered_expr) catch
            error("Error when calculating units: The rhs of equation $key is not dimensionally homogeneous.
The right hand side expression was $(val.Expr)") end
      end
      expected != infered && error("Error when calculating units: Expected and given units for $key did not coincide.
I infered the units $infered but you assign it $expected.
The right hand side expression was $(val.Expr)")
    end
  end
end

####################################################################################
####################################################################################
########################       GENERATE JACOBIAN      ##############################
####################################################################################
####################################################################################


# Calculation Jacobian matrix of the model
function generate_jacobian_matrix(compressed_model::dt.OdeSorted)
  names_states = collect(keys(compressed_model.States))
  Jacobian = Array(Union(Expr, Symbol, Number),(length(names_states), length(names_states)))
  cs = 1
    for j in names_states
      cd = 1
      for i in compressed_model.NamesDerivatives
        @inbounds Jacobian[cd,cs] =  differentiate(compressed_model.SortedEquations[2][i].Expr, parse(j))
        cd += 1
      end
    cs += 1
  end
  return Jacobian
end

####################################################################################
####################################################################################
###################       GENERATE EXTENDED SYSTEM      ############################
####################################################################################
####################################################################################

# Calculate array of sensitivities
function generate_sensitivity_array(compressed_model::dt.OdeSorted)
  names_states = collect(keys(compressed_model.States))
  names_derivatives = ["d_"*x*"_dt" for x in names_states]
  Sensitivity = Array{Union(Expr, Symbol, Number), 1}[]
  names_parameters = collect(keys(compressed_model.Parameters))
  for i in names_parameters
    sens = Array(Union(Expr, Symbol, Number), length(names_states))
    c = 1
    for j in names_derivatives
      sens[c] = differentiate(compressed_model.SortedEquations[2][j].Expr, parse(i))
      c += 1
    end
    push!(Sensitivity, sens)
  end
  return Sensitivity
end
