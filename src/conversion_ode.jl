
# Generate OdeModel from OdeSource
function OdeModel(model::dt.OdeSource)

  # Names of the derivatives ordered by associated states
  names_states = collect(keys(model.States))
  names_derivatives = ASCIIString[]
  for i in names_states
    for (key,val) in model.Derivatives
      val.State == i && (push!(names_derivatives, key))
    end
  end

  # We want to have easy access to the equations associated to derivatives and observed variables
  names_observed = collect(keys(model.Observed))

  # Create equations from observed, variable and derivatives
  equations = OrderedDict{ASCIIString, dt.Equation}()
  for (key,val) in model.Observed
    equations[key] = dt.Equation(parse(val.Expr), val.Units)
  end
  for (key,val) in model.Variables
    equations[key] = dt.Equation(parse(val.Expr), val.Units)
  end
  for (key,val) in model.Derivatives
    (length(val.Expr) == 0) && error("Error in derivative $key: There is no expression!!")
    equations[key] = dt.Equation(parse(val.Expr), val.Units)
  end

  return dt.OdeModel(model.DynamicType, model.Constants, model.Parameters,
                     model.Forcings, model.States, equations,
                     names_derivatives, names_observed)

end


# Generate OdeSorted from OdeModel
function OdeSorted(model::dt.OdeModel)
  # Before we sort the equations we need to check that there is a minimum consistency to it
  # which means that we need to run some tests and extract info for the sorting
  checked_ode_model, rhs_to_lhs, lhs_to_rhs = check_lhs_rhs(model)

  # Create an array to store the different levels of equations
  # and push in the first level for inputs
  SortedEquations = OrderedDict{ASCIIString, dt.Equation}[]
  push!(SortedEquations, generate_level0(checked_ode_model))

  # Now we apply the algorithm to sort the equations
    unsorted_equations = collect(keys(checked_ode_model.Equations))
    while length(unsorted_equations) > 0
        nsorted = length(unsorted_equations)
        unsorted_keys = unsorted_equations
        for j in unsorted_keys
          rhs = Set(lhs_to_rhs[j])
          level = 2
          for i in 1:length(SortedEquations)
            if length(rhs) == 0
              break
            end
            rhs = setdiff(rhs, intersect(rhs, Set(collect(keys(SortedEquations[i])))))
            level = i + 1
          end
          if length(rhs) == 0
            deleteat!(unsorted_equations, findfirst(unsorted_equations, string(j))[1])
            if level <= length(SortedEquations)
              SortedEquations[level][j] = checked_ode_model.Equations[j]
            else
              tmp = OrderedDict{ASCIIString, dt.Equation}()
              tmp[j] =  checked_ode_model.Equations[j]
              push!(SortedEquations, tmp)
            end
          end
        end
        if nsorted == length(unsorted_equations)
          error("There seem to be some circular dependencies")
        end
    end
  return dt.OdeSorted(model.DynamicType, model.Constants, model.Parameters,
                     model.Forcings, model.States, SortedEquations,
                     model.NamesDerivatives, model.NamesObserved)
end

# Generate the first level of sorted equations
# It assumes index starts at 1 (good for Julia, R and Fortran, needs modification in C, C++ and Python)
function generate_level0(model::dt.OdeModel)
  equations = OrderedDict{ASCIIString, dt.Equation}()
  c = 1
  for (key,val) in model.Parameters
    equations[key] = dt.Equation(parse(key * "= params[$c]"), val.Units)
    c += 1
  end
  c = 1
  for (key,val) in model.States
    equations[key] = dt.Equation(parse(key * "= states[$c]"), val.Units)
    c += 1
  end
  c = 1
  for (key,val) in model.Forcings
    equations[key] = dt.Equation(parse(key * "= forcs[$c]"), val.Units)
    c += 1
  end
  for (key,val) in model.Constants
    rhs = float(val.Value * val.Units.f)
    equations[key] = dt.Equation(parse(key * "= $(rhs)"), val.Units)
    c += 1
  end
  return equations
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


# Compress model to level 2 (return derivatives)
function compress_model(model::dt.OdeSorted; level = 2)
  equations = deepcopy(model.SortedEquations)
  if length(equations) < level + 1
    return model
  end
  for i in linrange(length(equations), level + 1, length(equations) - level)
    for (key,val) in equations[i]
      equations[i][key].Expr = expand_expression(val.Expr, equations)
      equations[i-1][key] = equations[i][key]
    end
  end
  new_model = deepcopy(model)
  new_model.SortedEquations = equations[1:level]
  lhs_names = collect(keys(new_model.SortedEquations[level]))
  names_derivatives = ["d_$(i)_dt" for i in collect(keys(model.States))]
  # Delete everything that is not a derivative
  for i in setdiff(lhs_names, names_derivatives)
      delete!(new_model.SortedEquations[level], i)
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
