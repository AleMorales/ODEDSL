########## Check units ###############

function calculate_units(exp::Expr, units::Dict{Symbol, Dimension}, c::Int64)
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

function check_units(sorted_model::OdeSorted)
  given_units = Dict{Symbol, Dimension}()
  for (key,val) in sorted_model.Constants  given_units[symbol(key)] = val.Units.d end
  for (key,val) in sorted_model.Parameters  given_units[symbol(key)] = val.Units.d end
  for (key,val) in sorted_model.States  given_units[symbol(key)] = val.Units.d end
  for (key,val) in sorted_model.Forcings  given_units[symbol(key)] = val.Units.d end
  for level in sorted_model.SortedEquations
    for (key,val) in level  given_units[symbol(key)] = val.Units.d end
  end
  # Add units of time (as input that is always available)
  given_units[symbol("time")] = ODEDSL.DataTypes.s.d
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

########## Generate Jacobian ###############

using Calculus

function generate_jacobian_function_Julia(model::OdeSorted, name)
  jacobian_matrix = generate_jacobian_matrix(model)
  string_assignments = [string(x[2].Expr) for x in model.SortedEquations[1]]
  code = string_assignments[1]*"\n"
  for i = 2:length(string_assignments)
    code *= string_assignments[i]*"\n"
  end
  for i = 1:size(jacobian_matrix)[1], j = 1:size(jacobian_matrix)[1]
    if isa(jacobian_matrix[i,j], Number) && eval(jacobian_matrix[i,j]) == 0
      continue
    else
      code *= "J[$i,$j] = $(jacobian_matrix[i,j])\n"
    end
  end
  return_line = "return nothing\n"
  function_text = paste("\n",
  "@inbounds function jacobian_$name(time::Float64, states::Array{Float64,1},
   params::Array{Float64,1}, forcs::Array{Float64,1}, J)\n", code, return_line,"end")
  return function_text
end

# Calculation Jacobian matrix of the model
function generate_jacobian_matrix(compressed_model::OdeSorted)
  names_states = collect(keys(compressed_model.States))
  names_derivatives = ["d_"*x*"_dt" for x in names_states]
  Jacobian = Array(Union(Expr, Symbol, Number),(length(names_states), length(names_states)))
  cd = 1
  for i in names_derivatives
    cs = 1
    for j in names_states
      Jacobian[cd,cs] =  differentiate(compressed_model.SortedEquations[2][i].Expr, parse(j))
      cs += 1
    end
    cd += 1
  end
  return Jacobian
end

########## Generate Extended System ###############

# Calculate extended system
function generate_extended_system(compressed_model::OdeSorted, name)
  sens_array = generate_sensitivity_array(compressed_model)
  # Create an extended ode system

  # Calculate the Jacobian of the extended system
  return "dummy_function", "dummy_jacobian"
end
# Calculate array of sensitivities
function generate_sensitivity_array(compressed_model::OdeSorted)
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


####################################################################################
####################################################################################
##########################   JULIA CODE GENERATION #################################
####################################################################################
####################################################################################

# Create the function in Julia on the Equations section of the model Dict
function create_derivatives_julia(model::OdeSorted, states, name)
  code = ""
  names_derivatives = ["d_"*i*"_dt" for i in states]
  names_forcings = collect(keys(model.Forcings))
  for level in 1:length(model.SortedEquations)
    for (lhs, rhs) in model.SortedEquations[level]
        if level == 1
          if in(lhs, names_forcings)
            code *= "@inbounds const " * string(rhs.Expr) * "[time]\n"
          else
            code *= "@inbounds const " * string(rhs.Expr) * "\n"
          end
        elseif in(lhs, names_derivatives)
          c = findin(names_derivatives, [lhs])
          code  *=  "@inbounds ydot$c = " * string(rhs.Expr) * "\n"
        else
          code *= "const " * lhs * " = " * string(rhs.Expr) * "\n"
        end
    end
  end
  # Construct the function
  paste("\n","function derivatives_$name(time::Float64, states::Vector{Float64},
                      params::Vector{Float64}, forcs::Vector{sim.Forcing}, ydot::Vector{Float64})\n",
                      code,"nothing \n end")
end

# Create the function in Julia on the Equations section of the model Dict
function create_observed_julia(model::OdeSorted, observed, name)
  code = ""
  names_forcings = collect(keys(model.Forcings))
  for level in 1:length(model.SortedEquations)
    for (lhs, rhs) in model.SortedEquations[level]
        if level == 1
          if in(lhs, names_forcings)
            code *= "@inbounds const " * string(rhs.Expr) * "[time]\n"
          else
            code *= "@inbounds const " * string(rhs.Expr) * "\n"
          end
        elseif in(lhs, observed)
          c = findin(observed, [lhs])
          code  *=  "const $lhs = " * string(rhs.Expr) * "\n"
          code  *=  "@inbounds obs$c = $lhs\n"
        else
          code *= lhs * " = " * string(rhs.Expr) * "\n"
        end
    end
  end
  paste("\n","function observed_$name(time::Float64, states::Vector{Float64},
                      params::Vector{Float64}, forcs::Vector{sim.Forcing}, obs::Vector{Float64})\n",
                      code, "nothing \n end")
end

function generate_code_Julia!(ode_model::OdeSource; unit_analysis = true, name = "autogenerated_model", file = "autogenerated_model", jacobian = false, sensitivities = false)
  # Generate the observed variables (everything that is exported but it is not a time derivative)
  observed = ASCIIString[]
  for (key,val) in ode_model.Equations
    val.Exported && push!(observed, key)
  end
  names_derivatives = collect(keys(ode_model.States))
  for i in 1:length(names_derivatives)
    names_derivatives[i] = "d_"*names_derivatives[i]*"_dt"
  end
  deleteat!(observed, findin(observed, names_derivatives))
  coef_observed = OrderedDict{ASCIIString, Float64}()
  for i in observed
    coef_observed[i] = ode_model.Equations[i].Units.f
  end

  # Sort the equations
  sorted_model = sort_equations(ode_model);

  # Create the default arguments
  named_states = OrderedDict{ASCIIString, Float64}()
  coef_states = OrderedDict{ASCIIString, Float64}()
  for (key,val) in sorted_model.States
    named_states[key] = val.Value * val.Units.f
    coef_states[key] = val.Units.f
  end
  named_parameters = OrderedDict{ASCIIString, Float64}()
  coef_parameters = OrderedDict{ASCIIString, Float64}()
  for (key,val) in sorted_model.Parameters
    named_parameters[key] = val.Value * val.Units.f
    coef_parameters[key] = val.Units.f
  end
  forcings = OrderedDict{ASCIIString,Any}()
  coef_forcings = OrderedDict{ASCIIString,Float64}()
  c = 0
  for (key,value) in sorted_model.Forcings
      c += 1
      forcings[key] = (float(value.Time), float(value.Value)*value.Units.f)
      coef_forcings[key] = value.Units.f
  end

  # Created compressed model (only if Jacobian or Sensitivities are required!)
  jacobian_function = "jacobian_$name() = nothing"
  sensitivity_function = "() -> ()"
  sensitivity_jacobian_function = "() -> ()"
  if jacobian || sensitivities
    compressed_model = compress_model(sorted_model, level = 2)
    jacobian && (jacobian_function = generate_jacobian_function(compressed_model, name))
    sensitivities && ((sensitivity_function, sensitivity_jacobian_function) = generate_extended_system(compressed_model, name))
  end

  # Check the units
  unit_analysis && check_units(sorted_model)

  # Generate the rhs function
  model_function = create_derivatives_julia(sorted_model,collect(keys(named_states)), name)
  # Generate function with the observed values
  observed_function = create_observed_julia(sorted_model,observed,name)

  write_model_Julia!(named_states,coef_states, named_parameters, coef_parameters,
                     forcings, coef_forcings, observed,coef_observed,
                      model_function,observed_function, jacobian_function, name, file)
  nothing
end

function write_model_Julia!(States::OrderedDict{ASCIIString, Float64},
                            Coef_states::OrderedDict{ASCIIString, Float64},
                            Parameters::OrderedDict{ASCIIString, Float64},
                            Coef_parameters::OrderedDict{ASCIIString, Float64},
                            Forcings::OrderedDict{ASCIIString, Any},
                            Coef_forcings::OrderedDict{ASCIIString, Float64},
                            Observed::Array{ASCIIString, 1},
                            Coef_observed::OrderedDict{ASCIIString, Float64},
                            Model::ASCIIString,
                            Observed_model::ASCIIString,
                            Jacobian::ASCIIString,
                            name::ASCIIString,
                            file::ASCIIString)
    f = open("$(file).jl","w")
    println(f, "import SimulationModels; using Dierckx; sim = SimulationModels")
    println(f, "function generate_$name()")
    states_values = Float64[]
    states_names = ASCIIString[]
    states_coefs = Float64[]
    for (key,val) in States
      push!(states_values, val)
      push!(states_names, key)
      push!(states_coefs, Coef_states[key])
    end
    println(f, "State_names = $(string(states_names)[12:end])")
    println(f, "States = sim.InputVector{Float64}($states_values, $states_coefs, $(string(states_names)[12:end]))")
    parameters_values = Float64[]
    parameters_names = ASCIIString[]
    parameters_coefs = Float64[]
    for (key,val) in Parameters
      push!(parameters_values, val)
      push!(parameters_names, key)
      push!(parameters_coefs, Coef_parameters[key])
    end
    println(f, "Parameters = sim.InputVector{Float64}($parameters_values, $parameters_coefs, $(string(parameters_names)[12:end]))")
    forcings_names = ASCIIString[]
    forcing_coefs = Float64[]
    for (key, val) in Forcings
      println(f, "$key = sim.Forcing(Spline1D($(val[1]), $(val[2]), k = 1))")
      push!(forcings_names, key)
      push!(forcing_coefs, Coef_forcings[key])
    end
    list_of_forcings = replace("$forcings_names", "\"", "")
    println(f, "Forcings = sim.InputVector{sim.Forcing}($(list_of_forcings[12:end]), $forcing_coefs, $(string(forcings_names)[12:end]))")
    println(f, "Observed = $(string(Observed))")
    observed_coef = Float64[]
    for (key,val) in Coef_observed
      push!(observed_coef, val)
    end
    println(f, "Observed_coefs = $observed_coef")
    println(f, "$Model")
    println(f, "$(Observed_model)")
    println(f, "$Jacobian")
    println(f, "sim.OdeModel(States, Parameters, Forcings, State_names, Observed, Observed_coefs, derivatives_$name, observed_$name, jacobian_$name)")
    println(f, "end")
    close(f)
    nothing
end

function generate_code_Julia!(source::String; unit_analysis = false, name = "autogenerated_model", file = "autogenerated_model", jacobian = false, sensitivities = false)
  parsed_model = process_file(source)
  reaction_model = convert_master_equation(parsed_model)
  ode_model = convert_reaction_model(reaction_model)
  generate_code_Julia!(ode_model, unit_analysis = unit_analysis, name = name, file = file, jacobian = jacobian, sensitivities = sensitivities)
end



####################################################################################
####################################################################################
############################   R CODE GENERATION ###################################
####################################################################################
####################################################################################

# Because Julia will use "pretty printing" for scalar product, we need to fool the parser
# by temporarily substituying by another binary operator and then doing string replacement
function sub_product(ex::Expr)
  for i in 1:length(ex.args)
  if ex.args[i] == :*
    ex.args[i] = :.*
  elseif isa(ex.args[i], Expr)
     ex.args[i] = sub_product(ex.args[i])
  end
end
return ex
end

sub_product(ex::Symbol) = ex

# Create a return line when the output is a list with two numeric vectors
function create_return_line_R(states, observed)
    return_line = "return(list(c("
    for i in states
      if i != states[end]
        return_line *= i * ","
      else
        return_line *= i
      end
    end
    return_line *= ")"
    length(observed) > 0 && (return_line *= ", c(")
    for i in observed
      if i != observed[end]
        return_line = return_line * i * ", "
      else
        return_line = return_line  * i
      end
    end
    length(observed) > 0 && (return_line *= ")")
    return_line = return_line * "))"
end

# Create the function in R on the Equations section of the model Dict
function create_function_R!(model::OdeSorted, observed)
  code = ""
  for level in 1:length(model.SortedEquations)
    for (lhs, rhs) in model.SortedEquations[level]
        if level == 1
            code *= string(rhs.Expr) * "\n"
        else
                code *= lhs * " = " * replace(string(rhs.Expr), ".*", "*") * "\n"
        end
    end
  end
  # Determine what the time derivatives are
  time_derivatives = String[]
  for i in collect(keys(model.States))
    push!(time_derivatives, "d_"*i*"_dt")
  end
  # Return line
  return_line = create_return_line_R(time_derivatives,observed)
  # Return the output
  code = replace(code, "1.0 *", "")
  return  paste("\n","function(time, states, params, forcs) { \n", code, return_line,"}")
end

function generate_code_R!(ode_model::OdeSource; unit_analysis = true, name = "autogenerated_model", file = "autogenerated_model", jacobian = false, sensitivities = false)
  # Generate the observed variables (everything that is exported but it is not a time derivative)
  observed = ASCIIString[]
  for (key,val) in ode_model.Equations
    val.Exported && push!(observed, key)
  end
  names_derivatives = collect(keys(ode_model.States))
  for i in 1:length(names_derivatives)
    names_derivatives[i] = "d_"*names_derivatives[i]*"_dt"
  end
  deleteat!(observed, findin(observed, names_derivatives))
  coef_observed = OrderedDict{ASCIIString, Float64}()
  for i in observed
    coef_observed[i] = ode_model.Equations[i].Units.f
  end

  # Sort the equations
  sorted_model = sort_equations(ode_model)

  # Created compressed model (only if Jacobian or Sensitivities are required!)
  jacobian_function = "() -> ()"
  sensitivity_function = "() -> ()"
  sensitivity_jacobian_function = "() -> ()"
  if jacobian || sensitivities
    compressed_model = compress_model(sorted_model, level = 2)
    jacobian && (jacobian_function = generate_jacobian_function(compressed_model, name))
    sensitivities && ((sensitivity_function, sensitivity_jacobian_function) = generate_extended_system(compressed_model, name))
  end

  # Check the units
  unit_analysis &&  check_units(sorted_model)

  # Go through the equations and substitute * by ×
  for i in 1:length(sorted_model.SortedEquations)
      for (key,val) in sorted_model.SortedEquations[i]
         sorted_model.SortedEquations[i][key].Expr = sub_product(sorted_model.SortedEquations[i][key].Expr)
      end
  end

  # Generate the rhs function
  model_function = create_function_R!(sorted_model,observed)

  # Create the default arguments
  named_states = OrderedDict{ASCIIString, Float64}()
  coef_states = OrderedDict{ASCIIString, Float64}()
  for (key,val) in sorted_model.States
    named_states[key] = val.Value * val.Units.f
    coef_states[key] = val.Units.f
  end
  named_parameters = OrderedDict{ASCIIString, Float64}()
  coef_parameters = OrderedDict{ASCIIString, Float64}()
  for (key,val) in sorted_model.Parameters
    named_parameters[key] = val.Value * val.Units.f
    coef_parameters[key] = val.Units.f
  end
  forcings = OrderedDict{ASCIIString,Any}()
  coef_forcings = OrderedDict{ASCIIString,Float64}()
  c = 0
  for (key,value) in sorted_model.Forcings
      c += 1
      forcings[key] = (float(value.Time), float(value.Value)*value.Units.f)
      coef_forcings[key] = value.Units.f
  end
  write_model_R!(named_states,coef_states, named_parameters, coef_parameters,
                 forcings, coef_forcings, observed,coef_observed,
                  model_function, name, file)
  nothing
end

function write_model_R!(States::OrderedDict{ASCIIString, Float64},
                        Coef_states::OrderedDict{ASCIIString, Float64},
                        Parameters::OrderedDict{ASCIIString, Float64},
                        Coef_parameters::OrderedDict{ASCIIString, Float64},
                        Forcings::OrderedDict{ASCIIString, Any},
                        Coef_forcings::OrderedDict{ASCIIString, Float64},
                        Observed::Array{ASCIIString, 1},
                        Coef_observed::OrderedDict{ASCIIString, Float64},
                        Model::ASCIIString,
                        name::ASCIIString,
                        file::ASCIIString)
    f = open("$(file).R","w")
    println(f, "library(SimulationModels); library(RcppSundials)")
    println(f, "$name <- ODEmodel\$new(")
    transformed_states = string(States)[2:(end-1)]
    transformed_states = replace(transformed_states, "=>", "=")
    units = replace(string(Coef_states)[2:(end-1)], "=>", "=")
    println(f, "States = list(Values = c($transformed_states), Coefs = c($units)),")
    transformed_parameters = string(Parameters)[2:(end-1)]
    transformed_parameters = replace(transformed_parameters, "=>", "=")
    units = replace(string(Coef_parameters)[2:(end-1)], "=>", "=")
    println(f, "Parameters = list(Values = c($transformed_parameters), Coefs = c($units)),")
    if length(Forcings) > 0
      println(f, "Forcings = list(Values = list(")
      forcs = ""
      for (key,val) in Forcings
        forcs *= "$key = cbind(c($(string(val[1])[2:(end-1)])),c($(string(val[2])[2:(end-1)]))),"
      end
      forcs = forcs[1:(end-1)]
      println(f, forcs)
      units = replace(string(Coef_forcings)[2:(end-1)], "=>", "=")
      println(f, "), Coefs = c($units)),")
    end
    println(f, "Time = 0:1,")
    names_observed = string(Observed)[13:(end-1)]
    units = replace(string(Coef_observed)[2:(end-1)], "=>", "=")
    println(f, "Observed = list(names = c($names_observed), Coefs = c($units)),")
    println(f, """
    Settings = list(rtol = 1e-6,atol = 1e-10, maxsteps = 1e5, maxord = 5, hini = 0,
                      hmin = 0, hmax = 0, maxerr = 12, maxnonlin = 12,
                      maxconvfail = 12, method = "bdf", maxtime = 0, jacobian = 0),
""")
    println(f, "model = $Model)")
    close(f)
    nothing
end


function generate_code_R!(source::String; unit_analysis = false,name = "autogenerated_model", file = "autogenerated_model", jacobian = false, sensitivities = false)
  parsed_model = process_file(source)
  reaction_model = convert_master_equation(parsed_model)
  ode_model = convert_reaction_model(reaction_model)
  generate_code_R!(ode_model, unit_analysis = unit_analysis, name = name, file = file, jacobian = jacobian, sensitivities = sensitivities)
end

####################################################################################
####################################################################################
############################   RC++ CODE GENERATION ################################
####################################################################################
####################################################################################

# Create a return line with an array containing the derivatives and observed functions (for the STL version of the model)
function create_return_line_Rcpp(states, observed)
    return_line = "vector<double> derivatives{"
    for i in states
      if i != states[end]
        return_line = return_line  * i * ","
      else
        return_line = return_line  * i
      end
    end
    return_line = return_line * "};\n vector<double> observed{"
    for i in observed
      if i != observed[end]
        return_line = return_line * i * ", "
      else
        return_line = return_line  * i
      end
    end
    return_line = return_line * "};\n array<vector<double>,2> output{derivatives, observed};\n return output;"
end

# Substitute the annoying x^y for pow(x,y)
function substitute_power(ex::Expr)
  new_ex = deepcopy(ex)
  for i in 1:length(new_ex.args)
    if isa(new_ex.args[i], Expr)
      new_ex.args[i] = substitute_power(new_ex.args[i])
    elseif new_ex.args[i] == :^
      new_ex.args[i] = :pow
    end
  end
  return new_ex
end

substitute_power(ex::Any) = ex

# Create the function in C++ using STL with the RcppSundials API
function create_function_Rcpp!(model::OdeSorted, observed, name)
  code = ""
  for level in 1:length(model.SortedEquations)
    for (lhs, rhs) in model.SortedEquations[level]
        if level == 1
          if(ismatch(r"params|states|forcs", string(rhs.Expr)))
           number = match(r"(?<=\[)[\d]+(?=\])",string(rhs.Expr))
           new_expr = replace(string(rhs.Expr), "[$(number.match)]", "[$(number.match)-1]")
           code *= "const double $(new_expr);\n"
         else
           code *= "const double $(rhs.Expr);\n"
         end
        else
          mod_expr = substitute_power(rhs.Expr)
          code *= "const double $lhs" * " = " * replace(string(mod_expr), ".*", "*") * ";\n"
        end
    end
  end
  # Determine what the time derivatives are
  time_derivatives = String[]
  for i in collect(keys(model.States))
    push!(time_derivatives, "d_"*i*"_dt")
  end
  # Return line
  return_line = create_return_line_Rcpp(time_derivatives,observed)
  # Return the output
  code = replace(code, "1.0 *", "")
  up_boiler_plate =
"""
#include <array>
#include <vector>
#include <math.h>
using namespace std;

extern "C" {
  array<vector<double>, 2> $(name)(const double& t, const vector<double>& states,
            const vector<double>& params, const vector<double>& forcs) { \n
"""
  low_boiler_plate =
"""
\n
}
};
"""
  return  paste("\n",up_boiler_plate, code, return_line,low_boiler_plate)
end

function generate_code_Rcpp!(ode_model::OdeSource; unit_analysis = true, name = "autogenerated_model", file = "autogenerated_model", jacobian = false, sensitivities = false)
  # Generate the observed variables (everything that is exported but it is not a time derivative)
  observed = ASCIIString[]
  for (key,val) in ode_model.Equations
    val.Exported && push!(observed, key)
  end
  names_derivatives = collect(keys(ode_model.States))
  for i in 1:length(names_derivatives)
    names_derivatives[i] = "d_"*names_derivatives[i]*"_dt"
  end
  deleteat!(observed, findin(observed, names_derivatives))
  coef_observed = OrderedDict{ASCIIString, Float64}()
  for i in observed
    coef_observed[i] = ode_model.Equations[i].Units.f
  end

  # Sort the equations
  sorted_model = sort_equations(ode_model)

  # Created compressed model (only if Jacobian or Sensitivities are required!)
  jacobian_function = "() -> ()"
  sensitivity_function = "() -> ()"
  sensitivity_jacobian_function = "() -> ()"
  if jacobian || sensitivities
    compressed_model = compress_model(sorted_model, level = 2)
    jacobian && (jacobian_function = generate_jacobian_function(compressed_model, name))
    sensitivities && ((sensitivity_function, sensitivity_jacobian_function) = generate_extended_system(compressed_model, name))
  end

  # Check the units
  unit_analysis &&  check_units(sorted_model)

  # Go through the equations and substitute * by ×
  for i in 1:length(sorted_model.SortedEquations)
      for (key,val) in sorted_model.SortedEquations[i]
         sorted_model.SortedEquations[i][key].Expr = sub_product(sorted_model.SortedEquations[i][key].Expr)
      end
  end

  # Generate the rhs function
  model_function = create_function_Rcpp!(sorted_model,observed, name)

  # Create the default arguments
  named_states = OrderedDict{ASCIIString, Float64}()
  coef_states = OrderedDict{ASCIIString, Float64}()
  for (key,val) in sorted_model.States
    named_states[key] = val.Value * val.Units.f
    coef_states[key] = val.Units.f
  end
  named_parameters = OrderedDict{ASCIIString, Float64}()
  coef_parameters = OrderedDict{ASCIIString, Float64}()
  for (key,val) in sorted_model.Parameters
    named_parameters[key] = val.Value * val.Units.f
    coef_parameters[key] = val.Units.f
  end
  forcings = OrderedDict{ASCIIString,Any}()
  coef_forcings = OrderedDict{ASCIIString,Float64}()
  c = 0
  for (key,value) in sorted_model.Forcings
      c += 1
      forcings[key] = (float(value.Time), float(value.Value)*value.Units.f)
      coef_forcings[key] = value.Units.f
  end
  write_model_Rcpp!(named_states,coef_states, named_parameters, coef_parameters,
                 forcings, coef_forcings, observed,coef_observed,
                  model_function, name, file)
  nothing
end


function write_model_Rcpp!(States::OrderedDict{ASCIIString, Float64},
                        Coef_states::OrderedDict{ASCIIString, Float64},
                        Parameters::OrderedDict{ASCIIString, Float64},
                        Coef_parameters::OrderedDict{ASCIIString, Float64},
                        Forcings::OrderedDict{ASCIIString, Any},
                        Coef_forcings::OrderedDict{ASCIIString, Float64},
                        Observed::Array{ASCIIString, 1},
                        Coef_observed::OrderedDict{ASCIIString, Float64},
                        Model::ASCIIString,
                        name::ASCIIString,
                        file::ASCIIString)
    # Create the C++ file
    f = open("$(file).cpp","w")
    println(f, Model)
    close(f)
    f = open("$(file).R","w")
    println(f,
"""library(SimulationModels); library(RcppSundials)

  # These are the pointers to the C++ functions
  # Note that in the current version the Jacobian is a dummy function (no content)
  system("R CMD SHLIB $(file).cpp -o $(file).so")
  dyn.load("$(file).so")
  # so simulations must always be run with the option Jacobian 0
  $(name)_pointer = getNativeSymbolInfo(name = "$name",PACKAGE = "$file")\$address
""")
    println(f, "$name <- ODEmodel\$new(")
    transformed_states = string(States)[2:(end-1)]
    transformed_states = replace(transformed_states, "=>", "=")
    units = replace(string(Coef_states)[2:(end-1)], "=>", "=")
    println(f, "States = list(Values = c($transformed_states), Coefs = c($units)),")
    transformed_parameters = string(Parameters)[2:(end-1)]
    transformed_parameters = replace(transformed_parameters, "=>", "=")
    units = replace(string(Coef_parameters)[2:(end-1)], "=>", "=")
    println(f, "Parameters = list(Values = c($transformed_parameters), Coefs = c($units)),")
    if length(Forcings) > 0
      println(f, "Forcings = list(Values = list(")
      forcs = ""
      for (key,val) in Forcings
        forcs *= "$key = cbind(c($(string(val[1])[2:(end-1)])),c($(string(val[2])[2:(end-1)]))),"
      end
      forcs = forcs[1:(end-1)]
      println(f, forcs)
      units = replace(string(Coef_forcings)[2:(end-1)], "=>", "=")
      println(f, "), Coefs = c($units)),")
    end
    println(f, "Time = 0:1,")
    names_observed = string(Observed)[13:(end-1)]
    units = replace(string(Coef_observed)[2:(end-1)], "=>", "=")
    println(f, "Observed = list(names = c($names_observed), Coefs = c($units)),")
    println(f, """
    Settings = list(rtol = 1e-6,atol = 1e-10, maxsteps = 1e5, maxord = 5, hini = 0,
                      hmin = 0, hmax = 0, maxerr = 12, maxnonlin = 12,
                      maxconvfail = 12, method = "bdf", maxtime = 0, jacobian = 0),
    """)
    println(f, "model = $(name)_pointer)")
    close(f)
    nothing
end


function generate_code_Rcpp!(source::String; unit_analysis = false,name = "autogenerated_model", file = "autogenerated_model", jacobian = false, sensitivities = false)
  parsed_model = process_file(source)
  reaction_model = convert_master_equation(parsed_model)
  ode_model = convert_reaction_model(reaction_model)
  generate_code_Rcpp!(ode_model, unit_analysis = unit_analysis, name = name, file = file, jacobian = jacobian, sensitivities = sensitivities)
end
