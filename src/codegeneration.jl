####################################################################################
####################################################################################
###########################       CHECK UNITS      #################################
####################################################################################
####################################################################################

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

####################################################################################
####################################################################################
########################       GENERATE JACOBIAN      ##############################
####################################################################################
####################################################################################


# Calculation Jacobian matrix of the model
function generate_jacobian_matrix(compressed_model::OdeSorted)
  names_states = collect(keys(compressed_model.States))
  names_derivatives = ["d_"*x*"_dt" for x in names_states]
  Jacobian = Array(Union(Expr, Symbol, Number),(length(names_states), length(names_states)))
  cd = 1
  for i in names_derivatives
    cs = 1
    for j in names_states
      #Jacobian[cd,cs] =  analytical_derivative(compressed_model.SortedEquations[2][i].Expr, parse(j), sympy_vars)
      Jacobian[cd,cs] =  differentiate(compressed_model.SortedEquations[2][i].Expr, parse(j))
      cs += 1
    end
    cd += 1
  end
  return Jacobian
end

####################################################################################
####################################################################################
###################       GENERATE EXTENDED SYSTEM      ############################
####################################################################################
####################################################################################

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
##########################   SOURCE CODE GENERATION ################################
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

# Calculate extended system
function generate_extended_system_julia(compressed_model::OdeSorted, name)
  sens_array = generate_sensitivity_array(compressed_model)
  # Create an extended ode system

  # Calculate the Jacobian of the extended system
  return "dummy_function", "dummy_jacobian"
end


####################################################################################
####################################################################################
###########################   FULL MODEL GENERATION ################################
####################################################################################
####################################################################################

"""
`generate_code_Julia!` will generate a Julia source file containing the code necessary
to load the model and simulate (using the SimulationModels packge). This function
will return nothing. As a side effect, it will write a file to the hard-drive.

Arguments
----------

`source`: Compulsory argument that contains a path to the ode file

`unit_analysis`: Keyword boolean argument to indicate whether we want to perform
analysis of physical dimensions and unit conversion

`name`: Keyword string argument with the name of the model we want to generate.
This name will be appended to the default names of the functions being generated
as well as the factory function from which we can generate new instances of the model.

`file`: Keyword string argument with the path to the file where we want to generate
the model. It does not have to coincide with the `name` argument given before. Note that
the extension should not be included.

`jacobian`: Keyword boolean argument to indicate whether the Jacobian of the model should
be generated or not.

`sensitivities`: Keyword boolean argument to indicate whether the sensitivity functions
associated to the model should be generated or not.
"""
function generate_code_Julia!(source::String; unit_analysis = false, name = "autogenerated_model", file = "autogenerated_model", jacobian = false, sensitivities = false)
  parsed_model = process_file(source)
  reaction_model = convert_master_equation(parsed_model)
  ode_model = convert_reaction_model(reaction_model)
  generate_code_Julia!(ode_model, unit_analysis = unit_analysis, name = name, file = file, jacobian = jacobian, sensitivities = sensitivities)
end

"""
`generate_code_Julia!` will generate a Julia source file containing the code necessary
to load the model and simulate (using the SimulationModels packge). This function
will return nothing. As a side effect, it will write a file to the hard-drive.

Arguments
----------

`source`: Compulsory argument that contains an OdeSource model

`unit_analysis`: Keyword boolean argument to indicate whether we want to perform
analysis of physical dimensions and unit conversion

`name`: Keyword string argument with the name of the model we want to generate.
This name will be appended to the default names of the functions being generated
as well as the factory function from which we can generate new instances of the model.

`file`: Keyword string argument with the path to the file where we want to generate
the model. It does not have to coincide with the `name` argument given before. Note that
the extension should not be included.

`jacobian`: Keyword boolean argument to indicate whether the Jacobian of the model should
be generated or not.

`sensitivities`: Keyword boolean argument to indicate whether the sensitivity functions
associated to the model should be generated or not.
"""

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
    jacobian && (jacobian_function = generate_jacobian_function_Julia(compressed_model, name))
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


"""
`write_model_Julia!` will write the source file required to run the models in julia.

It is generated according to the interface provided by SimulationModels which will take care
of wrapping the corresponding solvers, converting units, etc.
"""
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
