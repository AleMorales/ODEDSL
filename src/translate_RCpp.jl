
####################################################################################
####################################################################################
############################   RC++ CODE GENERATION ###############################
####################################################################################
####################################################################################
# Generates code for SimulationModels in R

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

function sub_minmax(ex::Expr)
  new_ex = deepcopy(ex)
  for i in 1:length(new_ex.args)
    if new_ex.args[i] == :min || new_ex.args[i] == :max
        new_ex.args[i] == :min ?  (new_ex.args[i] = :Min) : (new_ex.args[i] = :Max)
    elseif isa(new_ex.args[i], Expr)
        new_ex.args[i] = sub_minmax(new_ex.args[i])
  end
end
return new_ex
end

sub_minmax(ex::Symbol) = ex


function create_derivatives_RCpp(model::dt.OdeSorted, name)
  code = ""
  names_forcings = collect(keys(model.Forcings))
  for level in 1:length(model.SortedEquations)
    for (lhs, rhs) in model.SortedEquations[level]
        if level == 1
          if (ismatch(r"params|states", string(rhs.Expr)))
           number = match(r"(?<=\[)[\d]+(?=\])",string(rhs.Expr))
           new_expr = replace(string(rhs.Expr), "[$(number.match)]", "[$(int(number.match) - 1)]")
           new_expr = replace(new_expr, "params[", "inp->p[")
           code *= "const double $(new_expr);\n"
          elseif (ismatch(r"forcs", string(rhs.Expr)))
            number = match(r"(?<=\[)[\d]+(?=\])",string(rhs.Expr))
            new_expr = replace(string(rhs.Expr), "forcs[$(number.match)]", "linterp(inp->f[$(int(number.match) - 1)], time)")
            code *= "const double $(new_expr);\n"
          else
           code *= "const double $(rhs.Expr);\n"
          end
        elseif in(lhs, model.NamesDerivatives)
          mod_expr = substitute_power(rhs.Expr)
          mod_expr = sub_minmax(mod_expr)
          c = findin(model.NamesDerivatives, [lhs]) - 1
          code *= "const double $lhs" * " = " * replace(replace(string(mod_expr), ".*", "*"), ":", "") * ";\n"
          code  *=  "ydot$c = $lhs;\n"
        else
          mod_expr = substitute_power(rhs.Expr)
          mod_expr = sub_minmax(mod_expr)
          code *= "const double $lhs" * " = " * replace(replace(string(mod_expr), ".*", "*"), ":", "") * ";\n"
        end
    end
  end
  # Return line
  return_line = "return 0;"
  # Pretty printing
  code = replace(code, "1.0 *", "")
  # Some boiler plate code. Includes function definitions
  up_boiler_plate =
"""
#include <SimulationModels.hpp>
#include <limits>
using namespace std;

double ifelse(bool condition, const double& result1, const double& result2) {
  if(condition) {
    return result1;
  } else {
    return result2;
  }
}

inline double heaviside(const double& arg) {
  return arg <= 0.0 ? 0.0 : 1.0;
}

inline double dirac(const double& arg) {
  return abs(arg) <= numeric_limits<double>::epsilon()  ? 0 : numeric_limits<double>::infinity();
}

inline double Min(const double& arg1) {
  return arg1;
}

inline double Max(const double& arg1) {
  return arg1;
}

inline double Min(const double& arg1, const double& arg2) {
  return arg1 <= arg2 ? arg1 : arg2;
}

inline double Max(const double& arg1, const double& arg2) {
  return arg1 >= arg2 ? arg1 : arg2;
}

extern "C" {
  int $(name)_derivatives(double time, N_Vector states_, N_Vector ydot_, void* inp_) {
  Inputs* inp = (Inputs*) inp_;
  auto states = NV_DATA_S(states_);
  auto ydot = NV_DATA_S(ydot_);
"""
  low_boiler_plate =
"""
\n
}
};
"""
  return  paste("\n",up_boiler_plate, code, return_line,low_boiler_plate)
end

function create_observed_RCpp(model::dt.OdeSorted, name)
  code = ""
  names_forcings = collect(keys(model.Forcings))
  for level in 1:length(model.SortedEquations)
    for (lhs, rhs) in model.SortedEquations[level]
        if level == 1
          if (ismatch(r"params|states", string(rhs.Expr)))
           number = match(r"(?<=\[)[\d]+(?=\])",string(rhs.Expr))
           new_expr = replace(string(rhs.Expr), "[$(number.match)]", "[$(int(number.match) - 1)]")
           new_expr = replace(new_expr, "params[", "inp.p[")
           code *= "const double $(new_expr);\n"
         elseif (ismatch(r"forcs", string(rhs.Expr)))
           number = match(r"(?<=\[)[\d]+(?=\])",string(rhs.Expr))
           new_expr = replace(string(rhs.Expr), "forcs[$(number.match)]", "linterp(inp.f[$(int(number.match) - 1)], time)")
           code *= "const double $(new_expr);\n"
         else
           code *= "const double $(rhs.Expr);\n"
          end
        elseif in(lhs, model.NamesObserved)
          mod_expr = substitute_power(rhs.Expr)
          mod_expr = sub_minmax(mod_expr)
          c = findin(model.NamesObserved, [lhs]) - 1
          code *= "const double $lhs" * " = " * replace(replace(string(mod_expr), ".*", "*"), ":", "") * ";\n"
          code  *=  "observed$c = $lhs;\n"
        else
          mod_expr = substitute_power(rhs.Expr)
          mod_expr = sub_minmax(mod_expr)
          code *= "const double $lhs" * " = " * replace(replace(string(mod_expr), ".*", "*"), ":", "") * ";\n"
        end
    end
  end
  # Return line
  return_line = "return 0;"
  # Return the output
  code = replace(code, "1.0 *", "")
  # We cannot have includes again
  up_boiler_plate =
"""
extern "C" {
  int $(name)_observer(const double& time, N_Vector states_, vector<double>& observed, const Inputs& inp) {
  auto states = NV_DATA_S(states_);
"""
  low_boiler_plate = "\n}};"

  return  paste("\n",up_boiler_plate, code, return_line,low_boiler_plate)
end

function generate_jacobian_function_RCpp(model::dt.OdeSorted, names_derivatives, name::ASCIIString)
  jacobian_matrix = generate_jacobian_matrix(model, names_derivatives)
  for i = 1:size(jacobian_matrix, 1), j = 1:size(jacobian_matrix, 2)
      if isa(jacobian_matrix[i,j], Expr)
        jacobian_matrix[i,j] = sub_product(jacobian_matrix[i,j])
      end
  end
  code = ""
  for (lhs, rhs) in model.SortedEquations[1]
      if (ismatch(r"params|states", string(rhs.Expr)))
       number = match(r"(?<=\[)[\d]+(?=\])",string(rhs.Expr))
       new_expr = replace(string(rhs.Expr), "[$(number.match)]", "[$(int(number.match) - 1)]")
       new_expr = replace(new_expr, "params[", "inp->p[")
       code *= "const double $(new_expr);\n"
     elseif (ismatch(r"forcs", string(rhs.Expr)))
       number = match(r"(?<=\[)[\d]+(?=\])",string(rhs.Expr))
       new_expr = replace(string(rhs.Expr), "forcs[$(number.match)]", "linterp(inp->f[$(int(number.match) - 1)], time)")
       code *= "const double $(new_expr);\n"
     else
       code *= "const double $(rhs.Expr);\n"
     end
  end
  for i = 1:size(jacobian_matrix, 1), j = 1:size(jacobian_matrix, 2)
    if isa(jacobian_matrix[i,j], Number) && eval(jacobian_matrix[i,j]) == 0
      continue
    else
      rhs = substitute_power(jacobian_matrix[i,j])
      rhs = sub_minmax(rhs)
      rhs = replace(replace(string(rhs), ".*", "*"), ":", "")
      code *= "DENSE_ELEM(Jac,$(i-1),$(j - 1)) = $rhs;\n"
    end
  end
  function_text = paste("\n",
  """
  extern "C" {
  int $(name)_jacobian(long int N, double time, N_Vector states_, N_Vector fy_, DlsMat Jac,
                    void *inputs_, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  Inputs* inp = (Inputs*) inputs_;
  auto states = NV_DATA_S(states_);\n
  """,
  code,
  "return 0;}};")
  return function_text
end

function generate_code_RCpp!(ode_model::dt.OdeModel; unit_analysis = true, name = "autogenerated_model", file = "autogenerated_model", sensitivities = false)

  # Extract the conversion coefficients for observed variables
  coef_observed = OrderedDict{ASCIIString, Float64}()
  for i in ode_model.NamesObserved
    coef_observed[i] = ode_model.Equations[i].Units.f
  end

  # Sort the equations
  sorted_model = OdeSorted(ode_model)

  # Create the default arguments with proper conversions and default values
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

  # Check the physical dimensions of the equations in the model (if demanded)
  unit_analysis &&  check_units(sorted_model)

  # Deep copy is not working with an array of ordered dicts with equations
  # So I need to do a manual copy of this
  original_equations = OrderedDict{ASCIIString, dt.Equation}[]
  for i in 1:length(sorted_model.SortedEquations)
    push!(original_equations, OrderedDict{ASCIIString, dt.Equation}())
    for (key,val) in sorted_model.SortedEquations[i]
      original_equations[i][key] = dt.Equation(val.Expr, val.Units)
    end
  end

  # Generate a version of the model without intermediate variables
  compressed_model = compress_model(sorted_model, level = 2)

  sorted_model.SortedEquations = original_equations

  # Generate the function to calculate the Jacobian
  jacobian_function = generate_jacobian_function_RCpp(compressed_model, sorted_model.NamesDerivatives, name)

  # Go through the equations and substitute * by .* (we need to avoid fancy scalar products in Julia)
  # This must always go after the Jacobian!
  for i in 1:length(sorted_model.SortedEquations)
      for (key,val) in sorted_model.SortedEquations[i]
         sorted_model.SortedEquations[i][key].Expr = sub_product(sorted_model.SortedEquations[i][key].Expr)
      end
  end

  # Generate the rhs function
  model_function = create_derivatives_RCpp(sorted_model, name)
  # Generate function with the observed values
  observed_function = create_observed_RCpp(sorted_model,name)


  # Write down the code!
  write_ode_model_RCpp!(named_states,coef_states, named_parameters, coef_parameters,
                 forcings, coef_forcings, ode_model.NamesObserved,coef_observed,
                  model_function, observed_function, jacobian_function, name, file)
  nothing
end


function write_ode_model_RCpp!(States::OrderedDict{ASCIIString, Float64},
                        Coef_states::OrderedDict{ASCIIString, Float64},
                        Parameters::OrderedDict{ASCIIString, Float64},
                        Coef_parameters::OrderedDict{ASCIIString, Float64},
                        Forcings::OrderedDict{ASCIIString, Any},
                        Coef_forcings::OrderedDict{ASCIIString, Float64},
                        Observed::Array{ASCIIString, 1},
                        Coef_observed::OrderedDict{ASCIIString, Float64},
                        Model::ASCIIString,
                        Observed_model::ASCIIString,
                        jacobian_function::ASCIIString,
                        name::ASCIIString,
                        file::ASCIIString)
    # Create the C++ file
    f = open("$(file).cpp","w")
    println(f, Model)
    println(f, Observed_model)
    println(f, jacobian_function)
    close(f)
    f = open("$(file).R","w")
    println(f, """
generate_$(name) = function() {
  ODEmodel\$new(

""")
    transformed_states = string(States)[2:(end-1)]
    transformed_states = replace(transformed_states, "=>", "=")
    units = replace(string(Coef_states)[2:(end-1)], "=>", "=")
    println(f, replace("States = list(Values = c($transformed_states), Coefs = c($units)),", "\"",""))
    transformed_parameters = string(Parameters)[2:(end-1)]
    transformed_parameters = replace(transformed_parameters, "=>", "=")
    units = replace(string(Coef_parameters)[2:(end-1)], "=>", "=")
    println(f, replace("Parameters = list(Values = c($transformed_parameters), Coefs = c($units)),", "\"",""))
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
    println(f, "Observed = list(Names = c($names_observed), Coefs = c($units)),")
    println(f, """
    Settings = list(rtol = 1e-6,
                    atol = 1e-6, maxsteps = 1e3, maxord = 5, hini = 1e-3,
                    hmin = 0, hmax = 100, maxerr = 5, maxnonlin = 10,
                    maxconvfail = 10, method = "bdf", maxtime = 0, jacobian = 1,
                    observer = 1, nder = $(length(Observed))),
    """)
    println(f,
    """
                model = getNativeSymbolInfo(name = "$(name)_derivatives",PACKAGE = "$name")\$address,
                observer = getNativeSymbolInfo(name = "$(name)_observer",PACKAGE = "$name")\$address,
                jacobian = getNativeSymbolInfo(name = "$(name)_jacobian",PACKAGE = "$name")\$address,
                extended_function = getNativeSymbolInfo(name = "$(name)_derivatives",PACKAGE = "$name")\$address,
                extended_jacobian = getNativeSymbolInfo(name = "$(name)_jacobian",PACKAGE = "$name")\$address)
    }
    """)
    close(f)
    nothing
end


function generate_code_RCpp!(source::ASCIIString; unit_analysis = false,name = "autogenerated_model", file = "autogenerated_model", jacobian = false, sensitivities = false)
  ode_model =  OdeModel(OdeSource(process_file(source)));
  generate_code_RCpp!(ode_model, unit_analysis = unit_analysis, name = name, file = file, sensitivities = sensitivities)
end
