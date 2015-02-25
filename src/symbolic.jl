module symbolics

using ODEDSL
using ODEDSL.DataTypes
using SymPy

# Create SymPy symbols for variabels in a model
# This is creating the names in the global namespace of the module
function generate_symbols!{T<:String}(symbols::Array{T,1})
  for i in symbols
    isdefined(parse(i)) || (@eval global $(parse(i)) = SymPy.Sym($i))
  end
end

# Remove input will remove any input to the model by taking the limits of the equations.
function remove_input(model::OdeSource, input::ASCIIString, limit_to::Number)
  # No side effects
  new_model = deepcopy(model)
  # Get a list of all variables on the rhs of the equations, as well as the dicts mapping each side to the other
  rhs_vars_array, rhs_to_lhs, lhs_to_rhs = ODEDSL.get_rhs(model.Equations)
  # Very ugly...populate the global namespace with symbolic variables..need a macro!!
  generate_symbols!(rhs_vars_array)
  # Iterate over equations and take the limit of each equation
  for (key,val) in model.Equations
    if in(input, lhs_to_rhs[key])
      new_expr = convert(Expr, limit(eval(val.Expr), parse(input), limit_to))
      if isa(new_expr, Number)
        new_model.Constants["k_"*key] = Parameter(eval(new_expr),val.Units)
        new_model.Equations[key] = Equation(parse("k_"*key), val.Exported, val.Units)
      else
        new_model.Equations[key] = Equation(new_expr, val.Exported, val.Units)
      end

    end
  end
  # This will automatically remove all those inputs that are no longer needed
  checked_new_model, rhs_to_lhs, lhs_to_rhs =  ODEDSL.check_lhs_rhs(new_model);
  return checked_new_model
end


end
