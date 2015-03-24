
# Rcpp derivatives (RcppSundials)
function create_derivatives_rcpp(language, sorted_model, name)
  code = ""
  for level in 1:length(sorted_model.SortedEquations)
    for (lhs, rhs) in sorted_model.SortedEquations[level]
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
          mod_expr = sub_minmax(mod_expr)
          code *= "const double $lhs" * " = " * replace(replace(string(mod_expr), ".*", "*"), ":", "") * ";\n"
        end
    end
  end
  # Return line
  return_line = create_return_line_Rcpp(sorted_model.NamesDerivatives,sorted_model.NamesObserved)
  # Return the output
  code = replace(code, "1.0 *", "")
  up_boiler_plate =
"""
  array<vector<double>, 2> $(name)(const double& t, const vector<double>& states,
            const vector<double>& params, const vector<double>& forcs) { \n
"""
  low_boiler_plate =
"""
}
"""
  return paste("\n", up_boiler_plate, code, return_line,low_boiler_plate)
end


# Rcpp jacobian (RcppSundials)
function create_jacobian_rcpp(language, compressed_model, jacobian_matrix, names_derivatives, name)
  # Substitute * in AST to avoid pretty printing of scalar product
  for i = 1:size(jacobian_matrix, 1), j = 1:size(jacobian_matrix, 2)
      if isa(jacobian_matrix[i,j], Expr)
        jacobian_matrix[i,j] = sub_product(jacobian_matrix[i,j])
      end
  end
  code = ""
  for (lhs, rhs) in compressed_model.SortedEquations[1]
      if (ismatch(r"params|states|forcs", string(rhs.Expr)))
       number = match(r"(?<=\[)[\d]+(?=\])",string(rhs.Expr))
       new_expr = replace(string(rhs.Expr), "[$(number.match)]", "[$(int(number.match) - 1)]")
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
      code *= "output.at($(i-1),$(j - 1)) = $rhs;\n"
    end
  end

  up_boiler_plate =
"""
arma::mat $(name)_jacobian(const double& t, const vector<double>& states,
          const vector<double>& params, const vector<double>& forcs) {
  arma::mat output = arma::eye(states.size(), states.size());\n
"""
  low_boiler_plate =
"""
    return output;
  }
"""

  return paste("\n", up_boiler_plate, code,low_boiler_plate)

end

# Rcpp extended (RcppSundials)
function create_extended_rcpp(language, compressed_model, names_derivatives, name)

end

# Rcpp extended jacobian (RcppSundials)
function create_extended_jacobian_rcpp(language, extended_model, names_derivatives, name)

end

# Rcpp dae (RcppSundials)
function create_dae_rcpp(language, sorted_model, name)

end

# Write the file with all the model functions
function write_code_rcpp!(dynamic_type, model_function, jacobian_function, extended_function, extended_jacobian_function, name, file)

  up_boiler_plate =
"""
#include <array>
#include <vector>
#include <math.h>
using namespace std;

extern "C" {
"""
  low_boiler_plate =
"""
};
"""

  f = open("$(file).cpp","w")
  println(f, up_boiler_plate)
  println(f, model_function)
  println(f, jacobian_function)
  println(f, extended_function)
  println(f, extended_jacobian_function)
  println(f, low_boiler_plate)
  close(f)

  nothing

end


############################ AUXILLIARY FUNCTIONS ##############################

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
