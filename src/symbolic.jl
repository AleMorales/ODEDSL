
# Add support for analytical derivatives of min and max using Heaviside function (taken from SymPy)
function differentiate(::SymbolParameter{:min}, args, wrt)
   k = 1
   derivative = "("*string(differentiate(args[k],wrt))*")"
   derivative *= "*heaviside(min("*replace(replace(string(args[[1:(k-1), (k+1):end]]), r"\{|\}", ""), ":","")*") - ($(args[k]))) + "
   for k in 2:length(args)
       derivative *= string(differentiate(args[k],wrt))
       derivative *= "*heaviside(min("*replace(replace(string(args[[1:(k-1), (k+1):end]]), r"\{|\}", ""), ":","")*") - ($(args[k]))) + "
   end
   parse(derivative[1:(end - 2)])
end

function differentiate(::SymbolParameter{:max}, args, wrt)
   k = 1
   derivative = "("*string(differentiate(args[k],wrt))*")"
   derivative *= "*heaviside(-max("*replace(replace(string(args[[1:(k-1), (k+1):end]]), r"\{|\}", ""), ":","")*") + ($(args[k]))) + "
   for k in 2:length(args)
       derivative *= string(differentiate(args[k],wrt))
       derivative *= "*heaviside(-max("*replace(replace(string(args[[1:(k-1), (k+1):end]]), r"\{|\}", ""), ":","")*") + ($(args[k]))) + "
   end
   parse(derivative[1:(end - 2)])
end

function differentiate(::SymbolParameter{:heaviside}, args, wrt)
   derivative = "("*string(differentiate(args[1],wrt))*")"
   derivative *= "*dirac("*replace(replace(string(args), r"\{|\}", ""), ":","")*")"
   parse(derivative)
end

differentiate(::SymbolParameter{:dirac}, args, wrt) = 0
