function paste(sep, s...)
    output = ""
    for i in s
        if(isa(i,String))
            output *= i * sep
            elseif isa(i,Array{ASCIIString,1}) || isa(i,Array{UTF8String,1}) || isa(i,Array{Any,1})
            output *= paste(sep,i) * sep
        end
    end
    return output
end
function paste{T <: String}(sep, s::Array{T,1})
    output = ""
    for i in s
        output *= i * sep
    end
    return output
end

# These functions are reserved in order to have a basic mathematical lingo to write the models
const list_of_functions = ["abs", "abs2", "acos", "acosd", "acosh", "acot", "acotd",
                     "acoth", "acsc", "acscd", "acsch", "asec", "asecd", "asech",
                     "asin", "asind", "asinh", "atan", "atan2", "atand", "atanh",
                     "cbrt", "cos", "cosd", "cosh", "cot", "cotd", "coth", "csc",
                     "cscd", "csch", "exp", "exp2", "expm1", "log", "log", "log10",
                     "log1p", "log2", "sec", "secd", "sech", "sin", "sind", "sinh",
                     "sqrt", "tan", "tand", "tanh", "hypot","*","+","-","/","^","%","⋅",
                     "√","ifelse","min","max","<",">","=","<=",">="];