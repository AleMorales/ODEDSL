function paste(sep, s...)
    output = ""
    for i in s
        if(isa(i,String))
            output *= i * sep
        elseif isa(i,Array{ASCIIString,1})
            output *= paste(sep,i) * sep
        end
    end
    return output
end
function paste(sep, s::Array{ASCIIString,1})
    output = ""
    for i in s
        output *= i * sep
    end
    return output
end


heaviside(x) = 0.5 * (sign(x) + 1)
dirac(x) = abs(x) <= eps(1.0)  ? 0 : Inf
min(n::Number) = n
max(n::Number) = n

# These functions are reserved in order to have a basic mathematical lingo to write the models
# It also includes operators, reserved words such as time and control flow functions
const list_of_functions = ["abs", "abs2", "acos", "acosd", "acosh", "acot", "acotd",
                     "acoth", "acsc", "acscd", "acsch", "asec", "asecd", "asech",
                     "asin", "asind", "asinh", "atan", "atan2", "atand", "atanh",
                     "cbrt", "cos", "cosd", "cosh", "cot", "cotd", "coth", "csc",
                     "cscd", "csch", "exp", "exp2", "expm1", "log", "log", "log10",
                     "log1p", "log2", "sec", "secd", "sech", "sin", "sind", "sinh",
                     "sqrt", "tan", "tand", "tanh", "hypot","*","+","-","/","^","%","⋅",
                     "√","ifelse","min","max","<",">","=","<=",">=","time"];
