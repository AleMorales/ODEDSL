module ODEDSL

using DataStructures
using Docile
using Lexicon
using Calculus

import Calculus.differentiate, Base.min, Base.max


# In order to enable docstring-style documentation of Julia code

@document

# Include the source files
include("utils.jl")
include("datatypes.jl")
using .DataTypes
include("parser.jl")
include("symbolic.jl")
include("conversion.jl")
include("codegeneration.jl")
include("translate_Rcpp.jl")

end
