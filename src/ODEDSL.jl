module ODEDSL

export generate_code_Julia!, generate_code_R!, generate_code_Rcpp!
# In order to enable docstring-style documentation of Julia code
using Docile
using Lexicon
@document

# Include the source files
include("utils.jl")
include("datatypes.jl")
include("parser.jl")
include("symbolic.jl")
include("conversion.jl")
include("codegeneration.jl")

end
