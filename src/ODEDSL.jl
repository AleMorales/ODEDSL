module ODEDSL

using Base.Test
using DataStructures
using Docile
using Lexicon
using Calculus
import JSON

import Calculus.differentiate, Base.min, Base.max, Base.deepcopy


# In order to enable docstring-style documentation of Julia code

@document

# Include the source files
include("utils.jl")
include("datatypes.jl")
import .DataTypes
dt = DataTypes
include("parser.jl")
include("conversion_reaction.jl")
include("conversion_ode.jl")
include("symbolic.jl")
include("checks.jl")
include("translate.jl")
include("translate_julia.jl")
include("translate_rcpp.jl")

end
