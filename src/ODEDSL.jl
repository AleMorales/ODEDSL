module ODEDSL

using Base.Test
using DataStructures
using Docile
using Lexicon
using Calculus

import Calculus.differentiate, Base.min, Base.max


# In order to enable docstring-style documentation of Julia code

@document

# Include the source files
include("utils.jl")
include("datatypes2.jl")
import .DataTypes
dt = DataTypes
include("parser2.jl")
include("conversion_reaction.jl")
include("conversion_ode.jl")
include("symbolic.jl")
include("checks2.jl")
include("translate_RCpp2.jl")

end
