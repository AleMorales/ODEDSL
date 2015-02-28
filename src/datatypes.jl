module DataTypes

using DataStructures
import Base.convert, Base.show, Base.showerror, Base.log, Base.log10, Base.exp, Base.isless, Base.ifelse, Base.sqrt

export Parameter, Forcing, Equation
export Species, Compartment, Reactant, Component
export MEReaction, Reaction
export MESource, ReactionSource, OdeSource, OdeSorted, OdeModel
export Dimension, Unit

###################################################
###  Dimensional analysis and unit conversion  ####
###################################################

immutable Dimension
  K::Rational{Int64}
  s::Rational{Int64}
  m::Rational{Int64}
  g::Rational{Int64}
  cd::Rational{Int64}
  mol::Rational{Int64}
  A::Rational{Int64}
end

immutable Unit
  d::Dimension
  f::Number
end

# Unit error (when we try to add or substract different units)
type UnitError <: Exception
end
showerror(io::IO, e::UnitError) = print(io, "There was an error with the units");


# Algebra with dimensions
*(u1::Dimension, u2::Dimension) =
  Dimension(u1.K + u2.K, u1.s + u2.s, u1.m + u2.m,
            u1.g + u2.g, u1.cd + u2.cd, u1.mol + u2.mol, u1.A + u2.A)
/(u1::Dimension, u2::Dimension) =
  Dimension(u1.K - u2.K, u1.s - u2.s, u1.m - u2.m,
            u1.g - u2.g, u1.cd - u2.cd, u1.mol - u2.mol, u1.A - u2.A)
+(u1::Dimension, u2::Dimension) = u1 == u2 ? u1 : throw(UnitError())
-(u1::Dimension, u2::Dimension) = u1 == u2 ? u1 : throw(UnitError())
^(u1::Dimension, s::Integer) =
  Dimension(u1.K*s, u1.s*s, u1.m*s, u1.g*s, u1.cd*s, u1.mol*s, u1.A*s)
^(u1::Dimension, s::Rational) =
  Dimension(u1.K*s, u1.s*s, u1.m*s, u1.g*s, u1.cd*s, u1.mol*s, u1.A*s)
^(u1::Dimension, s::Float64) =
  Dimension(u1.K*s, u1.s*s, u1.m*s, u1.g*s, u1.cd*s, u1.mol*s, u1.A*s)
*(u1::Dimension, s::Number) = u1
*(s::Number, u1::Dimension) = u1
/(u1::Dimension, s::Number) = u1
/(s::Number, u1::Dimension) = u1^(-s)
function exp(d::Dimension)
  d != none.d && (throw(UnitError()))
  none.d
end
log(d::Dimension) = exp(d::Dimension)
log10(d::Dimension) = exp(d::Dimension)
sqrt(d::Dimension) = d^(1//2)
isless(d1::Dimension, d2::Dimension) = d1 != d2 && (throw(UnitError()))
+(d::Dimension, s::Number) = d != none.d ? (throw(UnitError())) : d
+(s::Number, d::Dimension) = d + s
-(d::Dimension, s::Number) = d != none.d ? (throw(UnitError())) : d
-(s::Number, d::Dimension) = d - s
^(d1::Dimension, d2::Dimension) = d2 != none.d ? (throw(UnitError())) : d1
isless(d1::Dimension, s::Number) = d1 != none.d ? (throw(UnitError())) : d1
isless(s::Number, d1::Dimension) = isless(d1, s)
ifelse(d1::Dimension, d2::Dimension, s::Number) = d1 != d2 ? (throw(UnitError())) : d1
-(d::Dimension) = d
^(d1::Number, d2::Dimension) = d2 != none.d ? (throw(UnitError())) : none.d

# Algebra with units
*(u1::Unit, u2::Unit) = Unit(u1.d * u2.d, u1.f*u2.f)
/(u1::Unit, u2::Unit) = Unit(u1.d / u2.d, u1.f/u2.f)
+(u1::Unit, u2::Unit) = u1.d == u2.d ? u1 : throw(UnitError())
-(u1::Unit, u2::Unit) = u1.d == u2.d ? u1 : throw(UnitError())
^(u1::Unit, s::Integer) = Unit(u1.d^s, u1.f^s)
^(u1::Unit, s::Rational) = Unit(u1.d^s, u1.f^s)
*(u1::Unit, s::Number) = Unit(u1.d, u1.f*s)
*(s::Number, u1::Unit) = u1*s
/(u1::Unit, s::Number) = Unit(u1.d, u1.f/s)
/(s::Number, u1::Unit) = u1^(-s)

# Show pretty representation (still need to improve more...)
function show(io::IO, u::Unit)
  unit = u.d.mol != 0 ? ( u.d.mol - 1 < eps(1.0) ? "*mol" : ( u.d.mol < 0 ? "*1/mol^$(-float(u.d.mol))" : "*mol^$(float(u.d.mol))") ) : ""
  unit *= u.d.g != 0 ? ( u.d.g == 1 ? "*g" : ( u.d.g < 0 ? "*1/g^$(-float(u.d.g))" : "*g^$(float(u.d.g))") ) : ""
  unit *= u.d.m != 0 ? ( u.d.m == 1 ? "*m" : ( u.d.m < 0 ? "*1/m^$(-float(u.d.m))" : "*m^$(float(u.d.m))") ) : ""
  unit *= u.d.s != 0 ? ( u.d.s == 1 ? "*s" : ( u.d.s < 0 ? "*1/s^$(-float(u.d.s))" : "*s^$(float(u.d.s))") ) : ""
  unit *= u.d.K != 0 ? ( u.d.K == 1 ? "*K" : ( u.d.K < 0 ? "*1/K^$(-float(u.d.K))" : "*K^$(float(u.d.K))") ) : ""
  unit *= u.d.cd != 0 ? ( u.d.cd == 1 ? "*cd" : ( u.d.cd < 0 ? "*1/cd^$(-float(u.d.cd))" : "*cd^$(float(u.d.cd))") ) : ""
  unit *= u.d.A != 0 ? ( u.d.A == 1 ? "*A" : ( u.d.A < 0 ? "*1/A^$(-float(u.d.A))" : "*A^$(float(u.d.A))") ) : ""
  # Remove trailing "* "
  unit = unit[2:end]
  print(io, "$(float(u.f)) $unit")
end

function show(io::IO, d::Dimension)
  unit = d.mol != 0 ? ( d.mol == 1 ? "*mol" : ( d.mol < 0 ? "*1/mol^$(-float(d.mol))" : "*mol^$(float(d.mol))") ) : ""
  unit *= d.g != 0 ? ( d.g == 1 ? "*g" : ( d.g < 0 ? "*1/g^$(-float(d.g))" : "*g^$(float(d.g))") ) : ""
  unit *= d.m != 0 ? ( d.m == 1 ? "*m" : ( d.m < 0 ? "*1/m^$(-float(d.m))" : "*m^$(float(d.m))") ) : ""
  unit *= d.s != 0 ? ( d.s == 1 ? "*s" : ( d.s < 0 ? "*1/s^$(-float(d.s))" : "*s^$(float(d.s))") ) : ""
  unit *= d.K != 0 ? ( d.K == 1 ? "*K" : ( d.K < 0 ? "*1/K^$(-float(d.K))" : "*K^$(float(d.K))") ) : ""
  unit *= d.cd != 0 ? ( d.cd == 1 ? "*cd" : ( d.cd < 0 ? "*1/cd^$(-float(d.cd))" : "*cd^$(float(d.cd))") ) : ""
  unit *= d.A != 0 ? ( d.A == 1 ? "*A" : ( d.A < 0 ? "*1/A^$(-float(d.A))" : "*A^$(float(d.A))") ) : ""
  # Remove trailing "* "
  unit = unit[2:end]
  print(io, "$unit")
end


# Prefixes
const prefixes = {"" => 1//1,
                  "E" => int(1e18)//1, "P" => int(1e15)//1,
                  "T" => int(1e12)//1, "G" => int(1e9)//1, "M" => int(1e6), "k" => int(1e3),
                  "h" => int(1e2)//1, "da" => int(1e1)//1, "d" => 1//int(1e1), "c" => 1//int(1e2),
                  "m" => 1//int(1e3), "u" => 1//int(1e6), "n" => 1//int(1e9), "p" => 1//int(1e12),
                  "f" => 1//int(1e15), "a" => 1//int(1e18)}

# SI units
const SI_units = {"K" => Unit(Dimension(1//1,0,0,0,0,0,0),1//1),
                  "s" => Unit(Dimension(0,1//1,0,0,0,0,0),1//1),
                  "m" => Unit(Dimension(0, 0, 1//1, 0, 0, 0, 0), 1//1),
                  "g" => Unit(Dimension(0,0,0,1//1,0,0,0),1//1),
                  "cd" => Unit(Dimension(0,0,0,0,1//1,0,0),1//1),
                  "mol" => Unit(Dimension(0,0,0,0,0,1//1,0),1//1),
                  "A" => Unit(Dimension(0,0,0,0,0,0,1//1),1//1),
                  "none" => Unit(Dimension(0,0,0,0,0,0,0),1//1)}

# Automagically generate the SI units
for (pr,prval) = prefixes
    for (key,val) = SI_units
        unit_name = parse(pr*key)
        unit = Unit(val.d, val.f*prval)
        @eval const $unit_name = $unit
    end
end

# Derived units
const derived_units = {"M" => mol/dm^3,
                       "L" => dm^3,
                      "ha" => 100//1*m^2,
                      "Hz" => 1/s,
                      "N" => m*kg/s^2,
                      "Pa" => kg/m/s^2,
                      "J" => m^2*kg/s^2,
                      "W" => m^2*kg/s^3,
                      "C" => s*A,
                      "V" => m^2*kg/s^3/A,
                      "F" => s^4*A^2/kg/m^2,
                      "O" => m^2*kg/s^3/A^2,
                      "S" => s^3*A^2/kg/m^2}

# Automagically generate the derived units
for (pr,prval) = prefixes
    for (key,val) = derived_units
        unit_name = parse(pr*key)
        unit = Unit(val.d, val.f*prval)
        @eval const $unit_name = $unit
    end
end

# Create a unit object from a string that is parsed with these constants
Unit(s::ASCIIString) = s == "" ? none : eval(parse(s))
Dimension(u::Unit) = Dimension(u.d.K, u.d.s, u.d.m, u.d.g, u.d.cd, u.d.mol, u.d.A)

######################################
###  Some convenient conversions  ####
######################################
# Facilitate translating dicts into ordered dicts
function convert{T <: Any}(::Type{OrderedDict{String, T}}, x::Dict{String, T})
    out =  OrderedDict{String, T}()
    for (key,val) in x
        out[key] = val
    end
end

function convert{T <: Any}(::Type{OrderedDict{String, T}}, x::Nothing)
    return OrderedDict{String, T}()
end

##################################################
###  Datatypes for type description of model  ####
##################################################

# Basic datatypes
type Parameter
    Value::Float64
    Units::Unit
end
type Forcing
    Time::Array{Float64, 1}
    Value::Array{Float64, 1}
    Units::Unit
end
type Equation
    Expr::Union(Expr, Symbol, Number)
    Exported::Bool
    Units::Unit
end
type Species
    Value::Float64
    Compartment::String
    Units::Unit
end
type Compartment
    InputType::String # Parameter | State
    Value::Float64
    Units::Unit
end
type Reactant
    Name::String
    Stoichiometry::Float64
end
type Component
    Species::String
    Values::Array{Float64,1}
    Forms::Array{String,1}
end
type MEReaction
    Species::String
    Component::Array{String, 1}
    From::Array{String, 1}
    To::Array{String, 1}
    Expr::Union(Expr, Symbol)
    Units::Unit
end
type Reaction
    Substrates::Array{Reactant, 1}
    Products::Array{Reactant, 1}
    Expr::Union(Expr, Symbol)
    Compartment::String
    Exported::Bool
    Units::Unit
end

# Container for the structured description of the simulation language
abstract Source
type MESource <: Source
    Constants::Dict{String, Parameter}
    Parameters::OrderedDict{String, Parameter}
    Forcings::OrderedDict{String, Forcing}
    States::OrderedDict{String, Parameter}
    Species::Dict{String, Species}
    Components::Dict{String, Component}
    Compartments::Dict{String, Compartment}
    Equations::Dict{String, Equation}
    Reactions::Dict{String, Reaction}
    MEReactions::Dict{String, MEReaction}
    MEEquations::Dict{String, Equation}
end
type ReactionSource <: Source
    Constants::Dict{String, Parameter}
    Parameters::OrderedDict{String, Parameter}
    Forcings::OrderedDict{String, Forcing}
    States::OrderedDict{String, Parameter}
    Species::Dict{String, Species}
    Compartments::Dict{String, Compartment}
    Equations::Dict{String, Equation}
    Reactions::Dict{String, Reaction}
end
type OdeSource <: Source
    Constants::Dict{String, Parameter}
    Parameters::OrderedDict{String, Parameter}
    Forcings::OrderedDict{String, Forcing}
    States::OrderedDict{String, Parameter}
    Equations::Dict{String, Equation}
end
# Container where all rules have been transformed into an ODE system and the equations have been sorted
type OdeSorted
    Constants::Dict{String, Parameter}
    Parameters::OrderedDict{String, Parameter}
    Forcings::OrderedDict{String, Forcing}
    States::OrderedDict{String, Parameter}
    SortedEquations::Array{Dict{String, Equation}, 1}
end


##################################################
#######  Methods to print the datatypes  ########
##################################################

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


function show(io::IO, model::OdeSource)
  print(io,
"""
Constants:
----------
$(paste("\n", ["$(i[1]) = $(i[2])" for i in model.Constants]))
Parameters:
-----------
$(paste("\n", ["$(i[1]) = $(i[2])" for i in model.Parameters]))
Forcings:
---------
$(paste("\n", ["$(i[1]) = $(i[2])" for i in model.Forcings]))
States:
-------
$(paste("\n", ["$(i[1]) = $(i[2])" for i in model.States]))
Equations:
----------
$(paste("\n\n", ["$(i[1]) ($(i[2].Units.d)) = $(i[2])" for i in model.Equations]))
""")
end


function show(io::IO, model::OdeSorted)
  print(io,
"""
Constants:
----------
$(paste("\n", ["$(i[1]) = $(i[2])" for i in model.Constants]))
Parameters:
-----------
$(paste("\n", ["$(i[1]) = $(i[2])" for i in model.Parameters]))
Forcings:
---------
$(paste("\n", ["$(i[1]) = $(i[2])" for i in model.Forcings]))
States:
-------
$(paste("\n", ["$(i[1]) = $(i[2])" for i in model.States]))
Equations:
----------
$(prod([paste("\n\n", "$(j[2]) ($(j[2].Units))") for j in model.SortedEquations[1]]))
$(prod([paste("\n\n", ["$(i[1]) ($(i[2].Units.d)) = $(i[2])" for i in j]) for j in model.SortedEquations[2:end]]))
""")
end
show(io::IO, par::Parameter) = print(io,"$(par.Value*par.Units)")
show(io::IO, eq::Equation) = print(io,"$(string(eq.Expr))")
show(io::IO, f::Forcing) = print(io,"$(f.Value*f.Units.f) $(f.Units.d) at $(f.Time)")


end
