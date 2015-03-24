module DataTypes

using DataStructures
import Base.convert, Base.show, Base.showerror, Base.log, Base.log10,
       Base.exp, Base.isless, Base.ifelse, Base.sqrt, Base.abs,
       Base.==, Base.min, Base.max

export Parameter, Forcing, State, Variable, Observed
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

# +
+(d::Dimension, s::Number) = d != none.d ? (throw(UnitError())) : d
+(s::Number, d::Dimension) = d + s
+(u1::Dimension, u2::Dimension) = u1 == u2 ? u1 : throw(UnitError())
+(d::Dimension) = d

# -
-(d::Dimension, s::Number) = d != none.d ? (throw(UnitError())) : d
-(s::Number, d::Dimension) = d - s
-(u1::Dimension, u2::Dimension) = u1 == u2 ? u1 : throw(UnitError())
-(d::Dimension) = d

# *
*(u1::Dimension, u2::Dimension) =
  Dimension(u1.K + u2.K, u1.s + u2.s, u1.m + u2.m,
            u1.g + u2.g, u1.cd + u2.cd, u1.mol + u2.mol, u1.A + u2.A)
*(u1::Dimension, s::Number) = u1
*(s::Number, u1::Dimension) = u1

# /
/(u1::Dimension, u2::Dimension) =
  Dimension(u1.K - u2.K, u1.s - u2.s, u1.m - u2.m,
            u1.g - u2.g, u1.cd - u2.cd, u1.mol - u2.mol, u1.A - u2.A)
/(u1::Dimension, s::Number) = u1
/(s::Number, u1::Dimension) = u1^(-s)
exp(d::Dimension) = d != none.d ? (throw(UnitError())) : none.d

# ^
^(d1::Dimension, d2::Dimension) = d2 != none.d ? (throw(UnitError())) : d1
^(u1::Dimension, s::Integer) =
  Dimension(u1.K*s, u1.s*s, u1.m*s, u1.g*s, u1.cd*s, u1.mol*s, u1.A*s)
^(u1::Dimension, s::Rational) =
  Dimension(u1.K*s, u1.s*s, u1.m*s, u1.g*s, u1.cd*s, u1.mol*s, u1.A*s)
^(u1::Dimension, s::Float64) =
  Dimension(u1.K*s, u1.s*s, u1.m*s, u1.g*s, u1.cd*s, u1.mol*s, u1.A*s)
^(d1::Number, d2::Dimension) = d2 != none.d ? (throw(UnitError())) : none.d

# exp
exp(d::Dimension) = d != none.d ? (throw(UnitError())) : d

# log
log(d::Dimension) = d != none.d ? (throw(UnitError())) : d

# log10
log10(d::Dimension) = d != none.d ? (throw(UnitError())) : d

# sqrt
sqrt(d::Dimension) = d^(1//2)

# isless
isless(d1::Dimension, d2::Dimension) = d1 != d2 && (throw(UnitError()))
isless(d1::Dimension, s::Number) = d1 != none.d ? (throw(UnitError())) : d1
isless(s::Number, d1::Dimension) = isless(d1, s)

# ifelse comparison
ifelse(d1::Dimension, d2::Dimension, s::Number) = exp(d2::Dimension)
ifelse(d1::Dimension, s::Number, d2::Dimension) = d2 != none.d ? (throw(UnitError())) : none.d
ifelse(d1::Dimension, d2::Dimension, d3::Dimension) = d2 != d3 ? (throw(UnitError())) : d2
ifelse(d1::Dimension, s1::Number, s2::Number) = none.d

# abs
abs(d::Dimension) = d

# <
<(d1::Dimension, d2::Dimension) = d1 != d2 ? (throw(UnitError())) : none.d
<(d1::Dimension, n::Number) = d1 != none.d ? (throw(UnitError())) : none.d
<(n::Number, d1::Dimension) = d1 < n

# >
>(d1::Dimension, d2::Dimension) = d1 != d2 ? (throw(UnitError())) : none.d
>(d1::Dimension, n::Number) = d1 != none.d ? (throw(UnitError())) : none.d
>(n::Number, d1::Dimension) = d1 < n

# max function
max(d1::Dimension, d2::Dimension) = d1 != d2 ? (throw(UnitError())) : d1
max(d1::Dimension, n::Number) = d1 != none.d ? (throw(UnitError())) : none.d
max(n::Number, d1::Dimension) = max(d1,n)

# min function
min(d1::Dimension, d2::Dimension) = d1 != d2 ? (throw(UnitError())) : d1
min(d1::Dimension, n::Number) = d1 != none.d ? (throw(UnitError())) : none.d
min(n::Number, d1::Dimension) = min(d1,n)





# Basic algebra with units (for composite and derived units)
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
==(p1::Parameter, p2::Parameter) = p1.Value == p2.Value &&
                                   p1.Units == p2.Units

type Forcing
    Time::Array{Float64, 1}
    Value::Array{Float64, 1}
    Units::Unit
end
==(f1::Forcing, f2::Forcing) = f1.Time == f2.Time &&
                               f1.Value == f2.Value &&
                               f1.Units == f2.Units

type Variable
    Expr::ASCIIString
    Units::Unit
end
==(v1::Variable, v2::Variable) = v1.Expr == v2.Expr &&
                                 v1.Units == v2.Units

type Derivative
    Expr::ASCIIString
    State::ASCIIString
    Units::Unit
end
==(d1::Derivative, d2::Derivative) = d1.Expr == d2.Expr &&
                                     d1.State == d2.State &&
                                     d1.Units == d2.Units


type Observed
    Expr::ASCIIString
    Units::Unit
end
==(o1::Observed, o2::Observed) = o1.Expr == o2.Expr &&
                                 o1.Units == o2.Units

type State
    Value::Float64
    Units::Unit
    QSS::Bool
end
==(s1::State, s2::State) = s1.Value == s2.Value &&
                           s1.Units == s2.Units &&
                           s1.QSS == s2.QSS

type Reactant
    Name::ASCIIString
    Stoichiometry::Float64
end
==(r1::Reactant, r2::Reactant) = r1.Name == r2.Name &&
                                 r1.Stoichiometry == r2.Stoichiometry

type Component
    Name::ASCIIString
    Species::ASCIIString
    Values::Array{Float64,1}
    Forms::Array{ASCIIString,1}
end
==(c1::Component, c2::Component) = c1.Species == c2.Species &&
                                   c1.Values == c2.Values &&
                                   c1.Forms == c2.Forms


type Species
    Value::Float64
    Compartment::ASCIIString
    Units::Unit
    Components::Vector{Component}
    QSS::Bool
end
==(s1::Species, s2::Species) = s1.Value == s2.Value &&
                               s1.Compartment == s2.Compartment &&
                               s1.Units == s2.Units &&
                               s1.Components == s2.Components &&
                               s1.QSS == s2.QSS

type Reaction
    Substrates::Array{Reactant, 1}
    Products::Array{Reactant, 1}
    Expr::ASCIIString
    Compartment::ASCIIString
    Units::Unit
end
==(r1::Reaction, r2::Reaction) = r1.Substrates == r2.Substrates &&
                                 r1.Products == r2.Products &&
                                 r1.Expr == r2.Expr &&
                                 r1.Compartment == r2.Compartment &&
                                 r1.Units == r2.Units

type Equation
    Expr::Union(Symbol, Float64, Expr)
    Units::Unit
end

# Container for the structured description of the simulation language
type ReactionSource
    DynamicType::ASCIIString
    Constants::OrderedDict{ASCIIString, Parameter}
    Parameters::OrderedDict{ASCIIString, Parameter}
    Forcings::OrderedDict{ASCIIString, Forcing}
    Species::OrderedDict{ASCIIString, Species}
    Reactions::OrderedDict{ASCIIString, Reaction}
    States::OrderedDict{ASCIIString, State}
    Derivatives::OrderedDict{ASCIIString, Derivative}
    Variables::OrderedDict{ASCIIString, Variable}
    Observed::OrderedDict{ASCIIString, Observed}
end
type OdeSource
    DynamicType::ASCIIString
    Constants::OrderedDict{ASCIIString, Parameter}
    Parameters::OrderedDict{ASCIIString, Parameter}
    Forcings::OrderedDict{ASCIIString, Forcing}
    States::OrderedDict{ASCIIString, State}
    Derivatives::OrderedDict{ASCIIString, Derivative}
    Variables::OrderedDict{ASCIIString, Variable}
    Observed::OrderedDict{ASCIIString, Observed}
end
# Container where the different equations are parsed and stored in a single dict (we still conserve their meaning)
type OdeModel
    DynamicType::ASCIIString
    Constants::OrderedDict{ASCIIString, Parameter}
    Parameters::OrderedDict{ASCIIString, Parameter}
    Forcings::OrderedDict{ASCIIString, Forcing}
    States::OrderedDict{ASCIIString, State}
    Equations::OrderedDict{ASCIIString, Equation}
    NamesDerivatives::Vector{ASCIIString}
    NamesObserved::Vector{ASCIIString}
end
# Container where all rules have been transformed into an ODE system and the equations have been sorted
type OdeSorted
    DynamicType::ASCIIString
    Constants::OrderedDict{ASCIIString, Parameter}
    Parameters::OrderedDict{ASCIIString, Parameter}
    Forcings::OrderedDict{ASCIIString, Forcing}
    States::OrderedDict{ASCIIString, State}
    SortedEquations::Vector{OrderedDict{ASCIIString, Equation}}
    NamesDerivatives::Vector{ASCIIString}
    NamesObserved::Vector{ASCIIString}
end

function deepcopy(model::OdeSorted)

  new_eqns = OrderedDict{ASCIIString, dt.Equation}[]
  for i in 1:length(sorted_model.SortedEquations)
    push!(new_eqns, OrderedDict{ASCIIString, dt.Equation}())
    for (key,val) in sorted_model.SortedEquations[i]
      new_eqns[i][key] = dt.Equation(val.Expr, val.Units)
    end
  end

  OdeSorted(model.DynamicType,
            copy(model.Constants),
            copy(model.Parameters),
            copy(model.Forcings),
            copy(model.States),
            neq_eqns,
            copy(NamesDerivatives),
            copy(NamesObserved))

end


##################################################
#######  Methods to print the datatypes  ########
##################################################

function paste(sep, s...)
    output = ""
    for i in s
        if(isa(i,ASCIIString))
            output *= i * sep
            elseif isa(i,Array{ASCIIString,1})
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
Variables:
----------
$(paste("\n\n", ["$(i[1]) = $(parse(replace(i[2].Expr,"1.0", ""))) " for i in model.Variables]))
Observed:
----------
$(paste("\n\n", ["$(i[1]) = $(parse(replace(i[2].Expr,"1.0", ""))) " for i in model.Observed]))
Derivatives:
----------
$(paste("\n\n", ["$(i[1]) = $(parse(replace(i[2].Expr,"1.0", ""))) " for i in model.Derivatives]))
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
show(io::IO, eq::Variable) = print(io,"$(string(eq.Expr)) $(string(eq.Units))")
show(io::IO, eq::Derivative) = print(io,"$(string(eq.Expr))")
show(io::IO, eq::Observed) = print(io,"$(string(eq.Expr))")
show(io::IO, f::Forcing) = print(io,"$(f.Value*f.Units.f) $(f.Units.d) at $(f.Time)")


end
