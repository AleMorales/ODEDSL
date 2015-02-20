module DataTypes

using DataStructures
import Base.convert

export convert, Parameter, Forcing, Equation
export Species, Compartment, Reactant, Component
export MEReaction, Reaction, Constraint, Unit
export MESource, ReactionSource, OdeSource, OdeSorted, OdeModel


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

# Basic datatypes
type Parameter
    Value::Float64
    Units::Union(Expr, Symbol)
end
type Forcing   
    Time::Array{Float64, 1}
    Value::Array{Float64, 1}
    Units::Union(Expr, Symbol)    
end
type Equation
    Expr::Union(Expr, Symbol, Number)
    Units::Union(Expr, Symbol) 
    Exported::Bool    
end 
type Species
    Value::Float64
    Units::Union(Expr, Symbol)
    Compartment::String
end
type Compartment
    InputType::String # Parameter | Constant | State | Forcing
    Units::Union(Expr, Symbol) 
    Value::Float64
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
    Units::Union(Expr, Symbol)
end
type Reaction
    Substrates::Array{Reactant, 1}
    Products::Array{Reactant, 1}
    Expr::Union(Expr, Symbol)
    Units::Union(Expr, Symbol)
    Compartment::String
    Exported::Bool
end
type Constraint
    Constraint::Union(Expr, Symbol, Number)
    Value::Union(Expr, Symbol, Number)
    Comparison::Symbol
end
type Unit 
    Expr::Union(Expr, Symbol, Number)
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
    Constraints::Dict{String, Constraint}
    MEConstraints::Dict{String, Constraint}
    Units::Dict{String, Unit}
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
    Constraints::Dict{String, Constraint}
    Units::Dict{String, Unit}    
end
type OdeSource <: Source
    Constants::Dict{String, Parameter}    
    Parameters::OrderedDict{String, Parameter}
    Forcings::OrderedDict{String, Forcing}
    States::OrderedDict{String, Parameter}
    Equations::Dict{String, Equation}
    Constraints::Dict{String, Constraint}
    Units::Dict{String, Unit}    
end
# Container where all rules have been transformed into an ODE system and the equations have been sorted
type OdeSorted
    Constants::Dict{String, Parameter}
    Parameters::OrderedDict{String, Parameter}
    Forcings::OrderedDict{String, Forcing}
    States::OrderedDict{String, Parameter}
    Constraints::Dict{String, Constraint}
    SortedEquations::Array{Dict{String, Equation}, 1}
    Units::Dict{String, Unit}    
end

# Container that contains the generated code in Julia. it 
type OdeModel
    States::OrderedDict{String, Any}
    Parameters::OrderedDict{String, Any}
    Forcings::OrderedDict{String, Any}
    Observed::Array{String, 1}
    Model::Function
    InputConstraints::Function
    OutputConstraints::Function
    Jacobian::Function
end

end