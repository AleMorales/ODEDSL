# Type that contains the data and functions associated to an ODE model in state-space form
type OdeModel
  # Initial values for the state variables
    States::NamedArray{Float64,1,Array{Float64,1},(Dict{ASCIIString,Int64},)}
  # Fixed parameters  
    Parameters::NamedArray{Float64,1,Array{Float64,1},(Dict{ASCIIString,Int64},)}
  # Dynamic inputs
    Inputs::NamedArray{Function,1,Array{Function,1},(Dict{ASCIIString,Int64},)}
  # Functions implementing time derivatives and observed variables
    Model::Function
    Observer::Function
end


# Example: NamedArray([1,2.0], (["A" => 1,"B" => 2],),("States",))

# Getters and Setters

#=# Pretty-printing of lists of names
function pretty!(names, element)
    if(length(names) == 1)
        return "There is no $element called \"$(names[1])\""
    else
        output = "There are no $(element)s called \"$(names[1])\""
        if(length(names) > 2)
            for i in names[2:(end - 1)]
                output *= ", \"$i\""
            end
        end
        output *= " or \"$(names[end])\""
        return output
    end
end

# Check if all names_required are contained in actual_names and, if not, raise an informative error
# pointing to the element of OdeModel where the mismatch occured.
function check_names!(names_required, actual_names, element)
  not_contained = filter(x -> !in(x, actual_names), names_required)
  if(length(not_contained) > 0)
        error(pretty!(not_contained, element))
  end
end
=#

