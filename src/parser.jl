
using .DataTypes
using DataStructures

# `read_file` will read the entire file as a single string. From this string we can
# start parsing the model description
function read_file(filename::String)
    file = open(filename)
    contents = readall(file)
    close(file)
    return contents
end

# Look a pattern consisting of a `#` followed by 0 or more character of any type except newline
# and ending in newline
strip_comments(source::String) = replace(source, r"\#.*\n"," ")

# Keywords that id the type of sentence.
# Note that the "export" version of a start keyword must go first in order for the
# store_logical_lines function to work
const start_keywords = ["export_equation","export_reaction","variable_compartment",
                        "export_meequation","export_mereaction","constant_compartment",
                        "mereaction", "reaction", "equation", "state", "species",
                        "component", "parameter", "constant", "forcing", 
                        "meequation"]

# Tokenize the string by whitespace, tab or newline
tokenize(source::String) = split(source, r"\s|\n|\t")

# Separate into logical lines knowing that a valid sentence must begin with a start keyword
# Whitespace is ignored. To make sure there are no problems, the reconstructed sentences are always
# appended surrounded by whitespace
function logical_lines(source::Array{SubString{ASCIIString},1})
  lines = ASCIIString[]
  c = 0
  for i in source
    if i == ""
      continue
    elseif in(i, start_keywords)
      c += 1
      push!(lines, i)
    else
      lines[c] *= " $i "
    end
  end
  return lines
end

# Store the lines in a dictionary where the key indicates the type of sentence
function store_logical_lines(source::Array{ASCIIString, 1})
  lines = Dict{ASCIIString, Array{ASCIIString, 1}}()
  [lines[i] = ASCIIString[] for i = start_keywords]
  for i in source
    for j in start_keywords
      if contains(i, j)
        push!(lines[j], i)
        break
      end
    end
  end
  return lines
end

#### Parsers for each type of sentence ####

## Parse equation
function parse_equation(s::ASCIIString, exported)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  expr = ""
  # Process non-blank tokens
  for i in sentence
    i == "" && continue
    # The first token should be the keyword
    if c == 1
      c = 2
    # The second token should be the name of the variable
    elseif c == 2
      name = i
      c = 3
    # The third token should be the equal sign, otherwise we have an error
    elseif c == 3
      i != "=" ? error("Error in equation $name: I was expecting an equal = sign.") : (c = 4)
    # The fourth (and posterior) token should be the expression
    elseif c == 4
      expr *= i
    end
  end
  # Make sure that we have a valid mathematical expression. If not, just raise an error
  eq = try Equation(parse(expr), exported) catch 
      error("Error in equation $name: I don't understand the equation.") end
  return name, eq
end

## Parse reaction
function parse_reaction(s::ASCIIString, exported)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  compartment = ""
  expr = ""
  substrates = ASCIIString[]
  stoich_substrate = Float64[]
  products = ASCIIString[]
  stoich_product = Float64[]
  # Process non-blank tokens
  for i in sentence
    i == "" && continue
    # The first token should be the keyword
    if c == 1
      c = 2
    # The second token should be the name of the variable
    elseif c == 2
      name = i
      c = 3
    # The third token should be the in keyword
    elseif c == 3
      i != "in" ? 
          error("Error in reaction $name: I was expecting the \"in\" keyword to specify the compartment.") : 
          (c = 4)
    # The fourth token should be the compartment
    elseif c == 4
      compartment = i
      c = 5
    # The 5th token should be the = sign
    elseif c == 5
      i != "=" ? 
          error("Error in reaction $name: I was expecting an equal = sign.") : 
          (c = 6)
    # The 6th and posterior tokens until "for" should be the expression
    elseif (c == 6) & (i != "for")
      expr *= i
    # After the mathematical expression the "for" is used to indicate the start of chemical stoichiometry
    elseif (c == 6) & (i == "for")
      c = 7
    # All tokens between "for" and "->" indicate substrates
    elseif (c == 7) & (i != "->") & (i != "+")
      try 
        isa(parse(i), Number) ? 
              push!(stoich_substrate, parse(i)) :
              push!(substrates, i)
      catch
        error("Error in reaction $name: I could not parse the substrates correctly.")
      end
    # The token "->" is used to separate substrates from products
    elseif (c == 7) & (i == "->")
      c = 8
    # All tokens after "->" represent products
    elseif c == 8 & (i != "+")
      try
        isa(parse(i), Number) ? 
              push!(stoich_product, parse(i)) :
              push!(products, i)
      catch
        error("Error in reaction $name: I could not parse the products correctly.")
      end
    end
  end
  # Try to catch situations where the stoichiometry was not given or when the was a parsing error...
  (c == 6) && 
          error("Error in reaction $name: I am missing the stoichiometry of the reaction.")
  (c == 7) && 
          error("Error in reaction $name: I could not parse the stoichiometry correctly.")
  # Make sure that we have a valid mathematical expression. If not raise an error
  rhs = try parse(expr) catch 
            error("Error in reaction $name: I don't understand the equation...") end
  # Create the arrays of substrates and products (first make sure that all coefficients are provided)
  length(stoich_substrate) != length(substrates) &&
      error("Error in reaction $name: I could not parse the substrates correctly.")
  length(stoich_product) != length(products) &&
      error("Error in reaction $name: I could not parse the products correctly.")
  Substrates = [Reactant(substrates[i], stoich_substrate[i]) for i = 1:length(substrates)] 
  Products = [Reactant(products[i], stoich_product[i]) for i = 1:length(products)]
  # Construct the Reaction object and return it
  return name, Reaction(Substrates, Products, rhs, compartment, exported)
end


## Parse parameter
function parse_parameter(s::ASCIIString)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  value = 0.0
  # Process non-blank tokens
  for i in sentence
    i == "" && continue
    # The first token should be the keyword
    if c == 1
      c = 2
    # The second token should be the name of the parameter
    elseif c == 2
      name = i
      c = 3
    # The third token should be the equal sign, otherwise we have an error
    elseif c == 3
      i != "=" ? error("Error in parameter $name: I was expecting an equal = sign.") : (c = 4)
    # The fourth (and posterior) token should be the expression
    elseif c == 4
      try value = parse(i) catch error("Error in parameter $name: You must supply a value") end
    end
  end
  return name, Parameter(value)
end

## Parse variable_compartment
function parse_variable_compartment(s::ASCIIString)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  value = 0.0
  # Process non-blank tokens
  for i in sentence
    i == "" && continue
    # The first token should be the keyword
    if c == 1
      c = 2
    # The second token should be the name of the parameter
    elseif c == 2
      name = i
      c = 3
    # The third token should be the equal sign, otherwise we have an error
    elseif c == 3
      i != "=" ? error("Error in compartment $name: I was expecting an equal = sign.") : (c = 4)
    # The fourth (and posterior) token should be the expression
    elseif c == 4
      try value = parse(i) catch error("Error in compartment $name: You must supply a value") end
    end
  end
  return name, Compartment("State", value)
end

## Parse constant_compartment
function parse_constant_compartment(s::ASCIIString)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  value = 0.0
  # Process non-blank tokens
  for i in sentence
    i == "" && continue
    # The first token should be the keyword
    if c == 1
      c = 2
    # The second token should be the name of the parameter
    elseif c == 2
      name = i
      c = 3
    # The third token should be the equal sign, otherwise we have an error
    elseif c == 3
      i != "=" ? error("Error in compartment $name: I was expecting an equal = sign.") : (c = 4)
    # The fourth (and posterior) token should be the expression
    elseif c == 4
      try value = parse(i) catch error("Error in compartment $name: You must supply a value") end
    end
  end
  return name, Compartment("Constant", value)
end

## Parse forcing
function parse_forcing(s::ASCIIString)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  values = ""
  times = ""
  # Process non-blank tokens
  for i in sentence
    i == "" && continue
    # The first token should be the keyword
    if c == 1
      c = 2
    # The second token should be the name of the parameter
    elseif c == 2
      name = i
      c = 3
    # The third token should be the equal sign, otherwise we have an error
    elseif c == 3
      i != "=" ? error("Error in forcing $name: I was expecting an equal = sign.") : (c = 4)
    # The tokens between the = sign and "at" keyword should be the values
    elseif (c == 4) & (i != "at")
      values *= i
    elseif (c == 4) & (i == "at")
      c = 5
    elseif c == 5
      times *= i
    end
  end
  c == 4 && error("Error in forcing $name: I could not parse the times correctly.") 
  values_vec = eval(parse(values))
  times_vec = eval(parse(times))
  return name, Forcing(times_vec, values_vec)
end

## Parse species
function parse_species(s::ASCIIString)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  value = ""
  compartment = ""
  # Process non-blank tokens
  for i in sentence
    i == "" && continue
    # The first token should be the keyword
    if c == 1
      c = 2
    # The second token should be the name of the parameter
    elseif c == 2
      name = i
      c = 3
    # The third token should be the in keyword
    elseif c == 3
      i != "in" ? 
          error("Error in species $name: I was expecting the \"in\" keyword to specify the compartment.") : 
          (c = 4)
    # The fourth token should be the compartment
    elseif c == 4
      compartment = i
      c = 5
    # The 5th token should be the = sign
    elseif c == 5
      i != "=" ? 
          error("Error in species $name: I was expecting an equal = sign.") : 
          (c = 6)
    # The sixth token should be the value
    elseif c == 6
      try value = parse(i) catch error("Error in species $name: You must supply a value") end
    end
  end
  return name, Species(value, compartment)
end


#### Create a dictionary that stores all the elements of the model description ####

# Parse the model and store all the ouputs in a MESource object
function process_file(file::ASCIIString)
    # Initiliaze
    docs = {"constant" => Dict{String, Parameter}(),
            "parameter" => OrderedDict{String, Parameter}(),
            "forcing" => OrderedDict{String, Forcing}(),
            "state" => OrderedDict{String, Parameter}(),
            "species" => Dict{String, Species}(),
            "component" => Dict{String, Component}(),
            "equation" => Dict{String, Equation}(),
            "reaction" => Dict{String, Reaction}(),
            "mereaction" => Dict{String, MEReaction}(),
            "meequation" => Dict{String, Equation}(),
            "constant_compartment" => Dict{String, Compartment}(),
            "export_meequation" => Dict{String, Equation}(),
            "export_mereaction" => Dict{String, MEReaction}(),
            "export_reaction" => Dict{String, Reaction}(),
            "export_equation" => Dict{String, Equation}(),
            "variable_compartment" => Dict{String, Compartment}()} 

    # Read the file and parse it into logical lines
    source_file = read_file(file);
    source_file = strip_comments(source_file)
    source_file = tokenize(source_file)
    source_file = logical_lines(source_file)
    source_lines = store_logical_lines(source_file)

    # Process each type of input and store in the dictionaries
    for i = source_lines["export_equation"]
      name, object = parse_equation(i, true)
      docs["export_equation"][name] = object
    end
    for i = source_lines["export_reaction"]
      name, object = parse_reaction(i, true)
      docs["export_reaction"][name] = object
    end
    for i = source_lines["equation"]
      name, object = parse_equation(i, false)
      docs["equation"][name] = object
    end
    for i = source_lines["reaction"]
      name, object = parse_reaction(i, false)
      docs["reaction"][name] = object
    end
    for i = source_lines["variable_compartment"]
      name, object = parse_variable_compartment(i)
      docs["variable_compartment"][name] = object
    end
    for i = source_lines["constant_compartment"]
      name, object = parse_constant_compartment(i)
      docs["constant_compartment"][name] = object
    end
    for i = source_lines["constant"]
      name, object = parse_parameter(i)
      docs["constant"][name] = object
    end
    for i = source_lines["state"]
      name, object = parse_parameter(i)
      docs["state"][name] = object
    end
    for i = source_lines["forcing"]
      name, object = parse_forcing(i)
      docs["forcing"][name] = object
    end
    for i = source_lines["species"]
      name, object = parse_species(i)
      docs["species"][name] = object
    end

    # Merge exported and non-exported variables
    merge!(docs["equation"], docs["export_equation"])
    merge!(docs["reaction"], docs["export_reaction"])
    merge!(docs["mereaction"], docs["export_mereaction"])
    merge!(docs["meequation"], docs["export_meequation"])
    # Merge the different types of compartments
    docs["compartment"] = merge(docs["constant_compartment"], docs["variable_compartment"])

  return MESource(docs["constant"], docs["parameter"], docs["forcing"], docs["state"],
                   docs["species"], docs["component"], docs["compartment"], docs["equation"],
                   docs["reaction"], docs["mereaction"], docs["meequation"])        

end

# Still need to write the parsers for:
#"export_meequation","export_mereaction",
#"mereaction", #"component", "meequation"
