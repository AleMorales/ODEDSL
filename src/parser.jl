
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
const start_keywords = ["qsspecies", "qsstate", "species", "state", "reaction",
                       "variable", "derivative", "observed", "parameter", "constant",
                       "forcing", "component", "function", "compartment"]

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


## Parse species
# species <name> in <comparment> = <value> in <units>
# species <name> in <comparment> = <value>
# species <name> in <compartment>
# species <name>
function parse_species(s::ASCIIString)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  value = 0.0
  compartment = "__default_comp__"
  unit = ""
  comp_given = false
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
    # The third token could be the keyword "in" if there is a compartment
    elseif (c == 3) & (i == "in")
      comp_given = true
    # The fourth token could be the compartment (only if preceded by "in")
    elseif (c == 3) & (comp_given == true) & (i != "=")
      compartment = i
    # The 3rd or 5th token should be the = sign
    elseif (c == 3) & (i == "=")
      c = 4
    # The 4th or 6th token should be the value
    elseif (c == 4) & (i != "in")
      try value = parse(i) catch error("Error in species $name: You must supply a value") end
    # If there is an "in" keyword after the value, it indicates the units
    elseif (c == 4) & (i == "in")
      c = 5
    else
      unit *= i
    end
  end
  return name, dt.Species(value, compartment, dt.Unit(unit), {}, false)
end


## Parse qsspecies
# qsspecies <name> in <comparment> = <value> in <units>
# qsspecies <name> in <comparment> = <value>
# qsspecies <name> in <compartment>
# qsspecies <name>
function parse_qsspecies(s::ASCIIString)
  name, species = parse_species(s)
  species.QSS = true
  return name, species
end


## Parse parameter
# parameter <name> = <value> in <unit>
function parse_parameter(s::ASCIIString)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  value = 0.0
  unit = ""
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
    elseif (c == 4) & (i != "in")
      try value = parse(i) catch error("Error in parameter $name: You must supply a value") end
    elseif (c == 4) & (i == "in")
      c = 5
    elseif (c == 5)
      unit *= i
    end
  end
  !isa(value, Number) && error("Error in parameter $name: You must supply a value")
  return name, dt.Parameter(value, dt.Unit(unit))
end


## Parse constant
# constant <name> = <value> in <unit>
function parse_constant(s::ASCIIString)
  return parse_parameter(s)
end

## Parse state
# state <name> = <value> in <units>
# state <name> = <value>
# state <name>
function parse_state(s::ASCIIString)
  name, state_as_param =  parse_parameter(s)
  state = dt.State(state_as_param.Value, state_as_param.Units, false)
  return name, state
end

## Parse qsstate
# qsstate <name> = <value> in <units>
# qsstate <name> = <value>
# qsstate <name>
function parse_qsstate(s::ASCIIString)
  name, state_as_param = parse_parameter(s)
  state = dt.State(state_as_param.Value, state_as_param.Units, true)
  return name, state
end

## Parse reaction
# reaction <name> in <compartment> = <expression> for <stoichiometry> in <unit>
# reaction <name> in <compartment> = <expression> for <stoichiometry>
# reaction <name> = <expression> for <stoichiometry> in <unit>
# reaction <name> = <expression> for <stoichiometry>
function parse_reaction(s::ASCIIString)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  compartment = "__default_comp__"
  comp_given = false
  expr = ""
  unit = ""
  substrates = ""
  products = ""
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
    # The third token may be the in keyword
    elseif (c == 3) && (i == "in") && (comp_given == false)
      comp_given = true
    # The fourth token may be the compartment
    elseif (c == 3) && (comp_given == true) && (i != "=")
      compartment = i
    # The 3rd/5th token should be the = sign
    elseif c == 3
      i != "=" ?
          error("Error in reaction $name: I was expecting an equal = sign.") :
          (c = 4)
    # The 4th/6th and posterior tokens until "for" should be the expression
    elseif (c == 4) & (i != "for")
      expr *= i
    # After the mathematical expression the "for" is used to indicate the start of chemical stoichiometry
    elseif (c == 4) & (i == "for")
      c = 5
    # All tokens between "for" and "->" indicate substrates
    elseif (c == 5) & (i != "->")
      substrates *= i
    # The token "->" is used to separate substrates from products
    elseif (c == 5) && (i == "->")
      c = 6
    # All tokens after "->" represent products until the "in" keyword
    elseif (c == 6) & (i != "in")
      products *= i
    elseif (c == 6) & (i == "in")
      c = 7
    elseif (c == 7)
      unit *= i
    end
  end
  # Try to catch situations where the stoichiometry was not given or when the was a parsing error...
  (c == 4) &&
          error("Error in reaction $name: I am missing the stoichiometry of the reaction.")
  (c == 5) &&
          error("Error in reaction $name: I could not parse the stoichiometry correctly.")
  # Make sure that we have a valid mathematical expression. If not raise an error.
  # We don't store the parsed expression as it makes manipulation more difficult, but it is good to check
  # now if we can later parse it
  rhs = try parse(expr) catch
            error("Error in reaction $name: I don't understand the equation: $expr") end
  # Create the arrays of substrates and products (first make sure that all coefficients are provided)
  Substrates = components_in_stoich(substrates)
  Products = components_in_stoich(products)
  # Construct the Reaction object and return it
  return name, dt.Reaction(Substrates, Products, expr, compartment, dt.Unit(unit))
end

# Parse a combination of chemical species with or without pattern-matching
# n species[[component[form]...]]
# n species
function components_in_stoich(reactants::ASCIIString)
    output = dt.Reactant[]
    # Strip by + sign
    chunks = split(reactants, "+")
    (length(chunks[1]) == 0) && return output
    # Iterate over chunks
    for i in chunks
      # temporary strings to store stoichiometry and species description
      stoich_string = ""
      species_string = ""
      parse_stoich = true
      # Iterate over characters
        for j in i
          # If we are parsing stoich and we hit a letter, we are finished with stoich
          parse_stoich && ismatch(r"[a-zA-Z]", string(j)) && (parse_stoich = false)
          # Put the character in the right susbtring depending on whether we are parsing stoich or not
          parse_stoich ? stoich_string *= string(j) : species_string *= string(j)
        end
      push!(output, dt.Reactant(species_string, float(stoich_string)))
    end
    return output
end

## Parse variable
# variable <name> = <expression> in <unit>
function parse_variable(s::ASCIIString)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  expr = ""
  unit = ""
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
    # The fourth until the in keyword should be the mathematical expression
    elseif (c == 4) & (i != "in")
      expr *= i
    elseif (c == 4) & (i == "in")
      c = 5
    else
      unit *= i
    end
  end
  c == 3 && error("Error in variable $name: You have to specify an expression bound to the variable.")
  # Make sure that we have a valid mathematical expression. If not raise an error.
  # We don't store the parsed expression as it makes manipulation more difficult, but it is good to check
  # now if we can later parse it
  rhs = try parse(expr) catch
            error("Error in variable $name: I don't understand the equation...") end

  return name, dt.Variable(expr, dt.Unit(unit))
end

## Parse observed
# observed <name> = <expression> in <unit>
function parse_observed(s::ASCIIString)
  name, variable = parse_variable(s)
  return name, dt.Observed(variable.Expr, variable.Units)
end

## Parse derivative
# observed <name> of <state> = <expression> in <unit>
# observed <name> of <state> = <expression>
function parse_derivative(s::ASCIIString)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  state = ""
  expr = ""
  unit = ""
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
    elseif c == 3
      i != "of" ? error("Error in equation $name: I was expecting the keyword of to indicate the state the derivative is associated to.") : (c = 4)
# The fifth token should be the equal sign, otherwise we have an error
    elseif c == 4
      state = i
      c = 5
    elseif c == 5
      i != "=" ? error("Error in equation $name: I was expecting an equal = sign.") : (c = 6)
# The sixth until the in keyword should be the mathematical expression
    elseif (c == 6) & (i != "in")
      expr *= i
    elseif (c == 6) & (i == "in")
      c = 7
    else
      unit *= i
    end
  end
  c == 5 && error("Error in derivative $name: You need to supply a function.")

  rhs = try parse(expr) catch
            error("Error in derivative $name: I don't understand the equation...") end

  return name, dt.Derivative(expr, state, dt.Unit(unit))
end


## Parse forcing
# forcing <name> = <values> at <times> in <unit>
# forcing <name> = <values> at <times>
# forcing <name>
function parse_forcing(s::ASCIIString)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  values = ""
  times = ""
  unit = ""
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
    elseif (c == 5) & (i != "in")
      times *= i
    elseif (c == 5) & (i == "in")
      c = 6
    else
      unit *= i
    end
  end
  c == 4 && error("Error in forcing $name: I could not parse the times correctly.")
  c == 3 && (values = "[0.0, 0.0]"; times = "[0.0,0.1]")
  values_vec = eval(parse(values))
  times_vec = eval(parse(times))
  length(values_vec) != length(times_vec) && error("Error in forcing $name: Default values and times must have the same length")
  return name, dt.Forcing(times_vec, values_vec, dt.Unit(unit))
end

# Parse component
# component <name> <labels> of <species> = <values>
function parse_component(s::ASCIIString)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  values = ""
  species = ""
  forms = ""
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
    # The tokens between second and the keyword "of" should be an array of forms
    elseif c == 3
      forms *= i
      ismatch(r"\]" ,i) && (c = 4)
    # The token after the array of forms should be the keyword "of"
    elseif c == 4
      i != "of" ?
          error("Error in component $name: I was expecting the \"of\" keyword to specify the species.") :
          (c = 5)
    # After than we get the species
    elseif c == 5
      species = i
      c = 6
    # Then an = sign
    elseif c == 6
      i != "=" ?
          error("Error in component $name: I was expecting an equal = sign.") :
          (c = 7)
    # The remaining tokens are the values
    elseif (c == 7)
      values *= i
    end
  end
  c == 3 && error("Error in component $name: I could not parse it correctly. Did you keep spaces between the brackets and contents?")
  values_vec = eval(parse(values))
  forms_vec = eval(parse(forms))
  forms_vec = ASCIIString[string(i) for i in forms_vec]
  abs(sum(values_vec) - 1.0) > eps(1.0)  && error("Error in component $name: The sum of default initial values for the different forms must be equal to 1")
  return name, dt.Component(name, species, values_vec, forms_vec)
end

## Parse compartment
# compartment <name> is <constant/variable> = <value> in <unit>
function parse_compartment(s::ASCIIString)
  sentence = tokenize(strip(s))
  c = 1
  name = ""
  value = 0.0
  unit = ""
  ctype = ""
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
    elseif c == 3
      i != "is" ? error("Error in compartment $name: I was expecting they keyword 'is' followed by const/var") : (c = 4)
    elseif c == 4
      ctype = i
      c = 5
    # The fifth token should be the equal sign, otherwise we have an error
    elseif c == 5
      i != "=" ? error("Error in compartment $name: I was expecting an equal = sign.") : (c = 6)
    # The fourth (and posterior) token should be the expression
    elseif (c == 6) & (i != "in")
      try value = parse(i) catch error("Error in compartment $name: You must supply a value") end
    elseif (c == 6) & (i == "in")
      c = 7
    elseif (c == 7)
      unit *= i
    end
  end
  if ctype == "const"
    return name, dt.Parameter(value, dt.Unit(unit))
  elseif ctype == "var"
    return name, dt.State(value, dt.Unit(unit), false)
  else
    error("Error in compartment $name: You must specify the type of compartment with the keywords const or var")
  end
end

#### Create a dictionary that stores all the elements of the model description ####
function process_file(file::ASCIIString)
    source_file = read_file(file);
    process_text(source_file)
end

# Parse the model and store all the ouputs in a ReactionSource object
function process_text(source_file::ASCIIString)
    # OrderedDicts where the sentences will be stored according to typed structures
    sentences = {"constant" => OrderedDict{ASCIIString, dt.Parameter}(),
                 "parameter" => OrderedDict{ASCIIString, dt.Parameter}(),
                 "forcing" => OrderedDict{ASCIIString, dt.Forcing}(),
                 "state" => OrderedDict{ASCIIString, dt.State}(),
                 "qsstate" => OrderedDict{ASCIIString, dt.State}(),
                 "species" => OrderedDict{ASCIIString, dt.Species}(),
                 "qsspecies" => OrderedDict{ASCIIString, dt.Species}(),
                 "component" => OrderedDict{ASCIIString, dt.Component}(),
                 "variable" => OrderedDict{ASCIIString, dt.Variable}(),
                 "reaction" => OrderedDict{ASCIIString, dt.Reaction}(),
                 "compartment" => OrderedDict{ASCIIString, Union(dt.Parameter, dt.State)}(),
                 "derivative" => OrderedDict{ASCIIString, dt.Derivative}(),
                 "observed" => OrderedDict{ASCIIString, dt.Observed}()}

    # Parsers for every type of line
    parsers = {"constant" => parse_constant,
               "parameter" => parse_parameter,
               "forcing" => parse_forcing,
               "state" => parse_state,
               "qsstate" => parse_qsstate,
               "species" => parse_species,
               "qsspecies" => parse_qsspecies,
               "component" => parse_component,
               "variable" => parse_variable,
               "reaction" => parse_reaction,
               "compartment" => parse_compartment,
               "derivative" => parse_derivative,
               "observed" => parse_observed}

    # Read the file and parse it into logical lines
    source_file = strip_comments(source_file)
    source_file = tokenize(source_file)
    source_file = logical_lines(source_file)
    source_lines = store_logical_lines(source_file)

    # Process each type of input and store in the dictionaries
    for i in collect(keys(sentences))
      for j in source_lines[i]
        name, object = parsers[i](j)
        sentences[i][name] = object
      end
    end

    # Components are actually part of species.
    # Iterate over components and add them to the corresponding species
    for (key,val) in sentences["component"]
      push!(sentences["species"][val.Species].Components, val)
    end


    # Merge quasi-steady and dynamic variables
    merge!(sentences["species"], sentences["qsspecies"])
    merge!(sentences["state"], sentences["qsstate"])

    # Put compartments in the right compartments
    for (key,val) in sentences["compartment"]
      isa(val, dt.Parameter) ?
        sentences["parameter"][key] = val :
        sentences["state"][key] = val
    end

    # Check if we have any quasi steady state chemical species or state variable
    if (length(sentences["qsspecies"]) > 1) | (length(sentences["qsstate"]) > 1)
      DynamicType = "DAE"
    else
      DynamicType = "ODE"
    end

    # Create the ReactionSource object
    dt.ReactionSource(DynamicType,
                  sentences["constant"],
                  sentences["parameter"],
                  sentences["forcing"],
                  sentences["species"],
                  sentences["reaction"],
                  sentences["state"],
                  sentences["derivative"],
                  sentences["variable"],
                  sentences["observed"])

end
