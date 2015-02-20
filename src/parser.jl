
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
                        "export_meequation","export_mereaction",
                        "mereaction", "reaction", "equation", "state", "species",
                        "component", "parameter", "constant", "forcing", "constant_compartment",
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
