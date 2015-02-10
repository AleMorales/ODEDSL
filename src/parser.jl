
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

# Keywords that id the type of sentence
const start_keywords = ["mereaction", "reaction", "equation", "state", "species",
                        "component", "parameter", "constant", "forcing", "constant_compartment",
                        "meequation","export_equation","export_reaction","variable_compartment",
                        "export_meequation","export_mereaction"]

# Tokenize the string by whitespace
tokenize(source::String) = split(source, r"\s|\n|\t", keep = true)

# Separate into logical lines
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
