
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
strip_comments(source::String) = replace(source, r"\#.*\n","")

# Keywords that id the type of sentence
const keywords1 = ["mereaction", "reaction", "equation", "state", "species",
                  "component", "parameter", "constant", "forcing", "compartment",
                  "meequation"]

# Substitute \n and \t
strip_tabs_newline(source::String) = replace(source, r"\n|\t", "")

# Tokenize the string by whitespace
tokenize(source::String) = split(source, r"\s+", keep = false)
