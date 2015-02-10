
function read_file(filename::String)
    file = open(filename)
    contents = readall(file)
    close(file)
    return contents
end
