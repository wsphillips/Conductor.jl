
using DelimitedFiles

function read_swc_data(filename)
    data = readdlm(filename, ' ', Any, '\n'; comments = true)
    size(data, 2) == 7 && return data
    if size(data,2) == 8 && all(isequal.(data[:,1],""))
        return data[:,2:8]
    else
        throw("Input file is not formatted correctly.")
    end
end

root(swc) = findfirst(isequal(-1), swc[:,7])


