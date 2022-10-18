
enabled_tests = lowercase.(ARGS)

function addtests(fname)
    key = lowercase(splitext(fname)[1])
    if isempty(enabled_tests) || key in enabled_tests
        include(fname)
    end
end

addtests("hodgkinhuxley.jl")
addtests("prinzneuron.jl")
addtests("simplesynapse.jl")
