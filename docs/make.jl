using Documenter, Conductor

makedocs(sitename="Conductor.jl",
         doctest = false,
         modules = [Conductor],
         pages = ["Home" => "index.md"])

deploydocs(
    repo = "github.com/wsphillips/Conductor.jl.git",
    devbranch = "main",
    push_preview=true
)
