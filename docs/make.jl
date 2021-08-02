using Documenter, Conductor

makedocs(sitename="Conductor.jl",
         pages = ["Home" => "index.md"])

deploydocs(
    repo = "github.com/wsphillips/Conductor.jl.git",
    devbranch = "main")
