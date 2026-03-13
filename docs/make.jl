using Documenter
using PenguinExtendedStefan

makedocs(
    modules = [PenguinExtendedStefan],
    authors = "PenguinxCutCell contributors",
    sitename = "PenguinExtendedStefan.jl",
    format = Documenter.HTML(
        canonical = "https://PenguinxCutCell.github.io/PenguinExtendedStefan.jl",
        repolink = "https://github.com/PenguinxCutCell/PenguinExtendedStefan.jl",
        collapselevel = 2,
    ),
    pages = [
        "Overview" => "index.md",
        "Mathematical Model" => "model.md",
        "API" => "api.md",
        "Level-Set Workflow" => "levelset.md",
        "Validation" => "validation.md",
        "Examples" => "examples.md",
        "Diagnostics" => "diagnostics.md",
        "Limitations" => "limitations.md",
    ],
    pagesonly = true,
    warnonly = true,
    remotes = nothing,
)

if get(ENV, "CI", "") == "true"
    deploydocs(
        repo = "github.com/PenguinxCutCell/PenguinExtendedStefan.jl",
        push_preview = true,
    )
end
