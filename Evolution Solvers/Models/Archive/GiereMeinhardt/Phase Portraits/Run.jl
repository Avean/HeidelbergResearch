include("phase_portrait.jl")

using REPL.TerminalMenus

project_dir = @__DIR__

model_dirs = sort(filter(name -> begin
        path = joinpath(project_dir, name)
        isdir(path) && isfile(joinpath(path, "model_data.jl"))
    end,
    readdir(project_dir)
))

if isempty(model_dirs)
    error("Nie znaleziono żadnych folderów modeli zawierających plik model_data.jl.")
end

function load_model_description(model_file::String)
    model_module = Module(Symbol("DescriptionModule_", abs(hash(model_file))))

    Base.include(model_module, model_file)

    if isdefined(model_module, :get_description)
        return Core.eval(model_module, :(get_description()))
    else
        return "Brak opisu modelu."
    end
end

function load_model(model_file::String)
    model_module = Module(Symbol("SelectedModel_", abs(hash(model_file))))

    Base.include(model_module, model_file)

    if !isdefined(model_module, :get_model)
        error("Plik $model_file musi definiować funkcję get_model().")
    end

    return Core.eval(model_module, :(get_model()))
end

function shorten_text(s::String, maxlen::Int)
    s = replace(strip(s), r"\s+" => " ")

    chars = collect(s)

    if maxlen <= 0
        return ""
    elseif length(chars) <= maxlen
        return s
    elseif maxlen == 1
        return "…"
    else
        return join(chars[1:maxlen-1]) * "…"
    end
end

function build_menu_options(model_dirs::Vector{String}, descriptions::Vector{String})
    term_width = try
        displaysize(stdout)[2]
    catch
        100
    end

    name_width = maximum(length.(model_dirs)) + 4
    name_width = max(name_width, 18)

    desc_width = max(20, term_width - name_width - 8)

    return [
        rpad(model_dirs[i], name_width) * shorten_text(descriptions[i], desc_width)
        for i in eachindex(model_dirs)
    ]
end

descriptions = String[]

for name in model_dirs
    model_file = joinpath(project_dir, name, "model_data.jl")
    push!(descriptions, load_model_description(model_file))
end

println()
println("Wybierz model strzałkami i zatwierdź Enter:")
println()

menu_options = build_menu_options(model_dirs, descriptions)

menu = RadioMenu(menu_options; pagesize = min(10, length(menu_options)))
choice = request("Model:", menu)

if choice == -1
    error("Wybór modelu został anulowany.")
end

selected_model_name = model_dirs[choice]
selected_model_dir = joinpath(project_dir, selected_model_name)
selected_model_file = joinpath(selected_model_dir, "model_data.jl")

M = load_model(selected_model_file)

fig = phase_portrait(
    M.f,
    M.g;
    u0 = M.u0,
    u1 = M.u1,
    v0 = M.v0,
    v1 = M.v1,
    Nu = M.Nu,
    Nv = M.Nv,
    arrow_length = M.arrow_length,
    h = M.h,
    equilibria = M.equilibria,
    title = M.title
)

output_path = joinpath(selected_model_dir, M.output_filename)
savefig(fig, output_path)

println()
println("Wybrany model: $selected_model_name")
println("Opis: $(descriptions[choice])")
println("Wykres zapisany w:")
println(output_path)