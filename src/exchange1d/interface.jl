"""
    exchange1d(filenames::Vector{String})

Interactive text-based interface for 1D chemical exchange analysis.

Guides the user through:
1. Model selection
2. Molecule mapping (if required by model)
3. Peak integration parameters
4. Parameter review and editing
5. Fitting and display of results
6. Saving results

# Arguments
- `filenames`: Vector of paths to NMR experiment directories
"""
function exchange1d(filenames::Vector{String})
    @info "Exchange 1D analysis"

    # ── 1. Model selection ──────────────────────────────────────────────
    model = _prompt_model()
    model === nothing && return nothing

    # ── 2. Build problem ────────────────────────────────────────────────
    prob = ExchangeProblem(filenames, model)
    @info "Loaded $(length(prob.experiments)) experiments:"
    for expt in prob.experiments
        @info "  $(short_expt_path(expt)) ($(typeof(expt).name.name), $(expt.field_teslas) T)"
    end

    # ── 3. Molecule mapping (if model requires it) ──────────────────────
    if nmolecules(model) > 1
        ok = _prompt_moleculemap!(model, prob)
        ok || return nothing
    end

    # ── 4. Integration ──────────────────────────────────────────────────
    _prompt_integration!(prob)

    # ── 5. Parameters + fit loop ────────────────────────────────────────
    p0 = default_params(prob)
    while true
        p0 = _prompt_params(p0)
        p0 === nothing && return nothing

        @info "Fitting..."
        result = fit(prob, p0)
        _print_result(result)
        display(plot(result.plots...; size=(1200, 800)))

        action = _prompt_after_fit()
        if action == :save
            _save_results(result, prob)
            return result
        elseif action == :adjust
            continue  # loop back to edit p0 (the initial parameters)
        else  # :quit
            return result
        end
    end
end

function exchange1d(directory::String="")
    filenames = select_expts(directory)
    isempty(filenames) && return nothing
    return exchange1d(filenames)
end

# ═══════════════════════════════════════════════════════════════════════════
# Helper: short path for display
# ═══════════════════════════════════════════════════════════════════════════

function short_expt_path(expt::AbstractExperiment)
    # spec[:filename] returns processed data path e.g. .../sophia_trypsin/34/pdata/1
    # strip pdata/N suffix to get the experiment directory, then show parent/number
    expt_dir = dirname(dirname(expt.spec[:filename]))
    parent_folder = basename(dirname(expt_dir))
    folder_num = basename(expt_dir)
    return joinpath(parent_folder, folder_num)
end

# ═══════════════════════════════════════════════════════════════════════════
# Step 1: Model selection
# ═══════════════════════════════════════════════════════════════════════════

function _available_models()
    models = Type[]
    _collect_concrete_subtypes!(AbstractModel, models)
    return models
end

function _collect_concrete_subtypes!(T, out)
    for S in subtypes(T)
        if isabstracttype(S)
            _collect_concrete_subtypes!(S, out)
        else
            push!(out, S)
        end
    end
end

function _prompt_model()
    model_types = _available_models()
    instances = [T() for T in model_types]
    names = [modelname(m) for m in instances]
    push!(names, "Cancel")

    menu = RadioMenu(names)
    choice = request("Select exchange model:", menu)
    (choice == -1 || choice == length(names)) && return nothing

    model = instances[choice]
    @info "Selected model: $(modelname(model))"
    return model
end

# ═══════════════════════════════════════════════════════════════════════════
# Step 3: Molecule mapping
# ═══════════════════════════════════════════════════════════════════════════

function _prompt_moleculemap!(model, prob::ExchangeProblem)
    roles = molecules(model)

    # collect all unique sample molecule names across experiments
    all_names = String[]
    for expt in prob.experiments
        for name in keys(sampleconcentrations(expt))
            name in all_names || push!(all_names, name)
        end
    end
    sort!(all_names)

    @info "Sample molecules: $(join(all_names, ", "))"

    for (role, description) in roles
        options = copy(all_names)
        push!(options, "Cancel")
        menu = RadioMenu(options)
        choice = request("Assign :$role ($description) to:", menu)
        (choice == -1 || choice == length(options)) && return false
        model.moleculemap[role] = all_names[choice]
        @info "  :$role → $(all_names[choice])"
    end
    return true
end

# ═══════════════════════════════════════════════════════════════════════════
# Step 4: Integration
# ═══════════════════════════════════════════════════════════════════════════

function _prompt_integration!(prob::ExchangeProblem)
    first_expt = prob.experiments[1]
    spec = first_expt.spec
    ppm_axis = dims(spec, F1Dim)
    ppm_min = round(minimum(ppm_axis); digits=2)
    ppm_max = round(maximum(ppm_axis); digits=2)
    default_peak = round(spec[1, :offsetppm]; digits=2)

    @info "Spectral range: $ppm_min to $ppm_max ppm"

    peakppm = _prompt_value("Peak position in ppm", default_peak) do v
        ppm_min ≤ v ≤ ppm_max || error("Peak position $v ppm is outside spectral range ($ppm_min to $ppm_max)")
    end

    ppmwidth = _prompt_value("Integration width in ppm", 0.1) do v
        v > 0 || error("Width must be positive")
    end

    noiseppm = _prompt_value("Noise position in ppm", nothing) do v
        lo, hi = v - ppmwidth / 2, v + ppmwidth / 2
        ppm_min ≤ lo && hi ≤ ppm_max ||
            error("Noise region $lo..$hi ppm is outside spectral range ($ppm_min to $ppm_max)")
    end

    @info "Integrating: peak=$peakppm ppm, noise=$noiseppm ppm, width=$ppmwidth ppm"
    integrate!(prob, peakppm, noiseppm, ppmwidth)
end

"""Prompt for a Float64 value with optional default and validation. Retries on error."""
function _prompt_value(validate::Function, prompt::String, default)
    while true
        if default === nothing
            print("$prompt: ")
        else
            print("$prompt [$default]: ")
        end
        input = strip(readline())
        value = if isempty(input) && default !== nothing
            default
        else
            tryparse(Float64, input)
        end
        if value === nothing
            @warn "Invalid number"
            continue
        end
        try
            validate(value)
            return value
        catch e
            @warn e.msg
        end
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# Step 5: Parameter editing
# ═══════════════════════════════════════════════════════════════════════════

"""
    _prompt_params(p0) -> ComponentArray or nothing

Interactive parameter editor. Returns edited parameters, or nothing if cancelled.
"""
function _prompt_params(p0::ComponentArray)
    while true
        items = _flatten_params_items(p0)
        maxlen = maximum(item -> length(item.label), items)
        menu_items = [rpad(item.label, maxlen + 2) * "= " * _format_value(item.value) for item in items]
        push!(menu_items, "▶ Continue to fit")
        push!(menu_items, "✕ Cancel")

        menu = RadioMenu(menu_items)
        choice = request("Review parameters (select to edit, or continue):", menu)

        # cancel
        (choice == -1 || choice == length(menu_items)) && return nothing
        # continue
        choice == length(menu_items) - 1 && break

        item = items[choice]
        print("  New value for $(item.label) [$(item.value)]: ")
        input = strip(readline())
        if !isempty(input)
            newval = tryparse(Float64, input)
            if newval !== nothing
                _set_param!(p0, item.flat_index, newval)
            else
                @warn "Could not parse value: $input"
            end
        end
    end
    return p0
end

"""
A flattened parameter item: individual scalar elements, with array elements expanded.
- `label`: display name (e.g. "spin.delta[1]")
- `flat_index`: index into the underlying flat data array of the ComponentArray
- `value`: the scalar value
- `section`: top-level section name (e.g. "model", "spin", "nuisance")
"""
struct _ParamItem
    label::String
    flat_index::Int
    value::Any
    section::String
end

"""Flatten a ComponentArray into individual scalar items, expanding arrays."""
function _flatten_params_items(ca::ComponentArray, prefix="", section="";
                               flat_offset::Ref{Int}=Ref(0))
    items = _ParamItem[]
    for key in keys(ca)
        val = ca[key]
        full = isempty(prefix) ? string(key) : prefix * "." * string(key)
        sec = isempty(section) ? string(key) : section
        if val isa ComponentArray
            append!(items, _flatten_params_items(val, full, sec; flat_offset))
        elseif val isa AbstractVector
            for i in eachindex(val)
                flat_offset[] += 1
                push!(items, _ParamItem(full * "[$i]", flat_offset[], val[i], sec))
            end
        else
            flat_offset[] += 1
            push!(items, _ParamItem(full, flat_offset[], val, sec))
        end
    end
    return items
end

"""Set a single scalar parameter by flat index into the underlying data."""
function _set_param!(ca::ComponentArray, flat_index::Int, value)
    getdata(ca)[flat_index] = value
    return ca
end

# ═══════════════════════════════════════════════════════════════════════════
# Step 6: Display results + post-fit prompt
# ═══════════════════════════════════════════════════════════════════════════

function _print_result(result)
    items = _flatten_params_items(result.params)
    maxlen = maximum(item -> length(item.label), items)
    w = max(maxlen + 4, 16)

    println()
    println("═" ^ 60)
    println("  Fit results")
    println("═" ^ 60)

    prev_section = ""
    for item in items
        if item.section != prev_section
            prev_section != "" && println("  ", "─" ^ 56)
            prev_section = item.section
        end
        println("  ", rpad(item.label, w), _format_value(item.value))
    end

    println("  ", "─" ^ 56)
    println("  ", rpad("chi2", w), round(result.chi2; digits=2))
    println("  ", rpad("reduced chi2", w), round(result.reduced_chi2; digits=4))
    println("  ", rpad("nobs", w), result.nobs)
    println("  ", rpad("nparams", w), result.nparams)
    println("  ", rpad("dof", w), result.dof)
    println("═" ^ 60)
    println()
end

"""Format a value for display — handles Measurement and plain numbers."""
_format_value(v::Measurement) = string(v)
_format_value(v::Float64) = string(v)
_format_value(v) = string(v)

"""Prompt user after fit: save, adjust parameters, or quit."""
function _prompt_after_fit()
    menu = RadioMenu(["Save results", "Adjust parameters and refit", "Quit without saving"])
    choice = request("What next?", menu)
    if choice == 1
        return :save
    elseif choice == 2
        return :adjust
    else
        return :quit
    end
end

# ═══════════════════════════════════════════════════════════════════════════
# Step 7: Save results
# ═══════════════════════════════════════════════════════════════════════════

function _save_results(result, prob::ExchangeProblem)
    print("Save results to folder: ")
    input = strip(readline())
    isempty(input) && return nothing

    outputfolder = input
    prepare_outputfolder(outputfolder)

    # save combined plot
    plt = plot(result.plots...; size=(1200, 800))
    savefig(plt, joinpath(outputfolder, "exchange1d_fit.pdf"))
    @info "Saved $(joinpath(outputfolder, "exchange1d_fit.pdf"))"

    # save individual experiment plots
    for (i, p) in enumerate(result.plots)
        savefig(p, joinpath(outputfolder, "exchange1d_expt_$i.pdf"))
        @info "Saved $(joinpath(outputfolder, "exchange1d_expt_$i.pdf"))"
    end

    # save parameters as text
    paramfile = joinpath(outputfolder, "exchange1d_params.txt")
    open(paramfile, "w") do io
        println(io, "Exchange 1D fit results")
        println(io, "=" ^ 60)
        println(io)

        println(io, "Model: $(modelname(prob.model))")
        println(io, "Experiments:")
        for expt in prob.experiments
            println(io, "  $(short_expt_path(expt)) ($(typeof(expt).name.name), $(expt.field_teslas) T)")
        end
        println(io)

        # initial parameters
        items0 = _flatten_params_items(result.params0)
        maxlen = maximum(item -> length(item.label), items0)
        items_fit = _flatten_params_items(result.params)
        maxlen = max(maxlen, maximum(item -> length(item.label), items_fit))
        w = max(maxlen + 4, 16)

        println(io, "Initial parameters:")
        prev_section = ""
        for item in items0
            if item.section != prev_section
                prev_section != "" && println(io, "  ", "-" ^ 56)
                prev_section = item.section
            end
            println(io, "  ", rpad(item.label, w), _format_value(item.value))
        end

        println(io)
        println(io, "Fitted parameters:")
        prev_section = ""
        for item in items_fit
            if item.section != prev_section
                prev_section != "" && println(io, "  ", "-" ^ 56)
                prev_section = item.section
            end
            println(io, "  ", rpad(item.label, w), _format_value(item.value))
        end

        println(io)
        println(io, "chi2:          $(result.chi2)")
        println(io, "reduced chi2:  $(result.reduced_chi2)")
        println(io, "nobs:          $(result.nobs)")
        println(io, "nparams:       $(result.nparams)")
        println(io, "dof:           $(result.dof)")
    end
    @info "Saved $paramfile"

    return nothing
end
