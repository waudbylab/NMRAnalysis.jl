# ═══════════════════════════════════════════════════════════════════════════
# Console styling
# ═══════════════════════════════════════════════════════════════════════════

"""Bold, cyan section heading for a stage of the interactive workflow —
consistent with the styling `FitResult`'s `show` method already uses for its
own section titles."""
sectionheader(msg::AbstractString) = printstyled(msg * "\n"; bold=true, color=:cyan)

"""A quieter status line for routine progress/confirmation messages (e.g.
per-file loading) that don't need `@info`'s logging prefix or full emphasis
— `@info`/`@warn` are reserved here for things that genuinely need to stand
out."""
detail(msg::AbstractString) = printstyled(msg * "\n"; color=:light_black)

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
    println()
    sectionheader("Exchange 1D analysis")
    println()

    # ── 1. Model selection ──────────────────────────────────────────────
    model = _prompt_model()
    model === nothing && return nothing
    println()

    # ── 2. Build problem ────────────────────────────────────────────────
    prob = ExchangeProblem(filenames, model)
    println()
    sectionheader("Loaded $(length(prob.experiments)) experiments:")
    for expt in prob.experiments
        println("  $(short_expt_path(expt)) ($(typeof(expt).name.name), $(expt.field_teslas) T)")
    end
    println()

    # ── 3. Molecule mapping (if model requires it) ──────────────────────
    if nmolecules(model) > 1
        ok = _prompt_moleculemap!(model, prob)
        ok || return nothing
        _prompt_concentrations!(model, prob)
        println()
    end

    # ── 4. Integration ──────────────────────────────────────────────────
    _prompt_integration!(prob)
    println()

    # ── 5. Parameters + fit loop ────────────────────────────────────────
    p0 = defaultparams(prob)
    fixed = Set{Int}()
    while true
        p0 = _prompt_params(p0, prob, fixed)
        p0 === nothing && return nothing

        println()
        sectionheader("Fitting...")
        result = fit(prob, p0; fixed=fixed)
        display(result)

        plots = plot(result)
        display(combineplots(plots))

        overlays = overlayplots(result)
        isempty(overlays) || display(combineplots(overlays))

        action = _prompt_after_fit()
        if action == :save
            _save_results(result)
            return result
        elseif action == :adjust
            println()
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

function exchange1d(exptno::Integer)
    return exchange1d(string(exptno))
end

function exchange1d(exptnos::AbstractVector{<:Integer})
    return exchange1d(string.(exptnos))
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
    instances = sort!([T() for T in model_types]; by=modelorder)
    names = [modelname(m) for m in instances]
    push!(names, "Cancel")

    menu = RadioMenu(names)
    choice = request("Select exchange model:", menu)
    (choice == -1 || choice == length(names)) && return nothing

    model = instances[choice]
    detail("Selected model: $(modelname(model))")
    return model
end

# ═══════════════════════════════════════════════════════════════════════════
# Step 3: Molecule mapping
# ═══════════════════════════════════════════════════════════════════════════

function _prompt_moleculemap!(model, prob::ExchangeProblem)
    roles = molecules(model)

    # collect all unique sample molecule names, and the concentration(s) at
    # which each appears, across experiments
    all_names = String[]
    concentrations = Dict{String,Vector{Float64}}()
    for expt in prob.experiments
        for (name, conc) in sampleconcentrations(expt)
            name in all_names || push!(all_names, name)
            push!(get!(() -> Float64[], concentrations, name), conc)
        end
    end
    sort!(all_names)

    if isempty(all_names)
        detail("No sample concentration metadata found in any experiment — " *
               "enter a name and concentration for each molecule directly.")
        return _prompt_moleculemap_manual!(model, prob, roles)
    end

    println()
    sectionheader("Sample molecules:")
    for name in all_names
        println("  $name: $(_format_concrange(concentrations[name]))")
    end
    println()

    for (role, description) in roles
        options = copy(all_names)
        push!(options, "Cancel")
        menu = RadioMenu(options)
        choice = request("Assign :$role ($description) to:", menu)
        (choice == -1 || choice == length(options)) && return false
        model.moleculemap[role] = all_names[choice]
        detail("  :$role → $(all_names[choice])")
    end
    return true
end

"""
    _prompt_moleculemap_manual!(model, prob, roles) -> Bool

Fallback used by `_prompt_moleculemap!` when no experiment carries sample
concentration metadata: since there is no list of molecule names to pick
from, ask directly for a name and concentration for each role and populate
both `model.moleculemap` and `model.concentrations`.
"""
function _prompt_moleculemap_manual!(model, prob::ExchangeProblem, roles)
    for (role, description) in roles
        print("Name for :$role ($description): ")
        name = strip(readline())
        isempty(name) && return false
        model.moleculemap[role] = name

        conc = _prompt_value("Concentration of $name", nothing) do v
            return v > 0 || error("Concentration must be positive")
        end
        model.concentrations[name] = conc
        detail("  :$role → \"$name\" ($conc)")
    end
    return true
end

"""
    _prompt_concentrations!(model, prob::ExchangeProblem)

For each molecule role required by `model`, check whether every experiment
carries sample concentration metadata for it. If some (or all) experiments
lack it — e.g. because they weren't matched to a sample containing this
molecule — prompt the user for a fallback concentration value and store it in
`model.concentrations`, which `moleculeconcentration` (and hence
`modelconcentrations`) falls back to for experiments missing metadata.

Warning about a partial mismatch here, rather than only about a total one,
is what turns a `KeyError` buried deep in the fit (see issue #36) into an
actionable prompt at setup time.
"""
function _prompt_concentrations!(model, prob::ExchangeProblem)
    for (role, description) in molecules(model)
        name = model.moleculemap[role]
        haskey(model.concentrations, name) && continue

        unmatched = filter(expt -> !haskey(sampleconcentrations(expt), name), prob.experiments)
        isempty(unmatched) && continue

        if length(unmatched) < length(prob.experiments)
            paths = join(short_expt_path.(unmatched), "\n  ")
            @warn "\"$name\" ($description) is missing sample concentration metadata in " *
                  "$(length(unmatched))/$(length(prob.experiments)) experiments — they are " *
                  "not matched to a sample containing \"$name\":\n  $paths\n" *
                  "Enter a fallback concentration to use for these experiments."
        else
            detail("No sample concentration found for \"$name\" ($description)")
        end
        conc = _prompt_value("Concentration of $name", nothing) do v
            return v > 0 || error("Concentration must be positive")
        end
        model.concentrations[name] = conc
    end
    return nothing
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

    detail("Spectral range: $ppm_min to $ppm_max ppm")

    peakppm = _prompt_value("Peak position in ppm", default_peak) do v
        return ppm_min ≤ v ≤ ppm_max ||
               error("Peak position $v ppm is outside spectral range ($ppm_min to $ppm_max)")
    end

    ppmwidth = _prompt_value("Integration width in ppm", 0.1) do v
        return v > 0 || error("Width must be positive")
    end

    noiseppm = _prompt_value("Noise position in ppm", nothing) do v
        lo, hi = v - ppmwidth / 2, v + ppmwidth / 2
        return ppm_min ≤ lo && hi ≤ ppm_max ||
               error("Noise region $lo..$hi ppm is outside spectral range ($ppm_min to $ppm_max)")
    end

    detail("Integrating: peak=$peakppm ppm, noise=$noiseppm ppm, width=$ppmwidth ppm")
    return integrate!(prob, peakppm, noiseppm, ppmwidth)
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
    _prompt_params(p0, prob, fixed) -> ComponentArray or nothing

Interactive parameter editor. Returns edited parameters, or nothing if cancelled.

`fixed` is a `Set{Int}` of flat parameter indices (mutated in place) that are
held constant during fitting. Enter `fix` when editing a parameter to add it
to this set (freezing it at its current value), or `free` to remove it.

Parameters stored internally in log space (see `_islogparam`) are displayed
and entered here in linear space.
"""
function _prompt_params(p0::ComponentArray, prob::ExchangeProblem, fixed::Set{Int})
    state_labels = states(prob.model)
    fields = _unique_fields(prob.experiments)

    sectionheader("Parameter editor")
    println("  - select a parameter to view or edit its value")
    println("  - type 'fix' instead of a value to freeze it at its current value (excluded from fitting)")
    println("  - type 'free' to unfreeze a parameter previously fixed")
    println("  - press Enter with no input to leave it unchanged")
    println()

    while true
        items = _flatten_params_items(p0)
        labels = [_pretty_label(item, state_labels, fields) for item in items]
        maxlen = maximum(length, labels)
        menu_items = [rpad(label, maxlen + 2) * "= " *
                      _format_value(_displayvalue(item, p0)) *
                      (item.flat_index in fixed ? "  [fixed]" : "")
                      for (label, item) in zip(labels, items)]
        push!(menu_items, "▶ Continue to fit")
        push!(menu_items, "✕ Cancel")

        menu = RadioMenu(menu_items)
        choice = request("Review parameters (select to edit, 'fix'/'free' to lock, or continue):",
                         menu)

        # cancel
        (choice == -1 || choice == length(menu_items)) && return nothing
        # continue
        choice == length(menu_items) - 1 && break

        item = items[choice]
        label = labels[choice]
        statustag = item.flat_index in fixed ? " [fixed]" : ""
        print("  New value for $(label)$statustag [$(_format_value(_displayvalue(item, p0)))] " *
              "(or 'fix'/'free'): ")
        input = strip(readline())
        if isempty(input)
            # keep current value and fixed status
        elseif lowercase(input) in ("fix", "fixed")
            push!(fixed, item.flat_index)
        elseif lowercase(input) in ("free", "unfix", "unfixed")
            delete!(fixed, item.flat_index)
        else
            newval = tryparse(Float64, input)
            if newval === nothing
                @warn "Could not parse value: $input"
            elseif _islogparam(item)
                if newval <= 0
                    @warn "Value must be positive: $input"
                else
                    _set_param!(p0, item.flat_index, log(newval))
                end
            elseif _isdgparam(item)
                if !(0 < newval < 1)
                    @warn "Population must be between 0 and 1: $input"
                else
                    newdg = _dgvalue_for_population(p0.model, _paramkey(item), newval)
                    _set_param!(p0, item.flat_index, newdg)
                end
            else
                _set_param!(p0, item.flat_index, newval)
            end
        end
    end
    return p0
end

"""
A flattened parameter item: individual scalar elements, with array elements expanded.
- `label`: internal name (e.g. "spin.delta[1]")
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

"""
    _islogparam(item::_ParamItem) -> Bool

Whether a parameter is stored internally as its natural logarithm. By
convention this is signalled by a `log` prefix on the parameter name (e.g.
`model.logKd`, `model.logkoff`) — any positive-only parameter can be
switched to log-space storage (to keep it from going negative during
fitting) simply by naming it this way in `defaultparams`/
`default_spin_params`/`default_nuisance_params` and exponentiating it where
it's used. Users always see and enter the linear (exponentiated) value —
the internal representation is invisible to them either way.
"""
function _islogparam(item::_ParamItem)
    return startswith(String(_paramkey(item)), "log")
end

"""Raw ComponentArray key for an item, i.e. its label with the section prefix
(and any trailing `[i]` array index) stripped."""
function _paramkey(item::_ParamItem)
    name = item.label
    prefix = item.section * "."
    if startswith(name, prefix)
        name = name[(length(prefix) + 1):end]
    end
    m = match(r"\[(\d+)\]$", name)
    m !== nothing && (name = name[1:(m.offset - 1)])
    return Symbol(name)
end

"""
    _isdgparam(item::_ParamItem) -> Bool

Whether a parameter follows the `dG<X>` convention (`model.dGB`, `model.dGC`,
...): an internal, unconstrained free energy difference (relative to a
reference state, e.g. `dGB = G_B - G_A`) standing in for a population
fraction. Unlike `_islogparam`, the user-facing value (`pB`, `pC`, ...) is
not a function of this one parameter alone — it also depends on the
*other* `dG*` values currently stored in the same `model` section, since
the underlying population fractions must sum to 1 (see `_dgpopulation` /
`_dgvalue_for_population`).
"""
function _isdgparam(item::_ParamItem)
    item.section == "model" || return false
    return match(r"^dG[A-Z]$", String(_paramkey(item))) !== nothing
end

"""`dG<X>` keys (e.g. `:dGB`, `:dGC`) present in a model parameter section."""
function _dgkeys(pmodel::ComponentArray)
    return [k for k in keys(pmodel)
            if match(r"^dG[A-Z]$", String(k)) !== nothing]
end

"""
    _dgpopulation(pmodel::ComponentArray, key::Symbol)

Population fraction for `key` (e.g. `:dGB` → `pB`), reconstructed from all
`dG*` values currently stored in `pmodel` (the model parameter section),
using the same reparametrization as `TwoStateModel`/`ThreeStateModel`'s
`populations` (population of the implicit reference state, e.g. A, is
`1 / (1 + Σ exp(-dG))`).
"""
function _dgpopulation(pmodel::ComponentArray, key::Symbol)
    denom = 1 + sum(exp(-pmodel[k]) for k in _dgkeys(pmodel))
    return exp(-pmodel[key]) / denom
end

"""
    _dgvalue_for_population(pmodel::ComponentArray, key::Symbol, pnew)

Inverse of `_dgpopulation`: the `dG` value for `key` that reproduces the
population fraction `pnew`, holding all *other* `dG*` values in `pmodel`
fixed at their current values.
"""
function _dgvalue_for_population(pmodel::ComponentArray, key::Symbol, pnew)
    Sother = sum(k == key ? 0.0 : exp(-pmodel[k]) for k in _dgkeys(pmodel))
    return -log(pnew * (1 + Sother) / (1 - pnew))
end

"""Convert a parameter item's stored value to linear (display/entry) space.
`params` is the enclosing top-level ComponentArray (with `.model`/`.spin`/
`.nuisance` sections), needed to reconstruct `dG`-parametrized populations."""
function _displayvalue(item::_ParamItem, params::ComponentArray)
    _islogparam(item) && return exp(item.value)
    _isdgparam(item) && return _dgpopulation(params.model, _paramkey(item))
    return item.value
end

# ═══════════════════════════════════════════════════════════════════════════
# Pretty parameter labels
# ═══════════════════════════════════════════════════════════════════════════

const _PARAM_DISPLAY_NAMES = Dict("delta" => "δ",
                                  "R2" => "R₂",
                                  "R1" => "R₁",
                                  "kex" => "kex",
                                  "koffB" => "koff,B",
                                  "koffC" => "koff,C",
                                  "koff" => "koff",
                                  "Kd" => "Kd",
                                  "R1_I0" => "I₀ (R₁)",
                                  "R1_inv_factor" => "Inversion factor (R₁)")

const _SECTION_TITLES = Dict("model" => "Exchange parameters",
                             "spin" => "Spin parameters",
                             "nuisance" => "Nuisance parameters")

"""Collect unique magnetic field strengths from experiments."""
function _unique_fields(experiments)
    fields = Float64[]
    for expt in experiments
        expt.field_teslas ∉ fields && push!(fields, expt.field_teslas)
    end
    return sort!(fields)
end

"""Format a field strength for display, e.g. 22.31 → \"22.31 T\"."""
function _format_field(field_teslas::Float64)
    s = string(round(field_teslas; digits=2))
    s = rstrip(rstrip(s, '0'), '.')
    return s * " T"
end

"""
    _format_concrange(values::Vector{Float64}) -> String

Format the concentration(s) at which a molecule appears across experiments:
a single value if consistent, otherwise the range observed (which usually
signals a titration series or inconsistently-entered sample metadata).
"""
function _format_concrange(values::Vector{Float64})
    lo, hi = extrema(values)
    return lo == hi ? string(lo) : "$lo to $hi (varies across experiments)"
end

"""
    _pretty_label(item, state_labels, unique_fields) -> String

Convert an internal parameter label like `spin.R2_22p31T[1]` to a
human-readable label like `R₂ (A)` or `R₂ (A) [22.31 T]`.
"""
function _pretty_label(item::_ParamItem, state_labels::Vector{String},
                       unique_fields::Vector{Float64})
    # strip section prefix
    name = item.label
    prefix = item.section * "."
    if startswith(name, prefix)
        name = name[(length(prefix) + 1):end]
    end

    # strip log marker — displayed/entered in linear space regardless of
    # internal storage, so the label should read the same either way
    startswith(name, "log") && (name = name[4:end])

    # extract array index
    state_idx = nothing
    m = match(r"\[(\d+)\]$", name)
    if m !== nothing
        state_idx = parse(Int, m[1])
        name = name[1:(m.offset - 1)]
    end

    # strip field label (e.g. _22p31T)
    show_field = length(unique_fields) > 1
    field_val = nothing
    for f in unique_fields
        fl = string(field_label(f))
        if occursin(fl, name)
            name = replace(name, fl => "")
            name = replace(name, "__" => "_")
            name = strip(name, '_')
            field_val = f
            break
        end
    end

    # map to display name — dG<X> (e.g. dGB) is a population reparametrized
    # for fitting (see `_isdgparam`) but is shown/entered as p<X> (e.g. pB)
    dgmatch = match(r"^dG([A-Z])$", name)
    pretty = if dgmatch !== nothing
        "p" * dgmatch[1]
    else
        get(_PARAM_DISPLAY_NAMES, name, name)
    end

    # add state label
    if state_idx !== nothing
        if state_idx <= length(state_labels)
            pretty *= " ($(state_labels[state_idx]))"
        else
            pretty *= " [$state_idx]"
        end
    end

    # add field if multiple fields present
    if show_field && field_val !== nothing
        pretty *= " [$(_format_field(field_val))]"
    end

    return pretty
end

"""Format a value for display — handles Measurement and plain numbers."""
_format_value(v::Measurement) = string(v)
_format_value(v::Float64) = string(v)
_format_value(v) = string(v)

# ═══════════════════════════════════════════════════════════════════════════
# Step 6: Plotting + post-fit prompt
# ═══════════════════════════════════════════════════════════════════════════

"""
    combineplots(plots) -> Plot

Create a combined figure from individual experiment plots, with scaled font sizes
and figure dimensions so that the result is legible even with many experiments.
"""
function combineplots(plots)
    n = length(plots)
    ncols = min(n, 4)
    nrows = ceil(Int, n / ncols)

    # scale figure: each experiment column ~350px wide, each row pair ~280px tall
    w = max(1200, ncols * 350)
    h = max(800, nrows * 280)

    plt = plot(plots...; size=(w, h))

    for sp in plt.subplots
        sp[:titlefontsize] = 8

        # font sizes must be set on each axis object directly
        for axis in (:xaxis, :yaxis)
            sp[axis].plotattributes[:guidefontsize] = 7
            sp[axis].plotattributes[:tickfontsize] = 6
        end

        # marker/line sizes are in absolute points, not scaled to subplot
        # size, so combining many full-size plots into a shrunk grid makes
        # points and lines look huge relative to the (now much smaller)
        # subplot — scale them down to match the reduced font sizes. The
        # stroke outline in particular needs to be well below the marker
        # radius or small markers turn into solid blobs.
        for series in sp.series_list
            series[:markersize] = min(series[:markersize], 1.5)
            series[:markerstrokewidth] = min(series[:markerstrokewidth], 0.25)
            series[:linewidth] = min(series[:linewidth], 1)
        end
    end

    return plt
end

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

function _save_results(result::FitResult)
    print("Save results to folder: ")
    input = strip(readline())
    isempty(input) && return nothing

    outputfolder = input
    prepare_outputfolder(outputfolder)
    saved = String[]

    # save per-experiment plots, combined into one grid figure
    plots = plot(result)
    plt = combineplots(plots)
    savefig(plt, joinpath(outputfolder, "exchange1d_fit.pdf"))
    push!(saved, "exchange1d_fit.pdf")

    for (i, p) in enumerate(plots)
        savefig(p, joinpath(outputfolder, "exchange1d_expt_$i.pdf"))
        push!(saved, "exchange1d_expt_$i.pdf")
    end

    # save overlays of similar experiments separately, so the per-experiment
    # grid above doesn't grow denser (and less legible) as they're added
    overlays = overlayplots(result)
    if !isempty(overlays)
        overlayplt = combineplots(overlays)
        savefig(overlayplt, joinpath(outputfolder, "exchange1d_overlay.pdf"))
        push!(saved, "exchange1d_overlay.pdf")

        for (i, p) in enumerate(overlays)
            savefig(p, joinpath(outputfolder, "exchange1d_overlay_$i.pdf"))
            push!(saved, "exchange1d_overlay_$i.pdf")
        end
    end

    # save parameters as text
    paramfile = joinpath(outputfolder, "exchange1d_params.txt")
    open(paramfile, "w") do io
        return show(io, MIME("text/plain"), result)
    end
    push!(saved, "exchange1d_params.txt")

    # save data filenames, experiment parameters, and sample information
    infofile = joinpath(outputfolder, "exchange1d_experiments.txt")
    open(infofile, "w") do io
        return writeexperimentsummary(io, result.prob)
    end
    push!(saved, "exchange1d_experiments.txt")

    println()
    sectionheader("Saved to $outputfolder:")
    for name in saved
        println("  $name")
    end
    println()

    return nothing
end

"""
    writeexperimentsummary(io::IO, prob::ExchangeProblem)

Write a plain-text summary of the experiments underlying `prob`: the model,
its molecule-to-sample mapping (if any), and — for each experiment — its
source filename, key acquisition parameters (via `experimentinfo`), and
sample concentrations. Saved alongside fit results by `_save_results` so an
analysis can always be traced back to its source data (issue #37).
"""
function writeexperimentsummary(io::IO, prob::ExchangeProblem)
    model = prob.model
    println(io, "Exchange 1D — experiment summary")
    println(io, "Model: $(modelname(model))")
    println(io)

    if nmolecules(model) > 1
        println(io, "Molecule mapping:")
        for (role, description) in molecules(model)
            name = model.moleculemap[role]
            fallback = get(model.concentrations, name, nothing)
            suffix = fallback === nothing ? "" : " (fallback concentration: $fallback)"
            println(io, "  :$role ($description) -> \"$name\"$suffix")
        end
        println(io)
    end

    if prob.integration !== nothing
        println(io, "Integration parameters:")
        println(io, "  Peak position: $(prob.integration.peakppm) ppm")
        println(io, "  Noise position: $(prob.integration.noiseppm) ppm")
        println(io, "  Integration width: $(prob.integration.ppmwidth) ppm")
        println(io)
    end

    for (i, expt) in enumerate(prob.experiments)
        println(io, "Experiment $i: $(short_expt_path(expt))")
        println(io, "  File: $(expt.spec[:filename])")
        for (key, value) in experimentinfo(expt)
            println(io, "  $key: $value")
        end
        sc = sampleconcentrations(expt)
        if !isempty(sc)
            println(io, "  Sample concentrations:")
            for name in sort!(collect(keys(sc)))
                println(io, "    $name: $(sc[name])")
            end
        end
        println(io)
    end
    return nothing
end
