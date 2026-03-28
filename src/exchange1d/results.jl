"""Compact one-line display for FitResult."""
function Base.show(io::IO, result::FitResult)
    name = modelname(result.prob.model)
    return print(io, "FitResult($name, χ² = $(round(result.chi2; digits=2)), dof = $(result.dof))")
end

"""Pretty multi-line display for FitResult with parameter tables and fit statistics."""
function Base.show(io::IO, ::MIME"text/plain", result::FitResult)
    prob = result.prob
    state_labels = states(prob.model)
    fields = _unique_fields(prob.experiments)

    items0 = _flatten_params_items(result.params0)
    items_fit = _flatten_params_items(result.params)

    println(io)
    printstyled(io, "═"^60 * "\n"; bold=true)
    printstyled(io, "  Exchange 1D Fit Results\n"; bold=true)
    printstyled(io, "  Model: $(modelname(prob.model))\n"; bold=true)
    printstyled(io, "═"^60 * "\n"; bold=true)

    # print tables grouped by section
    sections = unique(item.section for item in items_fit)
    for section in sections
        sec_items0 = filter(i -> i.section == section, items0)
        sec_items = filter(i -> i.section == section, items_fit)
        title = get(_SECTION_TITLES, section, section)

        println(io)
        printstyled(io, "  $title\n"; bold=true, color=:cyan)

        labels = [_pretty_label(item, state_labels, fields) for item in sec_items]
        initial = [_format_value(item.value) for item in sec_items0]
        fitted = [_format_value(item.value) for item in sec_items]
        tdata = hcat(labels, initial, fitted)

        pretty_table(io, tdata;
                     header=["Parameter", "Initial", "Fitted"],
                     alignment=[:l, :r, :r],
                     tf=tf_unicode_rounded,
                     crop=:none,
                     header_crayon=Crayon(; bold=true),)
    end

    # fit statistics
    println(io)
    printstyled(io, "  Fit statistics\n"; bold=true, color=:cyan)
    stats = ["χ²" string(round(result.chi2; digits=2));
             "Reduced χ²" string(round(result.reduced_chi2; digits=4));
             "Observations" string(result.nobs);
             "Parameters" string(result.nparams);
             "DOF" string(result.dof)]
    pretty_table(io, stats;
                 header=["Statistic", "Value"],
                 alignment=[:l, :r],
                 tf=tf_unicode_rounded,
                 crop=:none,
                 header_crayon=Crayon(; bold=true),)
    return println(io)
end

"""Plot all experiments in a FitResult, returning a vector of per-experiment plots."""
Plots.plot(result::FitResult) = plot_result(result.prob, result)
