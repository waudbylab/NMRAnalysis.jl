"""Compact one-line display for FitResult."""
function Base.show(io::IO, result::FitResult)
    name = modelname(result.prob.model)
    return print(io,
                 "FitResult($name, χ² = $(round(result.chi2; digits=2)), dof = $(result.dof))")
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
    printstyled(io, "  NMRAnalysis.jl v$(pkgversion(NMRAnalysis))\n"; bold=true)
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
        initial = [_format_value(_displayvalue(item, result.params0))
                   for item in sec_items0]
        fitted = [_format_value(_displayvalue(item, result.params)) *
                  (item.flat_index in result.fixed ? " (fixed)" :
                   item.flat_index in result.atbound ? " (at bound)" : "")
                  for item in sec_items]
        tdata = hcat(labels, initial, fitted)

        pretty_table(io, tdata;
                     header=["Parameter", "Initial", "Fitted"],
                     alignment=[:l, :r, :r],
                     tf=tf_unicode_rounded,
                     crop=:none,
                     header_crayon=Crayon(; bold=true))
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
                 header_crayon=Crayon(; bold=true))

    # correlation matrix of the fitted (non-fixed) parameters
    freeitems = items_fit[result.freeidx]
    if length(freeitems) > 1
        clabels = [_pretty_label(item, state_labels, fields) for item in freeitems]
        n = length(clabels)
        body = [string(round(result.cor[i, j]; digits=2)) for i in 1:n, j in 1:n]

        println(io)
        printstyled(io, "  Correlation matrix\n"; bold=true, color=:cyan)
        pretty_table(io, hcat(clabels, body);
                     header=vcat(["Parameter"], clabels),
                     alignment=vcat([:l], fill(:r, n)),
                     tf=tf_unicode_rounded,
                     crop=:none,
                     header_crayon=Crayon(; bold=true))
    end

    # strongly correlated parameter pairs, if any — see strongcorrelations
    correlations = strongcorrelations(result.cor)
    if !isempty(correlations)
        label1 = [_pretty_label(freeitems[i], state_labels, fields)
                  for (i, _, _) in correlations]
        label2 = [_pretty_label(freeitems[j], state_labels, fields)
                  for (_, j, _) in correlations]
        rvals = [string(round(r; digits=3)) for (_, _, r) in correlations]

        println(io)
        printstyled(io, "  Strongly correlated parameters (|r| ≥ 0.95)\n"; bold=true,
                    color=:yellow)
        pretty_table(io, hcat(label1, label2, rvals);
                     header=["Parameter", "Parameter", "r"],
                     alignment=[:l, :l, :r],
                     tf=tf_unicode_rounded,
                     crop=:none,
                     header_crayon=Crayon(; bold=true))
    end

    return println(io)
end
