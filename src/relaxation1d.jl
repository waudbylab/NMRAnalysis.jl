# Register analysis rule for 1D nutation calibration experiments
register_analysis!(["1d", "relaxation"],
                   ["R1"],
                   exp -> relaxation1d(exp.filename),
                   "1D R1 relaxation")
register_analysis!(["1d", "relaxation"],
                   ["R2"],
                   exp -> relaxation1d(exp.filename),
                   "1D R2 relaxation")

"""
    relaxation1d(filename)
    relaxation1d(spec, signalselector, noiseselector)

Analyze 1D NMR relaxation experiments (T1 or T2) by fitting signal decay/recovery curves.

Interactively select integration and noise regions, then fit relaxation data to extract
relaxation rates and times. Supports standard exponential decay and inversion-recovery experiments.

Required annotations:
- `relaxation.duration`: relaxation times used in the experiment
- `relaxation.model`   : model to fit, either "exponential_decay" or "inversion_recovery"

# Arguments
- `filename`: Path to Bruker experiment folder

# Example
```julia
relaxation1d("path/to/experiment")
```
"""
function relaxation1d(filename::String)
    ispath(filename) || throw(ArgumentError("No such file or directory: $filename"))

    spec = loadnmr(filename)
    hasannotations(spec) ||
        throw(ArgumentError("Experiment must have annotations for analysis"))

    ir = annotations(spec, :relaxation, :model) == "inversion_recovery"
    slice = ir ? spec[:, end] : spec[:, 1] # for IR, use last slice where signal is most recovered
    signal, noise = get1dregionandnoise(slice)

    return relaxation1d(spec, signal, noise)
end

function relaxation1d(spec::NMRData{T,2}, signalselector, noiseselector) where {T}
    # get annotations
    ir = annotations(spec, :relaxation, :model) == "inversion_recovery"
    tau = annotations(spec, :relaxation, :duration)
    type = annotations(spec, :relaxation, :type)
    nuc = annotations(spec, :relaxation, :channel)

    # integrate regions
    noise = vec(data(sum(spec[noiseselector, :]; dims=F1Dim)))
    noise = std(noise)

    integrals = vec(data(sum(spec[signalselector, :]; dims=F1Dim)))

    # normalise by max value
    noise /= maximum(integrals)
    integrals /= maximum(integrals)

    # calculate least-squares fit
    model(t, p) = ir ? p[1] * (1 .- p[3] * exp.(-t * p[2])) : p[1] * exp.(-t * p[2])
    p0 = ir ? [1.0, 2 / maximum(tau), 2.0] : [1.0, 2 / maximum(tau)]

    fit = curve_fit(model, tau, integrals, p0)
    pars = coef(fit) .± stderror(fit)
    I0 = coef(fit)[1]
    rate = pars[2]

    # make output plot
    x = LinRange(0, maximum(tau) * 1.1, 100)
    yfit = model(x, coef(fit))

    p1 = scatter(tau, (integrals .± noise) / I0; label="observed",
                 frame=:box,
                 xlabel="Relaxation time / s",
                 ylabel="Integrated signal",
                 title="Relaxation rate = $rate s⁻¹",
                 grid=nothing)
    plot!(p1, x, yfit / I0;
          label="fit",
          z_order=:back)
    # add invisible hline to force axis to go to zero
    hline!(p1, [0]; color=:black, lw=0, primary=false)

    @info "Relaxation analysis results for $(spec[:filename]):"
    @info " - $type relaxation of $nuc"
    @info " - Model: $(ir ? "Inversion recovery" : "Exponential decay")"
    @info " - Integration region: $signalselector ppm"
    @info " - Noise region: $noiseselector ppm"
    @info " - Fitted relaxation rate: $rate s⁻¹"
    @info " - Fitted relaxation time: $(1/rate) s"
    if ir
        @info " - Inversion-recovery amplitude: $(pars[3])"
    end

    display(p1)

    return (rate=rate, relaxation_time=1 / rate, type=type, nucleus=nuc, plt=p1)
end