# Register analysis rule for 1D nutation calibration experiments
register_analysis!(["1d", "calibration"],
                   ["nutation"],
                   exp -> analyse_1d_nutation(exp.filename),
                   "1D nutation calibration")

"""
    analyse_1d_nutation(filename)

Analyse a 1D nutation experiment stored in `filename` and return a plot
showing the nutation data and fit.

Required annotations:
- `calibration.channel` : channel to use for calibration
- `calibration.power`   : power level used in the experiment
- `calibration.duration`: pulse durations used in the experiment
- `calibration.model`   : model to fit, either "sine_modulated" or "cosine_modulated"

Optional annotations:
- `calibration.offset`  : offset in ppm to select the observed signal
"""
function analyse_1d_nutation(filename)
    expt = loadnmr(filename)
    hasannotations(expt) ||
        throw(ArgumentError("Experiment must have annotations for calibration"))

    ref_p, ref_pl = referencepulse(expt, annotations(expt, :calibration, :channel))
    cal_pl = annotations(expt, :calibration, :power)
    rough_hz = hz(cal_pl, ref_pl, ref_p, 90.0)

    offset = annotations(expt, :calibration, :offset)
    if isnothing(offset)
        offset = 0.0
    end
    δobs = ppm(offset, dims(expt, F1Dim))

    t = annotations(expt, :calibration, :duration)

    y = data(expt[Near(δobs), :]) / expt[:noise]

    model(t, p) =
        if annotations(expt, :calibration, :model) == "sine_modulated"
            @. p[1] * sin(2π * p[2] * t) * exp(-p[3] * t)
        elseif annotations(expt, :calibration, :model) == "cosine_modulated"
            @. p[1] * cos(2π * p[2] * t) * exp(-p[3] * t)
        else
            error("Unknown calibration model")
        end
    p0 = [maximum(y), rough_hz, 1.0]
    fit = curve_fit(model, t, y, p0)

    nut_hz = coef(fit)[2] ± stderror(fit)[2]
    R = coef(fit)[3] ± stderror(fit)[3]
    inhomogeneity = R / (2π * nut_hz)

    pulse90 = 1 / (4 * nut_hz) # in s
    @info "Nutation calibration results for $filename:"
    @info " - Power level: $(db(cal_pl)) dB"
    @info " - Nutation frequency ν₁: $nut_hz Hz"
    @info " - 90° pulse length: $(1e6*pulse90) µs"
    @info " - Decay rate: $R s⁻¹"
    @info " - B₁ inhomogeneity (R/2πν₁): $(inhomogeneity * 100) %"

    # plot results
    tfit = LinRange(0, 1.1 * maximum(t), 200)
    yfit = model(tfit, coef(fit))
    plt = scatter(t, y .± 1; label="Data",
                  xguide="Pulse duration / s",
                  yguide="SNR",
                  title="Nutation calibration ($(short_expt_path(filename)))",
                  titlefontsize=13,
                  frame=:box, grid=nothing)
    plot!(plt, tfit, yfit; label="Fit ($nut_hz Hz)")
    return (pulse90=pulse90, nut_hz=nut_hz, power_level=cal_pl, plt=plt)
end