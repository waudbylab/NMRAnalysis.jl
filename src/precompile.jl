# Run the 1D and 2D analyses on synthetic spectra, so that a user's first analysis doesn't
# compile them: the fitting, the window up to the point it would be drawn (which needs a
# display, unavailable while precompiling), and the output folder Save writes. The spectra
# are written as Bruker data and read back, since the example data would need a download.

@setup_workload begin
    lorentzian(u, u0, w) = w^2 / (w^2 + (u - u0)^2)
    axis(dim, δ, bf, nucleus) = dim(δ; metadata=Dict{Symbol,Any}(:bf => bf, :label => nucleus,
                                                                  :window => ExponentialWindow(5.0, 0.1)))
    H = collect(range(10.0, 6.0; length=512))
    N = collect(range(130.0, 105.0; length=64))
    T = [0.01, 0.05, 0.1, 0.2, 0.4]
    peaks = [("A12", 8.2, 120.0, 5.0), ("G40", 7.5, 115.0, 2.0)]   # label, ¹H, ¹⁵N, R

    series1d = [sum(1e4 * exp(-R * t) * lorentzian(h, h0, 0.01) for (_, h0, _, R) in peaks)
                for h in H, t in T]
    series2d = [sum(1e4 * exp(-R * t) * lorentzian(h, h0, 0.02) * lorentzian(n, n0, 0.3)
                    for (_, h0, n0, R) in peaks) for h in H, n in N, t in T]
    noise(A) = A .+ 0.01 .* sin.(reshape(1:length(A), size(A)))   # a noise level that is not zero

    @compile_workload begin
        # a failure here should cost compilation, not stop the package precompiling at all
        try
            with_logger(NullLogger()) do
                mktempdir() do dir
                    cd(dir) do
                        savenmr("1", NMRData(noise(series1d),
                                             (axis(F1Dim, H, 600.13e6, "1H"), X2Dim(1:length(T)))))
                        savenmr("2", NMRData(noise(series2d),
                                             (axis(F1Dim, H, 600.13e6, "1H"),
                                              axis(F2Dim, N, 60.82e6, "15N"), X3Dim(1:length(T)))))
                        write("peaks.txt", join(("$l $h $n" for (l, h, n, _) in peaks), "\n"))

                        screenshot(nothing; save=true) do
                            return relaxation1d("1"; T, model=:exponential, prompt=false)
                        end
                        screenshot(nothing; save=true) do
                            return relaxation2d("2"; T, peaklist="peaks.txt", prompt=false)
                        end
                        screenshot(nothing; save=true) do
                            return titration2d("2"; L0=[0, 10, 20, 40, 80], peaklist="peaks.txt",
                                               prompt=false)
                        end
                    end
                end
            end
        catch e
            @warn "NMRAnalysis precompile workload failed" exception = e
        end
    end
end
