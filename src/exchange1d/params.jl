"""
    defaultparams(prob::ExchangeProblem) -> ComponentArray

Assemble a default parameter ComponentArray from the model and experiments.

Each experiment contributes spin parameters (via `default_spin_params`) and
nuisance parameters (via `default_nuisance_params`). Both return flat
`Vector{Pair{Symbol,Any}}` with self-describing keys (e.g. `:R2_14p1T`,
`:R1_14p1T_I0`). Duplicate entries are merged automatically.

The returned ComponentArray has three sections:
- `model`: exchange model parameters (kex, pB, Kd, koff, etc.)
- `spin`: chemical shifts (delta, ppm) and field-dependent relaxation rates
- `nuisance`: per-experiment amplitude and fitting parameters
"""
function defaultparams(prob::ExchangeProblem)
    N = nstates(prob.model)

    spin_pairs = _collect_pairs(expt -> default_spin_params(expt, N), prob.experiments)
    nuisance_pairs = _collect_pairs(default_nuisance_params, prob.experiments)

    return ComponentArray(;
                          model=defaultparams(prob.model),
                          spin=ComponentArray(; spin_pairs...),
                          nuisance=isempty(nuisance_pairs) ? ComponentArray() :
                                   ComponentArray(; nuisance_pairs...),)
end

"""Collect pairs from all experiments, keeping the first occurrence of each key."""
function _collect_pairs(f, experiments)
    pairs = Pair{Symbol,Any}[]
    seen = Set{Symbol}()
    for expt in experiments
        for pair in f(expt)
            if pair.first ∉ seen
                push!(seen, pair.first)
                push!(pairs, pair)
            end
        end
    end
    return pairs
end
