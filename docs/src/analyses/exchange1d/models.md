# Exchange Models

Six exchange models are available, selected interactively at the start of the
[analysis workflow](overview.md). Three describe
intramolecular exchange between conformations of a single species; three describe
bimolecular binding, where populations are derived from a binding equation and the
total concentrations of the interacting molecules.

| Model | States | Parameters | Use case |
|---|---|---|---|
| **No exchange** | 1 | — | Null model (no exchange) |
| **Two-state exchange** | 2 | `kex`, `pB` | Intramolecular exchange between two conformations |
| **Three-state exchange** | 3 | `koffB`, `pB`, `koffC`, `pC` | Linear three-state exchange (A ⇌ B, A ⇌ C) |
| **Two-state binding** | 2 | `Kd`, `koff` | Bimolecular binding; requires sample concentrations |
| **Three-state binding** | 3 | `Kd1`, `koff1`, `Kd2`, `koff2` | Two independent, non-competing binding sites (A ⇌ B, A ⇌ C) |
| **Induced fit binding** | 3 | `Kd`, `koff`, `kclose`, `kopen` | Bimolecular binding followed by a conformational change (A + L ⇌ B ⇌ C) |

New models are added by defining a subtype of `AbstractModel` — see
[Extending Exchange1D](../../advanced/extending_exchange1d.md) if you need a mechanism not
covered here.

## Parameter fitting conventions

Model parameters are not fitted directly in the units you enter them in:

- Rate constants and dissociation constants (`kex`, `koff`, `Kd`, `kclose`, `kopen`, …) are
  fitted internally in log space (stored as `logkex`, `logKd`, etc.), which keeps them
  positive during optimisation without needing bounds. You still view and enter them as
  ordinary linear values (s⁻¹, or concentration units) in the parameter editor.
- Populations for the purely intramolecular models (`pB`, `pC`) are fitted internally via a
  free-energy-like parameter `dGB = ln(pA/pB)` (and `dGC` similarly), which keeps them
  between 0 and 1 for any real fitted value. Again, you only ever see and enter the
  population fraction itself.

Binding-model populations are not fitted parameters at all — they are *derived* at each
titration point from `Kd`/`koff` (and, for induced fit, `kclose`/`kopen`) together with the
total concentrations of the interacting molecules, via the equations below.

## Intramolecular exchange models

These apply to a single observed molecule exchanging between internal conformations or
states; no sample concentration information is needed.

### No exchange

A single non-exchanging state — the null model, useful as a fitting baseline or when no
exchange contribution is expected.

### Two-state exchange

Exchange between two conformations, `A ⇌ B`, parameterised by the total exchange rate `kex`
and the minor-state population `pB`:

```math
\mathbf{K} = \begin{pmatrix} -k_\text{ex} p_B & k_\text{ex} p_A \\ k_\text{ex} p_B & -k_\text{ex} p_A \end{pmatrix}
```

with ``p_A = 1 - p_B``.

### Three-state exchange

A linear three-state system with a shared ground state, `A ⇌ B` and `A ⇌ C` (no direct
`B ⇌ C` exchange), parameterised by the two off-rates `koffB`/`koffC` and populations
`pB`/`pC`:

```math
\mathbf{K} = \begin{pmatrix}
-k_{\text{off},B}\frac{p_B}{p_A} - k_{\text{off},C}\frac{p_C}{p_A} & k_{\text{off},B} & k_{\text{off},C} \\
k_{\text{off},B}\frac{p_B}{p_A} & -k_{\text{off},B} & 0 \\
k_{\text{off},C}\frac{p_C}{p_A} & 0 & -k_{\text{off},C}
\end{pmatrix}
```

with ``p_A = 1 - p_B - p_C``.

## Binding models

Binding models describe a two-molecule system: an **observed** species (role `:A`, e.g. the
protein being followed) and a **titrant** (role `:X`, e.g. the ligand). After selecting a
binding model you are prompted to assign these roles to the sample components found in the
experiment metadata — see the molecule mapping step of the [analysis workflow](overview.md).
Total concentrations for each titration point are read from the sample metadata attached to
each experiment (falling back to a manually entered value if that metadata is missing).

!!! tip
    All binding models require every experiment in the joint fit to carry, or be given, the
    total concentrations of both the observed molecule and the titrant for that particular
    sample — this is what lets populations vary correctly across a titration series.

### Two-state binding

Simple 1:1 binding, `A + X ⇌ AX`, parameterised by the dissociation constant `Kd` and
off-rate `koff`. Free titrant concentration ``X_\text{free}`` is obtained from the standard
quadratic binding equation, given total concentrations ``A_0`` (observed) and ``X_0``
(titrant):

```math
X_\text{free} = \tfrac{1}{2}\left(X_0 - A_0 - K_d + \sqrt{(K_d + A_0 + X_0)^2 - 4A_0 X_0}\right)
```

The bound population is then ``p_B = (X_0 - X_\text{free})/A_0``, and the pseudo-first-order
on-rate is ``k_\text{on,eff} = k_\text{off}\, p_B / p_A``.

### Three-state binding

Two independent, non-competing binding sites on the same observed molecule, `A + X ⇌ B` and
`A + X ⇌ C`, each with its own `Kd`/`koff` (`Kd1`/`koff1` and `Kd2`/`koff2`). There is no
direct `B ⇌ C` exchange. Free titrant concentration is obtained from the corresponding
quadratic equation for two parallel 1:1 sites, and populations follow from detailed balance
across each independent binding step. Use this when a single observed resonance reports on
two chemically distinct, mutually exclusive binding events (e.g. two ligands competing for
the same molecule, or two non-interacting sites) rather than a single two-step mechanism —
for the latter, see induced fit binding below.

### Induced fit binding

A two-step binding mechanism in which the initial bimolecular complex undergoes a
subsequent conformational change:

```math
A + L \underset{k_\text{off}}{\overset{k_\text{on}}{\rightleftharpoons}} B \underset{k_\text{open}}{\overset{k_\text{close}}{\rightleftharpoons}} C
```

Free ligand and equilibrium species concentrations are obtained analytically from the
coupled mass-action/mass-balance equations (ported from the TITAN `bmInducedFit` model),
given the dissociation constant `Kd`, off-rate `koff` for the bimolecular step, and
`kclose`/`kopen` for the conformational step. There is no direct `A ⇌ C` exchange:

```math
\mathbf{K} = \begin{pmatrix}
-k_{AB} & k_\text{off} & 0 \\
k_{AB} & -k_\text{off}-k_\text{close} & k_\text{open} \\
0 & k_\text{close} & -k_\text{open}
\end{pmatrix}, \qquad k_{AB} = \frac{k_\text{off}}{K_d}\, L_\text{free}
```

Use this model — rather than three-state binding — when states B and C represent the
**same** bound complex before and after a conformational change (e.g. an "open" encounter
complex isomerising to a "closed", conformationally selected complex), so that all ligand
binds through a single bimolecular step.
