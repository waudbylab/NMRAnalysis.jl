# Peak Tracking

**Peak tracking** analyses track how peak positions (and linewidths) change from plane to plane.
It is the common core for measurements where peaks move: fit a [titration](titration.md) binding
isotherm, or measure a coupling/RDC (see [Couplings and RDCs](rdc.md)). Use `peaktrack2d`
directly when you just want the per-plane positions and linewidths without a physical model.

```julia
using NMRAnalysis

# A series of 2D spectra in which peaks move
peaktrack2d(["11/pdata/1", "12/pdata/1", "13/pdata/1"])
```

- `inputfilenames`: a single pseudo-3D dataset, or a vector of 2D spectra (one per plane).

Each peak holds an **independent position, linewidth and amplitude in every plane**, and each
plane is fitted separately. This keeps the fit well-conditioned and decoupled: adjusting one
plane does not disturb the others.

![Screenshot of peak fitting](../../assets/peaktrack-demo.mov)

## Adding and positioning peaks

- **(A) — add a peak.** Adding walks through the planes so you mark the peak in each one:
  press `A` with the cursor on the peak in the current plane, then for each subsequent plane
  (stepping forwards and wrapping back to the start) press:
  - `a` to mark the peak's position in that plane and advance,
  - `space` to copy the current position into all remaining planes and finish,
  - `esc` to cancel.

  When at least one peak has already been fitted, new peaks are pre-seeded with the average
  displacement pattern of the fitted peaks, so they inherit the common motion.
- **(T) — add and track.** Drops a peak and follows the intensity maximum across the planes
  automatically (anchored at the current plane, propagated outwards). Best for well-resolved
  peaks; not available for RDC experiments.
- **Drag** a peak's handle to correct its position in the current plane. The handle enlarges
  when hovered. Dragging re-estimates that plane's amplitude.
- **(D)** deletes and **(R)** renames the selected peak, as elsewhere.

## Visualisation

- Each peak's **trajectory** across the planes is drawn as a coloured polyline (red = needs
  fitting, blue = fitted, green = selected).
- The **Show all** toggle overlays every plane's contours faintly, for context when a peak
  moves a long way.

## Output

**Save to folder** writes `results.csv` with each peak's per-plane positions (`x[1]`, `x[2]`,
…), linewidths and amplitudes, and any model-derived parameters. **Load peak list** restores
them, including the fitting radii. The default summary plot shows the combined chemical-shift
perturbation Δδ against residue number, one series per plane.
