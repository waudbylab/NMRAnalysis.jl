# Peak tracking

`peaktrack2d` follows peaks that move from plane to plane, and reports each peak's
position, linewidths and amplitude in every plane. Each plane is fitted separately. The
same window is used by [`titration2d`](titration.md) and [`rdc2d`](rdc.md), which fit a
model to the positions.

## Usage

```julia
using NMRAnalysis
peaktrack2d(11:13)
```

### Example

The phosphorylation of c-Jun, followed by ¹H,¹⁵N SOFAST-HMQC spectra:

```julia
peaktrack2d(example("kinetics2d_cjun"))
```

```@setup peaktrack2d
using NMRAnalysis
NMRAnalysis.screenshot("peaktrack2d.png") do
    peaktrack2d(example("kinetics2d_cjun"))
end
```

![The peaktrack2d window, with every plane's contours shown](peaktrack2d.png)

![Tracking peaks that move from plane to plane](../../assets/peaktrack-demo.mov)

## Add a peak

| Key | Action |
|---|---|
| `T` | Add a peak under the cursor, and follow the intensity maximum through the planes |
| `L` | Place a peak along a line: hold `L` at one end of its path, drag to the other, and release. Or press `L` at each end |
| `A` | Mark the peak by hand in each plane: `A` to mark and move on, `Space` to copy the position to the remaining planes, `Esc` to cancel |
| `S` | Show or hide the other planes' contours |

`T` suits well-resolved peaks, `L` peaks in fast exchange that move in a straight line, and
`A` crowded regions. After a peak is placed, drag it to correct its position in the current
plane. `T` and `L` aren't available in `rdc2d`.

Each peak's path is drawn across the planes, in the colours of the peak markers.

## Output

`series.csv` holds the position, linewidths and amplitude of each peak in each plane, and
the summary plot shows the combined shift change against residue number.
