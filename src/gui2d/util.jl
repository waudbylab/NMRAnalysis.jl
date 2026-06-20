nearest(A::AbstractArray, t) = findmin(abs.(A .- t))[1]
findnearest(A::AbstractArray, t) = findmin(abs.(A .- t))[2]

function choptitle(title, maxlength=30)
    if length(title) > maxlength
        title[1:maxlength] * "…"
    else
        title
    end
end

function maskellipse!(mask, x, y, x0, y0, xradius, yradius)
    # @debug "masking ellipse at $x0, $y0 with radii $xradius, $yradius" maxlog=10
    x = x .- x0
    y = y' .- y0
    fx = @. yradius^2 * x^2
    fy = @. xradius^2 * y^2
    f = @. fx + fy - xradius^2 * yradius^2
    return mask[f .≤ 0] .= true
end

function flatten_with_nan_separator(vectors::Vector{Vector{Point2f}})
    isempty(vectors) && return Point2f[]

    separator = Point2f(NaN, NaN)
    result = reduce(vectors[2:end]; init=vectors[1]) do acc, subvector
        return vcat(acc, [separator], subvector)
    end

    # Remove the trailing separator if it exists
    return length(result) > 0 ? result[1:(end - 1)] : result
end

"Standard amino acid one-letter codes."
const STANDARD_RESIDUES = Set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])

"""
    ResidueLabel

Parsed representation of a peak label.

# Fields
- `resnum`: residue number (negative for non-standard residue types, e.g. the
  default unassigned `X#` peaks; `0` if the label contains no digits)
- `onelettercode`: residue one-letter code (`'?'` if none could be identified)
- `atom`: atom name (uppercase, verbatim), or `""` for backbone amides
"""
struct ResidueLabel
    resnum::Int
    onelettercode::Char
    atom::String
end

"""
    parse_label(label)::ResidueLabel

Best-effort parse of a peak label into residue number, one-letter code, and atom
name. Handles, e.g.:
- Amides: `"G10"`, reversed `"7A"`
- Methyls/sidechains: `"I13CD1"`, `"L26CD1"`/`"L26CD2"`, `"V70CG1"`, `"M98CE"`
- Nucleic-acid atoms: `"A12C8"`, `"G5C1'"`
- Non-standard/unassigned: `"X99"` → resnum `-99`

This function is deliberately permissive: it **never throws**. Unrecognised
input simply yields a best guess (e.g. `ResidueLabel(0, '?', "")` for an empty
or number-less label). Parsing is metadata only and must never block fitting,
saving, or loading — users may use any labelling convention they like.
"""
function parse_label(label)::ResidueLabel
    isempty(label) && return ResidueLabel(0, '?', "")

    numbers = match(r"\d+", label)
    isnothing(numbers) && return ResidueLabel(0, first_letter_or_unknown(label), "")

    resnum = parse(Int, numbers.match)
    numstart = numbers.offset
    numend = numbers.offset + length(numbers.match) - 1

    letterpos = findfirst(isletter, label)
    code = isnothing(letterpos) ? '?' : label[letterpos]

    # Determine the atom token. If the first letter precedes the number
    # (code-first, e.g. "G5C1'"), the atom is whatever follows the number.
    # If the number comes first (e.g. reversed "7A"), the letter after the
    # number is the residue code and there is no atom.
    atom = if !isnothing(letterpos) && letterpos < numstart
        numend < lastindex(label) ? label[nextind(label, numend):end] : ""
    elseif !isnothing(letterpos)
        letterpos < lastindex(label) ? label[nextind(label, letterpos):end] : ""
    else
        ""
    end
    atom = uppercase(strip(atom))

    # Non-standard residue types (e.g. unassigned "X#") sort as negative
    if !isnothing(letterpos) && !in(code, STANDARD_RESIDUES)
        resnum = -resnum
    end

    return ResidueLabel(resnum, code, atom)
end

first_letter_or_unknown(label) = (p = findfirst(isletter, label);
                                  isnothing(p) ? '?' : label[p])

"""
    extract_residue_number(label)::Int

Residue number for `label` (see [`parse_label`](@ref)). Negative for
non-standard residue types; `0` if no digits are present.
"""
extract_residue_number(label) = parse_label(label).resnum

"""
    has_atom_labels(labels)::Bool

`true` if any label in the iterable carries an atom name (e.g. methyls such as
`"I13CD1"`), used to decide between scatter and bar summary plots.
"""
has_atom_labels(labels) = any(!isempty(parse_label(l).atom) for l in labels)