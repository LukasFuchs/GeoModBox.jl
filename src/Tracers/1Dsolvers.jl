"""
    Itp1D_Centers2Markers!(Tm, xm, Tc, xc, Δx)

Interpolate a one-dimensional cell-centered field to marker positions using
linear interpolation.

For each marker coordinate in `xm`, the function identifies the two
neighboring cell centers in `xc` and interpolates the corresponding values
of `Tc`. The interpolated marker values are written in place to `Tm`.

# Arguments

- `Tm`: Array receiving the interpolated marker values.
- `xm`: Marker coordinates.
- `Tc`: Field values defined at the cell centers.
- `xc`: Cell-center coordinates.
- `Δx`: Uniform spacing between adjacent cell centers.

# Notes

The arrays `Tm` and `xm` must have the same length, while `Tc` and `xc`
must contain the same number of entries. At least two cell centers are
required.

Marker positions outside the interval spanned by `xc` are linearly
extrapolated using the first or last pair of cell centers.

# Example

```julia
xc = [0.5, 1.5, 2.5, 3.5]
Tc = [10.0, 20.0, 30.0, 40.0]

xm = [1.0, 2.0, 3.0]
Tm = similar(xm)

Itp1D_Centers2Markers!(Tm, xm, Tc, xc, 1.0)
```
"""
function Itp1D_Centers2Markers!(Tm, xm, Tc, xc, Δx, xmin)
    length(Tm) == length(xm) ||
        throw(DimensionMismatch(
            "`Tm` and `xm` must have the same length."
        ))

    length(Tc) == length(xc) ||
        throw(DimensionMismatch(
            "`Tc` and `xc` must have the same length."
        ))

    length(Tc) >= 2 ||
        throw(ArgumentError(
            "At least two cell-center values are required."
        ))

    Δx > 0 ||
        throw(ArgumentError(
            "`Δx` must be positive."
        ))

    @inbounds for k in eachindex(xm, Tm)

        # Index of the cell center immediately west of the marker.
        iW = floor(Int, (xm[k] - xc[1]) / Δx) + 1
        iW = clamp(iW, 1, length(Tc) - 1)

        iE = iW + 1

        # Linear interpolation coordinate.
        ξ = (xm[k] - xc[iW]) / Δx

        Tm[k] = (1.0 - ξ) * Tc[iW] + ξ * Tc[iE]
    end

    return Tm
end

"""
    Itp1D_Markers2Centers!(Tc, xc, Tm, xm, Δx)

Interpolate values from one-dimensional markers to cell centers using
linear weighted averaging.

Each marker contributes to the two adjacent cell centers according to its
relative position between them. Marker contributions and interpolation
weights are accumulated separately, after which the value at each cell
center is obtained by dividing the weighted sum by the total weight.

# Arguments

- `Tc`: Cell-centered field receiving the interpolated marker values.
- `xc`: Coordinates of the cell centers.
- `Tm`: Values carried by the markers.
- `xm`: Marker coordinates.
- `Δx`: Uniform spacing between adjacent cell centers.

# Notes

The arrays `Tm` and `xm` must have the same length, and `Tc` and `xc`
must contain the same number of entries. At least two cell centers are
required.

Markers located outside the interval spanned by `xc` are assigned to the
nearest cell center. Cell centers receiving no marker contributions retain
their previous value in `Tc`.

# Example

```julia
xc = [0.5, 1.5, 2.5, 3.5]
xm = [1.0, 2.0, 3.0]
Tm = [10.0, 20.0, 30.0]

Tc = zeros(length(xc))
```
"""
function Itp1D_Markers2Centers!(Tc, xc, Tm, xm, Δx, xmin)
    length(Tm) == length(xm) ||
        throw(DimensionMismatch(
            "`Tm` and `xm` must have the same length."
        ))

    length(Tc) == length(xc) ||
        throw(DimensionMismatch(
            "`Tc` and `xc` must have the same length."
        ))

    length(Tc) >= 2 ||
        throw(ArgumentError(
            "At least two cell centers are required."
        ))

    Δx > 0 ||
        throw(ArgumentError(
            "`Δx` must be positive."
        ))

    # Preserve the old field so that empty cells can retain their value.
    Tc_old = copy(Tc)

    fill!(Tc, 0.0)
    Wc = zeros(eltype(Tc), size(Tc))

    @inbounds for k in eachindex(xm, Tm)

        if xm[k] <= xc[1]
            # Constant extrapolation at the western boundary.
            Tc[1] += Tm[k]
            Wc[1] += 1.0

        elseif xm[k] >= xc[end]
            # Constant extrapolation at the eastern boundary.
            Tc[end] += Tm[k]
            Wc[end] += 1.0

        else
            # Cell center immediately west of the marker.
            iW = floor(Int, (xm[k] - xc[1]) / Δx) + 1
            iW = clamp(iW, 1, length(Tc) - 1)

            iE = iW + 1

            # Relative marker position between the neighboring centers.
            ξ  = (xm[k] - xc[iW]) / Δx
            wW = 1.0 - ξ
            wE = ξ

            Tc[iW] += wW * Tm[k]
            Tc[iE] += wE * Tm[k]

            Wc[iW] += wW
            Wc[iE] += wE
        end
    end

    # Normalize only centers that received marker contributions.
    @inbounds for i in eachindex(Tc)
        if Wc[i] > 0.0
            Tc[i] /= Wc[i]
        else
            Tc[i] = Tc_old[i]
        end
    end

    return Tc
end