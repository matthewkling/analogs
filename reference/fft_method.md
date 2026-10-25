# FFT (binned-convolution) execution method

Configures the binned-convolution backend for
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
and its wrappers. Pass the result as `method = fft_method(...)`, or use
`method = "fft"` for the defaults. The FFT method computes the same
kernel-weighted statistics as the exact method, approximately, without
visiting every focal-analog pair, which can make it much faster for
large problems.

## Usage

``` r
fft_method(bin = 0.25, check = 0, max_cells = 1e+08)
```

## Arguments

- bin:

  Grid spacing, as a fraction of each kernel's scale: `theta` for
  Gaussian and inverse kernels, `max` for uniform kernels (see Details).
  Default `0.25`.

- check:

  Number of randomly chosen focals to recompute with the exact engine
  for an accuracy report, attached to the output as the `fft_check`
  attribute (see
  [`metadata()`](https://matthewkling.github.io/analogs/reference/metadata.md)).
  Default `0` (no check).

- max_cells:

  Upper limit on the number of grid nodes, as a guard against memory
  exhaustion (each node of each grid takes 8 bytes). Default `1e8`.

## Value

An object of class `"analog_fft_method"`.

## Details

### How it works

Every statistic the package computes for a focal is built from
kernel-weighted sums over the pool, such as the sum of weights or the
sum of weights times `y`. The exact method evaluates these sums focal by
focal. The FFT method (a multivariate version of binned kernel
estimation; Wand 1994) evaluates them for all focals at once using an
internal grid. This grid is not the input raster: it is a coarse regular
lattice spanning the analog space, with one axis per environmental
variable when the environmental kernel weights or applies a limit, and
two geographic axes when the geographic kernel does. Nodes are spaced at
a fraction (`bin`) of each kernel's bandwidth. For example, with a
geographic `theta` of 400 km over a 4000 km domain and `bin = 0.25`, the
geographic axes have nodes every 100 km, about 40 per axis. Each sum is
then computed in three steps:

1.  **Bin.** Each pool point's contribution is split among the grid
    nodes surrounding it, in proportion to its closeness to each (linear
    binning). The grid then holds how much pool "mass" lies near each
    node.

2.  **Convolve.** Each node is assigned the kernel-weighted sum of the
    mass at every node: the value a hypothetical focal sitting exactly
    on that node would get. Done directly this would cost (number of
    nodes)^2 operations; the fast Fourier transform does it in about
    (number of nodes) x log(number of nodes), by transforming the grid
    and the kernel, multiplying them, and transforming back. The grid is
    padded by one kernel width of empty nodes along each axis, so that
    mass near one edge does not wrap around to the opposite edge.

3.  **Read off.** Each focal's value is interpolated from the grid nodes
    surrounding its exact location and environment (multilinear
    interpolation).

Because the convolved surface is a sum of kernel-shaped bumps, it
changes little between nodes spaced well within the kernel's width, so
the coarse grid loses little accuracy. The output keeps the resolution
of the focal data: two neighboring raster cells with different climates
read off different parts of the environmental axes and can get quite
different results. Statistics that need several sums (for example, a
regression needs weighted cross-products of the covariates and `y`)
repeat the bin and convolve steps once per sum. For leave-one-out
(`exclude_self = TRUE`), each focal's own binned contribution is
subtracted after read-off.

The grid covers only the region where focals and pool points are within
each other's reach, so focals confined to part of the pool's extent cost
less. A family's nodes are spaced at `bin * theta` (`bin * max` for a
uniform kernel, or `bin * min` if a geographic `min` is smaller than
`theta`), shrunk slightly so that nodes fall exactly on the ends of the
range.

### Accuracy and cost

`bin` controls the trade-off between accuracy and cost. For the Gaussian
kernel, binning error shrinks with the square of `bin`; for the inverse
kernel (which has a sharp peak at zero distance) and for hard edges
(uniform kernels, `max`/`min` cutoffs where the weight is still
appreciable, `count`) it shrinks roughly in proportion to `bin`. Hard
limits are anti-aliased on the grid (each node counts in proportion to
the fraction of its cell inside the limit), which reduces but does not
remove this edge error. Halving `bin` multiplies the number of grid
nodes by `2^d` for `d` active axes. The default `bin = 0.25` typically
gives errors well under 1%, but error at a given `bin` is larger in
sparse neighborhoods and for ratio statistics such as regression
coefficients; use `check` to measure it.

### When to use it

The FFT method's cost is set by the size of its grid rather than by the
number of focal-analog pairs. Each axis has roughly
`range / (bin * theta)` nodes, where `range` is the data's extent along
that axis, and grid size is the product over axes. So what matters is
how many axes are active and how narrow the kernels are relative to the
data's range:

- **Two active axes** (e.g. two environmental variables, or a geographic
  kernel alone): the grid stays small even for narrow kernels, and the
  FFT method is typically one to three orders of magnitude faster than
  the exact method for large numbers of focals.

- **Three or more active axes**: grid size grows with the product of the
  range-to-bandwidth ratios across axes. The FFT method pays off when
  the kernels are wide relative to the data's range. When a kernel is
  narrow, analog neighborhoods are small, the exact method's search
  index is efficient, and the FFT grid may exceed `max_cells`.

Speedups also shrink with more threads, since the exact method
parallelizes over focals more completely than the FFT method. When in
doubt, time both on a subset of focals. If `max_cells` is exceeded, use
a coarser `bin` or the exact method; for raster inputs,
[`tiled_analog_search()`](https://matthewkling.github.io/analogs/reference/tiled_analog_search.md)
accepts `method = "fft"` and builds a separate, smaller grid for each
tile and its halo (spanning only the environments found there), which is
the most effective way to cut peak memory when geographic axes are
active.

### Requirements

Geographic kernels require projected coordinates (ideally an equal-area
projection), since on a longitude/latitude grid a fixed-distance kernel
changes shape with latitude. The FFT method requires `select = "all"`
and a raw (non-index) `pool`, and does not support `stat = "none"`,
`x_cov`, `downsample < 1`, or design-based standard errors for
`"regression"`. Gaussian kernels without a `max` are truncated at 7
`theta`.

## References

Wand, M. P. (1994). Fast computation of multivariate kernel estimators.
*Journal of Computational and Graphical Statistics*, 3(4), 433-445.

## See also

[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)

## Examples

``` r
fft_method()
#> <analog_fft_method> bin = 0.25, check = 0, max_cells = 1e+08
fft_method(bin = 0.1, check = 500)
#> <analog_fft_method> bin = 0.1, check = 500, max_cells = 1e+08
```
