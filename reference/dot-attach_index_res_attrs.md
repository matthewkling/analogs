# Internal helper: attach index resolution parameters as result attributes

Copies the lattice resolution knobs and realized bin layout from the
`analog_index` object onto a query result, so they surface via
[`metadata()`](https://matthewkling.github.io/analogs/reference/metadata.md).
These are properties of the index (set at build time), reused here
rather than recomputed. Missing fields (e.g. on legacy index objects)
are simply skipped.

## Usage

``` r
.attach_index_res_attrs(out, index)
```
