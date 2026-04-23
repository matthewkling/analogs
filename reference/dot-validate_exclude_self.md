# Validate exclude_self parameter and its compatibility with other args

Called from analog_search (and anywhere else that surfaces
exclude_self). Enforces identical(x, pool), disallows pre-built indices,
and disallows downsampling. Also disallows progress (chunking is
incompatible with the simple j==i self-exclusion check).

## Usage

``` r
.validate_exclude_self(exclude_self, x, pool, downsample, progress)
```
