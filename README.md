### Fortran-based methods for computing the spatial heterogeneity of 2-D period patterns

 - `normalizedSpatialHet` Computes the spatial heterogeneity of a pattern, normalized between 0 and 1. This method using a iterative looping approach which can impose large computational overhead with patterns with large dimensions.
 - `monteCarloSpatialHet` Estimates the normalized spatial heterogeneity of a pattern by sampling a subset of the sub-rectangles itererated over in `normalizedSpatialHet`. This method should be preferred for large patterns, as it offers considerable speedup.
