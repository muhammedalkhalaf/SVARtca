# SVARtca 1.0.0

* Initial release (renamed from TCA to SVARtca to avoid CRAN name conflict).
* Implements Transmission Channel Analysis (Wegner, Lieb, Smeekes 2025).
* Three decomposition modes: overlapping, exhaustive 3-way, exhaustive 4-way.
* Systems form construction from VAR coefficient matrices.
* Integration with the `vars` package via `tca_from_var()`.
* Plotting with `plot_tca()` (ggplot2).
* Binary additivity validation with `tca_validate_additivity()`.
