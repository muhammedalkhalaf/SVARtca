# SVARtca 1.0.3

No numerical results change in this release. `tca_systems_form()`,
`tca_analyze()`, `tca_decompose_binary()`, `tca_from_var()` and
`plot_tca()` return exactly the same values as in version 1.0.2. The
changes concern the additivity diagnostic, documentation and attribution.

* Author attribution corrected. The methodology paper is by Enrico Wegner,
  Lenard Lieb, Stephan Smeekes and Ines Wilms (arXiv:2405.18987,
  doi:10.48550/arXiv.2405.18987). Earlier versions cited "Emanuel Wegner"
  and omitted Ines Wilms. Corrected in DESCRIPTION (including the
  contributor entries: Enrico Wegner is credited as author of the MATLAB
  TCA toolbox the package is ported from and as co-author of the
  methodology; Lieb, Smeekes and Wilms as authors of the methodology), in
  all help pages, the vignette and the paper drafts; the CITATION version
  note was updated.
* `tca_validate_additivity()` now performs a real check of Theorem 2(ii)
  of the paper. Previously it computed `through = total - not_through`
  and then tested `total - (through + not_through)`, which is zero by
  construction and could not fail. The through effect is now recomputed
  independently with AND conditions (a first-passage partition of the
  paths through the variable) and `through + not_through` is compared
  with `total` at every response variable other than the conditioning
  variable. New argument `pair` additionally checks the inclusion-exclusion
  identity `through(v1 and v2) = through(v1) + through(v2) - through(v1 or v2)`
  used by the 4-way mode, with the left-hand side computed directly by
  AND conditions on both variables. New argument `tol` (default 1e-10).
  The function also reports whether `B` is strictly lower triangular.
  The return value is still an invisible logical, now with attributes
  `residuals`, `pair_residual`, `max_residual` and `lower_triangular`.
  On a valid systems form the residual is of the order of 1e-16; on a
  corrupted `B` (for example a backward edge above the diagonal) the
  residual is nonzero and the check fails.
* `print.tca_result()` documentation corrected: the default `target` is
  variable 1 (the first variable in the original ordering), as the code
  has always done, not the first intermediate.
* `tca_systems_form()` documentation: the `Psis` argument is now
  documented, with a note that the matrices must be the reduced-form MA
  coefficients on the reduced-form errors `u_t`, not the structural MA
  matrices on the structural shocks (`Psi_j = Phi0 Psi^s_j Phi0^{-1}`).
  Passing structural matrices gives wrong impulse responses at horizons
  1 and above.
* `tca_decompose_binary()` documentation now states that
  `total = through + not_through` holds by construction and points to
  `tca_validate_additivity()` for the independent check.
* Tests: the additivity tests can now fail, and the 3-way and 4-way
  tests compare every channel with a brute-force enumeration of the paths
  in the systems form graph (K = 4, h = 2) instead of checking a sum that
  holds by construction; new tests check that a
  corrupted `B` is reported, that the inclusion-exclusion identity holds,
  that `tca_decompose_binary()` matches an independent AND-based through
  effect and a hand-computed IRF recursion, and that the worked example
  of Section 2.2 of the paper, evaluated at the package's own parameter
  values alpha1 = -0.4, alpha2 = 0.3, alpha3 = 0.5, alpha4 = 1.5, gives
  TE = 0.6884057971, IE = 0.2811094453 and DE = 0.4072963518, equal to the
  paper's symbolic formulas evaluated at the same values.

# SVARtca 1.0.0

* Initial release (renamed from TCA to SVARtca to avoid CRAN name conflict).
* Implements Transmission Channel Analysis (Wegner, Lieb, Smeekes and Wilms
  2025) [citation corrected in 1.0.3].
* Three decomposition modes: overlapping, exhaustive 3-way, exhaustive 4-way.
* Systems form construction from VAR coefficient matrices.
* Integration with the `vars` package via `tca_from_var()`.
* Plotting with `plot_tca()` (ggplot2).
* Binary additivity validation with `tca_validate_additivity()`.
