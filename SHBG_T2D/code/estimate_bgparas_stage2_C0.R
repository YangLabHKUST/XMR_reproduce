# Two-stage LDSC with zero cross-population intercepts for the XMR analysis.
# The public XMR interface also supports unconstrained cross-population intercepts.
if (!"cross_intercept_zero" %in% names(formals(XMR::estimate_bgparas))) {
  stop("Install the XMR version specified in environment/README.md.")
}
estimate_bgparas_stage2_C0 <- function(..., cross_intercept_zero = TRUE) {
  XMR::estimate_bgparas(..., cross_intercept_zero = cross_intercept_zero)
}
