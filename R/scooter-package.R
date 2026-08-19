#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import dplyr
#' @import ggplot2
#' @import Seurat
#' @import tibble
#' @importFrom glue glue
## usethis namespace: end
NULL

.onLoad <- function(libname, pkgname) {
  # run_umap() always uses uwot (never umap.method = "umap-learn"), so Seurat's
  # one-time "default method changed from Python UMAP to uwot" warning is never
  # actionable here
  options(Seurat.warn.umap.uwot = FALSE)
}
