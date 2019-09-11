# #' @import RColorBrewer
# #' @import ggsci
# NULL

# environment that holds various global variables and settings
# scooter_global <- new.env(parent = emptyenv())
#
# scooter_global$test_var = "hello"
#
# assign(".sessionId", "xyz123", envir = parent.env(environment()))
#

#####

# old zzz.R
# .onLoad <- function(libname, pkgname) {
#   op <- options()
#   op.scooter <- list(
#     scooter.colors_samples =
#       c(RColorBrewer::brewer.pal(5, "Set1"), RColorBrewer::brewer.pal(8, "Dark2"), ggsci::pal_igv("default")(51)),
#     scooter.colors_clusters =
#       c(ggsci::pal_d3("category10")(10), ggsci::pal_d3("category20b")(20), ggsci::pal_igv("default")(51))
#   )
#   toset <- !(names(op.scooter) %in% names(op))
#   if (any(toset)) options(op.scooter[toset])
#
#   invisible()
# }

# more about .onLoad():
# https://stackoverflow.com/questions/28246952/global-variable-in-a-package-which-approach-is-more-recommended
# https://github.com/tidyverse/ggplot2/blob/1c09bae2aa24320bcb4891c664cccc16efc86a8a/R/zzz.r

# should try panglaodb + xcell
