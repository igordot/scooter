#' Get an example counts matrix.
#'
#' Get a small matrix of raw counts from the PBMC dataset.
#'
#' @return A matrix of raw counts.
#'
#' @examples
#' pbmc_mat <- get_test_counts_matrix()
#'
#' @importFrom utils read.table
#' @export
get_test_counts_matrix <- function() {
  pbmc_raw <- read.table(
    file = system.file("extdata", "pbmc.txt", package = "scooter"),
    as.is = TRUE,
    check.names = FALSE
  )
  as.matrix(pbmc_raw)
}


