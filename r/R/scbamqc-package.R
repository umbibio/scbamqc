#' @keywords internal
#' @section Workflow:
#' 1. [add_bam_qc()]: add the table from `scbamqc bam`.
#' 2. [basic_qc()]: remove empty barcodes and rare genes.
#' 3. [qc_score()]: score each barcode.
#' 4. [select_cells()]: select candidate cells.
#'
#' @importFrom mclust mclustBIC
"_PACKAGE"
