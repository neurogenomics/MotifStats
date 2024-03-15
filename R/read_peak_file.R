#' Import a narrowPeak file as a GRanges object
#'
#' \code{read_peak_file()} reads a narrowPeak file as a GRanges object and adds
#' an additional metadata column representing the actual position of the peak
#' summit.
#'
#' @importFrom rtracklayer import
#' @import GenomicRanges
#'
#' @param file_path Path to a narrowPeak file
#'
#' @examples
#' \dontrun{
#' peak_file <- system.file("extdata",
#'                          "rep1_peaks.narrowPeak",
#'                          package = "MotifStats"
#'                          )
#' peak_obj <- read_peak_file(file_path = peak_file)
#' }
#'
#' @export

read_peak_file <- function(file_path){
  obj <- rtracklayer::import(file_path, format = "narrowPeak")
  names(obj) <- GenomicRanges::mcols(obj)$name

  # Add metadata summit column
  GenomicRanges::mcols(obj)$summit <-
    GenomicRanges::start(obj) + GenomicRanges::mcols(obj)$peak

  return(obj)
}

