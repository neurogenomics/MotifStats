#' Import a SEACR BED file as a GRanges object
#'
#' \code{read_peak_file()} reads a SEACR BED file as a GRanges object and adds
#' an additional metadata column representing the position of the peak summit.
#'
#' @importFrom utils read.table
#' @import IRanges
#' @import GenomicRanges
#'
#' @param file_path Path to a \code{BED} file.
#'
#' @examples
#' \dontrun{
#' peak_file <- system.file("extdata",
#'                          "ctcf_seacr.stringent.bed",
#'                          package = "MotifStats"
#'                          )
#' peak_obj <- read_peak_file_seacr(file_path = peak_file)
#' }
#'
#' @export
read_peak_file_seacr <- function(file_path){
  # Read BED file as table
  obj <- utils::read.table(file_path, header = TRUE)
  names(obj) <- c("chr",
                  "start",
                  "end",
                  "total_signal",
                  "max_signal",
                  "max_signal_region")

  # Create summit column
  separated_cols <- strsplit(as.character(obj$max_signal_region), "[:|-]")
  obj$max_chr <- sapply(separated_cols, function(x) x[1])
  obj$max_start <- sapply(separated_cols, function(x) x[2])
  obj$max_end <- sapply(separated_cols, function(x) x[3])
  obj$summit <- floor((as.numeric(obj$max_start) + as.numeric(obj$max_end)) / 2)

  # Create GRanges object
  gr_ranges <- IRanges::IRanges(start = obj$start, end = obj$end)

  gr_obj <- GenomicRanges::GRanges(
    seqnames = obj$chr,
    ranges = gr_ranges,
    strand = "*"
  )

  GenomicRanges::mcols(gr_obj)$summit <- obj$summit
  GenomicRanges::mcols(gr_obj)$name <- paste0("peak_", 1:nrow(obj))
  names(gr_obj) <- GenomicRanges::mcols(gr_obj)$name

  return(gr_obj)
}
