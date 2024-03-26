#' Read motif file and convert to Position Weight Matrix (PWM)
#'
#' \code{read_motif_file()} reads a motif file and converts to a PWM. The
#' function supports multiple motif formats, including "jaspar", "meme",
#' "transfac" and "uniprobe".
#'
#' @importFrom TFBSTools PWMatrix
#' @import universalmotif
#'
#' @param motif_file Path to the motif file.
#' @param motif_id ID of the motif (e.g. "MA1930.1").
#' @param file_format Character string specifying the format of the motif file.
#' The options are "homer", "jaspar", "meme", "transfac" and "uniprobe"
#'
#' @returns A \code{universalmotif} motif object.
#'
#' @export
read_motif_file <- function(motif_file,
                            motif_id = "Unknown",
                            file_format) {
  read_functions <- list(
    homer = universalmotif::read_homer,
    jaspar = universalmotif::read_jaspar,
    meme = universalmotif::read_meme,
    transfac = universalmotif::read_transfac,
    uniprobe = universalmotif::read_uniprobe
  )

  if (!file_format %in% names(read_functions)) {
    stop(
      "Unsupported file format. The motif file must be one of",
      "homer, jaspar, meme, transfac or uniprobe."
    )
  }

  # Dynamically select and call the appropriate read function
  read_function <- read_functions[[file_format]]
  motif <- read_function(motif_file)

  return(motif)
}

