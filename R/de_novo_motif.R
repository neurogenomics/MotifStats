#' Perform de novo motif discovery
#'
#' \code{de_novo_motif()} is a wrapper of several de novo motif discovery
#' methods, including Streme (meme suite), meme and rGADEM.
#'
#' Streme is preferable over Meme, except when you have a small number of
#' input sequences (< 50).
#'
#' @param peak_file Character string specifying the path to the narrowPeak file.
#' @param method De novo motif discovery method to use. The choices are
#' "Streme", "Meme" and "rGADEM". Please note that Streme requires the MEME
#' suite to be installed on your machine. You can find installation instructions
#' here.
#' @param control The default is to use the top 0.1 peaks as input (ranked by
#' q-value) and the remaining 0.9 peaks as the background. If you wish to
#' supply your own control sequences these must be passed as a DNAStringSet
#' object.
#' @param ... Additional parameters to be passed to Streme, Meme or rGADEM. See
#' the respective documentation for more details.
#'
#' @returns The output of either Streme, Meme or rGADEM.
#'
#' @examples
#' \dontrun{
#' peak_file <- system.file("extdata",
#'                          "rep1_peaks.narrowPeak",
#'                          package = "MotifStats"
#'                          )
#'
#' data("creb_motif", package = "MotifStats")
#'
#' de_novo_motif(peak_file = peak_file,
#'               pwm = creb_motif,
#'               plot = TRUE,
#'               min_score = 0.8
#'               )
#' }
#'
#' @export
de_novo_motif <- function(peak_file,
                          control = "weak_peaks",
                          method,
                          ...) {
  valid_methods <- c("Streme", "Meme", "rGADEM")
  if (!method %in% valid_methods) {
    stopper(
      "Unsupported method. The de novo motif discovery method must be one",
      "of Streme, Meme or rGADEM."
    )
  }
  peaks <- read_peak_file(peak_file)

  if (control == "weak_peaks") {
    top_index <- round(length(peaks) * 0.1)
    ordered_peaks <- peaks[order(peaks$qValue, decreasing = TRUE),]
    top_peaks <- ordered_peaks[seq_len(top_index),]
    bottom_peaks <- ordered_peaks[!seq_len(top_index),]

    input_peak_sequences <-
      BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                       top_peaks)
    control_peak_sequences <-
      BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                       bottom_peaks)

  } else if (control == "shuffle") {
    input_peak_sequences <-
      BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                       peaks)
    control_peak_sequences <- "shuffle"

  } else {
    if (!inherits(control, "DNAStringSet")) {
      stopper("Control sequences must be provided as a DNAStringSet object.")
    }

    input_peak_sequences <-
      BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                       peaks)
    control_peak_sequences <-
      BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
                       control)

  }

  if (method == "Streme") {
    motifs <- memes::runStreme(input = input_peak_sequences,
                               control = control_peak_sequences,
                               outdir = "./streme_results",
                               ...)
  } else if (method == "Meme") {
    motifs <- memes::runMeme(input = input_peak_sequences,
                             control = control_peak_sequences,
                             outdir = "./meme_results",
                             ...)
  } else if (method == "rGADEM") {
    motifs <- rGADEM::GADEM(Sequences = input_peak_sequences,
                            ...)
  }

  return(motifs)
}

# de_novo_motif(peak_file = peak_file)
#
# peak_file <- "./rep1_peaks.narrowPeak"
# peaks <- read_peak_file(peak_file)
# peak_sequences <-
#    BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
#                     peaks)
# #
# #
# # fa <- system.file("extdata/fasta_ex/fa1.fa", package = "memes")
# # dreme_out <- memes::runDreme(fa, "shuffle", evalue = 39, outdir = tempdir())
# # ?memes::runStreme()
# #
# # BiocManager::install("rGADEM")
# # library(rGADEM)
# # peak_file <- "./rep1_peaks.narrowPeak"
# # peaks <- read_peak_file(peak_file)
# # peak_sequences <-
# #   BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
# #                    peaks)
# # gad <- rGADEM::GADEM(Sequences = peak_sequences)
# # rGADEM::getPWM(gad
# # )
#
# ?memes::runStreme()


