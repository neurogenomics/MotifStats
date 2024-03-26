#' Calculate the distance between peak summits and motifs
#'
#' \code{summit_to_motif()} calculates the distance between each motif and its
#' nearest peak summit.
#'
#' This function is designed to work with MACS2/3 narrowPeak files.
#'
#' @import GenomicRanges
#' @importFrom memes runFimo
#'
#' @inheritParams peak_proportion
#'
#' @returns A list containing an expanded GRanges peak object with metadata
#' columns relating to motif positions along with a vector of summit-to-motif
#' distances for each valid peak.
#'
#' @examples
#' \dontrun{
#' pwm <-
#'   read_motif_file(motif_file = "./MA0018.3.meme", file_format = "meme")
#'
#' summit_to_motif(
#'   peak_input = peak_file,
#'   pwm = creb_motif,
#'   fp_rate = 5e-02,
#'   genome_build = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
#')
#'}
#' @export
summit_to_motif <- function(peak_input,
                            motif,
                            fp_rate = 5e-02,
                            genome_build) {
  peaks_and_seqs <- check_peak_input(peak_input = peak_input,
                                     genome_build = genome_build)
  peaks <- peaks_and_seqs[[1]] # 1st position of peaks_and_seqs vector
  peak_sequences <- peaks_and_seqs[[2]] # 2nd position of peaks_and_seqs vector

  fimo_threshold <- fp_rate / (2 * mean(GenomicRanges::width(peaks)))
  messager("The p-value threshold for motif scanning with FIMO is",
           fimo_threshold)
  fimo_df <- memes::runFimo(peak_sequences,
                            motifs = motif,
                            thresh = fimo_threshold,
                            text = FALSE)
  index_to_repeat <- base::match(as.vector(GenomicRanges::seqnames(fimo_df)),
                                 names(peaks))
  expanded_peaks <- peaks[index_to_repeat]

  # Calculate motif start and end points
  mcols(expanded_peaks)$motif_start <-
    GenomicRanges::start(expanded_peaks) + GenomicRanges::start(fimo_df)
  mcols(expanded_peaks)$motif_end <-
    GenomicRanges::start(expanded_peaks) + GenomicRanges::end(fimo_df)

  # Calculate motif centre and distance from centre to summit
  motif_centre <-
    (mcols(expanded_peaks)$motif_start + mcols(expanded_peaks)$motif_end) / 2
  distance_to_summit <- abs(mcols(expanded_peaks)$summit - motif_centre)

  mcols(expanded_peaks)$distance_to_summit <- distance_to_summit

  return(list(peak_set = expanded_peaks,
              distance_to_summit = distance_to_summit))
}


read_peak_file("CTCF_peaks.narrowPeak")
