#' Calculate the distance between peak summits and motifs
#'
#' \code{summit_to_motif()} calculates the distance between each motif and its
#' nearest summit. A peak summit represents the exact nucleotide where signal
#' enrichment is highest.
#'
#' @import GenomicRanges
#' @import rtracklayer
#' @import TFBSTools
#'
#' @inheritParams peak_proportion
#' @param control_input Method for generating control sequences.
#' \code{"shuffle"} will shuffle the input sequences, whereas \code{"adjacent"}
#' selects real sequences immediately adjacent to each peak. The default is
#' \code{"shuffle"}. This argument is ignored if \code{optimal_min_score} is set
#' to \code{FALSE}.
#'
#' @returns GRanges peak object with metadata columns relating to motif
#' positions.
#'
#' @examples
#' \dontrun{
#' pwm <-
#'   read_motif_file(motif_file = "./MA0018.3.meme", file_format = "meme")
#'
#' summit_to_motif(peak_file = "read1_01_peaks.narrowPeak",
#'                 pwm = pwm,
#'                 min_score = NULL,
#'                 optimal_min_score = TRUE,
#'                 seed = 123
#'                 )
#'                 }
#' @export
summit_to_motif <- function(peak_input,
                            control_input = "shuffle",
                            pwm,
                            genome_build,
                            shuffle_k = 2,
                            min_score = 0.8,
                            optimal_min_score = FALSE,
                            num_samples = 4,
                            sample_size = 250,
                            num_threads = 1L,
                            deoptim_lower = 0.5,
                            deoptim_upper = 1,
                            deoptim_np = 20,
                            deoptim_itermax = 25,
                            deoptim_strategy = 6,
                            seed = NULL) {
  if(!is.null(seed)) set.seed(seed)

  peaks_and_seqs <- check_peak_input(peak_input = peak_input,
                                     genome_build = genome_build)
  peaks <- peaks_and_seqs[[1]] # 1st position of peaks_and_seqs vector
  peak_sequences <- peaks_and_seqs[[2]] # 2nd position of peaks_and_seqs vector

  if (optimal_min_score) {
    min_score <- optimal_min_score(peaks = peaks,
                                   control_input = control_input,
                                   pwm = pwm,
                                   genome_build = genome_build,
                                   num_samples = num_samples,
                                   sample_size = sample_size,
                                   num_threads = num_threads,
                                   deoptim_lower = deoptim_lower,
                                   deoptim_upper = deoptim_upper,
                                   deoptim_np = deoptim_np,
                                   deoptim_itermax = deoptim_itermax,
                                   deoptim_strategy = deoptim_strategy,
                                   seed = seed)
    messager("Optimal min_score threshold found:", min_score,
             "\n\nThis threshold will now be applied to the full dataset")

  } else if (!is.null(min_score)) {
    min_score <- min_score
  } else {
    stopper("Either optimal_min_score is set to TRUE or min_score must be",
            "supplied.")
  }

  # TO-DO:
  # if we want to use fimo for motif finding?
  fimo_df <- memes::runFimo(peak_sequences, motifs = pwm, thresh = 1e-4)
  index_to_repeat <- base::match(as.vector(GenomicRanges::seqnames(fimo_df)), names(peaks))
  expanded_peaks <- peaks[index_to_repeat]
  print(expanded_peaks)
  mcols(expanded_peaks)$motif_start <-
    GenomicRanges::start(expanded_peaks) + GenomicRanges::start(fimo_df)
  mcols(expanded_peaks)$motif_end <-
    GenomicRanges::start(expanded_peaks) + GenomicRanges::end(fimo_df)
  motif_center <-
    (mcols(expanded_peaks)$motif_start + mcols(expanded_peaks)$motif_end) / 2
  distance_to_summit <- abs(mcols(expanded_peaks)$summit - motif_center)
  mcols(expanded_peaks)$distance_to_summit <- distance_to_summit


  # hits <- get_motif_hits(
  #   peak_sequences = peak_sequences,
  #   pwm = pwm,
  #   min_score = min_score,
  #   genome_build = genome_build
  # )
  # # Write the motif hits to GFF3
  # motif_coord_df <- TFBSTools::writeGFF3(hits)
  #
  # # Find peaks with multiple motifs
  # index_to_repeat <- base::match(motif_coord_df$seqname, names(peaks))
  #
  # # Expand the peak object to repeat the peaks that have >1 motif
  # expanded_peaks <- peaks[index_to_repeat]
  #
  # # Add the motif start and end data for each peak
  # mcols(expanded_peaks)$motif_start <-
  #   GenomicRanges::start(expanded_peaks) + motif_coord_df$start
  # mcols(expanded_peaks)$motif_end <-
  #   GenomicRanges::start(expanded_peaks) + motif_coord_df$end
  #
  # # Calculate motif centres and distance to summit
  # motif_center <-
  #   (mcols(expanded_peaks)$motif_start + mcols(expanded_peaks)$motif_end) / 2
  # distance_to_summit <- abs(mcols(expanded_peaks)$summit - motif_center)
  #
  # # Adding distance to metadata
  # mcols(expanded_peaks)$distance_to_summit <- distance_to_summit

  return(
    list(
      peak_set = expanded_peaks,
      distance_to_summit = distance_to_summit
    )
  )
}

# CTCF motif
# pwm <- read_motif_file("MA1930.2.jaspar", file_format = "jaspar")
# one <- summit_to_motif("./both_01_peaks.narrowPeak",
#                         pwm = pwm)
# two <- summit_to_motif("./read1_01_peaks.narrowPeak",
#                         pwm = pwm)
# par(mfrow = c(1,2))
# hist(one$distance_to_summit, breaks = 40, xlim = c(0, 1000))
# hist(two$distance_to_summit, breaks = 40, xlim = c(0, 1000))
#
# median(one$distance_to_summit)
# median(two$distance_to_summit)
