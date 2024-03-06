#' Calculate the proportion of peaks with a motif
#'
#' \code{peak_proportion()} calculates the proportion of peaks containing a
#' motif. To generate a plot the peak file must contain a column that can be
#' used for ranking (e.g. p-values, q-values or signal values).
#'
#' @import rtracklayer
#' @import GenomicRanges
#' @import BSgenome
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import TFBSTools
#' @import ggplot2
#'
#' @param peak_file Character string specifying the path to the narrowPeak file.
#' @param pwm An object of class PWMatrix.
#' @param plot Logical specifying whether or not to create a plot. Note that the
#' peak file must contain a column suitable for ranking.
#' @param min_score Numeric specifying the minimum score a sequence must pass
#' to be considered a motif. The default is 0.8.
#' @param out_dir Path to the output directory where the filtered peak file
#' will be written.
#'
#' @examples
#' \dontrun{
#' pwm <-
#'   read_motif_file(motif_file = "./MA0018.3.meme", file_format = "meme")
#'
#' peak_proportion(peak_file = "SK032_1_SK035_1_SK035_11_R1_peaks.narrowPeak",
#'                 pwm = pwm,
#'                 plot = TRUE,
#'                 min_score = 0.8
#'                 )
#'                 }
#'
#' @export

# to-do: add a genome build option where a user specifies something like
# genome_build = BSgenome.Hsapiens.UCSC.hg38
peak_proportion <- function(peak_file,
                            pwm,
                            plot = FALSE,
                            min_score = 0.8
                            ) { # add outdir option for writing "true" peaks
  peaks <- read_peak_file(peak_file)
  hits <- get_motif_hits(peaks = peaks,
                         pwm = pwm,
                         min_score = min_score
                         )
  sig <- TFBSTools::relScore(hits)

  # 1 motif per peak
  onehit_peak_names <- names(sig)[sapply(sig, length) > 0]
  onehit_peaks <- peaks[names(peaks) %in% onehit_peak_names, ]

  # > 1 motif per peak
  mt1_peak_names <- names(sig)[sapply(sig, length) > 1]
  mt1_peaks <- peaks[names(peaks) %in% mt1_peak_names, ]

  if (plot) {
    ranked_peaks <- peaks[order(-GenomicRanges::mcols(peaks)$qValue)]
    motif_hits_df <-  data.frame(peak_name = names(ranked_peaks),
                                 peak_order = seq_along(ranked_peaks))

    motif_hits_df$contains_motif = motif_hits_df$peak_name %in% onehit_peak_names
    motif_hits_df = motif_hits_df[order(-motif_hits_df$peak_order),]

    motif_hits_df$perc_peaks <-
      with(motif_hits_df,
           (cumsum(contains_motif) / length(peak_order)))

    motif_hits_df$perc_peaks <- round(motif_hits_df$perc_peaks, 2)

    plot <- ggplot2::ggplot(motif_hits_df,
                            aes(x = peak_order, y = perc_peaks)) +
      geom_line(size = 1, color = "blue") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_blank()
      ) +
      xlab("Peak rank") +
      ylab("Cumulative proportion")

    print(plot)
  }

  return(
    list(
      one_hit_peak_prop = onehit_peaks,
      proportion_1_motif = length(onehit_peaks) / length(peaks),
      mt1_hit_peak_prop = mt1_peaks,
      proportion_mt1_motifs = length(mt1_peaks) / length(peaks)
    )
  )
}
