#' Find the optimal min_score threshold for \code{searchSeq()}
#'
#' \code{optimal_min_score()} finds the optimal min_score threshold for
#' TFBSTools::searchSeq by maximising the difference between the proportion
#' of motifs found in a subset of the peak data and a set of background control
#' sequences.
#'
#' @importFrom BSgenome getSeq
#' @import DEoptim
#' @import parallel
#'
#' @inheritParams peak_proportion
#'
#' @keywords internal

# Calculate the proportion of peaks (or ctrl sequences) with a motif
calculate_score <- function(peaks, pwm, param, genome_build) {
  hits <- get_motif_hits(
    peak_sequences = peaks,
    pwm = pwm,
    min_score = param,
    genome_build = genome_build
  )

  hit_peak_names <- names(hits)[sapply(hits, length) > 0]
  hit_peaks <- peaks[names(peaks) %in% hit_peak_names]

  return(length(hit_peaks) / length(peaks))
}

optimise_sample <- function(sample_data) {
  if (!is.null(sample_data$seed)) set.seed(sample_data$seed)

  loss_function <- function(param) {
    peak_score <- calculate_score(peaks = sample_data$peak_sequences,
                                  pwm = sample_data$pwm,
                                  param = param,
                                  genome_build = sample_data$genome_build)
    ctrl_score <- calculate_score(peaks = sample_data$ctrl_sequences,
                                  pwm = sample_data$pwm,
                                  param = param,
                                  genome_build = sample_data$genome_build)
    return(ctrl_score - peak_score)
  }

  control <- DEoptim::DEoptim.control(
    NP = sample_data$deoptim_np,
    itermax = sample_data$deoptim_itermax,
    strategy = sample_data$deoptim_strategy
  )
  result <- DEoptim::DEoptim(
    loss_function,
    lower = sample_data$deoptim_lower,
    upper = sample_data$deoptim_upper,
    control = control
  )

  return(result$optim$bestmem)
}

optimal_min_score <- function(peaks,
                              control_input = "shuffle",
                              pwm,
                              genome_build,
                              shuffle_k = NULL,
                              num_samples = 4,
                              sample_size = 250,
                              num_threads = 4L,
                              deoptim_lower = 0.5,
                              deoptim_upper = 1,
                              deoptim_np = 20,
                              deoptim_itermax = 25,
                              deoptim_strategy = 6,
                              seed = NULL) {
  total_peaks_needed <- num_samples * sample_size
  if (length(peaks) < total_peaks_needed) {
    messager(
      "Not enough peaks for the desired number of samples without repetition.",
      "Adjusting sample size."
    )
    sample_size <- floor(length(peaks) / num_samples)
  }

  shuffled_peaks <- sample(peaks, total_peaks_needed, replace = FALSE)
  samples <- split(shuffled_peaks, rep(1:num_samples, each = sample_size))

  samples_data <- lapply(samples, function(sample_peaks) {

    peak_sequences <- BSgenome::getSeq(genome_build, sample_peaks)
    if(control_input == "shuffle"){
      ctrl_sequences <- shuffle_sequences(peak_sequences, k = 2)
    } else if(control_input == "adjacent"){
      ctrl <- adjacent_sequences(sample_peaks, genome_build = genome_build)
      ctrl_sequences <- BSgenome::getSeq(genome_build, ctrl)
    }

    list(
      sample_peaks = sample_peaks,
      peak_sequences = peak_sequences,
      ctrl_sequences = ctrl_sequences,
      pwm = pwm,
      genome_build = genome_build,
      deoptim_lower = deoptim_lower,
      deoptim_upper = deoptim_upper,
      deoptim_np = deoptim_np,
      deoptim_itermax = deoptim_itermax,
      deoptim_strategy = deoptim_strategy,
      seed = seed
    )
  })

  if (num_threads > 1L) {
    cl <- parallel::makeCluster(min(num_threads, parallel::detectCores()))
    messager("Initialised cluster object:", cl)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl,
                            varlist = c("calculate_score",
                                        "adjacent_sequences",
                                        "optimise_sample"))
    message("Multithreaded optimisation starting...")
    optimal_params <- parallel::parLapply(cl, samples_data, optimise_sample)
  } else {
    message("Single-threaded optimisation starting...")
    optimal_params <- lapply(samples_data, optimise_sample)
  }

  messager(optimal_params)
  optimal_param <- mean(unlist(optimal_params))

  return(optimal_param)
}
