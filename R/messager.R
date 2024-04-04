#' Print messages
#'
#' Conditionally print messages.
#'
#' Allows developers to easily control verbosity of functions, and meet
#' Bioconductor requirements that dictate the message must first be stored in a
#' variable before passing to \link[base]{message}.
#'
#' @param v Whether to print messages or not.
#'
#' @return NULL
#' @keywords internal
messager <- function(..., v = TRUE) {
  if(v){
    msg <- paste(...)
    try({message(msg)})
  }
}
