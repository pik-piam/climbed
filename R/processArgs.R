#' Process and Standardize Input Arguments
#'
#' The function processes and standardizes the input arguments \code{tLim}, \code{std},
#' and \code{ssp}, ensuring they are in the required format for further processing in
#' the \code{start} function or other applications.
#'
#' @param tLim A numeric value, numeric vector, or named list/vector with names
#'  \code{"HDD"} and/or \code{"CDD"}.
#'   \itemize{
#'     \item If a single numeric value is provided, it is used for both \code{"HDD"}
#'     and \code{"CDD"} without creating a sequence.
#'     \item If a named list/vector is provided, each element is kept as is,
#'     without generating sequences.
#'   }
#'
#' @param std A numeric value or a named numeric vector with names \code{"tLim"} and \code{"tAmb"}.
#'   \itemize{
#'     \item If a single numeric value is provided, it is applied to both \code{"tLim"} and \code{"tAmb"}.
#'     \item If a named numeric vector is provided, it must have names \code{"tLim"} and \code{"tAmb"}.
#'   }
#'
#' @param ssp A character string or vector specifying the SSP (Shared Socioeconomic Pathways) scenarios.
#'
#' @return A list containing standardized arguments:
#'   \describe{
#'     \item{\code{tLim}}{A named list with \code{"HDD"} and \code{"CDD"} components, each being
#'     numeric values or sequences.}
#'     \item{\code{std}}{A named numeric vector with names \code{"tLim"} and \code{"tAmb"}.}
#'     \item{\code{ssp}}{A character vector of SSP scenarios.}
#'   }
#'
#' @author Hagen Tockhorn
#'
#' @export

processArgs <- function(tLim, std, ssp) {

  # Process tLim ---------------------------------------------------------------

  # Check and process tLim
  if (is.list(tLim) || is.vector(tLim)) {
    # Convert to list for uniform processing
    tLim <- as.list(tLim)

    # Check that names are "HDD" and "CDD"
    validNames <- c("HDD", "CDD")
    if (!all(names(tLim) %in% validNames)) {
      stop("tLim must have names 'HDD' and/or 'CDD'.")
    }

    # Process each element - keep values as is
    for (name in names(tLim)) {
      tLim_value <- tLim[[name]]
      if (!is.numeric(tLim_value)) {
        stop(paste("tLim['", name, "'] must be numeric or numeric sequences.", sep = ""))
      }
      # No spreading of values, keep as is
    }
  } else if (is.numeric(tLim) && length(tLim) == 1) {
    # Single numeric value, apply to both HDD and CDD without spreading
    tLim <- list(
      HDD = tLim,
      CDD = tLim
    )
  } else {
    stop("tLim must be a numeric value or a named list/vector with 'HDD' and 'CDD'.")
  }


  # Process std ----------------------------------------------------------------

  if (is.numeric(std)) {
    if (length(std) == 1) {
      # Single numeric value, apply to both tLim and tAmb
      std <- c("tLim" = std, "tAmb" = std)
    } else if (length(std) == 2) {
      # Check if names are "tLim" and "tAmb"
      if (is.null(names(std))) {
        names(std) <- c("tLim", "tAmb")
      } else if (!all(names(std) %in% c("tLim", "tAmb"))) {
        stop("std must have names 'tLim' and 'tAmb'.")
      }
    } else {
      stop("std must be a single numeric value or a named numeric vector of length 2 with 'tLim' and 'tAmb'.")
    }
  } else {
    stop("std must be numeric.")
  }


  # Process ssp ----------------------------------------------------------------

  if (is.character(ssp)) {
    ssp <- as.vector(ssp) %>%
      tolower()
  } else {
    stop("ssp must be a character string or vector.")
  }

  # Return standardized arguments
  return(list(tLim = tLim, std = std, ssp = ssp))
}
