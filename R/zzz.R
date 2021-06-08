## Register S3 methods for dplyr verbs to avoid having to export generics
.onLoad <- function(libname, pkgname) {
  register_s3_method("dplyr", "arrange", "ma_psychmeta")
  register_s3_method("dplyr", "arrange", "ma_table")

  register_s3_method("dplyr", "filter", "ma_psychmeta")
  register_s3_method("dplyr", "filter", "ma_table")

  register_s3_method("dplyr", "mutate", "ma_psychmeta")
  register_s3_method("dplyr", "mutate", "ma_table")

  register_s3_method("dplyr", "rename", "ma_psychmeta")
  register_s3_method("dplyr", "rename", "ma_table")

  register_s3_method("dplyr", "select", "ma_psychmeta")
  register_s3_method("dplyr", "select", "ma_table")
}

## Messages to be displayed when the user attaches psychmeta:
.onAttach <- function(libname, pkgname) {
  cli_available <- requireNamespace("cli", quietly = TRUE)
  if (requireNamespace("crayon", quietly = TRUE)) {
    green <- crayon::green
    white <- crayon::white
    red <- crayon::red
    bold <- crayon::bold
    italic <- crayon::italic
  } else {
    green <- white <- red <- bold <- italic <- identity
  }

  version <- read.dcf(file = system.file("DESCRIPTION", package = pkgname),
                      fields = "Version")
  packageStartupMessage(
    white("----------------------------------------------------- ",
          bold(paste(pkgname, "version", version)),
          " --")
    )
  packageStartupMessage(
    "\nPlease report any bugs to ",
    italic("github.com/psychmeta/psychmeta/issues"),
    "\nor ",
    italic("issues@psychmeta.com")
    )
  packageStartupMessage(
    "\nWe work hard to produce these open-source tools for the R community.",
    "\nPlease cite psychmeta when you use it in your research:",
    "\n  Dahlke, J. A., & Wiernik, B. M. (2019). psychmeta: An R package for",
    "\n    psychometric meta-analysis. ",
    italic("Applied Psychological Measurement, 43"), "(5), 415-416.",
    "\n    https://doi.org/10.1177/0146621618795933"
    )

  # Check if there is an internet connection. If there is, check whether the
  # local version of psychmeta is up to date compared to the CRAN version.
  packageStartupMessage(
    white(
      paste0("\n-----------------------------------------------------",
             paste0(rep_len("-",
                            nchar(paste(pkgname, "version", version)) - 13),
                    collapse = ""), " "), crayon::bold("Version check"), " --"))
  suppressWarnings(
    version_check <-
      try(unlist(strsplit(
        rawToChar(
          curl::curl_fetch_memory("https://CRAN.R-project.org/package=psychmeta")$content
        ), "\n")),
        silent = TRUE)
    )
  if (!inherits(version_check, "try-error")) {

    check_version <- function(cran_version, sys_version) {
      if (cran_version == sys_version) {
           out_of_date <- FALSE
           ahead_of_cran <- FALSE
           best_version <- "Current"
           development <- FALSE
      } else {
        cran_v_num <- as.numeric(strsplit(cran_v_char, "[.]")[[1]])
        sys_v_num <- as.numeric(strsplit(sys_version, "[.]")[[1]])

        if (length(cran_v_num) == 3) cran_v_num <- c(cran_v_num, "0")
        if (length(sys_v_num) == 3) sys_v_num <- c(sys_v_num, "0")

        cran_v_num <- cran_v_num
        sys_v_num <- sys_v_num
        desired_digits <- 4
        necessary_digits_cran <- desired_digits - nchar(cran_v_num)
        necessary_digits_sys <- desired_digits - nchar(sys_v_num)
        for(i in 1:length(sys_v_num)){
          cran_v_num[i] <- paste0(c(rep("0", necessary_digits_cran[i]),
                                    cran_v_num[i]),
                                  collapse = "")
          sys_v_num[i] <- paste0(c(rep("0", necessary_digits_sys[i]),
                                   sys_v_num[i]),
                                 collapse = "")
        }

        cran_v_num <- paste(cran_v_num, collapse = "")
        sys_v_num <- paste(sys_v_num, collapse = "")

        best_version <- sort(c(CRAN = cran_v_num, Local = sys_v_num),
                             decreasing = TRUE)[1]
        out_of_date <- names(best_version) == "CRAN"
        ahead_of_cran <- names(best_version) == "Local"
        if (out_of_date) best_version <- "CRAN"
        if (ahead_of_cran) best_version <- "Local"

        development <- zapsmall(as.numeric(sys_v_num[4])) > 0
      }

      data.frame(best_version = best_version,
                 cran_version = cran_version,
                 sys_version = sys_version,
                 out_of_date = out_of_date,
                 ahead_of_cran = ahead_of_cran,
                 development = development,
                 stringsAsFactors = FALSE)
    }

    cran_v_char <- version_check[grep("Version:", version_check) + 1]
    cran_v_char <- regmatches(cran_v_char,
                              regexpr("(?:\\d+\\.?)+", cran_v_char))
    if (length(cran_v_char) > 0) {
      vcheck <- check_version(cran_version = cran_v_char,
                              sys_version = version)

      if (vcheck$best_version == "CRAN") {
         packageStartupMessage(red(c(
           if (cli_available) cli::symbol$cross,
           " Oh no! It looks like your copy of psychmeta is out of date!"
         )))
         packageStartupMessage("\nTo update to the latest version, run:")
         packageStartupMessage('  install.packages("psychmeta")')
      } else if (vcheck$best_version == "Current") {
         packageStartupMessage(green(c(
           if (cli_available) cli::symbol$tick,
           " Yay! Your copy of psychmeta is up to date!"
         )))
      } else if(vcheck$best_version == "Local") {
         packageStartupMessage(green(c(
           if (cli_available) cli::symbol$tick,
           " Kudos! Your copy of psychmeta is more recent than the current CRAN release!"
         )))
      }

    } else {
      packageStartupMessage(c(
        if (cli_available) cli::symbol$cross,
        " CRAN version of psychmeta not found."
      ))
    }

  } else {
    packageStartupMessage(c(
      if (cli_available) cli::symbol$cross,
      " Version check not run."
    ))
  }

  # Report if version is development version
  sys_v_num <- as.numeric(strsplit(version, "[.]")[[1]])
  if (length(sys_v_num) >= 4) {
    if (sys_v_num[4] > 0) {
      packageStartupMessage(
        "NOTE: You are currently using an UNRELEASED development build (augmentation of release ",
        paste0("v", paste(sys_v_num[1:3], collapse = "."), ")")
      )
    }
  }
}

.support_unicode <- function(override = NULL) {
    if (isTRUE(override)) {
      TRUE
    } else {
      l10n_info()$`UTF-8` |
        isTRUE(.Options$cli.unicode) |
        nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))
    }
}


#' Retrieve the NEWS file for the psychmeta package
#'
#' @description
#' This function gives a shortcut to the \code{utils::news(package = "psychmeta")} function and displays psychmeta's NEWS file, which contains version information, outlines additions and changes to the package, and describes other updates.
#'
#' @export
#'
#' @examples
#' psychmeta_news()
psychmeta_news <- function(){
  if (requireNamespace("commonmark", quietly = TRUE)) {
    utils::news(package = "psychmeta")
  } else {
    message("Run install.packages('commonmark') to view psychmeta news.")
  }
}

register_s3_method <- function(pkg, generic, class, fun = NULL) {
  stopifnot(is.character(pkg), length(pkg) == 1)
  envir <- asNamespace(pkg)

  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)
  if (is.null(fun)) {
    fun <- get(paste0(generic, ".", class), envir = parent.frame())
  }
  stopifnot(is.function(fun))


  if (pkg %in% loadedNamespaces()) {
    registerS3method(generic, class, fun, envir = envir)
  }

  # Always register hook in case package is later unloaded & reloaded
  setHook(
    packageEvent(pkg, "onLoad"),
    function(...) {
      registerS3method(generic, class, fun, envir = envir)
    }
  )
}

# Supplemental imports
#' @importFrom rlang .data
#' @import mathjaxr
