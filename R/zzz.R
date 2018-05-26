## Messages to be displayed when the user loads psychmeta:
#' @importFrom rlang .data
.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
    packageStartupMessage("This is ", paste(pkgname, "version", version))
    packageStartupMessage("Please report any bugs to github.com/psychmeta/psychmeta/issues or issues@psychmeta.com")
    packageStartupMessage("\nWe work hard to produce these open-source tools for the R community, \nplease cite psychmeta when you use it in your research: \nDahlke, J. A., & Wiernik, B. M. (in press). psychmeta: An R package for \n   psychometric meta-analysis. Applied Psychological Measurement.")
    packageStartupMessage("\nFind info about psychmeta on the web at psychmeta.com and twitter.com/psychmetaR")

    # Check if there is an internet connection. If there is, check whether the local version of psychmeta is up to date compared to the CRAN version.
    if(try(is.character(RCurl::getURL("http://www.r-pkg.org/badges/version/psychmeta")), silent = TRUE) == TRUE){
         pkg_badge <- xml2::read_html("http://www.r-pkg.org/badges/version/psychmeta")
         cran_v_char <- gsub(x = stringr::str_split(as.character(pkg_badge), "\n")[[1]][9], pattern = " ", replacement = "")
         cran_v_num <- as.numeric(stringr::str_split(cran_v_char, "[.]")[[1]])
         sys_v_char <- stringr::str_split(version, "[.]")[[1]]
         sys_v_num <- as.numeric(sys_v_char)

         if(cran_v_char == version){
              out_of_date <- FALSE
         }else{
              best_version <- sort(c(CRAN = cran_v_char, Local = version), decreasing = TRUE)[1]
              out_of_date <- names(best_version) == "CRAN"
         }

         if(length(cran_v_num) == 3) cran_v_num <- c(cran_v_num, 0)
         if(length(sys_v_num) == 3) sys_v_num <- c(sys_v_num, 0)
         vcheck_equal <- sys_v_num == cran_v_num
         vcheck_greater <- zapsmall(sys_v_num) > zapsmall(cran_v_num)
         vcheck <- vcheck_equal | vcheck_greater
         development <- all(vcheck) & any(vcheck_greater)

         # use_symbols <- getOption("cli.unicode") | !(.Platform$OS.type == "windows")

         if(out_of_date){
              version_message <- "Oh no! It looks like your copy of psychmeta is out of date!"
              # if(use_symbols) version_message <- paste0(crayon::red(cli::symbol$cross), " ", version_message)
              packageStartupMessage(paste0("\n", version_message))
              packageStartupMessage("No worries, it's easy to obtain the latest version - just run the following command: \n")
              packageStartupMessage('                       update.packages("psychmeta")')
         }else{
              if(development){
                   version_message <- "Kudos! Your copy of psychmeta is more recent than the current CRAN release!"
              }else{
                   version_message <- "Yay! Your copy of psychmeta is up to date!"
              }
              # if(use_symbols) version_message <- paste0(crayon::green(cli::symbol$tick), " ", version_message)
              packageStartupMessage(paste0("\n", version_message))
         }

         if(sys_v_num[4] > 0)
              packageStartupMessage(paste0("NOTE: You are currently using an UNRELEASED development build of psychmeta (augmentation of release v", paste(sys_v_char[1:3], collapse = "."), ")"))
    }
}


#' Retrieve the NEWS file for the psychmeta package
#'
#' @description
#' This function gives a shortcut to the \code{utils::news(package = "psychmeta")} function and displays psychmeta's NEWS file, which contains version information, outlines additions and changes to the package, and describes other updates.
#'
#' @export
#'
#' @importFrom utils news
#'
#' @examples
#' psychmeta_news()
psychmeta_news <- function(){
     news(package = "psychmeta")
}

