## Global variables to be called from data frames with known variables names:
globalVariables(c(".", "Value",                                  ## Global variables defined for the create_ad function
                  "Analysis_ID",                                 ## Global variable defined for the organize_moderators function
                  "vi", "yi",                                    ## Global variables defined for the create_ad function
                  "mean_r", "sd_res", "rtpa", "rxpa", "rtya",    ## Global variables defined for the .ma_r_ic function
                  "rxy",                                         ## Global variables defined for the .ma_r_bb function
                  "d", "n_1", "n_2", "n_1_split", "n_2_split",   ## Global variables defined for the .ma_d_bb function
                  "difference", "construct_x", "construct_y"))



## Messages to be displayed when the user loads psychmeta:
.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), fields="Version")
    packageStartupMessage("This is ", paste(pkgname, version))
    packageStartupMessage("Please report any bugs to improve functionality. \n")
    packageStartupMessage("We work hard to produce these open-source tools for the R community - please cite psychmeta when you use it in your research: \nDahlke, J. A. & Wiernik, B. M. (2017). psychmeta: Psychometric meta-analysis toolkit. R package version ", version)
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

