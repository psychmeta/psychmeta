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
        
         check_version <- function(cran_version, sys_version){
              cran_v_char <- cran_version
              cran_v_num <- as.numeric(stringr::str_split(cran_v_char, "[.]")[[1]])
              sys_v_char <- stringr::str_split(sys_version, "[.]")[[1]]
              sys_v_num <- as.numeric(sys_v_char)
              
              if(length(cran_v_num) == 3) cran_v_num <- c(cran_v_num, "0")
              if(length(sys_v_num) == 3) sys_v_num <- c(sys_v_num, "0")
              
              .cran_v_num <- cran_v_num
              .sys_v_num <- sys_v_num
              desired_digits <- 4
              necessary_digits_cran <- desired_digits - nchar(cran_v_num)
              necessary_digits_sys <- desired_digits - nchar(sys_v_num)
              for(i in 1:length(sys_v_num)){
                   .cran_v_num[i] <- paste0(c(rep("0", necessary_digits_cran[i]), .cran_v_num[i]), collapse = "")
                   .sys_v_num[i] <- paste0(c(rep("0", necessary_digits_sys[i]), .sys_v_num[i]), collapse = "")
              }
              
              .cran_v_num <- paste(.cran_v_num, collapse = "")
              .sys_v_num <- paste(.sys_v_num, collapse = "")
              
              if(cran_v_char == sys_version){
                   out_of_date <- FALSE
                   ahead_of_cran <- FALSE
                   best_version <- c(Current = cran_version)
              }else{
                   best_version <- sort(c(CRAN = .cran_v_num, Local = .sys_v_num), decreasing = TRUE)[1]
                   out_of_date <- names(best_version) == "CRAN"
                   ahead_of_cran <- names(best_version) == "Local"
                   if(names(best_version) == "CRAN") best_version <- c(CRAN = cran_version)
                   if(names(best_version) == "Local") best_version <- c(Local = sys_version)
              }
              
              vcheck_devnum <- zapsmall(as.numeric(sys_v_num[4])) > 0
              
              as.data.frame(list(best_version = names(best_version), 
                                 cran_version = cran_version, 
                                 sys_version = sys_version, 
                                 out_of_date = out_of_date, 
                                 ahead_of_cran = ahead_of_cran, 
                                 development = vcheck_devnum))
         }
         
         pkg_badge <- xml2::read_html("http://www.r-pkg.org/badges/version/psychmeta")
         cran_v_char <- gsub(x = stringr::str_split(as.character(pkg_badge), "\n")[[1]][9], pattern = " ", replacement = "")
         vcheck <- check_version(cran_version = cran_v_char, sys_version = version)
         
         use_symbols <- .Platform$OS.type != "windows"
         
         if(vcheck$best_version == "CRAN"){
              version_message <- "Oh no! It looks like your copy of psychmeta is out of date!"
              if(use_symbols) version_message <- paste0(crayon::red(cli::symbol$cross), " ", version_message)
              packageStartupMessage(paste0("\n", version_message))
              packageStartupMessage("No worries, it's easy to obtain the latest version - just run the following command: \n")
              packageStartupMessage('                       update.packages("psychmeta")')
         }else if(vcheck$best_version == "Current"){
              version_message <- "Yay! Your copy of psychmeta is up to date!"
              if(use_symbols) version_message <- paste0(crayon::green(cli::symbol$tick), " ", version_message)
              packageStartupMessage(paste0("\n", version_message))
         }else if(vcheck$best_version == "Local"){
              version_message <- "Kudos! Your copy of psychmeta is more recent than the current CRAN release!"
              if(use_symbols) version_message <- paste0(crayon::green(cli::symbol$tick), " ", version_message)
              packageStartupMessage(paste0("\n", version_message))
         }

    }
    
    sys_v_char <- stringr::str_split(version, "[.]")[[1]]
    sys_v_num <- as.numeric(sys_v_char)
    if(length(sys_v_num) == 3) sys_v_num <- c(sys_v_num, 0)
    if(sys_v_num[4] > 0)
         packageStartupMessage(paste0("NOTE: You are currently using an UNRELEASED development build (augmentation of release v", paste(sys_v_char[1:3], collapse = "."), ")"))
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

