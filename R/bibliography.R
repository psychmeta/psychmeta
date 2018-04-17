#' Generate a list of references included in meta-analyses
#'
#' This function generates a list of studies contributing to a meta-analysis
#'
#' @param ma_obj A psychmeta meta-analysis object with \code{citekeys} supplied.
#' @param additional_citekeys Additional citekeys to include in the reference list.
#' @param bib A BibTeX file containing the citekeys for the meta-analyses.
#' @param analyses Which analyses to extract references for? See \code{\link{filter_ma}} for details.
#' @param match Match \code{all} or \code{any} of the filter criteria? See \code{\link{filter_ma}} for details.
#' @param style References should be formatted in what style? Can be the style ID for any \url{https://github.com/citation-style-language/styles}{CSL style} (formatted examples of styles are available from the \url{https://zotero.org/styles}{Zotero Style Repository}). Defaults to APA style. (Requires an internet connection to retrieve styles. If unavailable, references will be rendered in Chicago style.)
#' @param output_format The format of the output reference list. Available options are Word (default), HTML, PDF (requires the \code{tinytex} package), Rmarkdown, plain text, and BibLaTeX. Returning only the item citekeys is also possible.
#' @param file The filename or filepath for the output file. If \code{NULL}, file will be saved as \code{reference_list}. Set to \code{"console"} or \code{"print"} to output directly to the R console.
#' @param header A list of YAML header parameters to pass to \code{link{rmarkdown::render}}.
#'
#' @return A formatted reference list.
#' @export
#'
#' @importFrom rmarkdown render
#' @importFrom RefManageR PrintBibliography
#' @importFrom RefManageR WriteBib
#' @importFrom RefManageR ReadBib
#' @importFrom RCurl url.exists
#' @importFrom RefManageR BibOptions
#'
#' @examples
#' \dontrun{
#' ma_obj <- ma_r(ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'                construct_x = x_name, construct_y = y_name, sample_id = sample_id,
#'                moderators = moderator, citekey = citekey, data = data_r_meas_multi)
#'
#' generate_bib(ma_obj, analyses="all", match=c("all", "any"),
#' bib=system.file("sample_bibliography.bib", package="psychmeta"), style="apa",
#' output_format="word", header=list())
#' }
generate_bib <- function(ma_obj=NULL, additional_citekeys=NULL, bib=NULL, analyses="all", match=c("all", "any"),
                         style="apa", output_format=c("word", "html", "pdf", "text", "Rmd", "biblatex", "citekeys"),
                         file=NULL, header=list()){

     match <- match.arg(match)
     output_format <- match.arg(output_format)

     if(is.null(ma_obj) & is.null(additional_citekeys)) stop("Either ma_obj or additional_citekeys must be provided.")

     if(!is.null(ma_obj)) {
          # Get the requested meta-analyses
          ma_obj_filtered <- filter_ma(ma_obj=ma_obj, analyses=analyses, match=match)

          # Get the citekeys for requested metas
          citekeys <- if(is.null(ma_obj_filtered$construct_pairs)) {
               unique(c(additional_citekeys, unlist(strsplit(paste(unlist(lapply(ma_obj_filtered$barebones$escalc_list, function(y) as.character(y$citekey))), collapse=", "), ", "))))
          } else {
               unique(c(additional_citekeys, unlist(strsplit(paste(unlist(lapply(ma_obj_filtered$construct_pairs, function(x) lapply(x$barebones$escalc_list, function(y) as.character(y$citekey)))), collapse=", "), ", "))))
          }
     } else citekeys <- additional_citekeys

     # Set BibOptions to accept entries with missing data
     Bib_check.entries_original <- BibOptions()[["check.entries"]]
     BibOptions(check.entries = FALSE)

     # Read in .bib file
     bib <- ReadBib(bib)

     # Set the output file name
     if(is.null(file)) {
          name_check <- TRUE
     } else if(file != "console") {
          name_check <- TRUE
     } else name_check <- FALSE

     if(name_check) {
          if(output_format == "pdf") {
               if(is.null(file)) file <- "reference_list.pdf" else {
                    if(!grep("\\.pdf$", file)) file <- paste0(file, ".pdf")
               }
          } else if (output_format == "html") {
               if(is.null(file)) file <- "reference_list.html" else {
                    if(!grep("\\.html$", file)) file <- paste0(file, ".html")
               }
          } else if(output_format == "word") {
               if(is.null(file)) file <- "reference_list.docx" else {
                    if(!grep("\\.docx$", file)) file <- paste0(file, ".docx")
               }
          } else if(output_format == "Rmd") {
               if(is.null(file)) file <- "reference_list.Rmd" else {
                    if(!grep("\\.Rmd$", file)) file <- paste0(file, ".Rmd")
               }
          } else if(output_format == "biblatex") {
               if(is.null(file)) file <- "reference_list.bib" else {
                    if(!grep("\\.bib$", file)) file <- paste0(file, ".bib")
               }
          } else if(output_format == "text") {
               if(is.null(file)) file <- "reference_list.txt" else {
                    if(!grep("\\.txt$", file)) file <- paste0(file, ".txt")
               }
          } else if(output_format == "citekeys") {
               if(is.null(file)) file <- "citekeys.txt" else {
                    if(!grep("\\.txt$", file)) file <- paste0(file, ".txt")
               }
          }
     }

     # List the citations to include
     citations <- paste0("@", citekeys, collapse=", ")

     # Render the output

     if(file == "console") {
          if(output_format == "citekeys") {
               return(citations)
          } else if(output_format == "biblatex") {
               print(bib[citekeys], .opts = list(style = "Biblatex"))
          } else {
               print(bib[citekeys], .opts = list(style = "text", bib.style = "authoryear"))
          }
     } else {
          if(output_format %in% c("Rmd", "pdf", "html", "word")) {

               ## Fill in critical header slots
               if(is.null(header$title)) {
                    header$title <- "'Sources contributing to meta-analyses'"
               } else header$title <- paste0("'", header$title, "'")
               header$bibliography <- gsub(pattern="\\.(Rmd|pdf|docx|html)$",
                                           replacement="\\.bib", x=file)
               if( url.exists( paste0("https://zotero.org/styles/", gsub("\\.csl$", "", "apa") ) ) ) {
                    header$csl <- paste0("https://zotero.org/styles/", gsub("\\.csl$", "", "apa") )
               }

               ## Create the markdown header and document
               header <- paste(names(header), header, sep=": ", collapse="\n")
               document <- sprintf("---\n%s\nnocite: |\n  %s\n---\n", header, citations)

               ## Write the bibliography and Rmd file
               suppressMessages(WriteBib(bib[citekeys],
                                         file = gsub(pattern="\\.(Rmd|pdf|docx|html)$",
                                                     replacement="\\.bib", x=file)))
               write(document, file = gsub(pattern="\\.(pdf|docx|html)$",
                                           replacement="\\.Rmd", x=file))
               if(output_format != "Rmd") {
                    render(gsub(pattern="\\.(pdf|docx|html)$",
                                replacement="\\.Rmd", x=file),
                           output_file = file,
                           output_format = paste0(output_format, "_document") )
               }
          } else if(output_format == "biblatex") {
               suppressMessages(WriteBib(bib[citekeys], file = file))
          } else if(output_format == "text") {
               sink(file)
               print(bib[citekeys], .opts = list(style = "text", bib.style = "authoryear"))
               sink()
          } else if(output_format == "citekeys") {
               sink(file)
               print(citations)
               sink()
          }
     }

     # Reset BibOptions to original
     BibOptions(check.entries = Bib_check.entries_original)
}

