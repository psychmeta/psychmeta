#' Format numbers for presentation
#'
#' A function to format decimal digits, leading zeros, and sign characters.
#'
#' @encoding UTF-8
#'
#' @param x A vector, matrix, or data.frame of numbers to format
#' @param digits The number of decimal digits desired (used strictly; default: 2)
#' @param decimal.mark The character to use for the decimal point (defaults to locale default: \code{getOption("OutDec")})
#' @param leading0 How to print leading zeros on decimals. Can be logical to print (\code{TRUE}) or suppress (\code{FALSE}) leading zeros or a character string to subsitute for leading zeros. If \code{"figure"} (default), leading zeros are replaced with figure-spaces (\code{U+2007}: "<U+2007>") if a column contains any absolute values greater than 1 and suppressed otherwise.
#' @param neg.sign Character to use as negative sign. Defaults to minus-sign (\code{U+2212}: "<U+2212>").
#' @param pos.sign Character to use as positive sign. Set to \code{FALSE} to suppress. If \code{"figure"} (default), the positive sign is a figure-space (\code{U+2007}: "<U+2007>") if a column contains any negative numbers and suppressed otherwise.
#' @param drop0integer Logical. Should trailing decimal zeros be dropped for integers?
#' @param big.mark Character to mark between each \code{big.interval} digits \emph{before} the decimal point. Set to \code{FALSE} to suppress. Defaults to the SI/ISO 31-0 standard-recommened thin-spaces (\code{U+202F}: "<U+202F>").
#' @param big.interval See \code{big.mark} above; defaults to 3.
#' @param small.mark Character to mark between each \code{small.interval} digits \emph{after} the decimal point. Set to \code{FALSE} to suppress. Defaults to the SI/ISO 31-0 standard-recommened thin-spaces (\code{U+202F}: "<U+202F>").
#' @param small.interval See \code{small.mark} above; defaults to 3.
#'
#' @importFrom stringr str_replace
#'
#' @export
#' @examples
#' num_format(x = c(10000, 1000, 2.41, -1.20, 0.41, -0.20))
#'
#' # By default, num_format() uses your computer locale's default decimal mark as
#' # the decimal point. To force the usage of "." instead (e.g., for submission to
#' # an American journal), set decimal.mark = ".":
#' num_format(x = .41, decimal.mark = ".")
#'
#' # By default, num_format() separates groups of large digits using thin spaces.
#' # This is following the international standard for scientific communication (SI/ISO 31-0),
#' # which advises against using "." or "," to seprate digits because doing so can lead
#' # to confusion for human and computer readers because "." and "," are also used
#' # as decimal marks in various countries. If you prefer to use commmas to separate
#' # large digit groups, set big.mark = ",":
#' num_format(x = 10000, big.mark = ",")
num_format <- function(x, digits = 2L, decimal.mark = getOption("OutDec"), leading0 = "figure", neg.sign = "minus",
                       pos.sign = "figure", drop0integer = TRUE, big.mark = "thinspace",
                       big.interval = 3L, small.mark = "thinspace", small.interval = 3L) {

        is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {abs(x - round(x)) < tol}
        all_equal_vector <- function(x, tol = .Machine$double.eps^0.5) {diff(range(x)) < tol}

        if(is.null(dim(x))) {
                x_type <- "vector"
                x <- as.data.frame(x)
        } else if("tbl_df" %in% class(x)) {
                x_type <- "tibble"
                x <- as.data.frame(x)
        } else if("matrix" %in% class(x)) {
                x_type <- "matrix"
                x <- as.data.frame(x)
        } else x_type <- "other"

        which_integers <- is.wholenumber(x)

        if(neg.sign == "minus") neg.sign <- "\u2212"
        if(big.mark == "thinspace") big.mark <- "\u202F"
        if(small.mark == "thinspace") small.mark <- "\u202F"

        if(pos.sign == FALSE) flag <- "" else flag <- "+"
        out <- x

        # Initial formatting of numbers

        out[which_integers] <- x[which_integers] %>%
                formatC(digits = digits, format = "f", flag = flag,
                        decimal.mark = decimal.mark,
                        big.mark = big.mark, big.interval = big.interval,
                        small.mark = small.mark, small.interval = small.interval,
                        drop0trailing = drop0integer)

        out[!which_integers] <- x[!which_integers] %>%
                formatC(digits = digits, format = "f", flag = flag,
                        decimal.mark = decimal.mark,
                        big.mark = big.mark, big.interval = big.interval,
                        small.mark = small.mark, small.interval = small.interval,
                        drop0trailing = FALSE)

        # Clean up unicode big.mark and small.mark
        out[] <- mutate_all(out,
                            funs(str_replace(.data$.,
                                             paste0("(",paste(rev(strsplit(sub(" ", big.mark, " "),"")[[1]]), collapse=""),"|",sub(" ", big.mark, " "),")"),
                                             big.mark)))
        out[] <- mutate_all(out,
                            funs(str_replace(.data$.,
                                             paste0("(",paste(rev(strsplit(sub(" ", small.mark, " "),"")[[1]]), collapse=""),"|",sub(" ", small.mark, " "),")"),
                                             small.mark)))

        # Clean up leading zeros
        switch(leading0,
               "TRUE" = {},
               "FALSE" = out[] <- sapply(out, function(x) str_replace(x, paste0("^(-?)0", decimal.mark), paste0("\\1", decimal.mark))),
               figure = {
                       out <- mutate_if(out, apply(x, 2, function(i) {any(abs(i) >= 1)}),
                                        function(i) str_replace(i, paste0("^(\\+|-?)0", decimal.mark), paste0("\\1\u2007", decimal.mark)))
                       out <- mutate_if(out, apply(x, 2, function(i) {!any(abs(i) >= 1)}),
                                        function(i) str_replace(i, paste0("^(\\+|-?)0", decimal.mark), paste0("\\1", decimal.mark)))
               },
               # else =
               {out[] <- sapply(out, function(x) str_replace(x, paste0("^(\\+|-?)0", decimal.mark), paste0("\\1", leading0, decimal.mark)))}
        )

        # Clean up positive signs
        switch(pos.sign,
               "TRUE" = {},
               "FALSE" = {},
               figure = {
                       out <- mutate_if(out, apply(x, 2, function(i) {!any(i < 0)}),
                                        function(i) str_replace(i, "^\\+", ""))
                       out <- mutate_if(out, apply(x, 2, function(i) {any(i < 0)}),
                                        function(i) str_replace(i, "^\\+", "\u2007"))
               },
               # else =
               {out[] <- sapply(out, function(i) str_replace(i, "^\\+", pos.sign))}
        )

        # Clean up negative signs
        switch(neg.sign,
               "TRUE" = {},
               "FALSE" = {},
               "-" = {},
               # else =
               {out[] <- sapply(out, function(x) str_replace(x, "^-", neg.sign))}
        )

        if(x_type == "tibble") {
                out <- as_tibble(out, validate = FALSE)
        } else if(x_type == "vector") {
                out <- unlist(out)
        } else if(x_type == "matrix") {
                out <- as.matrix(out)
        }

        return(out)
}


#' Write a summary table of meta-analytic results
#'
#' @param ma_obj A psychmeta meta-analysis object.
#' @param file The filename or filepath for the output file. If \code{NULL}, file will be saved as \code{psychmeta_output}. Set to \code{"console"} or \code{"print"} to output directly to the R console.
#' @param show_msd Logical. Should means and standard deviations of effect sizes be shown (default \code{TRUE})
#' @param show_conf Logical. Should confidence intervals be shown (default: \code{TRUE})?
#' @param show_cred Logical. Should credibility intervals be shown (default: \code{TRUE})?
#' @param show_se Logical Should standard errors be shown (default: \code{FALSE})?
#' @param show_var Logical. Should variances be shown (default: \code{FALSE})?
#' @param analyses Which analyses to extract references for? See \code{\link{filter_ma}} for details.
#' @param match Match \code{all} or \code{any} of the filter criteria? See \code{\link{filter_ma}} for details.
#' @param case_sensitive Logical scalar that determines whether character values supplied in \code{analyses} should be treated as case sensitive (\code{TRUE}, default) or not (\code{FALSE}).
#' @param output_format The format of the output tables. Available options are Word (default), HTML, PDF (requires LaTeX, see the \code{tinytex} package), Rmarkdown, and plain text.
#' @param ma_method Meta-analytic methods to be included. Valid options are: \code{"ad"}, \code{"ic"}, and \code{"bb"}. Multiple methods are permitted. By default, results are given for one method with order of priority: 1. \code{"ad"}, 2. \code{"ic"}, 3. \code{"bb"}.
#' @param correction_type Type of meta-analytic corrections to be incldued. Valid options are: "ts" (default), "vgx", and "vgy". Multiple options are permitted.
#' @param bib A BibTeX file containing the citekeys for the meta-analyses. If not \code{NULL}, a bibliography will be included with the meta-analysis table. See \code{\link{generate_bib}} for additional arguments controlling the bibliography.
#' @param title.bib The title to give to the bibliography. If \code{NULL}, defaults to "Sources Contributing to Meta-Analyses"
#' @param additional_citekeys Additional citekeys to include in the reference list.
#' @param style What style should references be formatted in? Can be a file path or URL for a \url{https://github.com/citation-style-language/styles}{CSL citation style} or the style ID for any style available from the \url{https://zotero.org/styles}{Zotero Style Repository}). Defaults to APA style. (Retrieving a style by ID requires an internet connection. If unavailable, references will be rendered in Chicago style.).
#' @param header A list of YAML header parameters to pass to \code{link{rmarkdown::render}}.
#' @param digits The number of decimal digits desired (used strictly; default: 2)
#' @param decimal.mark The character to use for the decimal point (defaults to locale default: \code{getOption("OutDec")})
#' @param leading0 How to print leading zeros on decimals. See \code{\link{num_format}} for details.
#' @param neg.sign Character to use as negative sign. See \code{\link{num_format}} for details.
#' @param pos.sign Character to use as positive sign. See \code{\link{num_format}} for details.
#' @param drop0integer Logical. Should trailing decimal zeros be dropped for integers?
#' @param big.mark Character to separate groups of large digits. See \code{\link{num_format}} for details.
#' @param big.interval See \code{big.mark} above; defaults to 3.
#' @param small.mark Character to sparate groups of decimal digits. See \code{\link{num_format}} for details.
#' @param small.interval See \code{small.mark} above; defaults to 3.
#' @param conf_format How should confidence intervals be formatted? Options are:
#' \itemize{
#' \item{parentheses}{Bounds are enclosed in parentheses and separated by a comma: (LO, UP).}
#' \item{brackets}{Bounds are enclosed in square brackets and separated by a comma: [LO, UP].}
#' \item{columns}{Bounds are shown in individual columns.}
#' }
#' @param cred_format How should credility intervals be formatted? Options are the same as for \code{conf_format} above.
#' @param symbol_es For meta-analyses of generic (non-r, non-d) effect sizes, the symbol used for the effect sizes (default: \code{symbol_es = "ES"}).
#' @param caption Caption to print before tables. Either a character scalar or a named character vector with names corresponding to combinations of \code{ma_method} and \code{correction_type} (i.e., \code{bb}, \code{ic_ts}, \code{ad_vgx}, etc.).
#' @param verbose Logical. Should detailed SD and variance components be shown (default: \code{FALSE})?
#' @param save_build_files Should the RMarkdown and BibTeX (files) files used to generate the output be saved (default: \code{TRUE})?
#' @param ... Additional arguments (not used).
#'
#' @return Formatted tables of meta-analytic output.
#' @export
#'
#' @importFrom rmarkdown render
#' @importFrom dplyr case_when
#' @importFrom stringi stri_write_lines
#'
#' @examples
#' \dontrun{
#' ## Create output table for meta-analysis of correlations:
#' ma_r_obj <- ma_r(ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'                  construct_x = x_name, construct_y = y_name,
#'                  moderators = moderator, data = data_r_meas_multi)
#' ma_r_obj <- ma_r_ad(ma_obj = ma_r_obj, correct_rr_x = FALSE, correct_rr_y = FALSE)
#' metabulate(ma_obj = ma_r_obj, file = "meta tables correlations")
#'
#' ## Create output table for meta-analysis of d values:
#' ma_d_obj <- ma_d(ma_method = "ic", d = d, n1 = n1, n2 = n2, ryy = ryyi,
#'                  construct_y = construct, data = data_d_meas_multi)
#' ma_d_obj <- ma_d_ad(ma_obj = ma_d_obj, correct_rr_g = FALSE, correct_rr_y = FALSE)
#' metabulate(ma_obj = ma_d_obj, file = "meta tables d values")
#'
#' ## Create output table for meta-analysis of generic effect sizes:
#' dat <- data.frame(es = data_r_meas_multi$rxyi,
#'                   n = data_r_meas_multi$n,
#'                   var_e = (1 - data_r_meas_multi$rxyi^2)^2 / (data_r_meas_multi$n - 1))
#' ma_obj <- ma_generic(es = es, n = n, var_e = var_e, data = dat)
#' metabulate(ma_obj = ma_obj, file = "meta tables generic es")
#' }
metabulate <- function(ma_obj, file, show_msd = TRUE, show_conf = TRUE, show_cred = TRUE,
                       show_se = FALSE, show_var = FALSE, analyses="all",
                       match=c("all", "any"), case_sensitive = TRUE,
                       output_format=c("word", "html", "pdf", "text", "rmd"),
                       ma_method = "ad", correction_type = "ts",
                       bib = NULL, title.bib = NULL, additional_citekeys = NULL, style = "apa",
                       header = NULL, digits = 2L, decimal.mark = getOption("OutDec"),
                       leading0 = "figure", neg.sign = "minus", pos.sign = "figure",
                       drop0integer = TRUE, big.mark = "thinspace", big.interval = 3L,
                       small.mark = "thinspace", small.interval = 3L,  conf_format = "parentheses",
                       cred_format = "brackets", symbol_es = "ES", caption = "Results of meta-analyses",
                       save_build_files = TRUE, verbose = FALSE, ...){

        # Match arguments
        output_format <- tolower(output_format)
        output_format <- match.arg(output_format)
        ma_method <- match.arg(ma_method, c("bb", "ic", "ad"), several.ok = TRUE)
        correction_type <- match.arg(correction_type, c("ts", "vgx", "vgy"), several.ok = TRUE)
        conf_format <- match.arg(conf_format, c("parentheses", "brackets", "columns"))
        cred_format <- match.arg(cred_format, c("parentheses", "brackets", "columns"))

        # Get the requested meta-analyses
        ma_obj <- filter_ma(ma_obj = ma_obj,
                            analyses = analyses,
                            match = match,
                            case_sensitive = case_sensitive)

        # Get summary.ma_psychmeta object
        if(!"summary.ma_psychmeta" %in% class(ma_obj)) ma_obj <- summary(ma_obj)

        ma_metric <- ma_obj$ma_metric
        if(ma_metric %in% c("r_order2", "d_order2")) stop("metabulate does not currently support second-order meta-analyses.", call=FALSE)
        ma_methods <- ma_obj$ma_methods
        es_type <- case_when(
                ma_metric %in% c("r_as_r", "d_as_r") ~ "r",
                ma_metric %in% c("r_as_d", "d_as_d") ~ "d",
                TRUE ~ "generic"
        )
        ma_method <- ma_method[ma_method %in% ma_methods]
        if(length(ma_method) == 0) {
                ma_method <- case_when(
                        "ad" %in% ma_methods ~ "ad",
                        "ic" %in% ma_methods ~ "ic",
                        TRUE ~ "bb"
                )
        }
        ma_type <- c(ma_method[ma_method == "bb"],
                     paste(rep(ma_method[ma_method != "bb"], each = length(correction_type)),
                           correction_type, sep = "_"))

        meta_tables <- ma_obj$meta_tables
        conf_level <- attributes(ma_obj$ma_obj)$inputs$conf_level
        cred_level <- attributes(ma_obj$ma_obj)$inputs$cred_level

        # Generate formatted tables (includes generation of captions)
        meta_tables <- .metabulate(meta_tables = meta_tables, ma_type = ma_type,
                                   show_msd = show_msd, show_conf = show_conf,
                                   show_cred = show_cred, show_se = show_se,
                                   show_var = show_var, es_type = es_type, symbol_es = symbol_es,
                                   digits = digits, decimal.mark = decimal.mark,
                                   leading0 = leading0, neg.sign = neg.sign,
                                   pos.sign = pos.sign, drop0integer = drop0integer,
                                   big.mark = big.mark, big.interval = big.interval,
                                   small.mark = small.mark, small.interval = small.interval,
                                   conf_format = conf_format, cred_format = cred_format,
                                   verbose = verbose, conf_level = conf_level, cred_level = cred_level)

        # Set the output file name
        if(is.null(file)) file <- "psychmeta_output"
        if(file != "console") file <- .filename_check(file, output_format)

        # Assign values to citekeys and citations, convert bib to R bibliography
        if(!is.null(bib)) bib <- .generate_bib(ma_obj, bib, additional_citekeys)

        # Render the output
        .psychmeta_render(file = file, output_format = output_format, meta_tables = meta_tables,
                          ma_type = ma_type, es_type = es_type,
                          bib = bib$bib, citations = bib$citations, citekeys = bib$citekeys,
                          title.bib = title.bib, style = style, caption = caption,
                          save_build_files = save_build_files, header = header)

}



#' Generate a list of references included in meta-analyses
#'
#' This function generates a list of studies contributing to a meta-analysis
#'
#' @param ma_obj A psychmeta meta-analysis object with \code{citekeys} supplied.
#' @param bib A BibTeX file containing the citekeys for the meta-analyses.
#' @param additional_citekeys Additional citekeys to include in the reference list.
#' @param analyses Which analyses to extract references for? See \code{\link{filter_ma}} for details.
#' @param match Match \code{all} or \code{any} of the filter criteria? See \code{\link{filter_ma}} for details.
#' @param case_sensitive Logical scalar that determines whether character values supplied in \code{analyses} should be treated as case sensitive (\code{TRUE}, default) or not (\code{FALSE}).
#' @param style What style should references be formatted in? Can be a file path or URL for a \url{https://github.com/citation-style-language/styles}{CSL citation style} or the style ID for any style available from the \url{https://zotero.org/styles}{Zotero Style Repository}). Defaults to APA style. (Retrieving a style by ID requires an internet connection. If unavailable, references will be rendered in Chicago style.).
#' @param output_format The format of the output reference list. Available options are Word (default), HTML, PDF (requires LaTeX, see the \code{tinytex} package), or Rmarkdown, plain text, and BibLaTeX. Returning only the item citekeys is also possible.
#' @param file The filename or filepath for the output file. If \code{NULL}, file will be saved as \code{reference_list}. Set to \code{"console"} or \code{"print"} to output directly to the R console.
#' @param title.bib The title to give to the bibliography. If \code{NULL}, defaults to "Sources Contributing to Meta-Analyses"
#' @param save_build_files Should the BibTeX and RMarkdown files used to generate the bibliography be saved (default: \code{TRUE})?
#' @param header A list of YAML header parameters to pass to \code{link{rmarkdown::render}}.
#'
#' @return A formatted reference list.
#' @export
#'
#' @importFrom rmarkdown render
#' @importFrom dplyr case_when
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
generate_bib <- function(ma_obj=NULL, bib=NULL, additional_citekeys=NULL,
                         analyses="all", match=c("all", "any"), case_sensitive = TRUE,
                         style="apa", output_format=c("word", "html", "pdf", "text", "rmd", "biblatex", "citekeys"),
                         file=NULL, title.bib = NULL, save_build_files = TRUE, header=list()){



        output_format <- match.arg(tolower(output_format))

        if("summary.ma_psychmeta" %in% class(ma_obj)) ma_obj <- ma_obj$ma_obj

        if(class(ma_obj))

        # Get the requested meta-analyses
        ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive)

        # Assign values to citekeys and citations, convert bib to R bibliography
        bib <- .generate_bib(ma_obj, bib, additional_citekeys)
        citekeys <- bib$citekeys
        citations <- bib$citations
        bib <- bib$bib

        # Set the output file name
        if(is.null(file)) file <- "reference_list"
        if(file != "console") file <- .filename_check(file, output_format)

        # Render the output
        switch(output_format,

               citekeys = {
                       if(file == "console") {
                               return(citations)
                       } else {
                               writeLines(citations, con = file)
                       }
               },

               biblatex = {
                       if(file == "console") {
                               print(bib[citekeys], .opts = list(style = "Biblatex"))
                       } else {
                               suppressMessages(WriteBib(bib[citekeys], file = file))
                       }
               },
               # else =
                        {
                        .psychmeta_render(file = file, output_format = output_format,
                                          bib = bib, citations = citations,
                                          citekeys = citekeys, title.bib = title.bib, style = style,
                                          save_build_files = save_build_files, header = header)
                }
               )
}


.generate_bib <- function(ma_obj, bib, additional_citekeys){

        # TODO: Make RefManageR, RCurl, rmarkdown optional packages

        if(is.null(ma_obj) & is.null(additional_citekeys)) stop("Either ma_obj or additional_citekeys must be provided.")

        # Compile unique citekeys from meta-analyses and additionally supplied list
        citekeys <-
             unique(c(additional_citekeys,
                      unlist(map(get_metafor(ma_obj),
                                 ~ strsplit(as.character(.x$barebones$citekey), ", ")))
                      ))

        # Render citekeys as Markdown citations
        citations <- paste0("@", citekeys, collapse=", ")

        # Set BibOptions to accept entries with missing data
        Bib_check.entries_original <- BibOptions()[["check.entries"]]
        BibOptions(check.entries = FALSE)

        # Read in .bib file
        bib <- ReadBib(bib)

        # Reset BibOptions to original
        BibOptions(check.entries = Bib_check.entries_original)

        return(list(citekeys = citekeys, citations = citations, bib = bib))

}

.filename_check <- function(file, output_format){
        case_when(
                output_format == "pdf"      ~ if(!grepl("\\.pdf$",  file, ignore.case = TRUE)) paste0(file, ".pdf")  else file,
                output_format == "html"     ~ if(!grepl("\\.html$", file, ignore.case = TRUE)) paste0(file, ".html") else file,
                output_format == "word"     ~ if(!grepl("\\.docx$", file, ignore.case = TRUE)) paste0(file, ".docx") else file,
                output_format == "rmd"      ~ if(!grepl("\\.rmd$",  file, ignore.case = TRUE)) paste0(file, ".Rmd")  else file,
                output_format == "biblatex" ~ if(!grepl("\\.bib$",  file, ignore.case = TRUE)) paste0(file, ".bib")  else file,
                output_format == "text"     ~ if(!grepl("\\.txt$",  file, ignore.case = TRUE)) paste0(file, ".txt")  else file,
                output_format == "citekeys" ~ if(!grepl("\\.txt$",  file, ignore.case = TRUE)) paste0(file, ".txt")  else file,
                TRUE ~ file)
}

.clean_style_name <- function(style) {
        if(grepl("://", style)) {
                attributes(style) <- list(source = "url")
                return(style)
        } else if(grepl("(~|:/|:\\\\)", style)) {
                attributes(style) <- list(source = "local")
                return(style)
        } else {
                style <- paste0("https://zotero.org/styles/", str_replace(style, "\\.csl$", "") )
                attributes(style) <- list(source = "Zotero")
                return(style)
        }
}

.psychmeta_render <- function(file, output_format, meta_tables = NULL, ma_type = NULL, es_type = NULL,
                              bib = NULL, citations = NULL, citekeys = NULL, title.bib = NULL,
                              style = style, save_build_files = FALSE, header = list(), caption = NULL){

        if(output_format == "rmd") save_build_files <- TRUE
        if(!is.null(style)) style <- .clean_style_name(style)

        if(file == "console") {

                switch(output_format,

                       text = {if(!is.null(meta_tables)) {
                                 for(i in names(meta_tables)) {
                                   if(!is.null(caption)) {
                                     if(length(caption) > 1) {
                                       cat("\n\n", caption[[i]], "\n", rep("=", nchar(caption[[i]])))
                                     } else cat(caption, "\n", rep("=", nchar(caption)))
                                   }
                                   cat("\n")
                                   print(meta_tables[[i]])
                                   cat("\n", attr(meta_tables, "footnotes")[[i]], "\n")
                                 }
                               }

                               if(!is.null(bib)) {
                                       if(is.null(title.bib)) title.bib <- "Sources Contributing to Meta-Analyses"
                                       cat(rep("\n", 2*as.numeric(is.null(meta_tables))),
                                           title.bib, "\n",
                                           rep("=", nchar(title.bib)), "\n\n"
                                           )
                                       # TODO: Replace this with a call to citation.js to use CSL styles
                                       print(bib[citekeys], .opts = list(style = "text", bib.style = "authoryear"))
                               }
                       },

                       # else =
                              {# TODO: If the bug with `bibliography: ` YAML metadata gets fixed, move this line to
                               # the same metadata block as citations and style below.
                               if(!is.null(bib)) {
                                       # Write the bibliography file
                                       if(save_build_files) {
                                               bib_file <- str_replace(file, "\\.(Rmd|pdf|docx|html)$", "\\.bib")
                                       } else bib_file <- tempfile("psychmeta.bib")
                                       suppressMessages(WriteBib(bib[citekeys],
                                                                 file = bib_file))

                                       sprintf("---### This metadata line must be placed in your RMarkdown document main YAML header! ###\nbibliography: %s\n---\n",
                                               bib_file)
                               }

                                if(!is.null(meta_tables)) {
                                  for(i in names(meta_tables)) {
                                    cat("\n\n")
                                    knitr::kable(meta_tables[[i]], caption = if(length(caption) > 1) caption[[i]] else caption)
                                    cat("\n\n")
                                    cat("\n", attr(meta_tables, "footnotes")[[i]], "\n")
                                  }

                                }

                               if(!is.null(bib)) {
                                       if(is.null(title.bib)) title.bib <- "# Sources Contributing to Meta-Analyses"

                                       cat(rep("\n", 2*as.numeric(is.null(meta_tables))),
                                           title.bib,
                                           "\n\n---\n"
                                           )
                                       if(!is.null(style)) sprintf("csl: %s\n", style )
                                       sprintf("nocite: |\n  %s\n---\n", citations)
                               }
                       }
                )

        } else {

                switch(output_format,

                       text = {sink(file = file)
                               if(!is.null(meta_tables)) {
                                 for(i in names(meta_tables)) {
                                   if(!is.null(caption)) {
                                     if(length(caption) > 1) {
                                       cat("\n\n", caption[[i]], "\n", rep("=", nchar(caption[[i]])))
                                     } else cat(caption, "\n", rep("=", nchar(caption)))
                                   }
                                   cat("\n")
                                   print(meta_tables[[i]])
                                   cat("\n", attr(meta_tables, "footnotes")[[i]], "\n")
                                 }
                               }

                               if(!is.null(bib)) {
                                       if(is.null(title.bib)) title.bib <- "Sources Contributing to Meta-Analyses"
                                       cat(rep("\n", 2*as.numeric(is.null(meta_tables))),
                                           title.bib, "\n",
                                           rep("=", nchar(title.bib)), "\n\n"
                                       )
                                       # TODO: Replace this with a call to citation.js to use CSL styles
                                       print(bib[citekeys], .opts = list(style = "text", bib.style = "authoryear"))
                               }
                               sink()
                       },

                       # else =
                              {# Fill in critical header slots and write .bib file if necessary
                                if(is.null(header$title)) {
                                        if(is.null(meta_tables)) {
                                                if(is.null(title.bib)) {
                                                        header$title <- "Sources Contributing to Meta-Analyses"
                                                } else {
                                                        header$title <- title.bib
                                                        title.bib <- NULL
                                                }
                                        } else {
                                                header$title <- "Results of Meta-Analyses"
                                                if(is.null(title.bib)) {
                                                        title.bib <- "\n\n# Sources Contributing to Meta-Analyses\n\n"
                                                } else {
                                                        title.bib <- sprintf("\n\n# %s\n\n", title.bib)
                                                }
                                        }
                                }

                                header$output <- paste0(output_format, "_document")

                                if(!is.null(bib)) {
                                        # Write the bibliography file
                                        if(save_build_files) {
                                                bib_file <- str_replace(file, "\\.(Rmd|pdf|docx|html)$", "\\.bib")
                                        } else bib_file <- tempfile("psychmeta.bib")
                                        suppressMessages(WriteBib(bib[citekeys],
                                                                  file = bib_file))

                                        header$bibliography <- bib_file

                                        if(!is.null(style)) {
                                                if(attr(style, "source") %in% c("url", "Zotero")) {
                                                        if(url.exists(style)) {
                                                                header$csl <- style
                                                        } else {
                                                                message(sprintf("Caution: Style not found at %s\n         Check the %s or specify a local CSL style file.\n         References formatted using to the Chicago Manual of Style.",
                                                                                style,
                                                                                if(attr(style, "source") == "url") "URL" else "style name"))
                                                        }
                                                } else {
                                                        if(file.exists(style)) {
                                                                header$csl <- style
                                                        } else {
                                                                message(sprintf("Caution: Style not found at %s\n         Check the file path or specify a CSL style name from the Zotero Style Repository (https://zotero.org/styles).\n         References formatted using to the Chicago Manual of Style.",
                                                                                style))
                                                        }
                                                }
                                        }
                                }

                               ## Create the markdown header and document
                               header <- paste(names(header), header, sep=": ", collapse="\n")
                               mt_output <- character(0)
                               for(i in names(meta_tables)) {
                                       if(is.character(meta_tables[[i]])) enc2native(meta_tables[[i]])
                                       mt_output <-
                                               paste0(mt_output, "\n",
                                                      paste0(knitr::kable(meta_tables[[i]], caption = {if(length(caption) > 1) caption[[i]] else caption}), collapse="\n"),
                                                      "\n",
                                                      paste(attr(meta_tables, "footnotes")[[i]]), sep = "\n")
                                       }

                               document <- paste0("---\n", header, "\n---\n\n",
                                                  mt_output, "\n\n",
                                                  paste0("# ", title.bib, "\n\n")[!is.null(title.bib) & !is.null(bib)],
                                                  sprintf('---\nnocite: |\n  %s', citations)[!is.null(bib)]
                               )
                               document <- gsub("<U\\+(.{1,4})>", "\\\\u{\\1}", x = document)

                               # Create Rmd and output files
                               if(save_build_files) {
                                       rmd_document <- str_replace(file, "\\.(pdf|docx|html)$", "\\.rmd")

                               } else rmd_document <- tempfile("psychmeta.rmd")

                               stringi::stri_write_lines(document, rmd_document)

                               if(output_format != "rmd") {
                                       render(rmd_document,
                                              output_file = file,
                                              output_dir  = getwd(),
                                              encoding = "UTF-8")
                               }
                       }
                )
        }
}


#' Internal function for .metabulating results tables
#'
#' @keywords internal
.metabulate <- function(meta_tables, ma_type = "ad_ts", show_msd = TRUE, show_conf = TRUE, show_cred = TRUE,
                        show_se = FALSE, show_var = FALSE, es_type = NULL, symbol_es = "ES",
                        digits = 2L, decimal.mark = getOption("OutDec"), leading0 = "figure", neg.sign = "\u2212",
                        pos.sign = "figure", drop0integer = TRUE, big.mark = "\u202F",
                        big.interval = 3L, small.mark = "\u202F", small.interval = 3L,
                        conf_format = "parentheses", cred_format = "brackets",
                        verbose = FALSE, conf_level = .95, cred_level = .80) {

        if(es_type == "r") {
                meta_tables <- list(bb = meta_tables$barebones,
                                  ic_ts  = meta_tables$individual_correction$true_score,
                                  ic_vgx = meta_tables$individual_correction$validity_generalization_x,
                                  ic_vgy = meta_tables$individual_correction$validity_generalization_y,
                                  ad_ts  = meta_tables$artifact_distribution$true_score,
                                  ad_vgx = meta_tables$artifact_distribution$validity_generalization_x,
                                  ad_vgy = meta_tables$artifact_distribution$validity_generalization_y
                )[ma_type]
        } else if(es_type == "d") {
                meta_tables <- list(bb = meta_tables$barebones,
                                  ic_ts  = meta_tables$individual_correction$latentGroup_latentY,
                                  ic_vgx = meta_tables$individual_correction$observedGroup_latentY,
                                  ic_vgy = meta_tables$individual_correction$latentGroup_observedY,
                                  ad_ts  = meta_tables$artifact_distribution$latentGroup_latentY,
                                  ad_vgx = meta_tables$artifact_distribution$observedGroup_latentY,
                                  ad_vgy = meta_tables$artifact_distribution$latentGroup_observedY
                )[ma_type]
        } else {
                meta_tables <- list(bb = meta_tables$barebones)
        }

        length_intial <- max(sapply(meta_tables, function(x) length(colnames(x)[1:which(colnames(x) == "analysis_id")][colnames(x)[1:which(colnames(x) == "analysis_id")] %in% c("analysis_id", "pair_id", "analysis_type")])))
        length_moderators <- max(sapply(meta_tables, function(x) if(1 + which(colnames(x) == "analysis_type") != which(colnames(x) == "k")) length(colnames(x)[(1 + which(colnames(x) == "analysis_type")):(which(colnames(x) == "k") -1)]) else 0))

        # Select, rearrange, and format columns of meta_tables
        .arrange_format_columns <- function(ma_table) {
                ma_table <- mutate(ma_table,
                                   spacer1 = " ",
                                   spacer2 = " ",
                                   spacer3 = " ",
                                   spacer4 = " ",
                                   spacer5 = " ")
                x <- colnames(ma_table)

                # Select columns
                col_initial <- x[1:which(x == "analysis_type")]
                col_initial <- col_initial[!col_initial %in% c("analysis_id", "pair_id", "analysis_type")]

                col_moderators <- if(1 + which(x == "analysis_type") != which(x == "k")) x[(1 + which(x == "analysis_type")):(which(x == "k") - 1)] else NULL

                col_sampsize <- x[which(x %in% c("k", "N"))]

                if(show_msd == TRUE) {

                        col_m_bb   <- x[which(x %in% c("mean_r", "mean_d", "mean_es"))]
                        col_sd_bb  <- if(verbose == TRUE) {
                                x[which(x %in% c("sd_r", "sd_d", "sd_es", "sd_e", "sd_art", "sd_pre", "sd_res"))]
                        } else col_sd_bb <- x[which(x %in% c("sd_r", "sd_d", "sd_es", "sd_res"))]

                        col_m_cor  <- x[which(x %in% c("mean_rho", "mean_delta"))]
                        col_sd_cor <- if(verbose == TRUE) {
                                x[which(x %in% c("sd_r_c", "sd_d_c", "sd_e_c", "sd_art_c", "sd_pre_c", "sd_rho", "sd_delta"))]
                        } else col_sd_cor <- x[which(x %in% c("sd_r_c", "sd_d_c", "sd_rho", "sd_delta"))]

                } else {
                        col_m_bb   <- NULL
                        col_sd_bb  <- NULL
                        col_m_cor  <- NULL
                        col_sd_cor <- NULL
                }

                if(show_se == TRUE) {
                        col_se_bb  <- x[which(x %in% c("se_r", "se_d", "se_es"))]
                        col_se_cor <- x[which(x %in% c("se_r_c", "se_d_c"))]

                } else {
                        col_se_bb  <- NULL
                        col_se_cor <- NULL
                }
                if(show_var == TRUE) {

                        col_var_bb  <- if(verbose == TRUE) {
                                x[which(x %in% c("var_r", "var_d", "var_es", "var_e", "var_art", "var_pre", "var_res"))]
                        } else col_var_bb <- x[which(x %in% c("var_r", "var_d", "var_es", "var_res"))]
                        col_var_cor <- if(verbose == TRUE) {
                                x[which(x %in% c("var_r_c", "var_d_c", "var_e_c", "var_art_c", "var_pre_c", "var_rho", "var_delta"))]
                        } else col_var_cor <- x[which(x %in% c("var_r_c", "var_d_c", "var_rho", "var_delta"))]

                } else {
                        col_var_bb <- NULL
                        col_var_cor <- NULL
                }

                if(show_conf == TRUE) col_conf <- grep("^CI_.+", x, value = TRUE) else col_conf <- NULL
                if(show_cred == TRUE) col_cred <- grep("^CV_.+", x, value = TRUE) else col_cred <- NULL

                # Rearrange columns
                ma_table <- ma_table[,c(col_initial, col_moderators, col_sampsize,
                                        col_m_bb, col_se_bb, col_sd_bb,
                                        "spacer1"[show_var], col_var_bb,
                                        "spacer2"[!is.null(col_m_cor)], col_m_cor, col_se_cor, col_sd_cor,
                                        "spacer3"[!is.null(col_m_cor) & show_var], col_var_cor,
                                        "spacer4"[show_conf], col_conf,
                                        "spacer5"[show_cred], col_cred)]

                colnames(ma_table)[colnames(ma_table) %in% col_conf] <- c("ci_lower", "ci_upper")
                if(show_conf == TRUE) col_conf <- c("ci_lower", "ci_upper")
                colnames(ma_table)[colnames(ma_table) %in% col_cred] <- c("cv_lower", "cv_upper")
                if(show_cred == TRUE) col_cred <- c("cv_lower", "cv_upper")

                # Format columns
                ma_table[, col_sampsize] <- num_format(ma_table[, col_sampsize], digits = 0L, decimal.mark = decimal.mark,
                                                       leading0 = leading0, neg.sign = neg.sign, pos.sign = pos.sign,
                                                       drop0integer = drop0integer, big.mark = big.mark, big.interval = big.interval,
                                                       small.mark = small.mark, small.interval = small.interval)

                numeric_columns <- c(col_m_bb, col_se_bb, col_sd_bb, col_var_bb, col_m_cor, col_se_cor, col_sd_cor, col_var_cor, col_conf, col_cred)

                ma_table[, numeric_columns] <-
                        num_format(ma_table[, numeric_columns], digits = digits, decimal.mark = decimal.mark,
                                                       leading0 = leading0, neg.sign = neg.sign, pos.sign = pos.sign,
                                                       drop0integer = drop0integer, big.mark = big.mark, big.interval = big.interval,
                                                       small.mark = small.mark, small.interval = small.interval)

                # Format the interval columns
                if(show_conf == TRUE) {
                        switch(conf_format,
                               parentheses = {ma_table <- rename(select(mutate(ma_table, ci_lower = paste0("(", .data$ci_lower, ", ", .data$ci_upper, ")")),
                                                                        -.data$ci_upper), conf_int = .data$ci_lower)},
                               brackets = {ma_table <- rename(select(mutate(ma_table, ci_lower = paste0("[", .data$ci_lower, ", ", .data$ci_upper, "]")),
                                                                     -.data$ci_upper), conf_int = .data$ci_lower)},
                               )
                }
                if(show_cred == TRUE) {
                        switch(cred_format,
                               parentheses = {ma_table <- rename(select(mutate(ma_table, cv_lower = paste0("(", .data$cv_lower, ", ", .data$cv_upper, ")")),
                                                                        -.data$cv_upper), cred_int = .data$cv_lower)},
                               brackets = {ma_table <- rename(select(mutate(ma_table, cv_lower = paste0("[", .data$cv_lower, ", ", .data$cv_upper, "]")),
                                                                     -.data$cv_upper), cred_int = .data$cv_lower)},
                        )
                }

                return(ma_table)
        }

        # Rename columns
        .rename_columns <- function(ma_table, symbol_es, conf_level, cred_level) {
          name_map <- c(
            group_contrast    = "**Group Contrast**",
            construct_x       = "**Construct X**",
            construct_y       = "**Construct Y**",
            k                 = "**_k_**",
            N                 = "**_N_**",

            mean_r            = "**$\\overline{r}$**",
            var_r             = "**$\\sigma^{2}_{r}$**",
            sd_r              = "**$SD_{r}$**",
            se_r              = "**$SE_{\\overline{r}}$**",

            mean_d            = "**$\\overline{d}$",
            var_d             = "**$\\sigma^{2}_{d}$",
            sd_d              = "**$SD_{d}$",
            se_d              = "**$SE_{\\overline{d}}$",

            mean_es           = paste0("**$\\overline{", symbol_es, "}$**"),
            var_es            = paste0("**$\\sigma^{2}_{", symbol_es, "}$**"),
            sd_es             = paste0("**$SD{", symbol_es, "}$"),
            se_es             = paste0("**$SE_{\\overline{", symbol_es, "}}$**"),

            var_e             = "**$\\sigma^{2}_{e}$**",
            var_res           = "**$\\sigma^{2}_{res}$**",
            sd_e              = "**$SD_{e}$**",
            sd_res            = "**$SD_{res}$**",
            var_art           = "**$\\sigma^{2}_{art}$**",
            var_pre           = "**$\\sigma^{2}_{pre}$**",
            sd_art            = "**$SD_{art}$**",
            sd_pre            = "**$SD_{pre}$**",

            mean_rho          = "**$\\overline{\\mathrm{\\rho}}$**",
            var_r_c           = "**$\\sigma^{2}_{r_{c}}$**",

            var_rho           = "**$\\sigma^{2}_{\\mathrm{\\rho}}$**",
            sd_r_c            = "**$SD_{r_{c}}$**",
            se_r_c            = "**$SE_{\\overline{\\mathrm{\\rho}}}**$",
            sd_rho            = "**$SD_{\\mathrm{\\rho}}$**",

            mean_delta        = "**$\\overline{\\mathrm{\\delta}}$**",
            var_d_c           = "**$\\sigma^{2}_{d_{c}}$**",
            sd_d_c            = "**$SD_{d_{c}}$**",
            se_d_c            = "**$SE_{\\overline{\\mathrm{\\delta}}}$**",
            var_delta         = "**$\\sigma^{2}_{\\mathrm{\\delta}}$**",
            sd_delta          = "**$SD_{\\mathrm{\\delta}}$**",

            var_e_c           = "**$\\sigma^{2}_{e_{c}}$**",
            var_art_c         = "**$\\sigma^{2}_{art_{c}}$**",
            var_pre_c         = "**$\\sigma^{2}_{pre_{c}}$**",
            sd_e_c            = "**$SD_{e_{c}}$**",
            sd_art_c          = "**$SD_{art_{c}}$**",
            sd_pre_c          = "**$SD_{pre_{c}}$**",

            ci_lower          = paste0("**", conf_level*100, "% conf. int.**"),
            ci_upper          = " ",
            cv_lower          = paste0("**", cred_level*100, "% cred. int.**"),
            cv_upper          = " ",

            conf_int          = paste0("**", conf_level*100, "% conf. int.**"),
            cred_int          = paste0("**", cred_level*100, "% cred. int.**")
          )

          formatted_names <- names(name_map)
          names(formatted_names) <- name_map

          ma_table <- rename(ma_table, !!formatted_names[formatted_names %in% colnames(ma_table)])

          colnames(ma_table)[grep("spacer\\d+", colnames(ma_table))] <- " "

          ma_table

        }

        .format_meta_table <- function(ma_table, symbol_es, conf_level, cred_level) {
          ma_table <- ma_table %>%
            .arrange_format_columns() %>%
            .rename_columns(symbol_es, conf_level, cred_level)
          ma_table
        }

        meta_tables <- lapply(meta_tables, .format_meta_table, symbol_es = symbol_es, conf_level = conf_level, cred_level = cred_level)

        # Table footnotes
        footnote_list <-
          ### TODO: Add notes about actual corrections applied
          if(es_type == "r") {
            c(bb     = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{r}$\u00a0=\u00a0 mean observed correlation; "[show_msd],
                                 "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],
                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{r}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{r}$."[show_cred]),

                 ic_ts  = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{r}$\u00a0=\u00a0mean observed correlation; "[show_msd],
                                 "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{r}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{r}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $r$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],

                                 "$\\overline{\\mathrm{\\rho}}$\u00a0=\u00a0mean true-score correlation; "[show_msd],
                                 "$SE_{\\overline{\\mathrm{\\rho}}}$\u00a0=\u00a0standard error of $\\overline{\\mathrm{\\rho}}$; "[show_se],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{\\mathrm{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\rho}$; "[show_msd & !verbose],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{r_{c}}$; $SD_{\\mathrm{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\rho}$; "[show_msd & verbose],

                                 "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{\\mathrm{\\rho}}$\u00a0=\u00a0residual variance of $\\mathrm{\\rho}$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $r_{c}$; $\\sigma^{2}_{\\mathrm{\\rho}}$\u00a0=\u00a0residual variance of $\\mathrm{\\rho}$; "[show_var & verbose],

                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{\\mathrm{\\rho}}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{\\mathrm{\\rho}}$; "[show_cred],

                                 "correlations corrected individually."),

                 ic_vgx = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{r}$\u00a0=\u00a0mean observed correlation; "[show_msd],
                                 "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{r}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{r}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $r$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],

                                 "$\\overline{\\mathrm{\\rho}}$\u00a0=\u00a0mean operational validity (*X* measured with error); "[show_msd],
                                 "$SE_{\\overline{\\mathrm{\\rho}}}$\u00a0=\u00a0standard error of $\\overline{\\mathrm{\\rho}}$; "[show_se],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{\\mathrm{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\rho}$; "[show_msd & !verbose],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{r_{c}}$; $SD_{\\mathrm{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\rho}$; "[show_msd & verbose],

                                 "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{\\mathrm{\\rho}}$\u00a0=\u00a0residual variance of $\\mathrm{\\rho}$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $r_{c}$; $\\sigma^{2}_{\\mathrm{\\rho}}$\u00a0=\u00a0residual variance of $\\mathrm{\\rho}$; "[show_var & verbose],

                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{\\mathrm{\\rho}}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{\\mathrm{\\rho}}$; "[show_cred],

                                 "correlations corrected individually."),

                 ic_vgy = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{r}$\u00a0=\u00a0mean observed correlation; "[show_msd],
                                 "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{r}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{r}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $r$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],

                                 "$\\overline{\\mathrm{\\rho}}$\u00a0=\u00a0mean operational validity (*Y* measured with error); "[show_msd],
                                 "$SE_{\\overline{\\mathrm{\\rho}}}$\u00a0=\u00a0standard error of $\\overline{\\mathrm{\\rho}}$; "[show_se],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{\\mathrm{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\rho}$; "[show_msd & !verbose],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{r_{c}}$; $SD_{\\mathrm{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\rho}$; "[show_msd & verbose],

                                 "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{\\mathrm{\\rho}}$\u00a0=\u00a0residual variance of $\\mathrm{\\rho}$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $r_{c}$; $\\sigma^{2}_{\\mathrm{\\rho}}$\u00a0=\u00a0residual variance of $\\mathrm{\\rho}$; "[show_var & verbose],

                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{\\mathrm{\\rho}}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{\\mathrm{\\rho}}$; "[show_cred],

                                 "correlations corrected individually."),

                 ad_ts  = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{r}$\u00a0=\u00a0mean observed correlation; "[show_msd],
                                 "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{r}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{r}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $r$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],

                                 "$\\overline{\\mathrm{\\rho}}$\u00a0=\u00a0mean true-score correlation; "[show_msd],
                                 "$SE_{\\overline{\\mathrm{\\rho}}}$\u00a0=\u00a0standard error of $\\overline{\\mathrm{\\rho}}$; "[show_se],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{\\mathrm{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\rho}$; "[show_msd & !verbose],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{r_{c}}$; $SD_{\\mathrm{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\rho}$; "[show_msd & verbose],

                                 "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{\\mathrm{\\rho}}$\u00a0=\u00a0residual variance of $\\mathrm{\\rho}$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $r_{c}$; $\\sigma^{2}_{\\mathrm{\\rho}}$\u00a0=\u00a0residual variance of $\\mathrm{\\rho}$; "[show_var & verbose],

                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{\\mathrm{\\rho}}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{\\mathrm{\\rho}}$; "[show_cred],

                                 "correlations corrected using artifact distributions."),

                 ad_vgx = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{r}$\u00a0=\u00a0mean observed correlation; "[show_msd],
                                 "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{r}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{r}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $r$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],

                                 "$\\overline{\\mathrm{\\rho}}$\u00a0=\u00a0mean operational validity (*X* measured with error); "[show_msd],
                                 "$SE_{\\overline{\\mathrm{\\rho}}}$\u00a0=\u00a0standard error of $\\overline{\\mathrm{\\rho}}$; "[show_se],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{\\mathrm{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\rho}$; "[show_msd & !verbose],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{r_{c}}$; $SD_{\\mathrm{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\rho}$; "[show_msd & verbose],

                                 "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{\\mathrm{\\rho}}$\u00a0=\u00a0residual variance of $\\mathrm{\\rho}$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $r_{c}$; $\\sigma^{2}_{\\mathrm{\\rho}}$\u00a0=\u00a0residual variance of $\\mathrm{\\rho}$; "[show_var & verbose],

                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{\\mathrm{\\rho}}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{\\mathrm{\\rho}}$; "[show_cred],

                                 "correlations corrected using artifact distributions."),

                 ad_vgy = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{r}$\u00a0=\u00a0mean observed correlation; "[show_msd],
                                 "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                 "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{r}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{r}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $r$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],

                                 "$\\overline{\\mathrm{\\rho}}$\u00a0=\u00a0mean operational validity (*Y* measured with error); "[show_msd],
                                 "$SE_{\\overline{\\mathrm{\\rho}}}$\u00a0=\u00a0standard error of $\\overline{\\mathrm{\\rho}}$; "[show_se],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{\\mathrm{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\rho}$; "[show_msd & !verbose],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{r_{c}}$; $SD_{\\mathrm{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\rho}$; "[show_msd & verbose],

                                 "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{\\mathrm{\\rho}}$\u00a0=\u00a0residual variance of $\\mathrm{\\rho}$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $r_{c}$; $\\sigma^{2}_{\\mathrm{\\rho}}$\u00a0=\u00a0residual variance of $\\mathrm{\\rho}$; "[show_var & verbose],

                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{\\mathrm{\\rho}}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{\\mathrm{\\rho}}$; "[show_cred],

                                 "correlations corrected using artifact distributions.")
            )

          } else if(es_type == "d") {
            ### TODO: Don't refer to latent/observed groups if group membership reliability is not corrected.
            c(bb     = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{d}$\u00a0=\u00a0 mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                 "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],
                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{d}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{d}$."[show_cred]),

                 ic_ts  = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{d}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                 "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{d}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{d}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $d$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],

                                 "$\\overline{\\mathrm{\\delta}}$\u00a0=\u00a0mean true-score Cohen's $d$ (Hedges' $g$) between latent groups; "[show_msd],
                                 "$SE_{\\overline{\\mathrm{\\delta}}}$\u00a0=\u00a0standard error of $\\overline{\\mathrm{\\delta}}$; "[show_se],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{\\mathrm{\\delta}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\delta}$; "[show_msd & !verbose],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{d_{c}}$; $SD_{\\mathrm{\\delta}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\delta}$; "[show_msd & verbose],

                                 "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{\\mathrm{\\delta}}$\u00a0=\u00a0residual variance of $\\mathrm{\\delta}$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $d_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $d{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $d_{c}$; $\\sigma^{2}_{\\mathrm{\\delta}}$\u00a0=\u00a0residual variance of $\\mathrm{\\delta}$; "[show_var & verbose],

                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{\\mathrm{\\delta}}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{\\mathrm{\\delta}}$; "[show_cred],

                                 "effect sizes corrected individually."),

                 ic_vgx = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{d}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                 "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{d}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{d}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $d$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],

                                 "$\\overline{\\mathrm{\\delta}}$\u00a0=\u00a0mean true-score Cohen's $d$ (Hedges' $g$) between observed groups; "[show_msd],
                                 "$SE_{\\overline{\\mathrm{\\delta}}}$\u00a0=\u00a0standard error of $\\overline{\\mathrm{\\delta}}$; "[show_se],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{\\mathrm{\\delta}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\delta}$; "[show_msd & !verbose],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{d_{c}}$; $SD_{\\mathrm{\\delta}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\delta}$; "[show_msd & verbose],

                                 "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{\\mathrm{\\delta}}$\u00a0=\u00a0residual variance of $\\mathrm{\\delta}$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $d_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $d{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $d_{c}$; $\\sigma^{2}_{\\mathrm{\\delta}}$\u00a0=\u00a0residual variance of $\\mathrm{\\delta}$; "[show_var & verbose],

                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{\\mathrm{\\delta}}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{\\mathrm{\\delta}}$; "[show_cred],

                                 "effect sizes corrected individually."),

                 ic_vgy = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{d}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                 "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{d}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{d}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $d$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],

                                 "$\\overline{\\mathrm{\\delta}}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$) between latent groups; "[show_msd],
                                 "$SE_{\\overline{\\mathrm{\\delta}}}$\u00a0=\u00a0standard error of $\\overline{\\mathrm{\\delta}}$; "[show_se],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{\\mathrm{\\delta}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\delta}$; "[show_msd & !verbose],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{d_{c}}$; $SD_{\\mathrm{\\delta}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\delta}$; "[show_msd & verbose],

                                 "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{\\mathrm{\\delta}}$\u00a0=\u00a0residual variance of $\\mathrm{\\delta}$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $d_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $d{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $d_{c}$; $\\sigma^{2}_{\\mathrm{\\delta}}$\u00a0=\u00a0residual variance of $\\mathrm{\\delta}$; "[show_var & verbose],

                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{\\mathrm{\\delta}}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{\\mathrm{\\delta}}$; "[show_cred],

                                 "effect sizes corrected individually."),

                 ad_ts  = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{d}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                 "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{d}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{d}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $d$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],

                                 "$\\overline{\\mathrm{\\delta}}$\u00a0=\u00a0mean true-score Cohen's $d$ (Hedges' $g$) between latent groups; "[show_msd],
                                 "$SE_{\\overline{\\mathrm{\\delta}}}$\u00a0=\u00a0standard error of $\\overline{\\mathrm{\\delta}}$; "[show_se],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{\\mathrm{\\delta}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\delta}$; "[show_msd & !verbose],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{d_{c}}$; $SD_{\\mathrm{\\delta}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\delta}$; "[show_msd & verbose],

                                 "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{\\mathrm{\\delta}}$\u00a0=\u00a0residual variance of $\\mathrm{\\delta}$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $d_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $d{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $d_{c}$; $\\sigma^{2}_{\\mathrm{\\delta}}$\u00a0=\u00a0residual variance of $\\mathrm{\\delta}$; "[show_var & verbose],

                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{\\mathrm{\\delta}}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{\\mathrm{\\delta}}$; "[show_cred],

                                 "effect sizes corrected using artifact distributions."),

                 ad_vgx = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{d}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                 "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{d}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{d}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $d$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],

                                 "$\\overline{\\mathrm{\\delta}}$\u00a0=\u00a0mean true-score Cohen's $d$ (Hedges' $g$) between observed groups; "[show_msd],
                                 "$SE_{\\overline{\\mathrm{\\delta}}}$\u00a0=\u00a0standard error of $\\overline{\\mathrm{\\delta}}$; "[show_se],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{\\mathrm{\\delta}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\delta}$; "[show_msd & !verbose],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{d_{c}}$; $SD_{\\mathrm{\\delta}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\delta}$; "[show_msd & verbose],

                                 "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{\\mathrm{\\delta}}$\u00a0=\u00a0residual variance of $\\mathrm{\\delta}$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $d_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $d{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $d_{c}$; $\\sigma^{2}_{\\mathrm{\\delta}}$\u00a0=\u00a0residual variance of $\\mathrm{\\delta}$; "[show_var & verbose],

                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{\\mathrm{\\delta}}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{\\mathrm{\\delta}}$; "[show_cred],

                                 "effect sizes corrected using artifact distributions"),

                 ad_vgy = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                 "$\\overline{d}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                 "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                 "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{d}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{d}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $d$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],

                                 "$\\overline{\\mathrm{\\delta}}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$) between latent groups; "[show_msd],
                                 "$SE_{\\overline{\\mathrm{\\delta}}}$\u00a0=\u00a0standard error of $\\overline{\\mathrm{\\delta}}$; "[show_se],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{\\mathrm{\\delta}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\delta}$; "[show_msd & !verbose],
                                 "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{d_{c}}$; $SD_{\\mathrm{\\delta}}$\u00a0=\u00a0residual standard deviation of $\\mathrm{\\delta}$; "[show_msd & verbose],

                                 "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{\\mathrm{\\delta}}$\u00a0=\u00a0residual variance of $\\mathrm{\\delta}$ ;"[show_var & !verbose],
                                 "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $d_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $d{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $d_{c}$; $\\sigma^{2}_{\\mathrm{\\delta}}$\u00a0=\u00a0residual variance of $\\mathrm{\\delta}$; "[show_var & verbose],

                                 "conf. int.\u00a0=\u00a0confidence interval around $\\overline{\\mathrm{\\delta}}$; "[show_conf],
                                 "cred. int.\u00a0=\u00a0credibility interval around $\\overline{\\mathrm{\\delta}}$; "[show_cred],

                                 "effect sizes corrected using artifact distributions")
            )

          } else {
            c(bb = paste0("*Note:* *k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                          "$\\overline{", symbol_es, "}$\u00a0=\u00a0 mean observed effect size ($", symbol_es, "$); "[show_msd],
                          "$SE_{\\overline{", symbol_es, "}}$\u00a0=\u00a0standard error of $\\overline{", symbol_es, "}$; "[show_se],
                          "$SD_{", symbol_es, "}$\u00a0=\u00a0observed standard deviation of $", symbol_es, "$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $", symbol_es, "$; "[show_msd & !verbose],
                          "$SD_{", symbol_es, "}$\u00a0=\u00a0observed standard deviation of $", symbol_es, "$; $SD_{e}$\u00a0=\u00a0predicted $SD_{", symbol_es, "}$ due to sampling error; $SD_{res}$\u00a0=\u00a0residual standard deviation of $", symbol_es, "$; "[show_msd & verbose],
                          "$\\sigma^{2}_{", symbol_es, "}$\u00a0=\u00a0observed variance of $", symbol_es, "$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $", symbol_es, "$ ;"[show_var & !verbose],
                          "$\\sigma^{2}_{", symbol_es, "}$\u00a0=\u00a0observed variance of $", symbol_es, "$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $", symbol_es, "$ due to sampling error; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $", symbol_es, "$; "[show_var & verbose],
                          "conf. int.\u00a0=\u00a0confidence interval around $\\overline{", symbol_es, "}$; "[show_conf],
                          "cred. int.\u00a0=\u00a0credibility interval around $\\overline{", symbol_es, "}$."[show_cred])
            )
          }

        attr(meta_tables, "footnotes") <- footnote_list[ma_type]

      return(meta_tables)

}
