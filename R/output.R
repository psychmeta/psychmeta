#' Format numbers for presentation
#'
#' A function to format numbers and logical values as characters for display purposes.
#' Includes control over formatting of decimal digits, leading zeros, sign characters,
#' and characters to replace logical, NA, NaN, and Inf values. Factors are converted
#' to strings. Strings are returned verbatim.
#'
#' @encoding UTF-8
#'
#' @param x A vector, matrix, or data.frame of numbers to format
#' @param digits The number of decimal digits desired (used strictly; default: 2)
#' @param decimal.mark The character to use for the decimal point (defaults to locale default: \code{getOption("OutDec")})
#' @param leading0 How to print leading zeros on decimals. Can be logical to print (\code{TRUE}) or suppress (\code{FALSE}) leading zeros or a character string to subsitute for leading zeros. If \code{"conditional"} (default), leading zeros are shown if a column contains any absolute values greater than 1 and suppressed otherwise. If \code{"figure"}, leading zeros are replaced with a figure space (\href{https://unicode-table.com/en/2007/}{\code{U+2007}}) if a column contains any absolute values greater than 1 and suppressed otherwise.
#' @param drop0integer Logical. Should trailing decimal zeros be dropped for integers?
#' @param neg.sign Character to use as negative sign. Defaults to minus-sign (\href{https://unicode-table.com/en/2212/}{\code{U+2212}}).
#' @param pos.sign Character to use as positive sign. Set to \code{FALSE} to suppress. If \code{"figure"} (default), the positive sign is a figure-space (\href{https://unicode-table.com/en/2007/}{\code{U+2007}}) if a column contains any negative numbers and suppressed otherwise.
#' @param big.mark Character to mark between each \code{big.interval} digits \emph{before} the decimal point. Set to \code{FALSE} to suppress. Defaults to the SI/ISO 31-0 standard-recommened thin-spaces (\href{https://unicode-table.com/en/202f/}{\code{U+202F}}).
#' @param big.interval See \code{big.mark} above; defaults to 3.
#' @param small.mark Character to mark between each \code{small.interval} digits \emph{after} the decimal point. Set to \code{FALSE} to suppress. Defaults to the SI/ISO 31-0 standard-recommened thin-spaces (\href{https://unicode-table.com/en/202f/}{\code{U+202F}}).
#' @param small.interval See \code{small.mark} above; defaults to 3.
#' @param na.mark Character to replace \code{NA} and \code{NaN} values. Defaults to em-dash (\href{https://unicode-table.com/en/2014/}{\code{U+2014}}))
#' @param lgl.mark A length 2 vector containing characters to replace \code{TRUE} and \code{FALSE}. Defaults to c("+", "\href{https://unicode-table.com/en/2212/}{\code{U+2212}}").
#' @param inf.mark A length 2 vector containing characters to replace \code{Inf} and \code{-Inf}. Defaults to c("+\href{https://unicode-table.com/en/221e/}{\code{U+221e}}", "\href{https://unicode-table.com/en/2212/}{\code{U+2212}}\href{https://unicode-table.com/en/221e/}{\code{U+221e}}").
#'
#' @export
#' @examples
#' # format_num() converts numeric values to characters with the specified formatting options.
#' # By default, thousands digit groups are separated by thin spaces, negative signs are replaced
#' # with minus signs, and positive signs and leading zeros are replaced with figure spaces
#' # (which have the same width as numbers and minus signs). These options ensure that all
#' # results will align neatly in columns when tabled.
#' format_num(x = c(10000, 1000, 2.41, -1.20, 0.41, -0.20))
#'
#' # By default, format_num() uses your computer locale's default decimal mark as
#' # the decimal point. To force the usage of "." instead (e.g., for submission to
#' # a U.S. journal), set decimal.mark = ".":
#' format_num(x = .41, decimal.mark = ".")
#'
#' # By default, format_num() separates groups of large digits using thin spaces.
#' # This is following the international standard for scientific communication (SI/ISO 31-0),
#' # which advises against using "." or "," to seprate digits because doing so can lead
#' # to confusion for human and computer readers because "." and "," are also used
#' # as decimal marks in various countries. If you prefer to use commmas to separate
#' # large digit groups, set big.mark = ",":
#' format_num(x = 10000, big.mark = ",")
format_num <- function(x, digits = 2L, decimal.mark = getOption("OutDec"),
                       leading0 = "conditional", drop0integer = FALSE,
                       neg.sign = "minus", pos.sign = "figure",
                       big.mark = "thinspace", big.interval = 3L,
                       small.mark = "thinspace", small.interval = 3L,
                       na.mark = "\u2014", lgl.mark = c("+", "\u2212"),
                       inf.mark = c("+\u221e", "\u2212\u221e") ){

        is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {if(is.numeric(x)) abs(x - round(x)) < tol else FALSE}

        # Input checking
        if(length(lgl.mark) == 1) lgl.mark <- c(lgl.mark, lgl.mark)
        if(length(inf.mark) == 1) inf.mark <- c(inf.mark, inf.mark)

        if(is.null(dim(x))) {
                x_type <- "vector"
                x <- as.data.frame(x, stringsAsFactors = FALSE)
        } else if("tbl_df" %in% class(x)) {
                x_type <- "tibble"
                x <- as.data.frame(x, stringsAsFactors = FALSE)
        } else if("matrix" %in% class(x)) {
                x_type <- "matrix"
                x <- as.data.frame(x, stringsAsFactors = FALSE)
        } else x_type <- "other"

        # Classify inputs
        which_logical  <- purrr::modify(x, is.logical) %>% as.matrix()
        which_integers <- as.matrix(purrr::modify(x, is.wholenumber))
        which_integers[is.na(which_integers)] <- FALSE
        which_infinite <- purrr::modify(x, is.infinite) %>% as.matrix()
        which_numeric <-
                purrr::modify(x, ~ is.numeric(.x) & !is.na(.x)) &
                !which_infinite &
                !which_integers %>% as.matrix()

        if(neg.sign == "minus") neg.sign <- "\u2212"
        if(big.mark == "thinspace") big.mark <- "\u202F"
        if(small.mark == "thinspace") small.mark <- "\u202F"

        if(pos.sign == FALSE) flag <- "" else flag <- "+"

        # Initial formatting for each type of data
        out <- x
        out[is.na(out)] <- na.mark
        ### Convert characters, factors, date/times, and and other unsupported
        ### types to character
        out <- purrr::modify_if(out, ~ !is.numeric(.x) & !is.logical(.x) & !is.complex(.x), as.character)
        out[which_logical & x == TRUE] <- lgl.mark[1]
        out[which_logical & x == FALSE] <- lgl.mark[2]
        out[which_infinite & x == Inf] <- inf.mark[1]
        out[which_infinite & x == -Inf] <- inf.mark[2]

        out <- purrr::modify_if(out, is.complex,
                                ~ format(.x, trim = TRUE, digits = digits, nsmall = digits,
                                         scientific = FALSE, big.mark = big.mark,
                                         big.interval = big.interval, small.mark = small.mark,
                                         small.interval = small.interval,
                                         decimal.mark = decimal.mark, drop0trailing = FALSE))

        out[which_integers] <- as.integer(x[which_integers]) %>%
                formatC(digits = digits, format = "f", flag = flag,
                        decimal.mark = decimal.mark,
                        big.mark = big.mark, big.interval = big.interval,
                        small.mark = small.mark, small.interval = small.interval,
                        drop0trailing = drop0integer)

        out[which_numeric] <- as.double(x[which_numeric]) %>%
                formatC(digits = digits, format = "f", flag = flag,
                        decimal.mark = decimal.mark,
                        big.mark = big.mark, big.interval = big.interval,
                        small.mark = small.mark, small.interval = small.interval,
                        drop0trailing = FALSE)


        # Clean up unicode big.mark and small.mark
        out[] <- dplyr::mutate_all(out,
                                   funs(stringr::str_replace_all(.data$.,
                                                                 paste0("(",paste(rev(strsplit(sub(" ", big.mark, " "),"")[[1]]), collapse=""),"|",sub(" ", big.mark, " "),")"),
                                                                 big.mark)))
        out[] <- dplyr::mutate_all(out,
                                   funs(stringr::str_replace_all(.data$.,
                                                                 paste0("(",paste(rev(strsplit(sub(" ", small.mark, " "),"")[[1]]), collapse=""),"|",sub(" ", small.mark, " "),")"),
                                                                 small.mark)))

        # Clean up leading zeros
        switch(leading0,
               "TRUE" = {},
               "FALSE" = out[] <- purrr::map(out, ~ stringr::str_replace_all(.x, paste0("^(\\+|-?)0", decimal.mark), paste0("\\1", decimal.mark))),
               conditional = {
                       out <- dplyr::mutate_if(out,
                                               sapply(x, function(i) {is.numeric(i) & !any(if(is.numeric(i)) abs(i) >= 1, na.rm = TRUE)}),
                                               function(i) stringr::str_replace_all(i, paste0("^(\\+|-?)0", decimal.mark), paste0("\\1", decimal.mark)))
               },
               figure = {
                       out <- dplyr::mutate_if(out,
                                               sapply(x, function(i) {is.numeric(i) & any(if(is.numeric(i)) abs(i) >= 1, na.rm = TRUE)}),
                                               function(i) stringr::str_replace_all(i, paste0("^(\\+|-?)0", decimal.mark), paste0("\\1\u2007", decimal.mark)))
                       out <- dplyr::mutate_if(out,
                                               sapply(x, function(i) {is.numeric(i) & !any(if(is.numeric(i)) abs(i) >= 1, na.rm = TRUE)}),
                                               function(i) stringr::str_replace_all(i, paste0("^(\\+|-?)0", decimal.mark), paste0("\\1", decimal.mark)))
               },
               # else =
               {out[] <- purrr::map(out, ~ stringr::str_replace_all(.x, paste0("^(\\+|-?)0", decimal.mark), paste0("\\1", leading0, decimal.mark)))}
        )

        # Clean up positive signs
        switch(pos.sign,
               "TRUE" = {},
               "FALSE" = {},
               figure = {
                       out <- dplyr::mutate_if(out,
                                               sapply(x, function(i) {is.numeric(i) & !any(if(is.numeric(i)) i < 0, na.rm = TRUE)}),
                                               function(i) stringr::str_replace_all(i, "^\\+", ""))
                       out <- dplyr::mutate_if(out,
                                               sapply(x, function(i) {is.numeric(i) & any(if(is.numeric(i)) i < 0, na.rm = TRUE)}),
                                               function(i) stringr::str_replace_all(i, "^\\+", "\u2007"))
               },
               # else =
               {out[] <- purrr::map(out, ~ stringr::str_replace_all(.x, "^\\+", pos.sign))}
        )

        # Clean up negative signs
        switch(neg.sign,
               "TRUE" = {},
               "FALSE" = {},
               "-" = {},
               # else =
               {out[] <- purrr::map(out, ~ stringr::str_replace_all(.x, "^-", neg.sign))}
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


#' Add metabulate equation commands and LaTeX dependencies
#'
#' \code{\link{metabulate}} requires several lines of code to correctly render meta-analysis results table column headings and footnotes. If \code{metabulate} is used to render files directly, these are added to the internal RMarkdown document. If you use \code{metabulate} output in a larger RMarkdown document, use this function to automatically add the necessary lines of code based on your chosen output format.
#'
#' @param latex Should required commands be included when converting to PDF, LaTeX, and related formats?
#' @param html Should required commands be included when converting to HTML and related formats?
#' @param word_proc Should required commands be included when converting to Word, ODT, and related formats?
#'
#' @section PDF and LaTeX output:
#' If \code{latex} is \code{TRUE} and you render to PDF, LaTeX, or other output
#' formats requiring LaTeX (e.g., \code{beamer_presentation}, see \code{\link[knitr]{is_latex_output}}),
#' a YAML metadata block with a \code{header-includes} argument calling the required
#' \code{unicode-math} LaTeX package is printed.
#'
#' An RMarkdown file can only include one \code{header-includes} metadata entry. If
#' your document already has one, set \code{latex} to \code{FALSE} and manually add
#' add the \code{unicode-math} package to your LaTeX header instead.
#'
#' (Note that \code{header-includes} is generally discouraged in favor of adding
#' an \code{include} argument to specific output formats, see
#' \url{http://rmarkdown.rstudio.com/pdf_document_format.html#includes}.)
#'
#' @section HTML output:
#' If \code{html} is \code{TRUE} and you render to HTML (or related formats, see
#' \code{\link[knitr]{is_html_output}}, the following LaTeX math commands are defined:
#' \itemize{
#'   \item{\code{symit}}{}
#'   \item{\code{symup}}{}
#'   \item{\code{symbfit}}{}
#'   \item{\code{symbfup}}{}
#' }
#'
#' If you define your own LaTeX or MathJax macros for these commands, set
#' \code{html} to \code{FALSE}.
#'
#' @section Microsoft Office and LibreOffice output:
#' If \code{word_proc} is \code{TRUE} and you render to Word or ODT (or related
#' formats such as PowerPoint), the following LaTeX math commands are defined:
#' \itemize{
#'   \item{\code{symit}}{}
#'   \item{\code{symup}}{}
#'   \item{\code{symbfit}}{}
#'   \item{\code{symbfup}}{}
#' }
#'
#' If you define your own LaTeX, Office, OpenDocument macros for these commands,
#' set \code{word_proc} to \code{FALSE}.
#'
#'
#' @return Requested commands are printed to the console.
#' @export
#'
#' @family output functions
#'
#' @examples
#' ## Include this line as 'asis' output in your RMarkdown document:
#' metabulate_rmd_helper()
#'
#' ## If you've already included \\usepackage{unicode-math} in your RMarkdown header
#' ## for PDF (and related formats) header, set latex to FALSE:
#' metabulate_rmd_helper(latex = FALSE)
metabulate_rmd_helper <- function(latex = TRUE, html = TRUE,
                                  word_proc = TRUE) {

        if(knitr::is_latex_output() & latex == TRUE) {
                cat("\n\n---\nheader-includes:\n  - \\usepackage{unicode-math}\n---\n\n")
        }

        if(knitr::is_html_output() & html == TRUE) {
                cat('\\newcommand{\\symup}{\\mathrm}\\newcommand{\\symbfit}{\\boldsymbol}\\newcommand{\\symbfup}{\\boldsymbol}\\newcommand{\\symit}{}\n\n')
        }

        if(!knitr::is_latex_output() & !knitr::is_html_output() & word_proc == TRUE ) {
                cat('\\newcommand{\\symup}{\\mathrm}\\newcommand{\\symbfup}{\\mathbfup}\\newcommand{\\symbfit}{\\mathbfit}\\newcommand{\\symit}{}\n\n')
        }

}

#' Write a summary table of meta-analytic results
#'
#' @param ma_obj A psychmeta meta-analysis object.
#' @param file The filename (optionally with a subfolder path) for the output file. If \code{NULL}, the function will output directly to the R console (also useful if you want to include psychmeta results in a larger RMarkdown document).
#' @param output_dir The filepath for the output directory/folder. Defaults to the current working directory.
#' @param show_msd Logical. Should means and standard deviations of effect sizes be shown (default \code{TRUE})
#' @param show_conf Logical. Should confidence intervals be shown (default: \code{TRUE})?
#' @param show_cred Logical. Should credibility intervals be shown (default: \code{TRUE})?
#' @param show_se Logical Should standard errors be shown (default: \code{FALSE})?
#' @param show_var Logical. Should variances be shown (default: \code{FALSE})?
#' @param analyses Which analyses to extract references for? See \code{\link{filter_ma}} for details.
#' @param match Match \code{all} or \code{any} of the filter criteria? See \code{\link{filter_ma}} for details.
#' @param case_sensitive Logical scalar that determines whether character values supplied in \code{analyses} should be treated as case sensitive (\code{TRUE}, default) or not (\code{FALSE}).
#' @param output_format The format of the output tables. Available options are Word (default), HTML, PDF (requires LaTeX and the \code{unicode-math} LaTeX package to be installed), ODT, rmd (Rmarkdown), and text (plain text). You can also specify the full name of another RMarkdown \code{\link[rmarkdown]{output_format}}.
#' @param ma_method Meta-analytic methods to be included. Valid options are: \code{"ad"}, \code{"ic"}, and \code{"bb"}. Multiple methods are permitted. By default, results are given for one method with order of priority: 1. \code{"ad"}, 2. \code{"ic"}, 3. \code{"bb"}.
#' @param correction_type Type of meta-analytic corrections to be incldued. Valid options are: "ts" (default), "vgx", and "vgy". Multiple options are permitted.
#' @param collapse_construct_labels  Should the construct labels for construct pairs with multiple rows of results be simplified so that only the first occurence of each set of construct names is shown (\code{TRUE}; default) or should construct labels be shown for each row of the table (\code{FALSE}).
#' @param bold_headers Logical. Should column headers be bolded (default: \code{TRUE})?
#' @param digits,decimal.mark,leading0,drop0integer,neg.sign,pos.sign,big.mark,big.interval,small.mark,small.interval,na.mark,lgl.mark,inf.mark Number formatting arguments. See \code{\link{format_num}} for details.
#' @param conf_format How should confidence intervals be formatted? Options are:
#' \itemize{
#' \item{\code{parentheses}}{: Bounds are enclosed in parentheses and separated by a comma: (LO, UP).}
#' \item{\code{brackets}}{: Bounds are enclosed in square brackets and separated by a comma: [LO, UP].}
#' \item{\code{columns}}{: Bounds are shown in individual columns.}
#' }
#' @param cred_format How should credility intervals be formatted? Options are the same as for \code{conf_format} above.
#' @param symbol_es For meta-analyses of generic (non-r, non-d) effect sizes, the symbol used for the effect sizes (default: \code{symbol_es = "ES"}).
#' @param caption Caption to print before tables. Either a character scalar or a named character vector with names corresponding to combinations of \code{ma_method} and \code{correction_type} (i.e., \code{bb}, \code{ic_ts}, \code{ad_vgx}, etc.).
#' @param header A list of YAML header parameters to pass to \code{link[rmarkdown]{render}}.
#' @param verbose Logical. Should detailed SD and variance components be shown (default: \code{FALSE})?
#' @param unicode Logical. If \code{output_format} is "text", should UTF-8 characters be used (defaults to system default).
#' @param bib A BibTeX file containing the citekeys for the meta-analyses. If provided and file is not \code{NULL}, a bibliography will be included with the meta-analysis table. See \code{\link{generate_bib}} for additional arguments controlling the bibliography.
#' @param title.bib The title to give to the bibliography (see \code{bib} above). If \code{NULL}, defaults to "Sources Contributing to Meta-Analyses"
#' @param style What style should the bibliography (see \code{bib} above) be formatted in? Can be a file path or URL for a \href{https://github.com/citation-style-language/styles}{CSL citation style} or the style ID for any style available from the \href{https://zotero.org/styles}{Zotero Style Repository}). Defaults to APA style. (Retrieving a style by ID requires an internet connection. If unavailable, references will be rendered in Chicago style.).
#' @param additional_citekeys Additional citekeys to include in the reference list (see \code{bib} above).
#' @param save_build_files Should the RMarkdown and BibLaTeX (if any) files used to generate the output be saved (default: \code{FALSE})?
#' @param ... Additional arguments to pass to \code{\link[rmarkdown]{render}}.
#'
#' @return A list of meta-analysis results \code{\link[tibble]{tibble}}s with "caption" and "footnote" attributes.
#'
#' If \code{file} is specified, formatted tables and bibliographies are exported in the requested \code{output_format}.
#'
#' @return Formatted tables of meta-analytic output.
#' @export
#'
#' @importFrom rmarkdown render
#' @importFrom dplyr case_when
#' @importFrom stringi stri_write_lines
#'
#' @family output functions
#'
#' @examples
#' \dontrun{
#' ## Create a results table for meta-analysis of correlations and output to Word:
#' ma_r_obj <- ma_r(ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'                  construct_x = x_name, construct_y = y_name,
#'                  moderators = moderator, data = data_r_meas_multi)
#' metabulate(ma_obj = ma_r_obj, file = "meta tables correlations",
#'            output_format = "word", output_dir = tempdir())
#'
#' ## Output to PDF:
#' metabulate(ma_obj = ma_r_obj, file = "meta tables correlations",
#'            output_format = "pdf", output_dir = tempdir())
#'
#' ## Output to ODT (LibreOffice):
#' metabulate(ma_obj = ma_r_obj, file = "meta tables correlations",
#'            output_format = "odt", output_dir = tempdir())
#'
#'
#' ## To produce Markdown tables to include inline in an RMarkdown report,
#' ## leave file == NULL and output_format to anything but "text":
#' ma_table <- metabulate(ma_obj = ma_r_obj, file = NULL, output_format = "rmd")
#'
#' ## Use the metabulate_rmd_helper() function to ensure all symbols render properly:
#' metabulate_rmd_helper()
#'
#' ## Then, add the formatted table to your document using your preferred table
#' ## formatting functions:
#'
#' #### Using just the \link{knitr} package:
#' knitr::kable(ma_table[[1]], caption = attr(ma_table[[1]], "caption"))
#' cat("\n", attr(ma_table[[1]], "footnote"))
#'
#' #### Using \link{knitr} plus the \link{kableExtra} package:
#' knitr::kable(ma_table[[1]], "latex", booktabs = TRUE,
#'                   caption = attr(ma_table[[1]], "caption")) %>%
#'        kableExtra::kable_styling(latex_options = c("striped", "hold_position")) %>%
#'        kableExtra::footnote(general = attr(ma_table[[1]], "footnote")
#'
#' # !!! Note: On Windows, R currently can only handle Unicode characters if kables
#' # are printed at top-level (e.g., in an if() statement, in a for() loop,
#' # or in lapply() or map() ). To correctly print Unicode metabulate tables, call
#' # kable() as a top-level function (as above).
#'
#'
#' ## Create output table for meta-analysis of d values:
#' ma_d_obj <- ma_d(ma_method = "ic", d = d, n1 = n1, n2 = n2, ryy = ryyi,
#'                  construct_y = construct, data = data_d_meas_multi)
#' ma_d_obj <- ma_d_ad(ma_obj = ma_d_obj, correct_rr_g = FALSE, correct_rr_y = FALSE)
#' metabulate(ma_obj = ma_d_obj, file = "meta tables d values", output_dir = tempdir())
#'
#' ## Create output table for meta-analysis of generic effect sizes:
#' dat <- data.frame(es = data_r_meas_multi$rxyi,
#'                   n = data_r_meas_multi$n,
#'                   var_e = (1 - data_r_meas_multi$rxyi^2)^2 / (data_r_meas_multi$n - 1))
#' ma_obj <- ma_generic(es = es, n = n, var_e = var_e, data = dat)
#' metabulate(ma_obj = ma_obj, file = "meta tables generic es", output_dir = tempdir())
#' }
metabulate <- function(ma_obj, file = NULL, output_dir = getwd(),
                       output_format=c("word", "html", "pdf", "odt", "text", "rmd"),
                       show_msd = TRUE, show_conf = TRUE, show_cred = TRUE,
                       show_se = FALSE, show_var = FALSE,
                       analyses="all", match=c("all", "any"), case_sensitive = TRUE,
                       ma_method = "ad", correction_type = "ts",
                       collapse_construct_labels  = TRUE, bold_headers = TRUE,
                       digits = 2L, decimal.mark = getOption("OutDec"),
                       leading0 = "conditional", drop0integer = FALSE,
                       neg.sign = "minus", pos.sign = "figure",
                       big.mark = "thinspace", big.interval = 3L,
                       small.mark = "thinspace", small.interval = 3L,
                       na.mark = "\u2014", lgl.mark = c("+", "\u2212"),
                       inf.mark = c("+\u221e", "\u2212\u221e"),
                       conf_format = "parentheses", cred_format = "parentheses",
                       symbol_es = "ES", caption = "Results of meta-analyses",
                       header = NULL, verbose = FALSE, unicode = NULL,
                       bib = NULL, title.bib = NULL, style = "apa", additional_citekeys = NULL,
                       save_build_files = FALSE, ...){

        # Match arguments
        output_format <- tolower(output_format)
        output_format <- tryCatch(match.arg(output_format),
                                  error = function(e) {
                                          if(length(output_format) > 1) {
                                                  stop("Multiple output formats specified. Please specify a single output format.", call. = FALSE)
                                          } else output_format})

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
                     if(length(ma_method[ma_method != "bb"]) > 0) paste(rep(ma_method[ma_method != "bb"], each = length(correction_type)), correction_type, sep = "_"))

        meta_tables <- ma_obj$meta_tables
        conf_level <- attributes(ma_obj$ma_obj)$inputs$conf_level
        cred_level <- attributes(ma_obj$ma_obj)$inputs$cred_level

        # Generate formatted tables (includes generation of captions and footnotes)
        meta_tables <- .metabulate(meta_tables = meta_tables, ma_type = ma_type,
                                   output_format = output_format, caption = caption,
                                   show_msd = show_msd, show_conf = show_conf,
                                   show_cred = show_cred, show_se = show_se,
                                   show_var = show_var,
                                   collapse_construct_labels = collapse_construct_labels, bold_headers = bold_headers,
                                   es_type = es_type, symbol_es = symbol_es,
                                   digits = digits, decimal.mark = decimal.mark,
                                   leading0 = leading0, neg.sign = neg.sign,
                                   pos.sign = pos.sign, drop0integer = drop0integer,
                                   big.mark = big.mark, big.interval = big.interval,
                                   small.mark = small.mark, small.interval = small.interval,
                                   na.mark = na.mark, inf.mark = inf.mark, lgl.mark = lgl.mark,
                                   conf_format = conf_format, cred_format = cred_format,
                                   verbose = verbose, unicode = unicode,
                                   conf_level = conf_level, cred_level = cred_level)

        # Set the output file name
        if(!is.null(file)) file <- .filename_check(file, output_format)

        # Assign values to citekeys and citations, convert bib to R bibliography
        if(!is.null(file)) {
                if(!is.null(bib)) warning("Bibliography not generated when file == NULL.\nTry generate_bib() to include a bibliography in an RMarkdown document.")
                bib <- NULL
        } else if(!is.null(bib)) {
                bib <- .generate_bib(ma_obj, bib, additional_citekeys)
        }

        # Render the output
        .psychmeta_render(file = file, output_format = output_format, output_dir = output_dir,
                          meta_tables = meta_tables, ma_type = ma_type, es_type = es_type,
                          bib = bib$bib, citations = bib$citations, citekeys = bib$citekeys,
                          title.bib = title.bib, style = style,
                          save_build_files = save_build_files, header = header, ...)

}



#' Generate a list of references included in meta-analyses
#'
#' This function generates a list of studies contributing to a meta-analysis
#'
#' @param ma_obj A psychmeta meta-analysis object with \code{citekeys} supplied.
#' @param bib A BibTeX file containing the citekeys for the meta-analyses.
#' @param title.bib The title to give to the bibliography. If \code{NULL}, defaults to "Sources Contributing to Meta-Analyses"
#' @param style What style should references be formatted in? Can be a file path or URL for a \href{https://github.com/citation-style-language/styles}{CSL citation style} or the style ID for any style available from the \href{https://zotero.org/styles}{Zotero Style Repository}). Defaults to APA style. (Retrieving a style by ID requires an internet connection. If unavailable, references will be rendered in Chicago style.).
#' @param additional_citekeys Additional citekeys to include in the reference list.
#' @param file The filename or filepath for the output file. If \code{NULL}, function will output directly to the R console (if \code{output_format} is "text", the formatted references in \code{\link[RefManageR]{BibOptions}} "authoryear" style; if "citekeys", the citekeys for included sources; otherwise, code to generate the bibliography in an RMarkdown document).
#' @param output_dir The filepath for the output file. Defaults to the current working directory.
#' @param output_format The format of the output reference list. Available options are Word (default), HTML, PDF (requires LaTeX to be installed), ODT, or Rmarkdown, plain text, and BibLaTeX. Returning only the item citekeys is also possible. You can also specify the full name of another RMarkdown \code{\link[rmarkdown]{output_format}}.
#' @param analyses Which analyses to extract references for? See \code{\link{filter_ma}} for details.
#' @param match Match \code{all} or \code{any} of the filter criteria? See \code{\link{filter_ma}} for details.
#' @param case_sensitive Logical scalar that determines whether character values supplied in \code{analyses} should be treated as case sensitive (\code{TRUE}, default) or not (\code{FALSE}).
#' @param save_build_files Should the BibTeX and RMarkdown files used to generate the bibliography be saved (default: \code{FALSE}; always \code{TRUE} if file is \code{NULL})?
#' @param header A list of YAML header parameters to pass to \code{link{rmarkdown::render}}.
#' @param ... Additional arguments to pass to \code{\link[rmarkdown]{render}}.
#'
#' @return A list containing the RefManageR BibEntry object. Additionally, a reference list formatted in the requested style and output_format is exported (or printed if file is "console").
#'
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
#' @family output functions
#'
#' @examples
#' \dontrun{
#' ## Run a meta-analysis using ma_r() and include a citekey argument to provide
#' ## citation information for each source contributing to the meta-analyses.
#' ma_obj <- ma_r(ma_method = "ic", rxyi = rxyi, n = n, rxx = rxxi, ryy = ryyi,
#'                construct_x = x_name, construct_y = y_name, sample_id = sample_id,
#'                moderators = moderator, citekey = citekey, data = data_r_meas_multi)
#'
#' ## Next, use generate_bib() to generate the bibliography for the retained studies.
#' ## The bib argument is the BibTeX or BibLaTeX .bib file containing the full
#' ## reference information for each of the citekeys included in the meta-analysis database.
#' generate_bib(ma_obj, bib = system.file("templates/sample_bibliography.bib", package="psychmeta"),
#'              file = "sample bibliography", output_dir = tempdir(), output_format = "word")
#' }
generate_bib <- function(ma_obj=NULL, bib=NULL, title.bib = NULL, style="apa",
                         additional_citekeys=NULL,
                         file = NULL, output_dir = getwd(),
                         output_format=c("word", "html", "pdf", "text", "odt", "rmd", "biblatex", "citekeys"),
                         analyses="all", match=c("all", "any"), case_sensitive = TRUE,
                         save_build_files = FALSE, header=list(), ...){

        output_format <- tolower(output_format)
        output_format <- tryCatch(match.arg(output_format),
                                  error = function(e) {
                                          if(length(output_format) > 1) {
                                                  stop("Multiple output formats specified. Please specify a single output format.", call. = FALSE)
                                          } else output_format})

        if("summary.ma_psychmeta" %in% class(ma_obj)) ma_obj <- ma_obj$ma_obj

        # Get the requested meta-analyses
        ma_obj <- filter_ma(ma_obj = ma_obj, analyses = analyses, match = match, case_sensitive = case_sensitive)

        # Assign values to citekeys and citations, convert bib to R bibliography
        bib <- .generate_bib(ma_obj, bib, additional_citekeys)
        citekeys <- bib$citekeys
        citations <- bib$citations
        bib <- bib$bib

        # Set the output file name
        if(!is.null(file)) file <- .filename_check(file, output_format)

        # Render the output
        switch(output_format,

               citekeys = {
                       if(!is.null(file)) {
                               return(citations)
                       } else {
                               writeLines(citations, con = file)
                       }
               },

               biblatex = {
                       if(!is.null(file)) {
                               print(bib[citekeys], .opts = list(style = "Biblatex"))
                       } else {
                               suppressMessages(WriteBib(bib[citekeys], file = file))
                       }
               },
               # else =
               {
                       .psychmeta_render(file = file, output_format = output_format,
                                         output_dir = output_dir, bib = bib, citations = citations,
                                         citekeys = citekeys, title.bib = title.bib, style = style,
                                         save_build_files = save_build_files, header = header, ...)
                       invisible(bib[citekeys])
               }
        )
}


.generate_bib <- function(ma_obj, bib, additional_citekeys){

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

        return(list(bib = bib, citekeys = citekeys, citations = citations))

}

.filename_check <- function(file, output_format){
        case_when(
                output_format == "pdf"      ~ if(!grepl("\\.pdf$",  file, ignore.case = TRUE)) {paste0(file, ".pdf")}  else {file},
                output_format == "html"     ~ if(!grepl("\\.html$", file, ignore.case = TRUE)) {paste0(file, ".html")} else {file},
                output_format == "word"     ~ if(!grepl("\\.docx$", file, ignore.case = TRUE)) {paste0(file, ".docx")} else {file},
                output_format == "rmd"      ~ if(!grepl("\\.rmd$",  file, ignore.case = TRUE)) {paste0(file, ".Rmd")}  else {file},
                output_format == "biblatex" ~ if(!grepl("\\.bib$",  file, ignore.case = TRUE)) {paste0(file, ".bib")}  else {file},
                output_format == "text"     ~ if(!grepl("\\.txt$",  file, ignore.case = TRUE)) {paste0(file, ".txt")}  else {file},
                output_format == "citekeys" ~ if(!grepl("\\.txt$",  file, ignore.case = TRUE)) {paste0(file, ".txt")}  else {file},
                output_format == "odt"      ~ if(!grepl("\\.odt$",  file, ignore.case = TRUE)) {paste0(file, ".odt")}  else {file},
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
                style <- paste0("https://zotero.org/styles/", stringr::str_replace(style, "\\.csl$", "") )
                attributes(style) <- list(source = "Zotero")
                return(style)
        }
}


.psychmeta_render <- function(file = NULL, output_format, meta_tables = NULL, ma_type = NULL,
                              es_type = NULL, bib = NULL, citations = NULL,
                              citekeys = NULL, title.bib = NULL, style = style,
                              save_build_files = FALSE, output_dir = NULL,
                              header = list(), ...){

        if(!output_format %in% c("text", "rmd", "md_document") & !rmarkdown::pandoc_available("2")) {
                stop("Output to ", output_format, " requires Pandoc 2.0.0 or greater.\nInstall from http://pandoc.org/installing.html or\nupdate to RStudio 1.2.679-1 or greater. ")
        }
        if(!is.null(style)) style <- .clean_style_name(style)

        if(is.null(output_dir)) output_dir <- getwd()

        if(is.null(file)) {

                switch(output_format,

                       text = {if(!is.null(meta_tables)) {
                               print(meta_tables)
                       }

                               if(!is.null(bib)) {
                                       if(is.null(title.bib)) title.bib <- "Sources Contributing to Meta-Analyses"
                                       cat(rep("\n", 2*as.numeric(is.null(meta_tables))),
                                           title.bib, "\n",
                                           rep("=", nchar(title.bib)), "\n\n", sep = ""
                                       )
                                       # TODO: Replace this with a call to citation.js to use CSL styles
                                       print(bib[citekeys], .opts = list(style = "text", bib.style = "authoryear"))
                               }

                               invisible(meta_tables[!is.null(meta_tables)])

                       },

                       # else =
                       {     if(!is.null(meta_tables)) {
                               # Prevent LaTeX from removing figure space characters
                               meta_tables <- sapply(names(meta_tables),
                                                     function(x) {
                                                             ma_table <- mutate_all(meta_tables[[x]],
                                                                                    funs(stringr::str_replace_all(.data$.,
                                                                                                                  "\u2007",
                                                                                                                  "\u2007\\\\phantom{\u2212}")))
                                                             attributes(ma_table) <- attributes(meta_tables[[x]])
                                                             return(ma_table)
                                                     },
                                                     simplify = FALSE, USE.NAMES = TRUE)
                               class(meta_tables) <- "metabulate"
                               return(meta_tables[!is.null(meta_tables)])
                       }

                               # TODO: If the bug with `bibliography: ` YAML metadata gets fixed, move this line to
                               # the same metadata block as citations and style below.
                               if(!is.null(bib)) {
                                       # Write the bibliography file
                                       # Ignore save_build_files
                                       bib_file <- file.path(output_dir, stringr::str_replace(file, "\\.(Rmd|pdf|docx|html|odt)$", "\\.bib"))
                                       suppressMessages(WriteBib(bib[citekeys],
                                                                 file = bib_file))

                                       sprintf("---\n### These metadata lines must be placed in your RMarkdown document main YAML header! ###\nbibliography: %s",
                                               bib_file)

                                       if(is.null(title.bib)) title.bib <- "# Sources Contributing to Meta-Analyses"

                                       cat(rep("\n", 2*as.numeric(is.null(meta_tables))),
                                           title.bib,
                                           "\n\n---\n"
                                       )

                                       if(!is.null(style)) sprintf("csl: %s\n", style )
                                       sprintf('nocite: |\n  %s\n---\n', citations)

                                       invisible(bib[citekeys])
                               }

                       }
                )

        } else {
                if(output_format == "rmd") save_build_files <- TRUE
                # Prevent LaTeX from removing figure space characters
                if(output_format != "text" & !is.null(meta_tables)) {
                        meta_tables <- sapply(names(meta_tables),
                                              function(x) {
                                                      ma_table <- mutate_all(meta_tables[[x]],
                                                                             funs(stringr::str_replace_all(.data$.,
                                                                                                           "\u2007",
                                                                                                           "\u2007\\\\phantom{\u2212}")))
                                                      attributes(ma_table) <- attributes(meta_tables[[x]])
                                                      return(ma_table)
                                              },
                                              simplify = FALSE, USE.NAMES = TRUE)
                        class(meta_tables) <- "metabulate"
                }

                # Fill in critical header slots
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
                                if(is.null(title.bib)) title.bib <- "Sources Contributing to Meta-Analyses"
                        }
                }

                switch(output_format,

                       text = {if(!is.null(bib)) {
                               sink("NUL")
                               # TODO: Replace this with a call to citation.js to use CSL styles
                               bibliography <- print(bib[citekeys], .opts = list(style = "text", bib.style = "authoryear"))
                               sink()
                       } else bibliography <- NULL

                       document <- c(paste(stringr::str_to_title(names(header)), ": ", header, collapse = "\n"), "\n\n",

                                     paste0("##### ", attr(meta_tables[[ma_type[1]]], "caption"))[!is.null(meta_tables[[ma_type[1]]])],
                                     suppressWarnings(knitr::knit_print(knitr::kable(meta_tables[[ma_type[1]]], align = attr(meta_tables[[ma_type[1]]], "align")))[!is.null(meta_tables[[ma_type[1]]])]),
                                     paste0("Note: ", attr(meta_tables[[ma_type[1]]], "footnote"), "\n\n\n")[!is.null(meta_tables[[ma_type[1]]])],
                                     "\\newpage"[!is.null(meta_tables[[ma_type[1]]]) & length(ma_type) > 1],

                                     paste0("##### ", attr(meta_tables[[ma_type[2]]], "caption"))[!is.null(meta_tables[[ma_type[2]]])],
                                     suppressWarnings(knitr::knit_print(knitr::kable(meta_tables[[ma_type[2]]], align = attr(meta_tables[[ma_type[2]]], "align")))[!is.null(meta_tables[[ma_type[2]]])]),
                                     paste0("Note: ", attr(meta_tables[[ma_type[2]]], "footnote"), "\n\n\n")[!is.null(meta_tables[[ma_type[2]]])],
                                     "\\newpage"[!is.null(meta_tables[[ma_type[2]]]) & length(ma_type) > 2],

                                     paste0("##### ", attr(meta_tables[[ma_type[3]]], "caption"))[!is.null(meta_tables[[ma_type[3]]])],
                                     suppressWarnings(knitr::knit_print(knitr::kable(meta_tables[[ma_type[3]]], align = attr(meta_tables[[ma_type[3]]], "align")))[!is.null(meta_tables[[ma_type[3]]])]),
                                     paste0("Note: ", attr(meta_tables[[ma_type[3]]], "footnote"), "\n\n\n")[!is.null(meta_tables[[ma_type[3]]])],
                                     "\\newpage"[!is.null(meta_tables[[ma_type[3]]]) & length(ma_type) > 3],

                                     paste0("##### ", attr(meta_tables[[ma_type[4]]], "caption"))[!is.null(meta_tables[[ma_type[4]]])],
                                     suppressWarnings(knitr::knit_print(knitr::kable(meta_tables[[ma_type[4]]], align = attr(meta_tables[[ma_type[4]]], "align")))[!is.null(meta_tables[[ma_type[4]]])]),
                                     paste0("Note: ", attr(meta_tables[[ma_type[4]]], "footnote"), "\n\n\n")[!is.null(meta_tables[[ma_type[4]]])],
                                     "\\newpage"[!is.null(meta_tables[[ma_type[4]]]) & length(ma_type) > 4],

                                     paste0("##### ", attr(meta_tables[[ma_type[5]]], "caption"))[!is.null(meta_tables[[ma_type[5]]])],
                                     suppressWarnings(knitr::knit_print(knitr::kable(meta_tables[[ma_type[5]]], align = attr(meta_tables[[ma_type[5]]], "align")))[!is.null(meta_tables[[ma_type[5]]])]),
                                     paste0("Note: ", attr(meta_tables[[ma_type[5]]], "footnote"), "\n\n\n")[!is.null(meta_tables[[ma_type[5]]])],
                                     "\\newpage"[!is.null(meta_tables[[ma_type[5]]]) & length(ma_type) > 5],

                                     paste0("##### ", attr(meta_tables[[ma_type[6]]], "caption"))[!is.null(meta_tables[[ma_type[6]]])],
                                     suppressWarnings(knitr::knit_print(knitr::kable(meta_tables[[ma_type[6]]], align = attr(meta_tables[[ma_type[6]]], "align")))[!is.null(meta_tables[[ma_type[6]]])]),
                                     paste0("Note: ", attr(meta_tables[[ma_type[6]]], "footnote"), "\n\n\n")[!is.null(meta_tables[[ma_type[6]]])],
                                     "\\newpage"[!is.null(meta_tables[[ma_type[6]]]) & length(ma_type) > 6],

                                     paste0("##### ", attr(meta_tables[[ma_type[7]]], "caption"))[!is.null(meta_tables[[ma_type[7]]])],
                                     suppressWarnings(knitr::knit_print(knitr::kable(meta_tables[[ma_type[7]]], align = attr(meta_tables[[ma_type[7]]], "align")))[!is.null(meta_tables[[ma_type[7]]])]),
                                     paste0("Note: ", attr(meta_tables[[ma_type[7]]], "footnote"), "\n\n\n")[!is.null(meta_tables[[ma_type[7]]])],
                                     "\\newpage"[!is.null(meta_tables) & !is.null(bib)],

                                     paste0(title.bib, "\n\n")[!is.null(title.bib) & !is.null(bib)],
                                     bibliography[!is.null(bib)]
                       )

                       stringi::stri_write_lines(document, file, sep="")
                       },

                       # else =
                       {
                               # Fill in critical header slots and write .bib file if necessary
                               if (output_format %in% c("word", "html", "pdf", "odt")) {
                                       header$output <- paste0("\n  ", output_format, "_document",
                                                               if(output_format == "pdf") paste0(":\n    latex_engine: lualatex\n    includes:\n      in_header: ", system.file('templates/header.tex', package='psychmeta')),
                                                               if(output_format == "word" & is.null(meta_tables)) paste0(":\n    reference_docx: ", system.file('templates/reference_docx.docx', package='psychmeta')),
                                                               if(output_format == "word" & !is.null(meta_tables)) paste0(":\n    reference_docx: ", system.file('templates/reference_docx_landscape.docx', package='psychmeta')),
                                                               if(output_format == "odt"& is.null(meta_tables)) paste0(":\n    reference_odt: ", system.file('templates/reference_odt.odt', package='psychmeta')),
                                                               if(output_format == "odt"& !is.null(meta_tables)) paste0(":\n    reference_odt: ", system.file('templates/reference_odt_landscape.odt', package='psychmeta'))
                                       )
                               } else header$output <- output_format

                               if(!is.null(bib)) {
                                       # Write the bibliography file
                                       if(save_build_files) {
                                               bib_file <- file.path(output_dir, stringr::str_replace(file, "\\.(Rmd|pdf|docx|html|odt)$", "\\.bib"))
                                       } else bib_file <- file.path(tempdir(), stringr::str_replace(file, "\\.(Rmd|pdf|docx|html|odt)$", "\\.bib"))
                                       suppressMessages(WriteBib(bib[citekeys],
                                                                 file = bib_file))

                                       header$bibliography <- if(.Platform$`OS.type` == "windows") {
                                               gsub("\\\\", "/", bib_file)
                                       } else bib_file


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

                               # Save meta_tables to an RData workspace for later loading by Rmarkdown
                               if(!is.null(meta_tables)) {
                                       if(save_build_files) {
                                               rdata_file <- file.path(output_dir, stringr::str_replace(file, "\\.(pdf|docx|html|odt)$", "\\.Rdata"))
                                       } else rdata_file <- file.path(tempdir(), stringr::str_replace(file, "\\.(pdf|docx|html|odt)$", "\\.Rdata"))
                                       save(meta_tables, file = rdata_file)

                                       rdata_file <- if(.Platform$`OS.type` == "windows") {
                                               gsub("\\\\", "/", rdata_file)
                                       } else rdata_file

                                       tables_document <- c(

                                               "```{r eval=knitr::is_html_output(), results='asis', echo = FALSE}",
                                               "cat('\\\\newcommand{\\\\symup}{\\\\mathrm}\\\\newcommand{\\\\symbfit}{\\\\boldsymbol}\\\\newcommand{\\\\symbfup}{\\\\boldsymbol}\\\\newcommand{\\\\symit}{}')",
                                               "```\n",
                                               "```{r eval=!knitr::is_latex_output() & !knitr::is_html_output(), results='asis', echo = FALSE}",
                                               "cat('\\\\newcommand{\\\\symup}{\\\\mathrm}\\\\newcommand{\\\\symbfup}{\\\\mathbfup}\\\\newcommand{\\\\symbfit}{\\\\mathbfit}\\\\newcommand{\\\\symit}{}')",
                                               "```\n",

                                               "```{r, echo=FALSE}",
                                               paste0("load('", rdata_file, "')"),
                                               "```\n",
                                               "\\blandscape\n",

                                               '```{r, results = "asis", echo = FALSE}',
                                               'ma_type <- names(meta_tables)',

                                               'cat(paste0("##### ", attr(meta_tables[[ma_type[1]]], "caption"), "\\n")[!is.null(meta_tables[[ma_type[1]]])])',
                                               'knitr::knit_print(knitr::kable(meta_tables[[ma_type[1]]], align = attr(meta_tables[[ma_type[1]]], "align")))',
                                               'cat(paste0("*Note:* ", attr(meta_tables[[ma_type[1]]], "footnote"), "\\n\\n")[!is.null(meta_tables[[ma_type[1]]])])',
                                               'cat("\\\\newpage"[!is.null(meta_tables[[ma_type[1]]])])',

                                               'cat(paste0("##### ", attr(meta_tables[[ma_type[2]]], "caption"), "\\n")[!is.null(meta_tables[[ma_type[2]]])])',
                                               'knitr::knit_print(knitr::kable(meta_tables[[ma_type[2]]], align = attr(meta_tables[[ma_type[2]]], "align")))',
                                               'cat(paste0("*Note:* ", attr(meta_tables[[ma_type[2]]], "footnote"), "\\n\\n")[!is.null(meta_tables[[ma_type[2]]])])',
                                               'cat("\\\\newpage"[!is.null(meta_tables[[ma_type[2]]])])',

                                               'cat(paste0("##### ", attr(meta_tables[[ma_type[3]]], "caption"), "\\n")[!is.null(meta_tables[[ma_type[3]]])])',
                                               'knitr::knit_print(knitr::kable(meta_tables[[ma_type[3]]], align = attr(meta_tables[[ma_type[3]]], "align")))',
                                               'cat(paste0("*Note:* ", attr(meta_tables[[ma_type[3]]], "footnote"), "\\n\\n")[!is.null(meta_tables[[ma_type[3]]])])',
                                               'cat("\\\\newpage"[!is.null(meta_tables[[ma_type[3]]])])',

                                               'cat(paste0("##### ", attr(meta_tables[[ma_type[4]]], "caption"), "\\n")[!is.null(meta_tables[[ma_type[4]]])])',
                                               'knitr::knit_print(knitr::kable(meta_tables[[ma_type[4]]], align = attr(meta_tables[[ma_type[4]]], "align")))',
                                               'cat(paste0("*Note:* ", attr(meta_tables[[ma_type[4]]], "footnote"), "\\n\\n")[!is.null(meta_tables[[ma_type[4]]])])',
                                               'cat("\\\\newpage"[!is.null(meta_tables[[ma_type[4]]])])',

                                               'cat(paste0("##### ", attr(meta_tables[[ma_type[5]]], "caption"), "\\n")[!is.null(meta_tables[[ma_type[5]]])])',
                                               'knitr::knit_print(knitr::kable(meta_tables[[ma_type[5]]], align = attr(meta_tables[[ma_type[5]]], "align")))',
                                               'cat(paste0("*Note:* ", attr(meta_tables[[ma_type[5]]], "footnote"), "\\n\\n")[!is.null(meta_tables[[ma_type[5]]])])',
                                               'cat("\\\\newpage"[!is.null(meta_tables[[ma_type[5]]])])',

                                               'cat(paste0("##### ", attr(meta_tables[[ma_type[6]]], "caption"), "\\n")[!is.null(meta_tables[[ma_type[6]]])])',
                                               'knitr::knit_print(knitr::kable(meta_tables[[ma_type[6]]], align = attr(meta_tables[[ma_type[6]]], "align")))',
                                               'cat(paste0("*Note:* ", attr(meta_tables[[ma_type[6]]], "footnote"), "\\n\\n")[!is.null(meta_tables[[ma_type[6]]])])',
                                               'cat("\\\\newpage"[!is.null(meta_tables[[ma_type[6]]])])',

                                               'cat(paste0("##### ", attr(meta_tables[[ma_type[7]]], "caption"), "\\n")[!is.null(meta_tables[[ma_type[7]]])])',
                                               'knitr::knit_print(knitr::kable(meta_tables[[ma_type[7]]], align = attr(meta_tables[[ma_type[7]]], "align")))',
                                               'cat(paste0("*Note:* ", attr(meta_tables[[ma_type[7]]], "footnote"), "\\n\\n")[!is.null(meta_tables[[ma_type[7]]])])',

                                               '```',

                                               "\n\n\\elandscape\n\n"
                                       )
                               } else tables_document <- NULL

                               ## Create the markdown header and document
                               header <- paste(names(header), header, sep=": ", collapse="\n")

                               document <- c("---", header, "---\n\n",
                                             tables_document,
                                             "<br/>\n\n<br/>\n\n<br/>\n\n",
                                             paste0("# ", title.bib, "\n\n")[!is.null(title.bib) & !is.null(bib)],
                                             sprintf('---\nnocite: |\n  %s\n---\n', citations)[!is.null(bib)]
                               )

                               # Create Rmd and output files
                               if(save_build_files) {
                                       rmd_document <- file.path(output_dir, stringr::str_replace(file, "\\.(pdf|docx|html|odt)$", "\\.rmd"))
                               } else rmd_document <- file.path(tempdir(), stringr::str_replace(file, "\\.(pdf|docx|html|odt)$", "\\.rmd"))

                               stringi::stri_write_lines(document, rmd_document)

                               if(output_format != "rmd") {
                                       rmarkdown::render(rmd_document,
                                                         output_file = file,
                                                         output_dir  = output_dir,
                                                         encoding = "UTF-8",
                                                         ...)
                               }
                               invisible(c(list(meta_tables = meta_tables)[!is.null(meta_tables)], list(bib = bib[citekeys])[!is.null(bib)]))
                       }
                )
        }
}


#' Internal function for .metabulating results tables
#'
#' @keywords internal
.metabulate <- function(meta_tables, ma_type = "ad_ts", output_format = "word",
                        caption = caption, show_msd = TRUE, show_conf = TRUE,
                        show_cred = TRUE, show_se = FALSE, show_var = FALSE,
                        collapse_construct_labels = TRUE, bold_headers = TRUE,
                        es_type = NULL, symbol_es = "ES", digits = 2L,
                        decimal.mark = getOption("OutDec"), leading0 = "figure",
                        neg.sign = "\u2212", pos.sign = "figure", drop0integer = TRUE,
                        big.mark = "\u202F", big.interval = 3L, small.mark = "\u202F",
                        na.mark = "\u2014", lgl.mark = c("+", "\u2212"),
                        inf.mark = c("\u2007\u221e", "\u2212\u221e"),
                        small.interval = 3L, conf_format = "parentheses",
                        cred_format = "parentheses", verbose = FALSE,
                        unicode = unicode, conf_level = .95, cred_level = .80) {

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

        length_initial <- max(sapply(meta_tables, function(x) {
                initial_names <- colnames(x)[1:which(colnames(x) == "analysis_type")]
                initial_names <- initial_names[!initial_names %in% c("analysis_id", "pair_id", "analysis_type")]
                length(initial_names)
        }))
        length_moderators <- max(sapply(meta_tables, function(x) (which(colnames(x) == "k") -1) - which(colnames(x) == "analysis_type") ))

        # Select, rearrange, and format columns of meta_tables
        .arrange_format_columns <- function(ma_table, collapse_construct_labels) {
                x <- colnames(ma_table)

                if(collapse_construct_labels & nrow(ma_table) > 1 & "pair_id" %in% x){
                        pair_ids <- unlist(ma_table[["pair_id"]])
                        delete_id <- FALSE
                        for(i in 2:length(pair_ids)) delete_id[i] <- pair_ids[i] == pair_ids[i-1]
                        if("construct_x" %in% x)    ma_table[["construct_x"]][delete_id] <- ""
                        if("group_contrast" %in% x) ma_table[["group_contrast"]][delete_id] <- ""
                        if("construct_y" %in% x)    ma_table[["construct_y"]][delete_id] <- ""
                        rm(pair_ids, delete_id)
                }

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
                ma_table <- ma_table[1:nrow(ma_table), c(col_initial, col_moderators, col_sampsize,
                                                         col_m_bb, col_se_bb, col_sd_bb,
                                                         col_var_bb,
                                                         col_m_cor, col_se_cor, col_sd_cor,
                                                         col_var_cor,
                                                         col_conf,
                                                         col_cred)]

                colnames(ma_table)[colnames(ma_table) %in% col_conf] <- c("ci_lower", "ci_upper")
                if(show_conf == TRUE) col_conf <- c("ci_lower", "ci_upper")
                colnames(ma_table)[colnames(ma_table) %in% col_cred] <- c("cv_lower", "cv_upper")
                if(show_cred == TRUE) col_cred <- c("cv_lower", "cv_upper")

                # Format columns
                ma_table[1:nrow(ma_table), col_sampsize] <-
                        format_num(ma_table[1:nrow(ma_table), col_sampsize], digits = 0L, decimal.mark = decimal.mark,
                                   leading0 = leading0, neg.sign = neg.sign, pos.sign = pos.sign,
                                   drop0integer = drop0integer, big.mark = big.mark, big.interval = big.interval,
                                   small.mark = small.mark, small.interval = small.interval, na.mark = na.mark,
                                   inf.mark = inf.mark, lgl.mark = lgl.mark)

                numeric_columns <- c(col_m_bb, col_se_bb, col_sd_bb, col_var_bb, col_m_cor, col_se_cor, col_sd_cor, col_var_cor, col_conf, col_cred)

                ma_table[1:nrow(ma_table), numeric_columns] <-
                        format_num(ma_table[1:nrow(ma_table), numeric_columns], digits = digits, decimal.mark = decimal.mark,
                                   leading0 = leading0, neg.sign = neg.sign, pos.sign = pos.sign,
                                   drop0integer = drop0integer, big.mark = big.mark, big.interval = big.interval,
                                   small.mark = small.mark, small.interval = small.interval, na.mark = na.mark,
                                   inf.mark = inf.mark, lgl.mark = lgl.mark)

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
        .rename_columns <- function(ma_table, formatted_strings, output_format, length_initial, length_moderators, bold_headers) {
                formatted_names <- names(formatted_strings$col_names)
                names(formatted_names) <- formatted_strings$col_names

                ma_table <- rename(ma_table, !!formatted_names[formatted_names %in% colnames(ma_table)])

                if(bold_headers & output_format != "text" & length_moderators > 0) {
                        colnames(ma_table)[(length_initial + 1):(length_initial + length_moderators)] <-
                                paste0("**", colnames(ma_table)[(length_initial + 1):(length_initial + length_moderators)], "**")
                }

                ma_table
        }

        .format_meta_table <- function(ma_table_name, meta_tables, formatted_strings, caption, output_format, length_initial, length_moderators) {
                ma_table <- meta_tables[[ma_table_name]] %>%
                        .arrange_format_columns(collapse_construct_labels = collapse_construct_labels) %>%
                        .rename_columns(formatted_strings, output_format, length_initial, length_moderators, bold_headers = bold_headers)
                attr(ma_table, "footnote") <- formatted_strings$footnote[[ma_table_name]]
                if(length(caption) > 1) attr(ma_table, "caption") <- caption[[ma_table_name]] else attr(ma_table, "caption") <- caption
                attr(ma_table, "align") <- c(rep("l", length_initial + length_moderators), rep("r", 2), rep("c", ncol(ma_table) - length_initial - length_moderators - 2))
                ma_table
        }

        ### TODO: Add notes about actual corrections applied
        formatted_strings <- .formatted_strings(output_format = output_format, es_type = es_type, symbol_es = symbol_es,
                                                conf_level = conf_level, cred_level = cred_level, corrections = NULL,
                                                show_msd = show_msd, show_se = show_se, show_var = show_var,
                                                verbose = verbose, unicode = unicode,
                                                show_conf = show_conf, show_cred = show_cred, bold_headers = bold_headers)

        meta_tables <- sapply(names(meta_tables), .format_meta_table,
                              meta_tables = meta_tables, caption = caption,
                              formatted_strings = formatted_strings,
                              output_format = output_format, length_initial = length_initial,
                              length_moderators = length_moderators,
                              simplify = FALSE, USE.NAMES = TRUE)

        class(meta_tables) <- "metabulate"

        return(meta_tables)

}


#' Internal list of strings with formatted variable names for output
#'
#' @keywords Internal
#'
.formatted_strings <- function(output_format, es_type, symbol_es = "ES", conf_level = .95, cred_level = .80, corrections = NULL,
                               show_msd = TRUE, show_se = FALSE, show_var = FALSE, verbose = FALSE, show_conf = TRUE, show_cred = TRUE,
                               unicode = NULL, bold_headers = TRUE) {

        col_names <- if(output_format == "text") {
                if(.support_unicode(unicode)) {
                        c(
                                group_contrast    = "Group Contrast",
                                construct_x       = "Construct X",
                                construct_y       = "Construct Y",
                                k                 = "k",
                                N                 = "N",

                                mean_r            = "r\u0305",
                                var_r             = "\u03C3\u00B2_r",
                                sd_r              = "SD_r",
                                se_r              = "SE_r\u0305",

                                mean_d            = "d\u0305",
                                var_d             = "\u03C3\u00B2_d",
                                sd_d              = "SD_d",
                                se_d              = "SE_r\u0305",

                                mean_es           = paste0(strsplit(symbol_es, NULL)[[1]], "\u305", collapse = ""),
                                var_es            = paste0("\u03C3\u00B2_", symbol_es),
                                sd_es             = paste0("SD_", symbol_es),
                                se_es             = paste0("SE_", paste0(strsplit(symbol_es, NULL)[[1]], "\u305", collapse = "")),

                                var_e             = "\u03C3\u00B2_e",
                                var_res           = "\u03C3\u00B2_res",
                                sd_e              = "SD_e",
                                sd_res            = "SD_res",
                                var_art           = "\u03C3\u00B2_art",
                                var_pre           = "\u03C3\u00B2_pre",
                                sd_art            = "SD_art",
                                sd_pre            = "SD_pre",

                                mean_rho          = "\u03C1\u0305",
                                var_r_c           = "\u03C3\u00B2_r(c)",

                                var_rho           = "\u03C3\u00B2_\u03C1",
                                sd_r_c            = "SD_r(c)",
                                se_r_c            = "SE_\u03C1\u0305",
                                sd_rho            = "SD_\u03C1",

                                mean_delta        = "\u03B4\u0305",
                                var_d_c           = "\u03C3\u00B2_d(c)",
                                sd_d_c            = "SD_d(c)",
                                se_d_c            = "SE_\u03B4\u0305",
                                var_delta         = "\u03C3\u00B2_\u03B4",
                                sd_delta          = "SD_\u03B4",

                                var_e_c           = "\u03C3\u00B2_e(c)",
                                var_art_c         = "\u03C3\u00B2_art(c)",
                                var_pre_c         = "\u03C3\u00B2_pre(c)",
                                sd_e_c            = "SD_e(c)",
                                sd_art_c          = "SD_art(c)",
                                sd_pre_c          = "SD_pre(c)",

                                ci_lower          = paste0(conf_level*100, "% CI Lower"),
                                ci_upper          = paste0(conf_level*100, "% CI Upper"),
                                cv_lower          = paste0(cred_level*100, "% CV Lower"),
                                cv_upper          = paste0(cred_level*100, "% CV Upper"),

                                conf_int          = paste0(conf_level*100, "% CI"),
                                cred_int          = paste0(cred_level*100, "% CV")
                        )
                } else {
                        c(
                                group_contrast    = "Group Contrast",
                                construct_x       = "Construct X",
                                construct_y       = "Construct Y",
                                k                 = "k",
                                N                 = "N",

                                mean_r            = "Mean r",
                                var_r             = "Var_r",
                                sd_r              = "SD_r",
                                se_r              = "SE_[Mean r]",

                                mean_d            = "Mean d",
                                var_d             = "Var_d",
                                sd_d              = "SD_d",
                                se_d              = "SE_[Mean d]",

                                mean_es           = paste0("Mean ", symbol_es),
                                var_es            = paste0("Var_", symbol_es),
                                sd_es             = paste0("SD_", symbol_es),
                                se_es             = paste0("SE_[Mean ", symbol_es, "]"),

                                var_e             = "Var_e",
                                var_res           = "Var_res",
                                sd_e              = "SD_e",
                                sd_res            = "SD_res",
                                var_art           = "Var_art",
                                var_pre           = "Var_pre",
                                sd_art            = "SD_art",
                                sd_pre            = "SD_pre",

                                mean_rho          = "Mean rho",
                                var_r_c           = "Var_r(c)",

                                var_rho           = "Var_rho",
                                sd_r_c            = "SD_r(c)",
                                se_r_c            = "SE_[Mean rho]",
                                sd_rho            = "SD_rho",

                                mean_delta        = "Mean delta",
                                var_d_c           = "Var_d(c)",
                                sd_d_c            = "SD_d(c)",
                                se_d_c            = "SE_[Mean delta]",
                                var_delta         = "Var_delta",
                                sd_delta          = "SD_delta",

                                var_e_c           = "Var_e(c)",
                                var_art_c         = "Var_art(c)",
                                var_pre_c         = "Var_pre(c)",
                                sd_e_c            = "SD_e(c)",
                                sd_art_c          = "SD_art(c)",
                                sd_pre_c          = "SD_pre(c)",

                                ci_lower          = paste0(conf_level*100, "% CI Lower"),
                                ci_upper          = paste0(conf_level*100, "% CI Upper"),
                                cv_lower          = paste0(cred_level*100, "% CV Lower"),
                                cv_upper          = paste0(cred_level*100, "% CV Upper"),

                                conf_int          = paste0(conf_level*100, "% CI"),
                                cred_int          = paste0(cred_level*100, "% CV")
                        )
                }
        } else {
                c(
                        group_contrast    = "**Group Contrast**",
                        construct_x       = "**Construct X**",
                        construct_y       = "**Construct Y**",
                        k                 = "**_k_**",
                        N                 = "**_N_**",

                        mean_r            = "$\\symbfit{\\overline{r}}$",
                        var_r             = "$\\symbfit{\\sigma^{2}_{r}}$",
                        sd_r              = "$\\symbfit{SD_{r}}$",
                        se_r              = "$\\symbfit{SE_{\\overline{r}}}$",

                        mean_d            = "$\\symbfit{\\overline{d}}$",
                        var_d             = "$\\symbfit{\\sigma^{2}_{d}}$",
                        sd_d              = "$\\symbfit{SD_{d}}$",
                        se_d              = "$\\symbfit{SE_{\\overline{d}}}$",

                        mean_es           = paste0("$\\symbfit{\\overline{", symbol_es, "}}$"),
                        var_es            = paste0("$\\symbfit{\\sigma^{2}_{", symbol_es, "}}$"),
                        sd_es             = paste0("$\\symbfit{SD_{", symbol_es, "}}$"),
                        se_es             = paste0("$\\symbfit{SE_{\\overline{", symbol_es, "}}}$"),

                        var_e             = "$\\symbfit{\\sigma^{2}_{e}}$",
                        var_res           = "$\\symbfit{\\sigma^{2}_{res}}$",
                        sd_e              = "$\\symbfit{SD_{e}}$",
                        sd_res            = "$\\symbfit{SD_{res}}$",
                        var_art           = "$\\symbfit{\\sigma^{2}_{art}}$",
                        var_pre           = "$\\symbfit{\\sigma^{2}_{pre}}$",
                        sd_art            = "$\\symbfit{SD_{art}}$",
                        sd_pre            = "$\\symbfit{SD_{pre}}$",

                        mean_rho          = "$\\symbfup{\\overline{\\rho}}$",
                        var_r_c           = "$\\symbfit{\\sigma^{2}_{r_{c}}}$",

                        var_rho           = "$\\symbfit{\\sigma^{2}_{\\symbfup{\\rho}}}$",
                        sd_r_c            = "$\\symbfit{SD_{r_{c}}}$",
                        se_r_c            = "$\\symbfit{SE_{\\symbfup{\\overline{\\rho}}}}$",
                        sd_rho            = "$\\symbfit{SD_{\\symbfup{\\rho}}}$",

                        mean_delta        = "$\\symbfit{\\overline{\\delta}}$",
                        var_d_c           = "$\\symbfit{\\sigma^{2}_{d_{c}}}$",
                        sd_d_c            = "$\\symbfit{SD_{d_{c}}}$",
                        se_d_c            = "$\\symbfit{SE_{\\overline{\\delta}}}$",
                        var_delta         = "$\\symbfit{\\sigma^{2}_{\\delta}}$",
                        sd_delta          = "$\\symbfit{SD_{\\delta}}$",

                        var_e_c           = "$\\symbfit{\\sigma^{2}_{e_{c}}}$",
                        var_art_c         = "$\\symbfit{\\sigma^{2}_{art_{c}}}$",
                        var_pre_c         = "$\\symbfit{\\sigma^{2}_{pre_{c}}}$",
                        sd_e_c            = "$\\symbfit{SD_{e_{c}}}$",
                        sd_art_c          = "$\\symbfit{SD_{art_{c}}}$",
                        sd_pre_c          = "$\\symbfit{SD_{pre_{c}}}$",

                        ci_lower          = paste0("**", conf_level*100, "% CI Lower**"),
                        ci_upper          = paste0("**", conf_level*100, "% CI Upper**"),
                        cv_lower          = paste0("**", cred_level*100, "% CV Lower**"),
                        cv_upper          = paste0("**", cred_level*100, "% CV Upper**"),

                        conf_int          = paste0("**", conf_level*100, "% CI**"),
                        cred_int          = paste0("**", cred_level*100, "% CV**")
                )
        }

        if(!bold_headers){
                var_labels <- names(col_names)
                col_names <- stringr::str_replace_all(col_names, pattern = "symbfup", replacement = "symup")
                col_names <- stringr::str_replace_all(col_names, pattern = "symbfit", replacement = "symit")
                col_names <- stringr::str_replace_all(col_names, pattern = "[*][*]", replacement = "")
                names(col_names) <- var_labels
        }

        footnote <- if(output_format == "text") {
                ### TODO: Add notes about actual corrections applied
                if(.support_unicode(unicode)) {
                        if(es_type == "r") {
                                c(bb     = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "r\u0305 =  mean observed correlation; "[show_msd],
                                                  "SE_r\u0305 = standard error of r\u0305; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_res = residual variance of r ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_e = predicted variance of r due to sampling error; \u03C3\u00B2_res = residual variance of r; "[show_var & verbose],
                                                  "CI = confidence interval around r\u0305; "[show_conf],
                                                  "CV = credibility interval around r\u0305."[show_cred]),

                                  ic_ts  = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "r\u0305 = mean observed correlation; "[show_msd],
                                                  "SE_r\u0305 = standard error of r\u0305; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_art = predicted SD_r due to artifacts; SD_pre = total predicted SD_r; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_res = residual variance of r ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_e = predicted variance of r due to sampling error; \u03C3\u00B2_art = predicted variance of r due to artifacts; \u03C3\u00B2_pre = total predicted variance of r; \u03C3\u00B2_res = residual variance of r; "[show_var & verbose],

                                                  "\u03C1\u0305 = mean true-score correlation; "[show_msd],
                                                  "SE_\u03C1\u0305 = standard error of \u03C1\u0305; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_\u03C1 = residual standard deviation of \u03C1; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_e = predicted SD_r(c) due to sampling error; SD_art(c) = predicted SD_r(c) due to artifacts; SD_pre(c) = total predicted SD_r(c); SD_\u03C1 = residual standard deviation of \u03C1; "[show_msd & verbose],

                                                  "\u03C3\u00B2_r(c) = observed variance of r(c); \u03C3\u00B2_\u03C1 = residual variance of \u03C1 ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_r(c) = observed variance of r(c); \u03C3\u00B2_e(c) = predicted variance of r(c) due to sampling error; \u03C3\u00B2_art(c) = predicted variance of r(c) due to artifacts; \u03C3\u00B2_pre(c) = total predicted variance of r(c); \u03C3\u00B2_\u03C1 = residual variance of \u03C1; "[show_var & verbose],

                                                  "CI = confidence interval around \u03C1\u0305; "[show_conf],
                                                  "CV = credibility interval around \u03C1\u0305; "[show_cred],

                                                  "correlations corrected individually."),

                                  ic_vgx = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "r\u0305 = mean observed correlation; "[show_msd],
                                                  "SE_r\u0305 = standard error of r\u0305; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_art = predicted SD_r due to artifacts; SD_pre = total predicted SD_r; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_res = residual variance of r ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_e = predicted variance of r due to sampling error; \u03C3\u00B2_art = predicted variance of r due to artifacts; \u03C3\u00B2_pre = total predicted variance of r; \u03C3\u00B2_res = residual variance of r; "[show_var & verbose],

                                                  "\u03C1\u0305 = mean operational validity (X measured with error); "[show_msd],
                                                  "SE_\u03C1\u0305 = standard error of \u03C1\u0305; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_\u03C1 = residual standard deviation of \u03C1; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_e = predicted SD_r(c) due to sampling error; SD_art(c) = predicted SD_r(c) due to artifacts; SD_pre(c) = total predicted SD_r(c); SD_\u03C1 = residual standard deviation of \u03C1; "[show_msd & verbose],

                                                  "\u03C3\u00B2_r(c) = observed variance of r(c); \u03C3\u00B2_\u03C1 = residual variance of \u03C1 ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_r(c) = observed variance of r(c); \u03C3\u00B2_e(c) = predicted variance of r(c) due to sampling error; \u03C3\u00B2_art(c) = predicted variance of r(c) due to artifacts; \u03C3\u00B2_pre(c) = total predicted variance of r(c); \u03C3\u00B2_\u03C1 = residual variance of \u03C1; "[show_var & verbose],

                                                  "CI = confidence interval around \u03C1\u0305; "[show_conf],
                                                  "CV = credibility interval around \u03C1\u0305; "[show_cred],

                                                  "correlations corrected individually."),

                                  ic_vgy = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "r\u0305 = mean observed correlation; "[show_msd],
                                                  "SE_r\u0305 = standard error of r\u0305; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_art = predicted SD_r due to artifacts; SD_pre = total predicted SD_r; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_res = residual variance of r ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_e = predicted variance of r due to sampling error; \u03C3\u00B2_art = predicted variance of r due to artifacts; \u03C3\u00B2_pre = total predicted variance of r; \u03C3\u00B2_res = residual variance of r; "[show_var & verbose],

                                                  "\u03C1\u0305 = mean operational validity (Y measured with error); "[show_msd],
                                                  "SE_\u03C1\u0305 = standard error of \u03C1\u0305; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_\u03C1 = residual standard deviation of \u03C1; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_e = predicted SD_r(c) due to sampling error; SD_art(c) = predicted SD_r(c) due to artifacts; SD_pre(c) = total predicted SD_r(c); SD_\u03C1 = residual standard deviation of \u03C1; "[show_msd & verbose],

                                                  "\u03C3\u00B2_r(c) = observed variance of r(c); \u03C3\u00B2_\u03C1 = residual variance of \u03C1 ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_r(c) = observed variance of r(c); \u03C3\u00B2_e(c) = predicted variance of r(c) due to sampling error; \u03C3\u00B2_art(c) = predicted variance of r(c) due to artifacts; \u03C3\u00B2_pre(c) = total predicted variance of r(c); \u03C3\u00B2_\u03C1 = residual variance of \u03C1; "[show_var & verbose],

                                                  "CI = confidence interval around \u03C1\u0305; "[show_conf],
                                                  "CV = credibility interval around \u03C1\u0305; "[show_cred],

                                                  "correlations corrected individually."),

                                  ad_ts  = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "r\u0305 = mean observed correlation; "[show_msd],
                                                  "SE_r\u0305 = standard error of r\u0305; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_art = predicted SD_r due to artifacts; SD_pre = total predicted SD_r; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_res = residual variance of r ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_e = predicted variance of r due to sampling error; \u03C3\u00B2_art = predicted variance of r due to artifacts; \u03C3\u00B2_pre = total predicted variance of r; \u03C3\u00B2_res = residual variance of r; "[show_var & verbose],

                                                  "\u03C1\u0305 = mean true-score correlation; "[show_msd],
                                                  "SE_\u03C1\u0305 = standard error of \u03C1\u0305; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_\u03C1 = residual standard deviation of \u03C1; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_e = predicted SD_r(c) due to sampling error; SD_art(c) = predicted SD_r(c) due to artifacts; SD_pre(c) = total predicted SD_r(c); SD_\u03C1 = residual standard deviation of \u03C1; "[show_msd & verbose],

                                                  "\u03C3\u00B2_r(c) = observed variance of r(c); \u03C3\u00B2_\u03C1 = residual variance of \u03C1 ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_r(c) = observed variance of r(c); \u03C3\u00B2_e(c) = predicted variance of r(c) due to sampling error; \u03C3\u00B2_art(c) = predicted variance of r(c) due to artifacts; \u03C3\u00B2_pre(c) = total predicted variance of r(c); \u03C3\u00B2_\u03C1 = residual variance of \u03C1; "[show_var & verbose],

                                                  "CI = confidence interval around \u03C1\u0305; "[show_conf],
                                                  "CV = credibility interval around \u03C1\u0305; "[show_cred],

                                                  "correlations corrected using artifact distributions."),

                                  ad_vgx = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "r\u0305 = mean observed correlation; "[show_msd],
                                                  "SE_r\u0305 = standard error of r\u0305; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_art = predicted SD_r due to artifacts; SD_pre = total predicted SD_r; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_res = residual variance of r ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_e = predicted variance of r due to sampling error; \u03C3\u00B2_art = predicted variance of r due to artifacts; \u03C3\u00B2_pre = total predicted variance of r; \u03C3\u00B2_res = residual variance of r; "[show_var & verbose],

                                                  "\u03C1\u0305 = mean operational validity (X measured with error); "[show_msd],
                                                  "SE_\u03C1\u0305 = standard error of \u03C1\u0305; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_\u03C1 = residual standard deviation of \u03C1; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_e = predicted SD_r(c) due to sampling error; SD_art(c) = predicted SD_r(c) due to artifacts; SD_pre(c) = total predicted SD_r(c); SD_\u03C1 = residual standard deviation of \u03C1; "[show_msd & verbose],

                                                  "\u03C3\u00B2_r(c) = observed variance of r(c); \u03C3\u00B2_\u03C1 = residual variance of \u03C1 ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_r(c) = observed variance of r(c); \u03C3\u00B2_e(c) = predicted variance of r(c) due to sampling error; \u03C3\u00B2_art(c) = predicted variance of r(c) due to artifacts; \u03C3\u00B2_pre(c) = total predicted variance of r(c); \u03C3\u00B2_\u03C1 = residual variance of \u03C1; "[show_var & verbose],

                                                  "CI = confidence interval around \u03C1\u0305; "[show_conf],
                                                  "CV = credibility interval around \u03C1\u0305; "[show_cred],

                                                  "correlations corrected using artifact distributions."),

                                  ad_vgy = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "r\u0305 = mean observed correlation; "[show_msd],
                                                  "SE_r\u0305 = standard error of r\u0305; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_art = predicted SD_r due to artifacts; SD_pre = total predicted SD_r; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_res = residual variance of r ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_r = observed variance of r; \u03C3\u00B2_e = predicted variance of r due to sampling error; \u03C3\u00B2_art = predicted variance of r due to artifacts; \u03C3\u00B2_pre = total predicted variance of r; \u03C3\u00B2_res = residual variance of r; "[show_var & verbose],

                                                  "\u03C1\u0305 = mean operational validity (Y measured with error); "[show_msd],
                                                  "SE_\u03C1\u0305 = standard error of \u03C1\u0305; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_\u03C1 = residual standard deviation of \u03C1; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_e = predicted SD_r(c) due to sampling error; SD_art(c) = predicted SD_r(c) due to artifacts; SD_pre(c) = total predicted SD_r(c); SD_\u03C1 = residual standard deviation of \u03C1; "[show_msd & verbose],

                                                  "\u03C3\u00B2_r(c) = observed variance of r(c); \u03C3\u00B2_\u03C1 = residual variance of \u03C1 ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_r(c) = observed variance of r(c); \u03C3\u00B2_e(c) = predicted variance of r(c) due to sampling error; \u03C3\u00B2_art(c) = predicted variance of r(c) due to artifacts; \u03C3\u00B2_pre(c) = total predicted variance of r(c); \u03C3\u00B2_\u03C1 = residual variance of \u03C1; "[show_var & verbose],

                                                  "CI = confidence interval around \u03C1\u0305; "[show_conf],
                                                  "CV = credibility interval around \u03C1\u0305; "[show_cred],

                                                  "correlations corrected using artifact distributions.")
                                )

                        } else if(es_type == "d") {
                                ### TODO: Don't refer to latent/observed groups if group membership reliability is not corrected.
                                c(bb     = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "d\u0305 =  mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[d\u0305] = standard error of d\u0305; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_res = residual variance of d ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_e = predicted variance of d due to sampling error; \u03C3\u00B2_res = residual variance of d; "[show_var & verbose],
                                                  "CI = confidence interval around d\u0305; "[show_conf],
                                                  "CV = credibility interval around d\u0305."[show_cred]),

                                  ic_ts  = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "d\u0305 = mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[d\u0305] = standard error of d\u0305; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_art = predicted SD_d due to artifacts; SD_pre = total predicted SD_d; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_res = residual variance of d ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_e = predicted variance of d due to sampling error; \u03C3\u00B2_art = predicted variance of d due to artifacts; \u03C3\u00B2_pre = total predicted variance of d; \u03C3\u00B2_res = residual variance of d; "[show_var & verbose],

                                                  "\u03B4\u0305 = mean true-score Cohen's d (Hedges' g) between latent groups; "[show_msd],
                                                  "SE_\u03B4\u0305 = standard error of \u03B4\u0305; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_\u03B4 = residual standard deviation of \u03B4; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_e = predicted SD_d(c) due to sampling error; SD_art(c) = predicted SD_d(c) due to artifacts; SD_pre(c) = total predicted SD_d(c); SD_\u03B4 = residual standard deviation of \u03B4; "[show_msd & verbose],

                                                  "\u03C3\u00B2_d(c) = observed variance of d(c); \u03C3\u00B2_\u03B4 = residual variance of \u03B4 ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_d(c) = observed variance of d(c); \u03C3\u00B2_e(c) = predicted variance of d(c) due to sampling error; \u03C3\u00B2_art(c) = predicted variance of d(c) due to artifacts; \u03C3\u00B2_pre(c) = total predicted variance of d(c); \u03C3\u00B2_\u03B4 = residual variance of \u03B4; "[show_var & verbose],

                                                  "CI = confidence interval around \u03B4\u0305; "[show_conf],
                                                  "CV = credibility interval around \u03B4\u0305; "[show_cred],

                                                  "effect sizes corrected individually."),

                                  ic_vgx = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "d\u0305 = mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[d\u0305] = standard error of d\u0305; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_art = predicted SD_d due to artifacts; SD_pre = total predicted SD_d; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_res = residual variance of d ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_e = predicted variance of d due to sampling error; \u03C3\u00B2_art = predicted variance of d due to artifacts; \u03C3\u00B2_pre = total predicted variance of d; \u03C3\u00B2_res = residual variance of d; "[show_var & verbose],

                                                  "\u03B4\u0305 = mean true-score Cohen's d (Hedges' g) between observed groups; "[show_msd],
                                                  "SE_\u03B4\u0305 = standard error of \u03B4\u0305; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_\u03B4 = residual standard deviation of \u03B4; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_e = predicted SD_d(c) due to sampling error; SD_art(c) = predicted SD_d(c) due to artifacts; SD_pre(c) = total predicted SD_d(c); SD_\u03B4 = residual standard deviation of \u03B4; "[show_msd & verbose],

                                                  "\u03C3\u00B2_d(c) = observed variance of d(c); \u03C3\u00B2_\u03B4 = residual variance of \u03B4 ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_d(c) = observed variance of d(c); \u03C3\u00B2_e(c) = predicted variance of d(c) due to sampling error; \u03C3\u00B2_art(c) = predicted variance of d(c) due to artifacts; \u03C3\u00B2_pre(c) = total predicted variance of d(c); \u03C3\u00B2_\u03B4 = residual variance of \u03B4; "[show_var & verbose],

                                                  "CI = confidence interval around \u03B4\u0305; "[show_conf],
                                                  "CV = credibility interval around \u03B4\u0305; "[show_cred],

                                                  "effect sizes corrected individually."),

                                  ic_vgy = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "d\u0305 = mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[d\u0305] = standard error of d\u0305; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_art = predicted SD_d due to artifacts; SD_pre = total predicted SD_d; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_res = residual variance of d ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_e = predicted variance of d due to sampling error; \u03C3\u00B2_art = predicted variance of d due to artifacts; \u03C3\u00B2_pre = total predicted variance of d; \u03C3\u00B2_res = residual variance of d; "[show_var & verbose],

                                                  "\u03B4\u0305 = mean observed Cohen's d (Hedges' g) between latent groups; "[show_msd],
                                                  "SE_\u03B4\u0305 = standard error of \u03B4\u0305; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_\u03B4 = residual standard deviation of \u03B4; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_e = predicted SD_d(c) due to sampling error; SD_art(c) = predicted SD_d(c) due to artifacts; SD_pre(c) = total predicted SD_d(c); SD_\u03B4 = residual standard deviation of \u03B4; "[show_msd & verbose],

                                                  "\u03C3\u00B2_d(c) = observed variance of d(c); \u03C3\u00B2_\u03B4 = residual variance of \u03B4 ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_d(c) = observed variance of d(c); \u03C3\u00B2_e(c) = predicted variance of d(c) due to sampling error; \u03C3\u00B2_art(c) = predicted variance of d(c) due to artifacts; \u03C3\u00B2_pre(c) = total predicted variance of d(c); \u03C3\u00B2_\u03B4 = residual variance of \u03B4; "[show_var & verbose],

                                                  "CI = confidence interval around \u03B4\u0305; "[show_conf],
                                                  "CV = credibility interval around \u03B4\u0305; "[show_cred],

                                                  "effect sizes corrected individually."),

                                  ad_ts  = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "d\u0305 = mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[d\u0305] = standard error of d\u0305; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_art = predicted SD_d due to artifacts; SD_pre = total predicted SD_d; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_res = residual variance of d ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_e = predicted variance of d due to sampling error; \u03C3\u00B2_art = predicted variance of d due to artifacts; \u03C3\u00B2_pre = total predicted variance of d; \u03C3\u00B2_res = residual variance of d; "[show_var & verbose],

                                                  "\u03B4\u0305 = mean true-score Cohen's d (Hedges' g) between latent groups; "[show_msd],
                                                  "SE_\u03B4\u0305 = standard error of \u03B4\u0305; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_\u03B4 = residual standard deviation of \u03B4; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_e = predicted SD_d(c) due to sampling error; SD_art(c) = predicted SD_d(c) due to artifacts; SD_pre(c) = total predicted SD_d(c); SD_\u03B4 = residual standard deviation of \u03B4; "[show_msd & verbose],

                                                  "\u03C3\u00B2_d(c) = observed variance of d(c); \u03C3\u00B2_\u03B4 = residual variance of \u03B4 ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_d(c) = observed variance of d(c); \u03C3\u00B2_e(c) = predicted variance of d(c) due to sampling error; \u03C3\u00B2_art(c) = predicted variance of d(c) due to artifacts; \u03C3\u00B2_pre(c) = total predicted variance of d(c); \u03C3\u00B2_\u03B4 = residual variance of \u03B4; "[show_var & verbose],

                                                  "CI = confidence interval around \u03B4\u0305; "[show_conf],
                                                  "CV = credibility interval around \u03B4\u0305; "[show_cred],

                                                  "effect sizes corrected using artifact distributions."),

                                  ad_vgx = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "d\u0305 = mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[d\u0305] = standard error of d\u0305; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_art = predicted SD_d due to artifacts; SD_pre = total predicted SD_d; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_res = residual variance of d ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_e = predicted variance of d due to sampling error; \u03C3\u00B2_art = predicted variance of d due to artifacts; \u03C3\u00B2_pre = total predicted variance of d; \u03C3\u00B2_res = residual variance of d; "[show_var & verbose],

                                                  "\u03B4\u0305 = mean true-score Cohen's d (Hedges' g) between observed groups; "[show_msd],
                                                  "SE_\u03B4\u0305 = standard error of \u03B4\u0305; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_\u03B4 = residual standard deviation of \u03B4; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_e = predicted SD_d(c) due to sampling error; SD_art(c) = predicted SD_d(c) due to artifacts; SD_pre(c) = total predicted SD_d(c); SD_\u03B4 = residual standard deviation of \u03B4; "[show_msd & verbose],

                                                  "\u03C3\u00B2_d(c) = observed variance of d(c); \u03C3\u00B2_\u03B4 = residual variance of \u03B4 ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_d(c) = observed variance of d(c); \u03C3\u00B2_e(c) = predicted variance of d(c) due to sampling error; \u03C3\u00B2_art(c) = predicted variance of d(c) due to artifacts; \u03C3\u00B2_pre(c) = total predicted variance of d(c); \u03C3\u00B2_\u03B4 = residual variance of \u03B4; "[show_var & verbose],

                                                  "CI = confidence interval around \u03B4\u0305; "[show_conf],
                                                  "CV = credibility interval around \u03B4\u0305; "[show_cred],

                                                  "effect sizes corrected using artifact distributions."),

                                  ad_vgy = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "d\u0305 = mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[d\u0305] = standard error of d\u0305; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_art = predicted SD_d due to artifacts; SD_pre = total predicted SD_d; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_res = residual variance of d ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_d = observed variance of d; \u03C3\u00B2_e = predicted variance of d due to sampling error; \u03C3\u00B2_art = predicted variance of d due to artifacts; \u03C3\u00B2_pre = total predicted variance of d; \u03C3\u00B2_res = residual variance of d; "[show_var & verbose],

                                                  "\u03B4\u0305 = mean observed Cohen's d (Hedges' g) between latent groups; "[show_msd],
                                                  "SE_\u03B4\u0305 = standard error of \u03B4\u0305; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_\u03B4 = residual standard deviation of \u03B4; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_e = predicted SD_d(c) due to sampling error; SD_art(c) = predicted SD_d(c) due to artifacts; SD_pre(c) = total predicted SD_d(c); SD_\u03B4 = residual standard deviation of \u03B4; "[show_msd & verbose],

                                                  "\u03C3\u00B2_d(c) = observed variance of d(c); \u03C3\u00B2_\u03B4 = residual variance of \u03B4 ;"[show_var & !verbose],
                                                  "\u03C3\u00B2_d(c) = observed variance of d(c); \u03C3\u00B2_e(c) = predicted variance of d(c) due to sampling error; \u03C3\u00B2_art(c) = predicted variance of d(c) due to artifacts; \u03C3\u00B2_pre(c) = total predicted variance of d(c); \u03C3\u00B2_\u03B4 = residual variance of \u03B4; "[show_var & verbose],

                                                  "CI = confidence interval around \u03B4\u0305; "[show_conf],
                                                  "CV = credibility interval around \u03B4\u0305; "[show_cred],

                                                  "effect sizes corrected using artifact distributions.")
                                )

                        } else {
                                c(bb = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                              paste0(paste0(strsplit(symbol_es, NULL)[[1]], "\u305", collapse = ""), " =  mean observed effect size (", symbol_es, "); ")[show_msd],
                                              paste0("SE_", paste0(strsplit(symbol_es, NULL)[[1]], "\u305", collapse = ""), " = standard error of ", paste0(strsplit(symbol_es, NULL)[[1]], "\u305", collapse = ""), "; ")[show_se],
                                              paste0("SD_", symbol_es, " = observed standard deviation of ", symbol_es, "; SD_res = residual standard deviation of ", symbol_es, "; ")[show_msd & !verbose],
                                              paste0("SD_", symbol_es, " = observed standard deviation of ", symbol_es, "; SD_e = predicted SD_", symbol_es, " due to sampling error; SD_res = residual standard deviation of ", symbol_es, "; ")[show_msd & verbose],
                                              paste0("\u03C3\u00B2_", symbol_es, " = observed variance of ", symbol_es, "; \u03C3\u00B2_res = residual variance of ", symbol_es, " ;")[show_var & !verbose],
                                              paste0("\u03C3\u00B2_", symbol_es, " = observed variance of ", symbol_es, "; \u03C3\u00B2_e = predicted variance of ", symbol_es, " due to sampling error; \u03C3\u00B2_res = residual variance of ", symbol_es, "; ")[show_var & verbose],
                                              paste0("CI = confidence interval around ", paste0(strsplit(symbol_es, NULL)[[1]], "\u305", collapse = ""), "; ")[show_conf],
                                              paste0("CV = credibility interval around ", paste0(strsplit(symbol_es, NULL)[[1]], "\u305", collapse = ""), ".")[show_cred])
                                )
                        }

                } else {
                        if(es_type == "r") {
                                c(bb     = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean r =  mean observed correlation; "[show_msd],
                                                  "SE_[Mean r] = standard error of Mean r; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "Var_r = observed variance of r; Var_res = residual variance of r ;"[show_var & !verbose],
                                                  "Var_r = observed variance of r; Var_e = predicted variance of r due to sampling error; Var_res = residual variance of r; "[show_var & verbose],
                                                  "CI = confidence interval around Mean r; "[show_conf],
                                                  "CV = credibility interval around Mean r."[show_cred]),

                                  ic_ts  = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean r = mean observed correlation; "[show_msd],
                                                  "SE_[Mean r] = standard error of Mean r; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_art = predicted SD_r due to artifacts; SD_pre = total predicted SD_r; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "Var_r = observed variance of r; Var_res = residual variance of r ;"[show_var & !verbose],
                                                  "Var_r = observed variance of r; Var_e = predicted variance of r due to sampling error; Var_art = predicted variance of r due to artifacts; Var_pre = total predicted variance of r; Var_res = residual variance of r; "[show_var & verbose],

                                                  "Mean rho = mean true-score correlation; "[show_msd],
                                                  "SE_[Mean rho] = standard error of Mean rho; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_rho = residual standard deviation of rho; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_e = predicted SD_r(c) due to sampling error; SD_art(c) = predicted SD_r(c) due to artifacts; SD_pre(c) = total predicted SD_r(c); SD_rho = residual standard deviation of rho; "[show_msd & verbose],

                                                  "Var_r(c) = observed variance of r(c); Var_rho = residual variance of rho ;"[show_var & !verbose],
                                                  "Var_r(c) = observed variance of r(c); Var_e(c) = predicted variance of r(c) due to sampling error; Var_art(c) = predicted variance of r(c) due to artifacts; Var_pre(c) = total predicted variance of r(c); Var_rho = residual variance of rho; "[show_var & verbose],

                                                  "CI = confidence interval around Mean rho; "[show_conf],
                                                  "CV = credibility interval around Mean rho; "[show_cred],

                                                  "correlations corrected individually."),

                                  ic_vgx = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean r = mean observed correlation; "[show_msd],
                                                  "SE_[Mean r] = standard error of Mean r; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_art = predicted SD_r due to artifacts; SD_pre = total predicted SD_r; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "Var_r = observed variance of r; Var_res = residual variance of r ;"[show_var & !verbose],
                                                  "Var_r = observed variance of r; Var_e = predicted variance of r due to sampling error; Var_art = predicted variance of r due to artifacts; Var_pre = total predicted variance of r; Var_res = residual variance of r; "[show_var & verbose],

                                                  "Mean rho = mean operational validity (X measured with error); "[show_msd],
                                                  "SE_[Mean rho] = standard error of Mean rho; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_rho = residual standard deviation of rho; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_e = predicted SD_r(c) due to sampling error; SD_art(c) = predicted SD_r(c) due to artifacts; SD_pre(c) = total predicted SD_r(c); SD_rho = residual standard deviation of rho; "[show_msd & verbose],

                                                  "Var_r(c) = observed variance of r(c); Var_rho = residual variance of rho ;"[show_var & !verbose],
                                                  "Var_r(c) = observed variance of r(c); Var_e(c) = predicted variance of r(c) due to sampling error; Var_art(c) = predicted variance of r(c) due to artifacts; Var_pre(c) = total predicted variance of r(c); Var_rho = residual variance of rho; "[show_var & verbose],

                                                  "CI = confidence interval around Mean rho; "[show_conf],
                                                  "CV = credibility interval around Mean rho; "[show_cred],

                                                  "correlations corrected individually."),

                                  ic_vgy = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean r = mean observed correlation; "[show_msd],
                                                  "SE_[Mean r] = standard error of Mean r; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_art = predicted SD_r due to artifacts; SD_pre = total predicted SD_r; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "Var_r = observed variance of r; Var_res = residual variance of r ;"[show_var & !verbose],
                                                  "Var_r = observed variance of r; Var_e = predicted variance of r due to sampling error; Var_art = predicted variance of r due to artifacts; Var_pre = total predicted variance of r; Var_res = residual variance of r; "[show_var & verbose],

                                                  "Mean rho = mean operational validity (Y measured with error); "[show_msd],
                                                  "SE_[Mean rho] = standard error of Mean rho; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_rho = residual standard deviation of rho; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_e = predicted SD_r(c) due to sampling error; SD_art(c) = predicted SD_r(c) due to artifacts; SD_pre(c) = total predicted SD_r(c); SD_rho = residual standard deviation of rho; "[show_msd & verbose],

                                                  "Var_r(c) = observed variance of r(c); Var_rho = residual variance of rho ;"[show_var & !verbose],
                                                  "Var_r(c) = observed variance of r(c); Var_e(c) = predicted variance of r(c) due to sampling error; Var_art(c) = predicted variance of r(c) due to artifacts; Var_pre(c) = total predicted variance of r(c); Var_rho = residual variance of rho; "[show_var & verbose],

                                                  "CI = confidence interval around Mean rho; "[show_conf],
                                                  "CV = credibility interval around Mean rho; "[show_cred],

                                                  "correlations corrected individually."),

                                  ad_ts  = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean r = mean observed correlation; "[show_msd],
                                                  "SE_[Mean r] = standard error of Mean r; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_art = predicted SD_r due to artifacts; SD_pre = total predicted SD_r; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "Var_r = observed variance of r; Var_res = residual variance of r ;"[show_var & !verbose],
                                                  "Var_r = observed variance of r; Var_e = predicted variance of r due to sampling error; Var_art = predicted variance of r due to artifacts; Var_pre = total predicted variance of r; Var_res = residual variance of r; "[show_var & verbose],

                                                  "Mean rho = mean true-score correlation; "[show_msd],
                                                  "SE_[Mean rho] = standard error of Mean rho; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_rho = residual standard deviation of rho; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_e = predicted SD_r(c) due to sampling error; SD_art(c) = predicted SD_r(c) due to artifacts; SD_pre(c) = total predicted SD_r(c); SD_rho = residual standard deviation of rho; "[show_msd & verbose],

                                                  "Var_r(c) = observed variance of r(c); Var_rho = residual variance of rho ;"[show_var & !verbose],
                                                  "Var_r(c) = observed variance of r(c); Var_e(c) = predicted variance of r(c) due to sampling error; Var_art(c) = predicted variance of r(c) due to artifacts; Var_pre(c) = total predicted variance of r(c); Var_rho = residual variance of rho; "[show_var & verbose],

                                                  "CI = confidence interval around Mean rho; "[show_conf],
                                                  "CV = credibility interval around Mean rho; "[show_cred],

                                                  "correlations corrected using artifact distributions."),

                                  ad_vgx = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean r = mean observed correlation; "[show_msd],
                                                  "SE_[Mean r] = standard error of Mean r; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_art = predicted SD_r due to artifacts; SD_pre = total predicted SD_r; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "Var_r = observed variance of r; Var_res = residual variance of r ;"[show_var & !verbose],
                                                  "Var_r = observed variance of r; Var_e = predicted variance of r due to sampling error; Var_art = predicted variance of r due to artifacts; Var_pre = total predicted variance of r; Var_res = residual variance of r; "[show_var & verbose],

                                                  "Mean rho = mean operational validity (X measured with error); "[show_msd],
                                                  "SE_[Mean rho] = standard error of Mean rho; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_rho = residual standard deviation of rho; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_e = predicted SD_r(c) due to sampling error; SD_art(c) = predicted SD_r(c) due to artifacts; SD_pre(c) = total predicted SD_r(c); SD_rho = residual standard deviation of rho; "[show_msd & verbose],

                                                  "Var_r(c) = observed variance of r(c); Var_rho = residual variance of rho ;"[show_var & !verbose],
                                                  "Var_r(c) = observed variance of r(c); Var_e(c) = predicted variance of r(c) due to sampling error; Var_art(c) = predicted variance of r(c) due to artifacts; Var_pre(c) = total predicted variance of r(c); Var_rho = residual variance of rho; "[show_var & verbose],

                                                  "CI = confidence interval around Mean rho; "[show_conf],
                                                  "CV = credibility interval around Mean rho; "[show_cred],

                                                  "correlations corrected using artifact distributions."),

                                  ad_vgy = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean r = mean observed correlation; "[show_msd],
                                                  "SE_[Mean r] = standard error of Mean r; "[show_se],
                                                  "SD_r = observed standard deviation of r; SD_res = residual standard deviation of r; "[show_msd & !verbose],
                                                  "SD_r = observed standard deviation of r; SD_e = predicted SD_r due to sampling error; SD_art = predicted SD_r due to artifacts; SD_pre = total predicted SD_r; SD_res = residual standard deviation of r; "[show_msd & verbose],
                                                  "Var_r = observed variance of r; Var_res = residual variance of r ;"[show_var & !verbose],
                                                  "Var_r = observed variance of r; Var_e = predicted variance of r due to sampling error; Var_art = predicted variance of r due to artifacts; Var_pre = total predicted variance of r; Var_res = residual variance of r; "[show_var & verbose],

                                                  "Mean rho = mean operational validity (Y measured with error); "[show_msd],
                                                  "SE_[Mean rho] = standard error of Mean rho; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_rho = residual standard deviation of rho; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected correlations (r(c)); SD_e = predicted SD_r(c) due to sampling error; SD_art(c) = predicted SD_r(c) due to artifacts; SD_pre(c) = total predicted SD_r(c); SD_rho = residual standard deviation of rho; "[show_msd & verbose],

                                                  "Var_r(c) = observed variance of r(c); Var_rho = residual variance of rho ;"[show_var & !verbose],
                                                  "Var_r(c) = observed variance of r(c); Var_e(c) = predicted variance of r(c) due to sampling error; Var_art(c) = predicted variance of r(c) due to artifacts; Var_pre(c) = total predicted variance of r(c); Var_rho = residual variance of rho; "[show_var & verbose],

                                                  "CI = confidence interval around Mean rho; "[show_conf],
                                                  "CV = credibility interval around Mean rho; "[show_cred],

                                                  "correlations corrected using artifact distributions.")
                                )

                        } else if(es_type == "d") {
                                ### TODO: Don't refer to latent/observed groups if group membership reliability is not corrected.
                                c(bb     = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean d =  mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[Mean d] = standard error of Mean d; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "Var_d = observed variance of d; Var_res = residual variance of d ;"[show_var & !verbose],
                                                  "Var_d = observed variance of d; Var_e = predicted variance of d due to sampling error; Var_res = residual variance of d; "[show_var & verbose],
                                                  "CI = confidence interval around Mean d; "[show_conf],
                                                  "CV = credibility interval around Mean d."[show_cred]),

                                  ic_ts  = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean d = mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[Mean d] = standard error of Mean d; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_art = predicted SD_d due to artifacts; SD_pre = total predicted SD_d; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "Var_d = observed variance of d; Var_res = residual variance of d ;"[show_var & !verbose],
                                                  "Var_d = observed variance of d; Var_e = predicted variance of d due to sampling error; Var_art = predicted variance of d due to artifacts; Var_pre = total predicted variance of d; Var_res = residual variance of d; "[show_var & verbose],

                                                  "Mean delta = mean true-score Cohen's d (Hedges' g) between latent groups; "[show_msd],
                                                  "SE_[Mean delta] = standard error of Mean delta; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_delta = residual standard deviation of delta; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_e = predicted SD_d(c) due to sampling error; SD_art(c) = predicted SD_d(c) due to artifacts; SD_pre(c) = total predicted SD_d(c); SD_delta = residual standard deviation of delta; "[show_msd & verbose],

                                                  "Var_d(c) = observed variance of d(c); Var_delta = residual variance of delta ;"[show_var & !verbose],
                                                  "Var_d(c) = observed variance of d(c); Var_e(c) = predicted variance of d(c) due to sampling error; Var_art(c) = predicted variance of d(c) due to artifacts; Var_pre(c) = total predicted variance of d(c); Var_delta = residual variance of delta; "[show_var & verbose],

                                                  "CI = confidence interval around Mean delta; "[show_conf],
                                                  "CV = credibility interval around Mean delta; "[show_cred],

                                                  "effect sizes corrected individually."),

                                  ic_vgx = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean d = mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[Mean d] = standard error of Mean d; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_art = predicted SD_d due to artifacts; SD_pre = total predicted SD_d; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "Var_d = observed variance of d; Var_res = residual variance of d ;"[show_var & !verbose],
                                                  "Var_d = observed variance of d; Var_e = predicted variance of d due to sampling error; Var_art = predicted variance of d due to artifacts; Var_pre = total predicted variance of d; Var_res = residual variance of d; "[show_var & verbose],

                                                  "Mean delta = mean true-score Cohen's d (Hedges' g) between observed groups; "[show_msd],
                                                  "SE_[Mean delta] = standard error of Mean delta; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_delta = residual standard deviation of delta; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_e = predicted SD_d(c) due to sampling error; SD_art(c) = predicted SD_d(c) due to artifacts; SD_pre(c) = total predicted SD_d(c); SD_delta = residual standard deviation of delta; "[show_msd & verbose],

                                                  "Var_d(c) = observed variance of d(c); Var_delta = residual variance of delta ;"[show_var & !verbose],
                                                  "Var_d(c) = observed variance of d(c); Var_e(c) = predicted variance of d(c) due to sampling error; Var_art(c) = predicted variance of d(c) due to artifacts; Var_pre(c) = total predicted variance of d(c); Var_delta = residual variance of delta; "[show_var & verbose],

                                                  "CI = confidence interval around Mean delta; "[show_conf],
                                                  "CV = credibility interval around Mean delta; "[show_cred],

                                                  "effect sizes corrected individually."),

                                  ic_vgy = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean d = mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[Mean d] = standard error of Mean d; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_art = predicted SD_d due to artifacts; SD_pre = total predicted SD_d; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "Var_d = observed variance of d; Var_res = residual variance of d ;"[show_var & !verbose],
                                                  "Var_d = observed variance of d; Var_e = predicted variance of d due to sampling error; Var_art = predicted variance of d due to artifacts; Var_pre = total predicted variance of d; Var_res = residual variance of d; "[show_var & verbose],

                                                  "Mean delta = mean observed Cohen's d (Hedges' g) between latent groups; "[show_msd],
                                                  "SE_[Mean delta] = standard error of Mean delta; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_delta = residual standard deviation of delta; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_e = predicted SD_d(c) due to sampling error; SD_art(c) = predicted SD_d(c) due to artifacts; SD_pre(c) = total predicted SD_d(c); SD_delta = residual standard deviation of delta; "[show_msd & verbose],

                                                  "Var_d(c) = observed variance of d(c); Var_delta = residual variance of delta ;"[show_var & !verbose],
                                                  "Var_d(c) = observed variance of d(c); Var_e(c) = predicted variance of d(c) due to sampling error; Var_art(c) = predicted variance of d(c) due to artifacts; Var_pre(c) = total predicted variance of d(c); Var_delta = residual variance of delta; "[show_var & verbose],

                                                  "CI = confidence interval around Mean delta; "[show_conf],
                                                  "CV = credibility interval around Mean delta; "[show_cred],

                                                  "effect sizes corrected individually."),

                                  ad_ts  = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean d = mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[Mean d] = standard error of Mean d; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_art = predicted SD_d due to artifacts; SD_pre = total predicted SD_d; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "Var_d = observed variance of d; Var_res = residual variance of d ;"[show_var & !verbose],
                                                  "Var_d = observed variance of d; Var_e = predicted variance of d due to sampling error; Var_art = predicted variance of d due to artifacts; Var_pre = total predicted variance of d; Var_res = residual variance of d; "[show_var & verbose],

                                                  "Mean delta = mean true-score Cohen's d (Hedges' g) between latent groups; "[show_msd],
                                                  "SE_[Mean delta] = standard error of Mean delta; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_delta = residual standard deviation of delta; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_e = predicted SD_d(c) due to sampling error; SD_art(c) = predicted SD_d(c) due to artifacts; SD_pre(c) = total predicted SD_d(c); SD_delta = residual standard deviation of delta; "[show_msd & verbose],

                                                  "Var_d(c) = observed variance of d(c); Var_delta = residual variance of delta ;"[show_var & !verbose],
                                                  "Var_d(c) = observed variance of d(c); Var_e(c) = predicted variance of d(c) due to sampling error; Var_art(c) = predicted variance of d(c) due to artifacts; Var_pre(c) = total predicted variance of d(c); Var_delta = residual variance of delta; "[show_var & verbose],

                                                  "CI = confidence interval around Mean delta; "[show_conf],
                                                  "CV = credibility interval around Mean delta; "[show_cred],

                                                  "effect sizes corrected using artifact distributions."),

                                  ad_vgx = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean d = mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[Mean d] = standard error of Mean d; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_art = predicted SD_d due to artifacts; SD_pre = total predicted SD_d; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "Var_d = observed variance of d; Var_res = residual variance of d ;"[show_var & !verbose],
                                                  "Var_d = observed variance of d; Var_e = predicted variance of d due to sampling error; Var_art = predicted variance of d due to artifacts; Var_pre = total predicted variance of d; Var_res = residual variance of d; "[show_var & verbose],

                                                  "Mean delta = mean true-score Cohen's d (Hedges' g) between observed groups; "[show_msd],
                                                  "SE_[Mean delta] = standard error of Mean delta; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_delta = residual standard deviation of delta; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_e = predicted SD_d(c) due to sampling error; SD_art(c) = predicted SD_d(c) due to artifacts; SD_pre(c) = total predicted SD_d(c); SD_delta = residual standard deviation of delta; "[show_msd & verbose],

                                                  "Var_d(c) = observed variance of d(c); Var_delta = residual variance of delta ;"[show_var & !verbose],
                                                  "Var_d(c) = observed variance of d(c); Var_e(c) = predicted variance of d(c) due to sampling error; Var_art(c) = predicted variance of d(c) due to artifacts; Var_pre(c) = total predicted variance of d(c); Var_delta = residual variance of delta; "[show_var & verbose],

                                                  "CI = confidence interval around Mean delta; "[show_conf],
                                                  "CV = credibility interval around Mean delta; "[show_cred],

                                                  "effect sizes corrected using artifact distributions."),

                                  ad_vgy = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                                  "Mean d = mean observed Cohen's d (Hedges' g); "[show_msd],
                                                  "SE_[Mean d] = standard error of Mean d; "[show_se],
                                                  "SD_d = observed standard deviation of d; SD_res = residual standard deviation of d; "[show_msd & !verbose],
                                                  "SD_d = observed standard deviation of d; SD_e = predicted SD_d due to sampling error; SD_art = predicted SD_d due to artifacts; SD_pre = total predicted SD_d; SD_res = residual standard deviation of d; "[show_msd & verbose],
                                                  "Var_d = observed variance of d; Var_res = residual variance of d ;"[show_var & !verbose],
                                                  "Var_d = observed variance of d; Var_e = predicted variance of d due to sampling error; Var_art = predicted variance of d due to artifacts; Var_pre = total predicted variance of d; Var_res = residual variance of d; "[show_var & verbose],

                                                  "Mean delta = mean observed Cohen's d (Hedges' g) between latent groups; "[show_msd],
                                                  "SE_[Mean delta] = standard error of Mean delta; "[show_se],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_delta = residual standard deviation of delta; "[show_msd & !verbose],
                                                  "SD_r(c) = observed standard deviation of corrected d values (d(c)); SD_e = predicted SD_d(c) due to sampling error; SD_art(c) = predicted SD_d(c) due to artifacts; SD_pre(c) = total predicted SD_d(c); SD_delta = residual standard deviation of delta; "[show_msd & verbose],

                                                  "Var_d(c) = observed variance of d(c); Var_delta = residual variance of delta ;"[show_var & !verbose],
                                                  "Var_d(c) = observed variance of d(c); Var_e(c) = predicted variance of d(c) due to sampling error; Var_art(c) = predicted variance of d(c) due to artifacts; Var_pre(c) = total predicted variance of d(c); Var_delta = residual variance of delta; "[show_var & verbose],

                                                  "CI = confidence interval around Mean delta; "[show_conf],
                                                  "CV = credibility interval around Mean delta; "[show_cred],

                                                  "effect sizes corrected using artifact distributions.")
                                )

                        } else {
                                c(bb = paste0("k = number of studies contributing to meta-analysis; N = total sample size; ",
                                              "Mean ", symbol_es, " =  mean observed effect size (", symbol_es, "); "[show_msd],
                                              "SE_[Mean ", symbol_es, "] = standard error of [Mean ", symbol_es, "]; "[show_se],
                                              "SD_", symbol_es, " = observed standard deviation of ", symbol_es, "; SD_res = residual standard deviation of ", symbol_es, "; "[show_msd & !verbose],
                                              "SD_", symbol_es, " = observed standard deviation of ", symbol_es, "; SD_e = predicted SD_", symbol_es, " due to sampling error; SD_res = residual standard deviation of ", symbol_es, "; "[show_msd & verbose],
                                              "Var_", symbol_es, " = observed variance of ", symbol_es, "; Var_res = residual variance of ", symbol_es, " ;"[show_var & !verbose],
                                              "Var_", symbol_es, " = observed variance of ", symbol_es, "; Var_e = predicted variance of ", symbol_es, " due to sampling error; Var_res = residual variance of ", symbol_es, "; "[show_var & verbose],
                                              "CI = confidence interval around Mean [", symbol_es, "]; "[show_conf],
                                              "CV = credibility interval around Mean [", symbol_es, "]."[show_cred])
                                )
                        }
                }
        } else {
                if(es_type == "r") {
                        c(bb     = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{r}$\u00a0=\u00a0 mean observed correlation; "[show_msd],
                                          "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],
                                          "CI\u00a0=\u00a0confidence interval around $\\overline{r}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\overline{r}$."[show_cred]),

                          ic_ts  = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{r}$\u00a0=\u00a0mean observed correlation; "[show_msd],
                                          "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{r}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{r}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $r$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],

                                          "$\\symup{\\overline{\\rho}}$\u00a0=\u00a0mean true-score correlation; "[show_msd],
                                          "$SE_{\\symup{\\overline{\\rho}}}$\u00a0=\u00a0standard error of $\\symup{\\overline{\\rho}}$; "[show_se],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{\\symup{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\symup{\\rho}$; "[show_msd & !verbose],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{r_{c}}$; $SD_{\\symup{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\symup{\\rho}$; "[show_msd & verbose],

                                          "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{\\symup{\\rho}}$\u00a0=\u00a0residual variance of $\\symup{\\rho}$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $r_{c}$; $\\sigma^{2}_{\\symup{\\rho}}$\u00a0=\u00a0residual variance of $\\symup{\\rho}$; "[show_var & verbose],

                                          "CI\u00a0=\u00a0confidence interval around $\\symup{\\overline{\\rho}}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\symup{\\overline{\\rho}}$; "[show_cred],

                                          "correlations corrected individually."),

                          ic_vgx = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{r}$\u00a0=\u00a0mean observed correlation; "[show_msd],
                                          "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{r}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{r}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $r$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],

                                          "$\\symup{\\overline{\\rho}}$\u00a0=\u00a0mean operational validity (*X* measured with error); "[show_msd],
                                          "$SE_{\\symup{\\overline{\\rho}}}$\u00a0=\u00a0standard error of $\\symup{\\overline{\\rho}}$; "[show_se],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{\\symup{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\symup{\\rho}$; "[show_msd & !verbose],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{r_{c}}$; $SD_{\\symup{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\symup{\\rho}$; "[show_msd & verbose],

                                          "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{\\symup{\\rho}}$\u00a0=\u00a0residual variance of $\\symup{\\rho}$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $r_{c}$; $\\sigma^{2}_{\\symup{\\rho}}$\u00a0=\u00a0residual variance of $\\symup{\\rho}$; "[show_var & verbose],

                                          "CI\u00a0=\u00a0confidence interval around $\\symup{\\overline{\\rho}}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\symup{\\overline{\\rho}}$; "[show_cred],

                                          "correlations corrected individually."),

                          ic_vgy = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{r}$\u00a0=\u00a0mean observed correlation; "[show_msd],
                                          "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{r}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{r}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $r$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],

                                          "$\\symup{\\overline{\\rho}}$\u00a0=\u00a0mean operational validity (*Y* measured with error); "[show_msd],
                                          "$SE_{\\symup{\\overline{\\rho}}}$\u00a0=\u00a0standard error of $\\symup{\\overline{\\rho}}$; "[show_se],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{\\symup{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\symup{\\rho}$; "[show_msd & !verbose],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{r_{c}}$; $SD_{\\symup{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\symup{\\rho}$; "[show_msd & verbose],

                                          "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{\\symup{\\rho}}$\u00a0=\u00a0residual variance of $\\symup{\\rho}$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $r_{c}$; $\\sigma^{2}_{\\symup{\\rho}}$\u00a0=\u00a0residual variance of $\\symup{\\rho}$; "[show_var & verbose],

                                          "CI\u00a0=\u00a0confidence interval around $\\symup{\\overline{\\rho}}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\symup{\\overline{\\rho}}$; "[show_cred],

                                          "correlations corrected individually."),

                          ad_ts  = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{r}$\u00a0=\u00a0mean observed correlation; "[show_msd],
                                          "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{r}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{r}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $r$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],

                                          "$\\symup{\\overline{\\rho}}$\u00a0=\u00a0mean true-score correlation; "[show_msd],
                                          "$SE_{\\symup{\\overline{\\rho}}}$\u00a0=\u00a0standard error of $\\symup{\\overline{\\rho}}$; "[show_se],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{\\symup{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\symup{\\rho}$; "[show_msd & !verbose],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{r_{c}}$; $SD_{\\symup{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\symup{\\rho}$; "[show_msd & verbose],

                                          "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{\\symup{\\rho}}$\u00a0=\u00a0residual variance of $\\symup{\\rho}$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $r_{c}$; $\\sigma^{2}_{\\symup{\\rho}}$\u00a0=\u00a0residual variance of $\\symup{\\rho}$; "[show_var & verbose],

                                          "CI\u00a0=\u00a0confidence interval around $\\symup{\\overline{\\rho}}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\symup{\\overline{\\rho}}$; "[show_cred],

                                          "correlations corrected using artifact distributions."),

                          ad_vgx = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{r}$\u00a0=\u00a0mean observed correlation; "[show_msd],
                                          "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{r}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{r}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $r$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],

                                          "$\\symup{\\overline{\\rho}}$\u00a0=\u00a0mean operational validity (*X* measured with error); "[show_msd],
                                          "$SE_{\\symup{\\overline{\\rho}}}$\u00a0=\u00a0standard error of $\\symup{\\overline{\\rho}}$; "[show_se],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{\\symup{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\symup{\\rho}$; "[show_msd & !verbose],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{r_{c}}$; $SD_{\\symup{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\symup{\\rho}$; "[show_msd & verbose],

                                          "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{\\symup{\\rho}}$\u00a0=\u00a0residual variance of $\\symup{\\rho}$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $r_{c}$; $\\sigma^{2}_{\\symup{\\rho}}$\u00a0=\u00a0residual variance of $\\symup{\\rho}$; "[show_var & verbose],

                                          "CI\u00a0=\u00a0confidence interval around $\\symup{\\overline{\\rho}}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\symup{\\overline{\\rho}}$; "[show_cred],

                                          "correlations corrected using artifact distributions."),

                          ad_vgy = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{r}$\u00a0=\u00a0mean observed correlation; "[show_msd],
                                          "$SE_{\\overline{r}}$\u00a0=\u00a0standard error of $\\overline{r}$; "[show_se],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & !verbose],
                                          "$SD_{r}$\u00a0=\u00a0observed standard deviation of $r$; $SD_{e}$\u00a0=\u00a0predicted $SD_{r}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{r}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{r}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $r$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{r}$\u00a0=\u00a0observed variance of $r$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $r$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $r$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $r$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $r$; "[show_var & verbose],

                                          "$\\symup{\\overline{\\rho}}$\u00a0=\u00a0mean operational validity (*Y* measured with error); "[show_msd],
                                          "$SE_{\\symup{\\overline{\\rho}}}$\u00a0=\u00a0standard error of $\\symup{\\overline{\\rho}}$; "[show_se],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{\\symup{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\symup{\\rho}$; "[show_msd & !verbose],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected correlations ($r_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{r_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{r_{c}}$; $SD_{\\symup{\\rho}}$\u00a0=\u00a0residual standard deviation of $\\symup{\\rho}$; "[show_msd & verbose],

                                          "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{\\symup{\\rho}}$\u00a0=\u00a0residual variance of $\\symup{\\rho}$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{r_{c}}$\u00a0=\u00a0observed variance of $r_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $r_{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $r_{c}$; $\\sigma^{2}_{\\symup{\\rho}}$\u00a0=\u00a0residual variance of $\\symup{\\rho}$; "[show_var & verbose],

                                          "CI\u00a0=\u00a0confidence interval around $\\symup{\\overline{\\rho}}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\symup{\\overline{\\rho}}$; "[show_cred],

                                          "correlations corrected using artifact distributions.")
                        )

                } else if(es_type == "d") {
                        ### f: Don't refer to latent/observed groups if group membership reliability is not corrected.
                        c(bb     = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{d}$\u00a0=\u00a0 mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                          "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],
                                          "CI\u00a0=\u00a0confidence interval around $\\overline{d}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\overline{d}$."[show_cred]),

                          ic_ts  = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{d}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                          "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{d}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{d}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $d$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],

                                          "$\\overline{\\delta}$\u00a0=\u00a0mean true-score Cohen's $d$ (Hedges' $g$) between latent groups; "[show_msd],
                                          "$SE_{\\overline{\\delta}}$\u00a0=\u00a0standard error of $\\overline{\\delta}$; "[show_se],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{\\delta}$\u00a0=\u00a0residual standard deviation of $\\delta$; "[show_msd & !verbose],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{d_{c}}$; $SD_{\\delta}$\u00a0=\u00a0residual standard deviation of $\\delta$; "[show_msd & verbose],

                                          "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{\\delta}$\u00a0=\u00a0residual variance of $\\delta$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $d_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $d{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $d_{c}$; $\\sigma^{2}_{\\delta}$\u00a0=\u00a0residual variance of $\\delta$; "[show_var & verbose],

                                          "CI\u00a0=\u00a0confidence interval around $\\overline{\\delta}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\overline{\\delta}$; "[show_cred],

                                          "effect sizes corrected individually."),

                          ic_vgx = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{d}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                          "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{d}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{d}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $d$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],

                                          "$\\overline{\\delta}$\u00a0=\u00a0mean true-score Cohen's $d$ (Hedges' $g$) between observed groups; "[show_msd],
                                          "$SE_{\\overline{\\delta}}$\u00a0=\u00a0standard error of $\\overline{\\delta}$; "[show_se],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{\\delta}$\u00a0=\u00a0residual standard deviation of $\\delta$; "[show_msd & !verbose],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{d_{c}}$; $SD_{\\delta}$\u00a0=\u00a0residual standard deviation of $\\delta$; "[show_msd & verbose],

                                          "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{\\delta}$\u00a0=\u00a0residual variance of $\\delta$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $d_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $d{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $d_{c}$; $\\sigma^{2}_{\\delta}$\u00a0=\u00a0residual variance of $\\delta$; "[show_var & verbose],

                                          "CI\u00a0=\u00a0confidence interval around $\\overline{\\delta}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\overline{\\delta}$; "[show_cred],

                                          "effect sizes corrected individually."),

                          ic_vgy = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{d}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                          "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{d}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{d}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $d$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],

                                          "$\\overline{\\delta}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$) between latent groups; "[show_msd],
                                          "$SE_{\\overline{\\delta}}$\u00a0=\u00a0standard error of $\\overline{\\delta}$; "[show_se],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{\\delta}$\u00a0=\u00a0residual standard deviation of $\\delta$; "[show_msd & !verbose],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{d_{c}}$; $SD_{\\delta}$\u00a0=\u00a0residual standard deviation of $\\delta$; "[show_msd & verbose],

                                          "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{\\delta}$\u00a0=\u00a0residual variance of $\\delta$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $d_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $d{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $d_{c}$; $\\sigma^{2}_{\\delta}$\u00a0=\u00a0residual variance of $\\delta$; "[show_var & verbose],

                                          "CI\u00a0=\u00a0confidence interval around $\\overline{\\delta}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\overline{\\delta}$; "[show_cred],

                                          "effect sizes corrected individually."),

                          ad_ts  = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{d}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                          "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{d}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{d}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $d$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],

                                          "$\\overline{\\delta}$\u00a0=\u00a0mean true-score Cohen's $d$ (Hedges' $g$) between latent groups; "[show_msd],
                                          "$SE_{\\overline{\\delta}}$\u00a0=\u00a0standard error of $\\overline{\\delta}$; "[show_se],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{\\delta}$\u00a0=\u00a0residual standard deviation of $\\delta$; "[show_msd & !verbose],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{d_{c}}$; $SD_{\\delta}$\u00a0=\u00a0residual standard deviation of $\\delta$; "[show_msd & verbose],

                                          "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{\\delta}$\u00a0=\u00a0residual variance of $\\delta$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $d_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $d{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $d_{c}$; $\\sigma^{2}_{\\delta}$\u00a0=\u00a0residual variance of $\\delta$; "[show_var & verbose],

                                          "CI\u00a0=\u00a0confidence interval around $\\overline{\\delta}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\overline{\\delta}$; "[show_cred],

                                          "effect sizes corrected using artifact distributions."),

                          ad_vgx = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{d}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                          "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{d}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{d}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $d$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],

                                          "$\\overline{\\delta}$\u00a0=\u00a0mean true-score Cohen's $d$ (Hedges' $g$) between observed groups; "[show_msd],
                                          "$SE_{\\overline{\\delta}}$\u00a0=\u00a0standard error of $\\overline{\\delta}$; "[show_se],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{\\delta}$\u00a0=\u00a0residual standard deviation of $\\delta$; "[show_msd & !verbose],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{d_{c}}$; $SD_{\\delta}$\u00a0=\u00a0residual standard deviation of $\\delta$; "[show_msd & verbose],

                                          "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{\\delta}$\u00a0=\u00a0residual variance of $\\delta$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $d_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $d{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $d_{c}$; $\\sigma^{2}_{\\delta}$\u00a0=\u00a0residual variance of $\\delta$; "[show_var & verbose],

                                          "CI\u00a0=\u00a0confidence interval around $\\overline{\\delta}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\overline{\\delta}$; "[show_cred],

                                          "effect sizes corrected using artifact distributions"),

                          ad_vgy = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                          "$\\overline{d}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$); "[show_msd],
                                          "$SE_{\\overline{d}}$\u00a0=\u00a0standard error of $\\overline{d}$; "[show_se],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & !verbose],
                                          "$SD_{d}$\u00a0=\u00a0observed standard deviation of $d$; $SD_{e}$\u00a0=\u00a0predicted $SD_{d}$ due to sampling error; $SD_{art}$\u00a0=\u00a0predicted $SD_{d}$ due to artifacts; $SD_{pre}$\u00a0=\u00a0total predicted $SD_{d}$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $d$; "[show_msd & verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{d}$\u00a0=\u00a0observed variance of $d$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $d$ due to sampling error; $\\sigma^{2}_{art}$\u00a0=\u00a0predicted variance of $d$ due to artifacts; $\\sigma^{2}_{pre}$\u00a0=\u00a0total predicted variance of $d$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $d$; "[show_var & verbose],

                                          "$\\overline{\\delta}$\u00a0=\u00a0mean observed Cohen's $d$ (Hedges' $g$) between latent groups; "[show_msd],
                                          "$SE_{\\overline{\\delta}}$\u00a0=\u00a0standard error of $\\overline{\\delta}$; "[show_se],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{\\delta}$\u00a0=\u00a0residual standard deviation of $\\delta$; "[show_msd & !verbose],
                                          "$SD_{r_{c}}$\u00a0=\u00a0observed standard deviation of corrected $d$ values ($d_{c}$); $SD_{e}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to sampling error; $SD_{art_{c}}$\u00a0=\u00a0predicted $SD_{d_{c}}$ due to artifacts; $SD_{pre_{c}}$\u00a0=\u00a0total predicted $SD_{d_{c}}$; $SD_{\\delta}$\u00a0=\u00a0residual standard deviation of $\\delta$; "[show_msd & verbose],

                                          "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{\\delta}$\u00a0=\u00a0residual variance of $\\delta$ ;"[show_var & !verbose],
                                          "$\\sigma^{2}_{d_{c}}$\u00a0=\u00a0observed variance of $d_{c}$; $\\sigma^{2}_{e_{c}}$\u00a0=\u00a0predicted variance of $d_{c}$ due to sampling error; $\\sigma^{2}_{art_{c}}$\u00a0=\u00a0predicted variance of $d{c}$ due to artifacts; $\\sigma^{2}_{pre_{c}}$\u00a0=\u00a0total predicted variance of $d_{c}$; $\\sigma^{2}_{\\delta}$\u00a0=\u00a0residual variance of $\\delta$; "[show_var & verbose],

                                          "CI\u00a0=\u00a0confidence interval around $\\overline{\\delta}$; "[show_conf],
                                          "CV\u00a0=\u00a0credibility interval around $\\overline{\\delta}$; "[show_cred],

                                          "effect sizes corrected using artifact distributions")
                        )

                } else {
                        c(bb = paste0("*k*\u00a0=\u00a0number of studies contributing to meta-analysis; *N*\u00a0=\u00a0total sample size; ",
                                      paste0("$\\overline{", symbol_es, "}$\u00a0=\u00a0 mean observed effect size ($", symbol_es, "$); ")[show_msd],
                                      paste0("$SE_{\\overline{", symbol_es, "}}$\u00a0=\u00a0standard error of $\\overline{", symbol_es, "}$; ")[show_se],
                                      paste0("$SD_{", symbol_es, "}$\u00a0=\u00a0observed standard deviation of $", symbol_es, "$; $SD_{res}$\u00a0=\u00a0residual standard deviation of $", symbol_es, "$; ")[show_msd & !verbose],
                                      paste0("$SD_{", symbol_es, "}$\u00a0=\u00a0observed standard deviation of $", symbol_es, "$; $SD_{e}$\u00a0=\u00a0predicted $SD_{", symbol_es, "}$ due to sampling error; $SD_{res}$\u00a0=\u00a0residual standard deviation of $", symbol_es, "$; ")[show_msd & verbose],
                                      paste0("$\\sigma^{2}_{", symbol_es, "}$\u00a0=\u00a0observed variance of $", symbol_es, "$; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $", symbol_es, "$ ;")[show_var & !verbose],
                                      paste0("$\\sigma^{2}_{", symbol_es, "}$\u00a0=\u00a0observed variance of $", symbol_es, "$; $\\sigma^{2}_{e}$\u00a0=\u00a0predicted variance of $", symbol_es, "$ due to sampling error; $\\sigma^{2}_{res}$\u00a0=\u00a0residual variance of $", symbol_es, "$; ")[show_var & verbose],
                                      paste0("CI\u00a0=\u00a0confidence interval around $\\overline{", symbol_es, "}$; ")[show_conf],
                                      paste0("CV\u00a0=\u00a0credibility interval around $\\overline{", symbol_es, "}$.")[show_cred])
                        )
                }
        }

        return(list(col_names = col_names,
                    footnote = footnote))
}
