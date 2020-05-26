#' Creates an .RData file for the convert_es function to be loaded in for testing.
#' Format of variables: actual__arguments__function

library(psychmeta)

# output_es = "r" ---------------------------------------------------------

actual__d_r__convert_es <- convert_es(es = 1, input_es = "d", output_es = "r", n1 = 100)

actual__d_rn2__convert_es <- convert_es(es = 1, input_es = "d", output_es = "r", n1 = 50, n2 = 50)

actual__F_r__convert_es <- convert_es(es = 10.3, input_es = "F", output_es = "r", n1 = 100, n2 = 150)

actual__t_r__convert_es <- convert_es(es = -1.3, input_es = "t", output_es = "r", n1 = 100, n2 = 140)

actual__chisq_r__convert_es <- convert_es(es = 1.3, input_es = "chisq", output_es = "r", n1 = 100, n2 = 100)

actual__p.chisq_r__convert_es <- convert_es(es = .021, input_es = "p.chisq", output_es = "r", n1 = 100, n2 = 100)

actual__or_r__convert_es <- convert_es(es = 4.37, input_es = "or", output_es = "r", n1 = 100, n2 = 100)

actual__lor_r__convert_es <- convert_es(es = 1.47, input_es = "lor", output_es = "r", n1 = 100, n2 = 100)

actual__r_r__convert_es <- convert_es(es = .3, input_es = "r", output_es = "r", n1 = 100)

# output_es = "d" ---------------------------------------------------------

actual__r_d__convert_es <- convert_es(es = .2, input_es = "r", output_es = "d", n1 = 100, n2 = 150)

actual__t_d__convert_es <- convert_es(es = -1.3, input_es = "t", output_es = "d", n1 = 100, n2 = 140)

actual__F_d__convert_es <- convert_es(es = 10.3, input_es = "F", output_es = "d", n1 = 100, n2 = 150)

actual__chisq_d__convert_es <- convert_es(es = 1.3, input_es = "chisq", output_es = "d", n1 = 100, n2 = 100)

actual__p.chisq_d__convert_es <- convert_es(es = .021, input_es = "p.chisq", output_es = "d", n1 = 100, n2 = 100)

actual__or_d__convert_es <- convert_es(es = 4.37, input_es = "or", output_es = "d", n1 = 100, n2 = 100)

actual__lor_d__convert_es <- convert_es(es = 1.47, input_es = "lor", output_es = "d", n1 = 100, n2 = 100)

actual__d_d__convert_es <- convert_es(es = .8, input_es = "d", output_es = "d", n1 = 64, n2 = 36)

# output_es = "A" ---------------------------------------------------------

actual__A_A__convert_es <- convert_es(es = .8, input_es = "A", output_es = "A", n1 = 64, n2 = 36)


#save.image(file = "test-convert_es-actual.RData")
