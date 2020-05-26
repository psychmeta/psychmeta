#Installs my version of psychmeta
devtools::install_github("wesley4546/psychmeta", ref = "r4")

#retarts R
.rs.restartR()

#Loads it
library(psychmeta)



# Testing Sandbox ---------------------------------------------------------


t <- convert_es(es = 1, input_es = "d", output_es = "r", n1 = 50, n2 = 50)
t


