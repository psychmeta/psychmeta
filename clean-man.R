## This script removes extraneous files after roxygen generates documentation and also cleans up Unicode usage in .Rd files.
## Source this script immediately after running devtools::document("psychmeta").

## Note: The base directory of this script is the folder containing the psychmeta cloned repo, for consistency with how the
## psychmeta working directory is defined during package building.

system("rm psychmeta/man/dot-*.*")

system("sed -i -e 's/lgl.mark = c('+', '-')/lgl.mark = c('+', '\\\\u2212')/g' psychmeta/man/format_num.Rd")
system("sed -i -e 's/inf.mark = c('+8', '-8')/inf.mark = c('+\\\\u221e', '\\\\u2212\\\\u221e')/g' psychmeta/man/format_num.Rd")
system("sed -i -e 's/neg.sign = '-'/neg.sign = '\\\\u2212'/g' psychmeta/man/format_num.Rd")
system("sed -i -e 's/<U+202F>/\\\\u202F/g' psychmeta/man/format_num.Rd")
system("sed -i -e 's/∞/\\\\u221e/g' psychmeta/man/format_num.Rd")
system("sed -i -e 's/−/\\\\u2212/g' psychmeta/man/format_num.Rd")
system("sed -i -e 's/ /\\\\u202f/g' psychmeta/man/format_num.Rd")
system("sed -i -e 's/lgl.mark = c('+', '-')/lgl.mark = c('+', '\\\\u2212')/g' psychmeta/man/metabulate.Rd")
system("sed -i -e 's/inf.mark = c('+8', '-8')/inf.mark = c('+\\\\u221e', '\\\\u2212\\\\u221e')/g' psychmeta/man/metabulate.Rd")
system("sed -i -e 's/—/\\\\u2014/g' psychmeta/man/format_num.Rd")

system("sed -i -e 's/neg.sign = '-'/neg.sign = '\\\\u2212'/g' psychmeta/man/metabulate.Rd")
system("sed -i -e 's/<U+202F>/\\\\u202F/g' psychmeta/man/metabulate.Rd")
system("sed -i -e 's/∞/\\\\u221e/g' psychmeta/man/metabulate.Rd")
system("sed -i -e 's/−/\\\\u2212/g' psychmeta/man/metabulate.Rd")
system("sed -i -e 's/ /\\\\u202f/g' psychmeta/man/metabulate.Rd")
system("sed -i -e 's/—/\\\\u2014/g' psychmeta/man/metabulate.Rd")

system("rm psychmeta/man/format_num.Rd-e")
system("rm psychmeta/man/metabulate.Rd-e")