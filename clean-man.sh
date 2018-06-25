#!/bin/sh
rm man/dot-*.*
sed -i -e 's/lgl.mark = c("+", "-")/lgl.mark = c("+", "\\u2212")/g' man/format_num.Rd
sed -i -e 's/inf.mark = c("+8", "-8")/inf.mark = c("+\\u221e", "\\u2212\\u221e")/g' man/format_num.Rd
sed -i -e 's/neg.sign = "-"/neg.sign = "\\u2212"/g' man/format_num.Rd
sed -i -e 's/<U+202F>/\\u202F/g' man/format_num.Rd
sed -i -e 's/∞/\\u221e/g' man/format_num.Rd
sed -i -e 's/−/\\u2212/g' man/format_num.Rd
sed -i -e 's/ /\\u202f/g' man/format_num.Rd
sed -i -e 's/lgl.mark = c("+", "-")/lgl.mark = c("+", "\\u2212")/g' man/metabulate.Rd
sed -i -e 's/inf.mark = c("+8", "-8")/inf.mark = c("+\\u221e", "\\u2212\\u221e")/g' man/metabulate.Rd
sed -i -e 's/neg.sign = "-"/neg.sign = "\\u2212"/g' man/metabulate.Rd
sed -i -e 's/<U+202F>/\\u202F/g' man/metabulate.Rd
sed -i -e 's/∞/\\u221e/g' man/metabulate.Rd
sed -i -e 's/−/\\u2212/g' man/metabulate.Rd
sed -i -e 's/ /\\u202f/g' man/metabulate.Rd
