#!/bin/sh
rm man/dot-*.*
sed -i -e 's/lgl.mark = c("+", "-")/lgl.mark = c("+", "−")/g' man/format_num.Rd
sed -i -e 's/inf.mark = c("+8", "-8")/inf.mark = c("+∞", "−∞")/g' man/format_num.Rd
sed -i -e 's/lgl.mark = c("+", "-")/lgl.mark = c("+", "−")/g' man/metabulate.Rd
sed -i -e 's/inf.mark = c("+8", "-8")/inf.mark = c("+∞", "−∞")/g' man/metabulate.Rd
