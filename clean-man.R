## This script removes extraneous files after roxygen generates documentation.
## Source this script immediately after running devtools::document("psychmeta").

## Note: The base directory of this script is the folder containing the psychmeta cloned repo, for consistency with how the
## psychmeta working directory is defined during package building.

system("rm psychmeta/man/dot-*.*")
