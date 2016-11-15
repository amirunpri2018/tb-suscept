.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Package contains source code for Blischak et al.")
}

# Utility function to force R to output a specific number of decimal places.
# Avoids Git thinking results have changed simply because of slight changes
# in insignificant digits.
# http://stackoverflow.com/a/12135122
specify_decimal <- function(x, k) format(round(x, k), nsmall = k)
