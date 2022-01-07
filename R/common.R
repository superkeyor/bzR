###**************************************************.
###*functions shared by different genetic data analysis modalities
###**************************************************.

#' Install a package from bioconductor
#' @description internally source/install biocLite.R, wrapper of \code{\link[BiocInstaller]{biocLite}}
#' @param pkg eg, 'biomaRt'
#' @param suppressUpdates if TRUE, do not ask to update other old packages and do not update them; if FALSE, will ask to choose
#' @return returns nothing
#' @export
bz.install = function(pkg,suppressUpdates=TRUE,...) {
    source("https://bioconductor.org/biocLite.R")
    biocLite(pkg,suppressUpdates=suppressUpdates,...)
}

