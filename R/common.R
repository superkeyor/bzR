###**************************************************.
###*functions shared by different genetic data analysis modalities
###**************************************************.

#' Install a package from bioconductor
#' @description internally source/install biocLite.R, wrapper of biocLite(pkg)
#' @param pkg eg, 'biomaRt'
#' @return returns nothing
#' @export
bz.install = function(pkg) {
    source("https://bioconductor.org/biocLite.R")
    biocLite(pkg)
}

