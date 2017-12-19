#' recode snps
#' @description recode snps
#' @param x df or vector (char or factor of char)
#' @param cols evaluated as cols=ez.selcol(x,cols), ignored if x is a vector
#' @param recodes vector with three elements corresponding to (minor-minor, minor-major or major-minor, major-major). Could be c(0,1,2), c(-1,0,1), c(1,2,3), c(0,1,1), c(1,0,0), c('AA','AB','BB')
#' @return returns a new df or vector. Regardless of input data type, if recodes are number, then the returned is also numeric; if recodes are char, then always factor of char.
#' @note assume biallelic, otherwise error. minor/major bases are calculated based on the actual input data
#' @export
snp.recode = function(x, cols=NULL, recodes=c(0,1,2)) {
    if (!is.data.frame(x)){
        bases = strsplit(paste(as.character(na.omit(x)), collapse=""),"")[[1]]
        if (length(unique(bases))>2) {
          if (is.null(cols)) stop(sprintf('input not biallelic: %s', toString(unique(bases),width = 300)))
          if (!is.null(cols)) stop(sprintf('col %s not biallelic: %s', cols, toString(unique(bases),width = 300)))
        }
        freqs = table(bases)
        minor = names(freqs)[which(freqs==min(freqs))][1]  # [1] in case 50%, 50%
        major = setdiff(names(freqs),minor)
        AA = paste0(minor,minor)
        AB = paste0(minor,major)
        BA = paste0(major,minor)
        BB = paste0(major,major)
        x[which(x==AA)] <- recodes[1] 
        x[which(x==AB)] <- recodes[2]
        x[which(x==BA)] <- recodes[2] 
        x[which(x==BB)] <- recodes[3] 
        x <- utils::type.convert(x, as.is = TRUE)
        if (is.character(x)) x=factor(x)
        result=x
    }
    if (is.data.frame(x) & !is.null(cols)) {
        cols=ez.selcol(x,cols)
        for (col in cols) {
            x[[col]]=snp.recode(x[[col]],cols=col,recodes=recodes)
        }
        result=x
    }
    return(result)
}