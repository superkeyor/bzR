###**************************************************.
###*wrapper of bioMart
###**************************************************.

#' Lists available archived versions of Ensembl (hosts), marts, snp attributes/filters, gene attributes/filters
#' @param what one of 'hosts', 'marts', 'snp', 'gene'
#' @param host default 'www.ensembl.org'. Other eg, 'grch37.ensembl.org', 'May2017.archive.ensembl.org'. See all, run \code{\link{mart.list}}
#' @export
mart.list <- function(what='snp',host=NULL,...) {
    if (is.null(host)) {host='www.ensembl.org'}

    if (what=='hosts') {
        # modified codes from listEnsemblArchives from bioMart v2.34 (not available in v2.30)
        html <- XML::htmlParse("http://www.ensembl.org/info/website/archives/index.html")
        archive_box <- XML::getNodeSet(html, path = "//div[@class='plain-box float-right archive-box']")[[1]]
        archive_box_string <- XML::toString.XMLNode(archive_box)
        archives <- strsplit(archive_box_string, split = "<li>")[[1]][-1]
        extracted <- stringr::str_extract_all(string = archives, 
                                   pattern = "Ensembl [A-Za-z0-9 ]{2,6}|http://.*ensembl\\.org|[A-Z][a-z]{2} [0-9]{4}")
        hosts <- do.call("rbind", extracted)
        colnames(hosts) <- c("url", "version", "date")
        hosts <- hosts[,c(2,3,1)]
        View(hosts)
    }

    if (what=='marts') {
        marts = biomaRt::listEnsembl(host=host, ...)
        View(marts)
    }

    if (what=='snp') {
        SNP_Attributes = biomaRt::listAttributes(biomaRt::useEnsembl(biomart="snp", dataset="hsapiens_snp", host=host, ...))
        SNP_Filters = biomaRt::listFilters(biomaRt::useEnsembl(biomart="snp", dataset="hsapiens_snp", host=host, ...))
        View(SNP_Attributes)
        View(SNP_Filters)
    }

    if (what=='gene') {
        Gene_Attributes = biomaRt::listAttributes(biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host=host, ...))
        Gene_Filters = biomaRt::listFilters(biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host=host, ...))
        View(Gene_Attributes)
        View(Gene_Filters)
    }

}

#' return a mart object representing snp:hsapiens_snp
#' @param host default 'www.ensembl.org'. Other eg, 'grch37.ensembl.org', 'May2017.archive.ensembl.org'. See all, run \code{\link{mart.list}}
#' @export
mart.snp = function(host=NULL, biomart="snp", dataset="hsapiens_snp", ...) {
    if (is.null(host)) {host='www.ensembl.org'}
    return(biomaRt::useEnsembl(biomart=biomart, dataset=dataset, host=host, ...))
}

#' retrieve human snp info from ensembl
#' @param values  the actual input data values, rs number, eg, 'rs2075507', c('rs2075507', 'rs547420070', 'rs77274555'); to search archived synonymous rs number, use mart.snpinfo2
#' @param filters the kind/type of your input data
#' @param attributes what to return, eg, 'synonym_name'
#' @param host default 'www.ensembl.org'. Other eg, 'grch37.ensembl.org', 'May2017.archive.ensembl.org'. See all, run \code{\link{mart.list}}
#' @return returns a data frame
#' @note
#' In case multple filters are in use, the values argument requires a list of values where each position in the list corresponds to the position of the filters in the filters argument. eg,
#' \cr filters=c("chr_name","start","end")
#' \cr values=list(8,148350, 158612)
#' @export
mart.snpinfo = function(values=c('rs2075507', 'rs547420070', 'rs77274555'),filters='snp_filter',attributes=c('refsnp_id','chr_name','chrom_start','chrom_end','allele',
        'allele_1','minor_allele','minor_allele_freq','ensembl_gene_stable_id'),host=NULL) {
    if (is.null(host)) {host='www.ensembl.org'}
    rs <- biomaRt::getBM(attributes=attributes, filters=filters, values=values, mart=mart.snp(host=host))
    
    # retrieve gene names with ensembl_gene_stable_id
    if ('ensembl_gene_stable_id' %in% colnames(rs)) {
        # thanks to getBM(uniqueRows = TRUE),
        # unfounded id or empty id '' in rs[['ensembl_gene_stable_id']] will be auto filtered out in rs2
        # but if rs[['ensembl_gene_stable_id']] is all empty '', rs2 fails
        # https://stackoverflow.com/a/46895403/2292993
        allempty=function(x){all(is.na(x) || is.null(x) || x == "" || x == 0)}
        values=rs[['ensembl_gene_stable_id']]
        attributes=c("ensembl_gene_id","hgnc_symbol","description")
        if (!allempty(values)) {
            rs2=biomaRt::getBM(filters='ensembl_gene_id', attributes=attributes, values=values, mart=mart.gene(host=host))
            # combine results
            rs=dplyr::left_join(rs,rs2,by=c('ensembl_gene_stable_id'='ensembl_gene_id'))
        } else {
            # simply create new columns
            rs[attributes]=''
        }
    }

    return(rs)
}

#' @rdname mart.snpinfo
#' @export
mart.snpinfo2 = function(values=c('rs2075507', 'rs547420070', 'rs77274555'),filters='snp_synonym_filter',attributes=c('refsnp_id','chr_name','chrom_start','chrom_end','allele',
        'allele_1','minor_allele','minor_allele_freq','ensembl_gene_stable_id'),host=NULL) {
    if (is.null(host)) {host='www.ensembl.org'}
    rs <- biomaRt::getBM(attributes=attributes, filters=filters, values=values, mart=mart.snp(host=host))
    
    # retrieve gene names with ensembl_gene_stable_id
    if ('ensembl_gene_stable_id' %in% colnames(rs)) {
        # thanks to getBM(uniqueRows = TRUE),
        # unfounded id or empty id '' in rs[['ensembl_gene_stable_id']] will be auto filtered out in rs2
        # but if rs[['ensembl_gene_stable_id']] is all empty '', rs2 fails
        # https://stackoverflow.com/a/46895403/2292993
        allempty=function(x){all(is.na(x) || is.null(x) || x == "" || x == 0)}
        values=rs[['ensembl_gene_stable_id']]
        attributes=c("ensembl_gene_id","hgnc_symbol","description")
        if (!allempty(values)) {
            rs2=biomaRt::getBM(filters='ensembl_gene_id', attributes=attributes, values=values, mart=mart.gene(host=host))
            # combine results
            rs=dplyr::left_join(rs,rs2,by=c('ensembl_gene_stable_id'='ensembl_gene_id'))
        } else {
            # simply create new columns
            rs[attributes]=''
        }
    }

    return(rs)
}

#' return a mart object representing ensembl:hsapiens_gene_ensembl
#' @param host default 'www.ensembl.org'. Other eg, 'grch37.ensembl.org', 'May2017.archive.ensembl.org'. See all, run \code{\link{mart.list}}
#' @export
mart.gene = function(host=NULL, biomart="ensembl", dataset="hsapiens_gene_ensembl", ...) {
    if (is.null(host)) {host='www.ensembl.org'}
    return(biomaRt::useEnsembl(biomart=biomart, dataset=dataset, host=host, ...))
}

#' retrieve human gene info from ensembl
#' @param values  the actual input data values, ensembl_gene_id, ENSG00000118473, c('ENSG00000118473', 'ENSG00000162426')
#' @param filters the kind/type of your input data
#' @param attributes what to return, max eternal=3(?), otherwise, Too many attributes selected for External References
#' \cr "hgnc_id"(HGNC:25412),"entrezgene" (84251),"kegg_enzyme","go_id"(GO:0030122),"ucsc"(uc057hhx.1) 
#' @param host default 'www.ensembl.org'. Other eg, 'grch37.ensembl.org', 'May2017.archive.ensembl.org'. See all, run \code{\link{mart.list}}
#' @return returns a data frame
#' @note
#' In case multple filters are in use, the values argument requires a list of values where each position in the list corresponds to the position of the filters in the filters argument. eg,
#' \cr filters=c("chr_name","start","end")
#' \cr values=list(8,148350, 158612)
#' @export
mart.geneinfo = function(values=c('ENSG00000118473', 'ENSG00000162426'),filters='ensembl_gene_id',attributes=c("entrezgene","ucsc","hgnc_id","ensembl_gene_id","hgnc_symbol","description","chromosome_name","start_position","end_position"),host=NULL) {
  if (is.null(host)) {host='www.ensembl.org'}
  rs <- biomaRt::getBM(attributes=attributes, filters=filters, values=values, mart=mart.gene(host=host))
  return(rs)
}
