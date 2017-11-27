###**************************************************.
###*wrapper of bioMart
###**************************************************.
# overview: 
# snpinfo() accepts rsxxxx and returns relevant info
# geneinfo(), geneconv(), genesnps() accepts gene id, gene name, and returns relevant info

#' Lists available archived versions of Ensembl (hosts), marts, snp attributes/filters, gene attributes/filters
#' @param what one of 'hosts', 'marts', 'snp', 'gene'
#' @param host default 'www.ensembl.org'. Other eg, 'grch37.ensembl.org', 'May2017.archive.ensembl.org'. See all, run \code{\link{mart.list}}
#' @export
mart.list <- function(what='snp',host=NULL,...) {
    if (!(what %in% c('hosts', 'marts', 'snp', 'gene'))) {
        ez.print("double check 'what', should be one of 'hosts', 'marts', 'snp', 'gene'")
        return(invisible(NULL))
    }

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
#' @param values  the actual input data values, cannot be all empty, but partial empty/NA fine, rs number, eg, 'rs2075507', c('rs2075507', 'rs547420070', 'rs77274555'); to search archived synonymous rs number, use mart.snpinfo2
#' @param filters row filters in the db, the kind/type of your input data
#' @param attributes column attributes, what to return, eg, 'synonym_name'
#' @param host default 'www.ensembl.org'. Other eg, 'grch37.ensembl.org', 'May2017.archive.ensembl.org'. See all, run \code{\link{mart.list}}
#' @return returns a data frame
#' @note
#' In case multple filters are in use, the values argument requires a list of values where each position in the list corresponds to the position of the filters in the filters argument. eg,
#' \cr filters=c("chr_name","start","end")
#' \cr values=list(8,148350, 158612)
#' \cr 
#' \cr for multiple regions, use: 
#' \cr filters="chromosomal_region" 
#' \cr values=c('17:43750137:43750137','17:43750172:43750172')
#' \cr
#' \cr values=list(chr_name=17,start=c(43750137,43750172),end=c(43750137,43750172)) 
#' \cr should retrieve 43750137:43750137, 43750172:43750172
#' \cr but it seems to run for 43750137:43750172, which is not desired
#' @export
mart.snpinfo = function(values=c('rs2075507', 'rs547420070', 'rs77274555'),filters='snp_filter',attributes=c('refsnp_id','validated','allele','allele_1','minor_allele','minor_allele_freq','clinical_significance','chr_name','chrom_start','chrom_end','ensembl_gene_stable_id'),host=NULL) {
    rs <- biomaRt::getBM(attributes=attributes, filters=filters, values=values, mart=mart.snp(host=host))
    
    # retrieve gene names with ensembl_gene_stable_id
    if ('ensembl_gene_stable_id' %in% colnames(rs)) {
        # getBM() values cannot be all empty (but partial empty/NA is fine, e.g., c('rs123','','rs456', NA), 
        # but it can give out all empty results
        # (or nothing—but all desired column names returned as a data frame, i.e. empty data frame)
        # also by default, it removes duplicated rows when returning
        # biomaRt doesn't return anything if it cannot find it. 
        # It is designed to work just like the BioMart web services at www.biomart.org, which behave the same.
        # https://stackoverflow.com/a/46895403/2292993
        allempty=function(x){all(is.na(x) || is.null(x) || x == "" || x == 0)}
        values=rs[['ensembl_gene_stable_id']]
        attributes=c("ensembl_gene_id","hgnc_symbol","description")
        if (!allempty(values)) {
            rs2=biomaRt::getBM(filters='ensembl_gene_id', attributes=attributes, values=values, mart=mart.gene(host=host))
            # combine results
            # even if rs2 is empty, join works fine
            rs=dplyr::left_join(rs,rs2,by=c('ensembl_gene_stable_id'='ensembl_gene_id'))
            rs=ez.recol(rs,'hgnc_symbol, description first')
        } else {
            # simply create new columns
            rs[attributes]=''
        }
    }

    return(rs)
}

#' @rdname mart.snpinfo
#' @export
mart.snpinfo2 = function(values='rs2097603',filters='snp_synonym_filter',attributes=c('refsnp_id','chr_name','chrom_start','chrom_end','allele',
        'allele_1','minor_allele','minor_allele_freq','ensembl_gene_stable_id'),host=NULL) {
    rs=mart.snpinfo(values=values,filters=filters,attributes=attributes,host=host)
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
#' @param values  the actual input data values, cannot be all empty, but partial empty/NA fine, ensembl_gene_id, ENSG00000118473, c('ENSG00000118473', 'ENSG00000162426'). If vector, should be the same id type
#' @param filters row filters in the db, the kind/type of your input data, 'hgnc_symbol', "ensembl_gene_id"
#' \cr "ensembl_gene_id"(ENSG00000118473),"hgnc_id"(HGNC:25412),"entrezgene" (84251),"kegg_enzyme"(00010+1.1.1.1),"go_id"(GO:0030122),"ucsc"(uc057hhx.1) 
#' @param attributes column attributes, what to return
#' @param host default 'www.ensembl.org'. Other eg, 'grch37.ensembl.org', 'May2017.archive.ensembl.org'. See all, run \code{\link{mart.list}}
#' @return returns a data frame
#' @note
#' In case multple filters are in use, the values argument requires a list of values where each position in the list corresponds to the position of the filters in the filters argument. eg,
#' \cr filters=c("chr_name","start","end")
#' \cr values=list(8,148350, 158612)
#' \cr 
#' \cr for multiple regions, use: 
#' \cr filters="chromosomal_region" 
#' \cr values=c('17:43750137:43750137','17:43750172:43750172')
#' \cr
#' \cr values=list(chr_name=17,start=c(43750137,43750172),end=c(43750137,43750172)) 
#' \cr should retrieve 43750137:43750137, 43750172:43750172
#' \cr but it seems to run for 43750137:43750172, which is not desired
#' @export
mart.geneinfo = function(values=c('COMT', 'CRHR1'),filters="hgnc_symbol",attributes=c("ensembl_gene_id","hgnc_symbol","description","chromosome_name","start_position","end_position"),host=NULL) {    
    rs <- biomaRt::getBM(attributes=attributes, filters=filters, values=values, mart=mart.gene(host=host))
    return(rs)
}

#' convert from gene name to id, id to name, id to id
#' @param values  the actual input data values, cannot be all empty, but partial empty/NA fine, ensembl_gene_id, ENSG00000118473, c('ENSG00000118473', 'ENSG00000162426'). If vector, should be the same id type
#' @param filters row filters in the db, from what id, 'hgnc_symbol'
#' \cr "ensembl_gene_id"(ENSG00000118473),"hgnc_id"(HGNC:25412),"entrezgene" (84251),"kegg_enzyme"(00010+1.1.1.1),"go_id"(GO:0030122),"ucsc"(uc057hhx.1) 
#' @param attributes column attributes, to what id(s), external max 3
#' @param host default 'www.ensembl.org'. Other eg, 'grch37.ensembl.org', 'May2017.archive.ensembl.org'. See all, run \code{\link{mart.list}}
#' @return returns a data frame
#' @note
#' To allow BioMart to return results in a reasonable amound of time, 
#' Ensembl has restricted the number of External references that you can select to 3. 
#' ensembl_gene_id-->ensembl_gene_id still works, but if the id does not exist, returns empty data frame
#' @export
mart.geneconv = function(values=c('ENSG00000118473', 'ENSG00000162426'),filters="ensembl_gene_id",attributes=c("ensembl_gene_id","hgnc_symbol","hgnc_id","entrezgene"),host=NULL) {    
    rs <- biomaRt::getBM(attributes=attributes, filters=filters, values=values, mart=mart.gene(host=host))
    return(rs)
}

#' retrieves a gene's all snps
#' @param values  the actual input data values, cannot be all empty, but partial empty/NA fine, ensembl_gene_id, ENSG00000118473, c('ENSG00000118473', 'ENSG00000162426'). If vector, should be the same id type
#' @param filters row filters in the db, input data type, 'hgnc_symbol'
#' \cr "ensembl_gene_id"(ENSG00000118473),"hgnc_id"(HGNC:25412),"entrezgene" (84251),"kegg_enzyme"(00010+1.1.1.1),"go_id"(GO:0030122),"ucsc"(uc057hhx.1) 
#' @param attributes column attributes, what to return
#' @param host default 'www.ensembl.org'. Other eg, 'grch37.ensembl.org', 'May2017.archive.ensembl.org'. See all, run \code{\link{mart.list}}
#' @return returns a data frame
#' @export
mart.genesnps = function(values='CRHR1',filters='hgnc_symbol',attributes=c("variation_name","validated","allele","minor_allele","minor_allele_freq","transcript_count","clinical_significance","ensembl_gene_id","chromosome_name","start_position","end_position","band","external_gene_name","description"),host=NULL) {
    rs <- biomaRt::getBM(attributes=attributes, filters='hgnc_symbol', values=values, mart=mart.gene())
    return(rs)
}
