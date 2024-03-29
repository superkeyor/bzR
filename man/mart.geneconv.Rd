% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomart.R
\name{mart.geneconv}
\alias{mart.geneconv}
\title{convert from gene name to id, id to name, id to id}
\usage{
mart.geneconv(
  values = c("ENSG00000118473", "ENSG00000162426"),
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol", "hgnc_id", "entrezgene"),
  host = NULL
)
}
\arguments{
\item{values}{the actual input data values, cannot be all empty, but partial empty/NA fine, ensembl_gene_id, ENSG00000118473, c('ENSG00000118473', 'ENSG00000162426'). If vector, should be the same id type}

\item{filters}{row filters in the db, from what id, 'hgnc_symbol'
\cr "ensembl_gene_id"(ENSG00000118473),"hgnc_id"(HGNC:25412),"entrezgene" (84251),"kegg_enzyme"(00010+1.1.1.1),"go_id"(GO:0030122),"ucsc"(uc057hhx.1)}

\item{attributes}{column attributes, to what id(s), external max 3}

\item{host}{default 'www.ensembl.org'. Other eg, 'grch37.ensembl.org', 'aug2017.archive.ensembl.org'. See all, run \code{\link{mart.list}}}
}
\value{
returns a data frame
}
\description{
convert from gene name to id, id to name, id to id
}
\note{
To allow BioMart to return results in a reasonable amound of time, 
Ensembl has restricted the number of External references that you can select to 3. 
ensembl_gene_id-->ensembl_gene_id still works, but if the id does not exist, returns empty data frame
}
