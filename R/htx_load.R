## freshen_genes = function (endpoint = URL_hsds(), svrtype = "hsds", dsetname = "/counts") 
## {
## #
## # for rhdf5client >= 1.5.3 we have to use HSDSArray instead of older class components
## #
##     ds = HSDSArray(endpoint = endpoint, svrtype = svrtype, domain = "/shared/bioconductor/htxcomp_genes.h5",
##         dsetname = dsetname)
##     ds
## }
## 
## freshen_txlevel = function (endpoint = URL_hsds(), svrtype = "hsds", dsetname = "/counts") 
## {
## #
## # for rhdf5client >= 1.5.3 we have to use HSDSArray instead of older class components
## #
##     ds = HSDSArray(endpoint = endpoint, svrtype = svrtype, domain = "/shared/bioconductor/basecounts.h5",
##         dsetname = dsetname)
##     ds
## }

# this had to change for 3.17 with some changes to S4Vectors
htx_check_cache = function (cache = BiocFileCache::BiocFileCache(), genesOnly=TRUE) 
{
# return 2-vector: c(action, rid if available)
    if (!requireNamespace("BiocFileCache")) stop("install BiocFileCache to use this function")
    if (genesOnly) {
       if (BiocManager::version() < "3.17") qans = BiocFileCache::bfcquery(cache, "rangedHtxGeneSE.rds")
       else qans = BiocFileCache::bfcquery(cache, "rangedHtxGeneSE_317.rds")
       }
    else qans = BiocFileCache::bfcquery(cache, "htxcompSE.rds") # transcript-level
    if (nrow(qans)<1) return(c("install", NA))
    chkupdate = BiocFileCache::bfcneedsupdate(cache, qans$rid[1])
    if (chkupdate) return(c("update", rid=qans$rid[1]))
    c("ok", qans$rid[1])
}


# WOULD importFrom BiocFileCache bfcinfo BiocFileCache bfcrpath
#' load a SummarizedExperiment shell for the Human Transcriptome Compendium
#' @importFrom S4Vectors mcols
#' @importFrom SummarizedExperiment rowRanges rowData<-
#' @param remotePath defaults to NULL, will be removed in future version
#' @param cache a BiocFileCache instance, defaulting to value of BiocFileCache()
#' @param genesOnly logical(1) if TRUE return reference to 
#' SummarizedExperiment with gene-level quantifications; in this 
#' case the remotePath value is
#' set to `https://s3.amazonaws.com/bcfound-bigrna/rangedHtxGeneSE.rds`.
#' As of Bioc 3.17, changes to S4Vectors necessitate different remote images
#' or redesign of handling HSDS instances.  Thus remotePath parameter is dropped
#' and version-specific code is added to the loader to select an appropriate path.
#' @return a RangedSummarizedExperiment instance
#' @examples
#' if (interactive()) {
#'   h = htx_load()
#'   h
#'   SummarizedExperiment::assay(h)
#'   }
#' @export
htx_load = function (remotePath = NULL,
    cache = BiocFileCache::BiocFileCache(), genesOnly=TRUE)
{
    if (!genesOnly) stop("transcript-level quantifications not available currently")
    if (!requireNamespace("BiocFileCache")) stop("install BiocFileCache to use this function")
    remotePath = "https://bcfound-bigrna.s3.amazonaws.com/rangedHtxGeneSE.rds"
    if (BiocManager::version() >= "3.17") remotePath = "https://bcfound-bigrna.s3.amazonaws.com/rangedHtxGeneSE_317.rds"
    chkans = htx_check_cache(cache)
    if (chkans[1] == "install") {
        message("adding RDS to local cache, future invocations will use local image")
        }
    else if (chkans[1] == "update") {
        message("updating local cache")
        tmp = BiocFileCache::bfcupdate(cache, rid=chkans[2], fpath=remotePath, download=TRUE)
        }
    else if (chkans[1] != "ok") stop("unintended response from cache check")
    path = BiocFileCache::bfcrpath(cache, remotePath)
    ans = readRDS(path)
#    if (BiocManager::version() >= "3.17") { # new reference needed
#           harr = rhdf5client::HSDSArray(rhdf5client::URL_hsds(), 
#              "hsds", "/shared/bioconductor/htxcomp_genes.h5", "/counts")
#           requireNamespace("S4Vectors")
#           sl = new("SimpleAssays", data=SimpleList(counts_lstpm=harr)) 
#           slot(ans, "assays") = sl
#           }
     ans
}


#' add gene-level rowData derived from transcript level rowRanges
#' @param x result of htx_load()
#' @return RangedSummarizedExperiment with enhanced rowData
#' @examples
#' # this function operates on a SummarizedExperiment that has
#' # transcript-level rowRanges but gene-level quantifications
#' addRD
#' @export
addRD = function(x) {
 txl = unlist(rowRanges(x), use.names=FALSE)
 drp = which(duplicated(txl$gene_id))
 txl = txl[-drp]
 rowData(x) = mcols(txl)[, c("gene_type", "gene_id", "gene_name", "havana_gene")]
 x
}
 
