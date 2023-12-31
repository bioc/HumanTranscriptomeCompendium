% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache_works.R
\name{get_ds4841}
\alias{get_ds4841}
\title{return instance of ssrch::DocSet with metadata on 4841 human
transcriptome studies in NCBI SRA}
\usage{
get_ds4841(
  cache = BiocFileCache::BiocFileCache(),
  csv_zip_path = path_doc4842()
)
}
\arguments{
\item{cache}{instance of `BiocFileCache`, defaults to `BiocFileCache::BiocFileCache()`}

\item{csv_zip_path}{a path leading to the zip file of CSV for metadata in the DocSet instance}
}
\value{
instance of DocSet as defined in ssrch package
}
\description{
return instance of ssrch::DocSet with metadata on 4841 human
transcriptome studies in NCBI SRA
}
\note{
will bind the correct value of `zipf` in `environment(ds4841@doc_retriever)`, which depends on details of installation
}
\examples{
get_ds4841()
}
