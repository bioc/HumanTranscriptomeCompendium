htx_clean_cache = function( cache = BiocFileCache::BiocFileCache() ) {
  ans = BiocFileCache::bfcquery(cache, "rangedHtxGene")
  if (length(ans$rid)>0) try( BiocFileCache::bfcremove(cache, ans$rid) )
}
 
