################################################################################
# Given the output of an eQTL mapping study, including a matrix of LOD scores
# for all genes at all markers, and a marker map, harvest the local and distal
# peaks above the user-specified threshold. Return a data.frame with the results.
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2025-07-21
################################################################################

require(qtl2)

# Arguments:
# lod: numeric matrix containing lod scores from gene eQTL mapping with 
#      qtl2:scan1(). Markers in rows and genes in columns. Must have rownames
#      and colnames. 
# map: qtl2-style marker map. List containing marker positions for each chromsome. 
#      Length equals number of chromosomes and must be named with chromosome IDs. 
#      Each element is a named numeric vector containing the marker position in Mb.
# annot: data.frame containing gene annotation, including gene ID, chromosome and
#        position. Must contain columns named gene_id, symbol, chr, and start. 
#        The genes in lod must appear in annot.
# thr: floating point number that is the LOD threshold above which QTL will be
#      harvested.
# local_radius: floating point number that is the number of Mb in each direction
#               around a gene in which to look for a local eQTL. 
#               e.g. a local radius of 5 Mb means a 10 Mb window.
#
# Returns: data.frame containing the following columns:
#          gene_id: gene identifier provided in lod and annot,
#          lodindex: LOD row in the lod object,
#          chr_qtl: chromosome where QTL is located,
#          pos_qtl: position of QTL in Mb,
#          lod: LOD of the QTL peak,
#          ci_lo: lower bound of QTL interval in Mb,
#          ci_hi: upper bound of QTL interval in Mb,
#          symbol: gene symbol from annot,
#          chr_gene: chromosome where gene is located,
#          pos_qtl: position of gene in Mb,
#          local_distal: either "local" or "distal", depending on the positions
#                        of the gene and QTL peaks.
harvest_eqtl = function(lod, map, annot, thr, local_radius) {

  # Verify that we have the required arguments.
  if(is.null(lod)) {
    stop('lod argument is missing. usage: harvest_eqtl(lod, map, thr)')
  } # if(is.null(lod))

  if(is.null(map)) {
    stop('lod argument is missing. usage: harvest_eqtl(lod, map, thr)')
  } # if(is.null(map))

  if(is.null(thr)) {
    stop('lod argument is missing. usage: harvest_eqtl(lod, map, thr)')
  } # if(is.null(thr))

  # Verify that the marker map and lod objects have the same markers.
  if(nrow(lod) == sum(sapply(map, length))) {

    map_markers = unlist(sapply(map, names))
    lod_markers = rownames(lod)

    if(!all(lod_markers == map_markers)) {
      stop('lod and map objects do not contain the same markers.')
    } # if(!all(lod_markers == map_markers))

  } else {
    stop('lod and map objects do not contain the same number of markers.')
  } # else
  
  # Change the gene annotation gene positions to Mb, if needed.
  annot = annot[,c('gene_id', 'symbol', 'chr', 'start')]
  if(max(annot$start) > 200) {
  
    annot$start = annot$start * 1e-6
    
  } # if(max(annot$start) > 200)
  
  # Get all of the QTL at or above the threshold, with support intervals.
  peaks = find_peaks(lod, map, threshold = thr, prob = 0.95)
  
  # Change the column name to make the join work.
  colnames(peaks)[colnames(peaks) == 'lodcolumn'] = 'gene_id'
  
  # Join the eQTL and annotation data.
  peaks = merge(peaks, annot, by = 'gene_id', suffixes = c('_qtl', '_gene'),
                all.x = TRUE, sort = FALSE)

  # Determine whether the eQTL peaks are local or distal.
  peaks$local_distal = ifelse(peaks$chr_gene == peaks$chr_qtl & 
                              abs(peaks$pos - peaks$start) <= local_radius,
                              'local',
                              'distal')
  
  # Fix some column names to clarify gene versus QTL positions.
  colnames(peaks)[colnames(peaks) == 'pos']   = 'pos_qtl'
  colnames(peaks)[colnames(peaks) == 'start'] = 'pos_gene'
  
  return(peaks)

} # harvest_eqtl()

