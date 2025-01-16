################################################################################
# This script contains code to estimate haplotype blocks in the DO using the
# 36-state genoprobs. It does this by crawling along the probabilities,
# looking for markers which are highly correlated with the current marker.
# The first marker to show low correlation with the current marker is taken
# to be the end of one haplotype block.
# The first function estimates haplotype blocks for ONE SAMPLE on ONE
# CHROMOSOME.
# The second function accepts a qtl2-style genoprobs object and at qtl2-style
# marker map and estimates the haplotype blocks for all samples on all 
# chromosomes.
# The functions are currently SLOW.
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2025-01-16

##########
# Estimate the haplotype blocks for ONE SAMPLE on ONE CHROMOSOME.
# This function crawls along the probability correlation matrix and
# searches for haplotype blocks by finding markers which are highly
# correlated with the first marker in the block. 
# Arguments:
# probs: numeric matrix containing the 36-state genoprobs for
#        one sample on one chromosome. Genotype states in rows
#        and markers in columns. Dimensions must have dimnames.
# mkrs: numeric vector containing the marker positions for the
#       markers in probs. Must be named. Markers must already 
#       be in the same order as probs.
get_blocks_1sample = function(probs, mkrs) {


  states  = rownames(probs)
  markers = colnames(probs)

  blocks = data.frame(prox_marker = rep('', 1000),
                      dist_marker = rep('', 1000),
                      gt          = rep('', 1000),
                      mean_prob   = rep(0, 1000),
                      prox_idx    = rep(0, 1000),
                      dist_idx    = rep(0, 1000))
  
  block_idx = 1
  prox_idx  = 1
  dist_idx  = 1
  
  # While the proximal index is less than the number of markers.
  while(prox_idx < length(markers)) {
  
    # Set the search range from the current marker to the end of the
    # chromosome.
    search_range = prox_idx:length(markers)
    
    # Get the correlation of the current marker with all distal markers.
    probs_cor = cor(probs[,prox_idx], probs[,search_range])[1,]
    
    # Find the distal end of the haplotype block by searching for the
    # first marker that has low correlation with the proximal marker.
    wh = which(probs_cor < 0.5)
    
    # This may happen at the end of the chromosome.
    if(length(wh) == 0) {
    
      dist_idx = length(markers)

    } else {
    
      dist_idx = min(wh) + prox_idx - 2

    } # else

    # If we have a one marker block, skip over it.  
    if(prox_idx != dist_idx) {
  
      # Set the proximal and distal marker names.
      blocks$prox_marker[block_idx] = markers[prox_idx]
      blocks$dist_marker[block_idx] = markers[dist_idx]
    
      # Get the mean probability for each genotype state.
      rm_probs = rowMeans(probs[,prox_idx:dist_idx])
      
      # Set the called genotype to the state with the highest probability.
      blocks$gt[block_idx]        = names(which.max(rm_probs))
      
      # Get the mean probability for the maximal state.
      blocks$mean_prob[block_idx] = max(rm_probs)
    
      # Set the proximal and distal maker indices.
      blocks$prox_idx[block_idx] = prox_idx
      blocks$dist_idx[block_idx] = dist_idx

      # Increment block index.
      block_idx = block_idx + 1
    
    } # if(prox_idx != dist_idx)
    
    # Increment proximal index for the next block.
    prox_idx  = dist_idx  + 1

  } # while(prox_idx < ncol(probs))

  # Trim down the blocks since we started with 1000.
  blocks = subset(blocks, prox_marker != '')
  
  return(blocks)

} # get_blocks_1sample()


##########
# This function estimates the haplotype blocks in the DO for all samples
# on all chromosomes.
# Arguments:
# probs: list which is a qtl2-style genoprobs object. Each list element
#        is a 3-dimensional array with samples in rows, 36 genotype states
#        in columns, and markers in slices. All dimensions must have dimnames.
#        The list element names are chromosomes. Marker names must match 
#        those in "map".
# map: list which is a qtl2-style marker map. Each list element is a numeric
#      vector containing marker positions in Mb. The vector names must be
#      marker names. The list element names are chromosomes. Marker names
#      must match those in "probs".
# NOTE: I tried to parallelize this with foreach and doParallel, but the
# workers couldn't find the get_blocks_1sample() function even though it
# was in the environment when I called %dopar%.
get_hap_blocks = function(probs, map) {

  # Verify that probs and map have the same length.
  if(length(probs) != length(map)) {
    stop(paste0('get_hap_blocks: probs (', length(probs), ') and map (',
         length(map), ') must have the same length.'))
  } # if(length(probs) != length(map)

  # Verify that probs and map have the same chromsome names.
  if(!all(names(probs) == names(map))) {
    stop(paste0('get_hap_blocks: probs and map must have the same chromosome ',
         'names in the same order.'))
  } # if(!all(names(probs) == names(map)))

  # Verify that we have a 36-state genoprobs object.
  if(ncol(probs[[1]]) != 36) {
    stop('get_hap_blocks: This functtion only works with the 36-state genoprobs.')
  } # if(ncol(probs[[1]]) != 36)

  # Verify that the markers are the same between probs and map.
  for(chr in seq_along(probs)) {

    if(!all(names(map[[chr]]) == dimnames(probs[[chr]])[[3]])) {
      stop(paste0('get_hap_blocks: Markers do not match between probs and map',
                  ' on chromosome', names(map)[chr]))
    } # if(!all(names(map[[chr]]) == dimnames(probs[[chr]])[[3]]))

  } # for(chr)

  # Estimate haplotype blocks for all samples on all chromosomes.
  hap_blocks = NULL

  for(chr in seq_along(probs)) {

    t1 = proc.time()
    message(paste('CHR:', chr))

    for(sample in 1:nrow(probs[[chr]])) {

      sample_blocks = get_blocks_1sample(probs[[chr]][sample,,], map[[chr]])
      sample_blocks = cbind(id = rownames(probs[[chr]][sample]), chr = chr, 
                            sample_blocks)

      hap_blocks = rbind(hap_blocks, sample_blocks)

    } # for(sample)

    message(proc.time()[3] - t1[3])

  } # for(chr)

  return(hap_blocks)

} # get_hap_blocks()


