################################################################################
# Given a qtl2-style genoprobs or allele probs object with markers at a 
# high-density, smooth the probs down to a lower density marker set.
# An example would be moving from 1 million markers to 200,000 markers.
# Note that the smoothed probs at the end of chromosomes may not include the
# full window size.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2025-04-24
################################################################################

library(IRanges)
library(qtl2)
library(qtl2convert)
library(tidyverse)

# Interpolate the probs object down to a lower marker density, smoothing within
# a given window.
# Arguments:
# probs: qtl2-style probs object (genoprobs or allele probs). Named list in
#        which each element is a 3D array containing the diplotype probs for 
#        one chromosome.
# map1: qtl2-style map object containing the marker positions in the probs
#       object. Named list in which each element is a named numeric vector
#       of marker positions in Mb.
# map2: qtl2-style map object containing the marker positions to output.
#        Named list in which each element is a named numeric vector
#       of marker positions in Mb.
# window: integer that is the width of the smoothing window. Should be an
#         odd number. i.e. 5 means that we will smooth +/- 2 markers around
#         the current marker.
# Returns: qtl2-style probs object at the marker resolution of map2.
smooth_genoprobs = function(probs, map1, map2, window = 5) {

  # Verify that the lengths of all objects is equal.
  if(length(probs) != length(map1)) {

    stop('probs and map1 do not have the same length. They must have the same length.')

  } # if(length(probs) != length(map1))

  if(length(probs) != length(map2)) {

    stop('probs and map1 do not have the same length. They must have the same length.')

  } # if(length(probs) != length(map2))

  # Verify that map1 is larger than map2.
  map1_len = sum(sapply(map1, length))
  map2_len = sum(sapply(map2, length))
  if(map1_len <= map2_len) {

    stop('map1 has fewer markers than map2. map1 must have more markers than map2.')

  } # if(map1_len <= map2_len)

  # Verify that the objects are the correct class.
  if(!'calc_genoprob' %in% class(probs)) {
  
    stop('probs must be of class calc_genoprob')
  
  } # if(!'calc_genoprob' %in% class(probs))

  # Order the objects in the same chromosome order.
  map1 = map1[names(probs)]
  map2 = map2[names(probs)]
  
  # Smooth each chromosome.
  new_probs = setNames(as.list(names(probs)), names(probs))

  for(i in seq_along(probs)) {
  
    print(paste('CHR:', names(probs)[i]))
  
    new_probs[[i]] = smooth_one_chr(probs[[i]], map1[[i]], map2[[i]], window)
  
  } # for(i)
  
  attributes(new_probs) = attributes(probs)
  
  return(probs)

} # smooth_genoprobs()


# Smooth one chromosome.
smooth_one_chr = function(pr, m1, m2, win) {

  # Get the +/- window radius.
  radius = floor(win / 2)

  # Convert the marker maps to IRanges objects..
  m1 = IRanges(start = m1 * 1e6, width = 1, names = names(m1))
  m2 = IRanges(start = m2 * 1e6, width = 1, names = names(m2))
  
  # Get the markers in m1 which are nearest the markers in m2.
  m1_near_m2 = nearest(m2, m1)
  m1_near_m2 = cbind(m1_near_m2 - radius, m1_near_m2 + radius)

  # Make sure that we don't go off the ends of the marker map.
  m1_near_m2[m1_near_m2[,1] < 0, 1]          = 1
  m1_near_m2[m1_near_m2[,2] > length(m1), 2] = length(m1)
  
  # Get the marker ranges in m1 to aggregate.
  m1_ranges = apply(m1_near_m2, 1, function(z) { z[1]:z[2] })
  
  new_pr = array(0, dim = c(nrow(pr), ncol(pr), length(m2)), 
                 dimnames = list(rownames(pr), colnames(pr), names(m2)))
  for(i in 1:dim(new_pr)[3]) {
  
    new_pr[,,i] = apply(pr[,,m1_ranges[[i]], drop = FALSE], 1:2, mean)
  
  } # for(i)
  
  return(new_pr)

} # smooth_one_chr()


############
# Test code.

base_dir = '/projects/korstanje-lab/Pureplex/TumorStudy_combined/results/quilt/20250421_tumorstudy_combined/2000/geno_probs'
probs = readRDS(file.path(base_dir, 'complete_alleleprobs.rds'))
map1  = readRDS(file.path(base_dir, 'complete_pmap.rds'))
map2  = readRDS(file.path(base_dir, 'grid_pmap.rds'))

np = smooth_genoprobs(probs, map1, map2, window = 5)

