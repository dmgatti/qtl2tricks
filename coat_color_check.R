################################################################################
# Test inferred genotype at coat color loci versus reported coat color for 
# Diversity Outbred mice.
# Users must report the coat color as 'agouti', 'albino', or 'black'. 
# Indeterminate coat colors may be coded as NA.
# We only use the allele probs at the agouti locus on Chr 2 and the tyrosinse
# locus on Chrr 7. Coordinates in GRCm39.
# Albino is epistatic (i.e. masks) to black. So test the albino locus
# first and the agouti locus second.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2025-05-23
################################################################################


# Coat color loci.
albino_locus = list(chr = '7', pos = 87.108308)
black_locus  = list(chr = '2', pos = 154.842726)

# Attempt to call the more probable marginal diplotype at a single marker.
call_diplotypes = function(pr) {

  # Round probs to get het and homozygous states.
  pr2 = round(2 * pr)

  # Get the letters of the diplotypes.
  wh = apply(pr2 > 0, 1, which)
  wh = lapply(wh, function(z) { sort(names(z)) })
  wh_len = sapply(wh, length) 

  dt = setNames(rep('', length(wh)), names(wh))

  # Set homozygotes.
  wh1 = which(wh_len == 1)
  dt[wh1] = paste0(wh[wh1], wh[wh1])

  # Set heterozygotes.
  wh2 = which(wh_len == 2)
  dt[wh2] = sapply(wh[wh2], paste0, collapse = '')

  # Rows that do not sum to 2 have some issue with the diplotypes.
  # Sometimes, more than 2 founders have high values, i.e. 0.33, 0.33, 0.33.
  # Or one founder is > 0.5 and the rest of the probs are spread out.
  # In this case, select the two highest founder probs.
  rs  = rowSums(pr2)
  wh0 = which(rs < 2)

  for(i in seq_along(wh0)) {

    tmp        = sort(pr[wh0[i],])[7:8]
    dt[wh0[i]] = paste0(sort(names(tmp)), collapse = '')

  } # for(i)

  dt = data.frame(id = names(dt),dt, round(pr, digits = 3))

  return(dt)

} # call_diplotypes()


# Compare the genotype at the albino and black loci with the reported coat 
# color.
# Arguments:
# pheno: data.frame with mouse IDs in the 'id' column and coat color in the
#        'coat' column. Coat color must be encoded as 'agouti', 'albino', 
#        'black' , or NA.
# probs: qtl2-style allele probs. A named list with one element per chromosome.
#        Each element is a three-dimensional numeric array with samples in 
#        rows, founders in columns, and markers in slices.
# map: qtl2-style marker map. A names list with one element per chromosome.
#      Each element is a named, numeric vector containing the position of
#      each marker in probs in Mb. The chromosomes and marker names must match
#      between probs and map.
# Returns: data.frame containing the mouse id, reported coat color, two-letter
#          diplotype at the marker nearest to the albino locus, inferred coat
#          color based on the diplotype, and a boolean column indicating whether
#          the reported and inferrred coat colors match.
compare_coat = function(pheno, probs, map) {

  # Verify that samples are aligned between pheno and probs.
  stopifnot(pheno$id == rownames(probs[[1]]))

  # Verify that probs and map have the same chromosomes in the same 
  # order.
  stopifnot(names(map) == names(probs))

  # Verify that the marker names are aligned between map and probs
  # on Chr 2 & 7.
  stopifnot(names(map[[black_locus$chr]]) == 
            dimnames(probs[[black_locus$chr]])[[3]])
  stopifnot(names(map[[albino_locus$chr]]) == 
            dimnames(probs[[albino_locus$chr]])[[3]])

  # Get albino diplotypes.
  # Find the nearest marker to the albino locus.
  outer_diff    = abs(albino_locus$pos -map[[albino_locus$chr]])
  albino_marker = names(outer_diff)[which.min(outer_diff)]

  # Get the allele probs at the albino marker.
  albino_probs = probs[[albino_locus$chr]][,,albino_marker]

  # Call diplotypes.
  albino_dt = call_diplotypes(albino_probs)
  colnames(albino_dt)[-1] = paste0('albino_', colnames(albino_dt)[-1])

  # Get black diplotypes.
  # Find the nearest marker to the black locus.
  outer_diff   = abs(black_locus$pos -map[[black_locus$chr]])
  black_marker = names(outer_diff)[which.min(outer_diff)]

  # Get the allele probs at the black marker.
  black_probs = probs[[black_locus$chr]][,,black_marker]

  # Call diplotypes.
  black_dt = call_diplotypes(black_probs)
  colnames(black_dt)[-1] = paste0('black_', colnames(black_dt)[-1])

  # Create return value data.frame.
  retval = data.frame(id   = pheno$id,
                      coat = pheno$coat)
  retval = merge(retval, albino_dt, by = 'id')
  retval = merge(retval, black_dt,  by = 'id')

  # Infer the coat color from the diplotype.
  # A/J (A) and NOD (D) contribute the albino allele.
  retval$geno_coat = ifelse(retval$albino_dt %in% c('AA', 'AD', 'DD'),
                            'albino', 'agouti')
  # A/J (A) and BL6 (B) contribute that black allele. Albino is epistatic
  # to black.
  retval$geno_coat = ifelse(retval$black_dt   %in% c('AA', 'AB', 'BB') &
                            !retval$albino_dt %in% c('AA', 'AD', 'DD'),
                            'black', retval$geno_coat)

  # If the phenotype coat color and the diplotype coat color match,
  # then "match" = TRUE, otherwise FALSE.
  retval$match = retval$coat == retval$geno_coat

  # Reorder the columns and return the results.
  return(retval[,c('id', 'coat', 'albino_dt', 'black_dt', 'geno_coat', 
                   'match', paste0('albino_', LETTERS[1:8]),
                   paste0('black_', LETTERS[1:8]))])

} # compare_coat()






# Prepare test data.
#library(qtl2convert)

#pheno = read.csv('C:/Users/c-dgatti/Documents/data/coat_color_demo.csv',
#                 colClasses = rep('character', 2))
#probs = readRDS('C:/Users/c-dgatti/Documents/TB/haplo_reconstr/tufts_do_alleleprobs_202409.rds')
#markers = read.csv('C:/Users/c-dgatti/Documents/muga/gm_uwisc_v4.csv')
#markers = markers[markers$chr %in% c(1:19, 'X'),]
#markers$pos = markers$bp_grcm39 * 1e-6
#map     = map_df_to_list(markers, pos_column = 'pos')

# Fix samples and markers.
#rn = sapply(strsplit(rownames(probs[[1]]), ';'), '[', 2)

#common_samples = rn[rn %in% pheno$id]
#pheno = pheno[pheno$id %in% common_samples,]
#pheno = pheno[match(common_samples, pheno$id),]

#for(i in seq_along(probs)) {

#  common_markers = intersect(names(map[[i]]), dimnames(probs[[i]])[[3]])

#  rownames(probs[[i]]) = rn
#  probs[[i]] = probs[[i]][common_samples,,common_markers]
#  map[[i]]   = map[[i]][common_markers]
#  stopifnot(names(map[[i]]) == dimnames(probs[[i]])[[3]])

#} # for(i)

#stopifnot(pheno$id == rownames(probs[[1]]))

#coat_comp = compare_coat(pheno, probs, map)

