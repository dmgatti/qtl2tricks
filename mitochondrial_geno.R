################################################################################
# Given the CC/DO founder and sample genotypes, assign each sample to one of the
# founder mitochondrial genotypes.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-02-23
################################################################################

# Assign mitochondrial genotypes.
# Genotypes can be 1-letter or 2-letter.
# Arguments:
# founders: character matrix containing genotypes. Markers in rows, samples
#           in columns, named as A:H. First column should be 'marker'.
# samples:  character matrix containing genotypes. Markers in rows, samples
#           in columns. First column should be 'marker'.
# We expect five groups:
# ABCD
# E
# F
# G
# H
# Returns: data.frame containing three columns.
#   id: the sample ID provided in the samples argument.
#   mito: character string indicating the most probable mitochondrial 
#         genotype.
#   cor: the correlation of the most probable mitochondrial genotype with
#        the corresponding founder.
mitochondrial_geno = function(founders, sample) {

  # Verify that the data contains the same markers. We'll trust that
  # they are mitochondrial, but perhaps adding markers to the arguments 
  # would be a good idea.
  stopifnot(nrow(founders) == nrow(samples))
  stopifnot(founders$marker == samples$marker)

  # Make substitution vector for mitochondrial groups.
  mito_groups = setNames(rep(c('ABCD', 'E', 'F', 'G', 'H'), c(4, 1, 1, 1, 1)),
                         LETTERS[1:8])

  # Convert founder data into a long numeric matrix with 8 columns,
  # one for each founder.
  f = founders
  rownames(f) = f$marker
  f = data.frame(t(f[,-1]))
  dn = dimnames(f)
  f = lapply(f, factor)
  f = matrix(as.numeric(unlist(f)), ncol = length(f[[1]]), byrow = TRUE, 
             dimnames = list(NULL, dn[[1]]))

  s = samples
  rownames(s) = s$marker
  s = data.frame(t(s[,-1]))
  dn = dimnames(s)
  s = lapply(s, factor)
  s = matrix(as.numeric(unlist(s)), ncol = length(s[[1]]), byrow = TRUE,
             dimnames = list(NULL, dn[[1]]))

  fg_cor = cor(f, s, use = 'pairwise')
  max_cor = apply(fg_cor, 2, max)
  max_grp = apply(fg_cor, 2, which.max)
  max_grp = setNames(colnames(f)[max_grp], names(max_grp))
  max_grp = setNames(mito_groups[max_grp], names(max_grp))

  stopifnot(names(max_cor) == names(max_grp))

  return(data.frame(id   = names(max_cor),
                    mito = max_grp,
                    cor  = max_cor))

} # mitochondrial_geno()

