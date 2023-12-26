################################################################################
# Plot a karyogram of one mouse using the allele probabilities.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-05-11
################################################################################

require(qtl2)

# Given a set of allele probs, plot the karyogram using the allele probabilities
# rather than converting to two alleles at each marker. This gives us some
# idea of model uncertainty in the plots. 
# Arguments:
# probs: qtl2-style allele probs object, created from qtl2::genoprobs_to_alleleprobs().
#        Named list containing allele probs for each chromosome. Samples in rows,
#        founders in columns, markers in slices (dim[3]).
# map:  qtl2-style marker map. Named numeric vector. Marker names in names.
#       Positions in vector.
# ind: integer indicating the individual in probs to plot.
# col: Named color vector containing colors to use when plotting.
plot_onegeno_probs = function(probs, map, ind = 1, col = qtl2::CCcolors) {
  
  par(plt = c(0.08, 0.95, 0.08, 0.95), mgp = c(1, 1, 0))
  
  plot(1, 1, col = 'white', xlim = c(0,21), ylim = c(200, 0), las = 1,
       xaxt = 'n', ann = F)
  axis(side = 1, at = 1:20, labels = names(probs))
  mtext(text = rownames(probs[[1]])[ind], side = 3, line = 0.5, cex = 2)
  abline(h = c(0, 50, 100, 150, 200), col = 'grey80')
  abline(v = 1:20, col = 'grey80')
  
  pr = probs[ind,]
  
  for(chr in names(pr)) {
    
    common_markers = intersect(names(map[[chr]]), dimnames(pr[[chr]])[[3]])
    map[[chr]]     = map[[chr]][common_markers]
    pr[[chr]]      = pr[[chr]][,,common_markers]
    pr[[chr]]      = apply(pr[[chr]], 2, cumsum)
    pr[[chr]]      = cbind(0, t(pr[[chr]]))
    
  } # for(chr)
  
  for(chr in seq_along(pr)) {
    
    chr_map = c(map[[chr]], rev(map[[chr]]))

    for(i in 1:8) {
      
      polygon(x = chr - 0.35 + c(pr[[chr]][,i], rev(pr[[chr]][,i + 1])) * 0.7, 
              y = chr_map, density = NULL, col = col[i], border = NA)
      
    } # for(i)
    
  } # for(pr)

  legend(x = 15, y = 140, legend = names(qtl2::CCcolors), pch = 15,
         col = qtl2::CCcolors, bg = 'white')
  
} # plot_onegeno_probs()


# Test code.
# markers = read.csv('C:/Users/c-dgatti/Documents/muga/gm_uwisc_v4.csv')
# markers$pos = markers$bp_grcm39 * 1e-6
# map = qtl2convert::map_df_to_list(markers, pos_column = 'pos')
# 
# for(s in 1:nrow(probs[[1]])) {
#   png(file.path('C:/Users/c-dgatti/Documents/projects/low_pass_DNAseq/figures', paste0(rownames(probs[[1]])[s], '_karyogram.png')), 
#       width = 1000, height = 1000, res = 128)
#   plot_onegeno_probs(probs, map, ind = s, col = qtl2::CCcolors)
#   legend(15, 137, legend = names(qtl2::CCcolors), fill = qtl2::CCcolors)
#   dev.off()
# }
