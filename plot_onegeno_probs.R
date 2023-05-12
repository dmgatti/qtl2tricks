################################################################################
# Plot a karyogram of one mouse using the allele probabilities.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-05-11
################################################################################

require(qtl2)

plot_onegeno_probs = function(probs, map, ind = 1, col = qtl2::CCcolors) {
  
  par(plt = c(0.08, 0.95, 0.08, 0.95), mgp = c(1, 1, 0))
  
  plot(1, 1, col = 'white', xlim = c(0,21), ylim = c(200, 0), las = 1,
       xaxt = 'n')
  axis(side = 1, at = 1:20, labels = names(probs))
  mtext(text = rownames(probs[[1]])[ind], side = 3, line = 0.5, ces = 2)
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
  
} # plot_onegeno_probs()


markers = read.csv('C:/Users/c-dgatti/Documents/muga/gm_uwisc_v4.csv')
markers$pos = markers$bp_grcm39 * 1e-6
map = qtl2convert::map_df_to_list(markers, pos_column = 'pos')
