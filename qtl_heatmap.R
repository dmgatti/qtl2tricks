################################################################################
# Given a qtl2-style LOD object and a corresponding marker map, make a heatmap
# of the LOD curves. Hierarchically cluster the phenotypes and truncate the
# maximum LOD to about 8 so that high LOD scores don't dominate the color 
# scale. We also condense the markers to about 1300 pixels.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2024-09-17
################################################################################

##### LIBRARIES #####

require(tidyverse)

##### FUNCTIONS #####

# Arguments:
# scan1_output: qtl-style scan1 LOD scores. Numeric matrix with marker names in
#               rownames and phenotypes in column names.
# map: qtl2-style marker map. Named list in which each element is a named numeric
#      vector with marker names in names. Names of list must be chromosomes.
# max_threshold: float which is the maximum truncation threshold. All LOD scores 
#                above this value will be truncated to this value so that large 
#                LOD scores don't dominate the color scale.
# min_threshold: float which is the minimum truncation threshold. All LOD scores
#                above this value will be truncated to this value to reduce clutter.
# viridis_scale: string containing one of the valid viridis color scales:
#                "magma", "inferno", "plasma", "viridis", "cividis", "rocket", 
#                "mako", "turbo"
# Returns: QTL heatmap as a ggplot tile-plot.
qtl_heatmap = function(scan1_output, map, max_threshold = 8, 
                       min_threshold = 4, viridis_scale = 'magma') {

  # Check for required arguments.
  if(is.null(scan1_output)) {
    stop('scan1_output is a required argument.')
  } # if(is.null(scan1_output))
  
  if(is.null(map)) {
    stop('map is a required argument.')
  } # if(is.null(map))
  
  # Verify that the markers in the lod and map are the same.
  stopifnot(rownames(scan1_output) == unlist(sapply(map, names)))
  
  # Split the LOD by chromosome.
  map_chr = factor(rep(names(map), sapply(map, length)))
  
  # Figure out how many markers to condense on each chromosome.
  # The mouse genome is about 2600 Mb long, so we will have one pixel every
  # 2 Mb.
  chr_len = sapply(X = map, FUN = max)
  num_pxl = floor(chr_len / 2)
  
  # Get the new marker positions.
  new_map = map
  for(i in seq_along(map)) {
    
    new_map[[i]]        = seq(from       = min(map[[i]]), 
                              to         = ceiling(chr_len[i]), 
                              length.out = num_pxl[i])
  
    names(new_map[[i]]) = paste(names(map)[i], 
                                round(new_map[[i]] * 1e6), 
                                sep = '_')
    
  } # for(i)
  
  # Fill in the condensed LOD matrix. We will get the maximum LOD in each
  # interval.
  lod = setNames(vector('list', length(new_map)), names(new_map))
  
  for(i in seq_along(map)) {
    
    # Get a factor with the marker bins that we'll be condensing.
    # I'm subtracting half of the gap between markers and adding one more marker
    # at the end, which is gap / 2 beyond the last marker.
    half_gap = mean(diff(new_map[[i]])) * 0.5
    brks     = cut(map[[i]], c(new_map[[i]] - half_gap, chr_len[i] + half_gap))
    
    # Get the maximum LOD in each marker bin for each phenotype.
    tmp      = data.frame(scan1_output[names(map[[i]]),])
    tmp      = split(x = tmp, f = brks)
    lod[[i]] = suppressWarnings(t(sapply(tmp, function(z) { apply(z, MARGIN = 2, FUN = max) })))
    
    # This happens when a new marker bin is empty, i.e. there are no markers.
    lod[[i]][is.infinite(lod[[i]])] = 0
    
    # Set the marker names on the new LOD.
    rownames(lod[[i]]) = names(new_map[[i]])
    
  } # for(i)
  
  # Combine the new condensed LOD scores.
  lod = do.call(rbind, lod)
  
  # Truncate the LOD scores.
  lod[lod < min_threshold] = min_threshold
  lod[lod > max_threshold] = max_threshold

  # Cluster the LOD scores.
  lod_cor = cor(lod)
  cl      = hclust(as.dist(1.0 - lod_cor), method = 'average')
  lod     = lod[,cl$order]

  # Make the heatmap.
  new_markers = data.frame(mkr = unlist(sapply(new_map, names)),
                           chr = rep(names(new_map), sapply(new_map, length)),
                           pos = unlist(new_map))
  lod = data.frame(lod) |>
          rownames_to_column(var = 'mkr') |>
          left_join(new_markers) |>
          group_by(chr) |>
          mutate(width = mean(diff(pos))) |>
          ungroup() |>
          pivot_longer(cols      = -c(mkr, chr, pos, width), 
                       names_to  = 'Phenotype',
                       values_to = 'LOD') |>
          mutate(chr       = factor(chr, levels = c(1:19, 'X')),
                 Phenotype = factor(Phenotype, levels = cl$labels[cl$order]))
  
  p = lod |>
        ggplot(aes(pos, Phenotype, fill = LOD, width = width)) +
          geom_tile(color = NA) +
          scale_fill_viridis_c(option = viridis_scale, direction = -1) +
          facet_wrap(~chr, nrow = 1, scales = 'free_x') +
          labs(x = 'Position (Mb)', y = '') +
          theme(axis.text.x   = element_text(angle = 90, hjust = 1, vjust = 0.5),
                panel.spacing = unit(0, "lines"))

  return(p)
        
} # qtl_heatmap()
