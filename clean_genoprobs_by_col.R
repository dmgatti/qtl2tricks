################################################################################
# Given a qtl2-style genoprobs object, 
################################################################################

clean_one_chr = function(pr, thr) {

      # Get the diplotype sums at each marker.
      cs = apply(pr, 2:3, sum)
    
      # Convert the colSums to a logical matrix in which TRUE means that 
      # the colSum is less than the threshold.
      lt_thr = cs < thr

      # Find markers which have at least one column with a sum below
      # the threshold.
      wh = which(colSums(lt_thr) > 0)

      # Loop through the markers that have low diplotype probs,
      # zero out the problem columns, and renormalize the probs at 
      # each marker.
      for(mkr in wh) {

        # Get the columns which need to be set to 0.
        cols2zero = which(lt_thr[,mkr])
      
        # Set the columns to zero.
        pr[,cols2zero,mkr] = 0
      
        # Renormalize to make the sum of each row at this marker equal to 1.
        pr[,,mkr] = pr[,,mkr] / rowSums(pr[,,mkr])

      } # for(mkr)
      
      return(pr)

} # clean_one_chr()



# probs: qtl2-style genoprobs object. Names list containing 3 dimensional
#        numeric arrays wiht samples in rows, diplotypes in columns, and 
#        makers in slices.
# thr: floating point nubmer that is the minimum sum of diplotype probs
#      in a singal column and one marker. If the diplotype sum is less 
#      than this value, all of the diplotype probs for that column will
#      be set to zero and the diplotype probabilities will be renormalized
#      to sum to 1.0.
clean_genoprobs_by_col = function(probs, thr = 1.0) {

  if(!'calc_genoprob' %in% class(probs)) {
  
    stop('probs must be of class calc_genoprob. Please pass in a qtl2-style genoprobs object.')
  
  } # if(!'calc_genoprob' %in% class(probs))

  # Loop through each chromosome.
  for(chr in seq_along(probs)) {

    print(paste('CHR:', names(probs)[chr]))

    if(attr(probs, 'is_x_chr')[chr] == TRUE) {
    
      # Figure out which samples are male and female.
      y_cols = grep('Y$', colnames(probs[[chr]]))
      
      # If there are not columns ending in 'Y', then process all samples together.
      if(length(y_cols) == 0 ) {
      
        probs[[chr]] = clean_one_chr(probs[[chr]], thr)       
      
      } else {

        # Process males and females separately.
        y_sums = rowMeans(probs[[chr]][,y_cols,])
      
        # This needs to be tested more.
        females = which(y_sums == 0)
        
        # If all or none of the samples are female, then process all samples together.
        if(length(females) == 0 | length(females) == nrow(probs[[chr]])) {
        
          probs[[chr]] = clean_one_chr(probs[[chr]], thr)
        
        } else {
        
          # Process each sex separately.
          males = setdiff(1:nrow(probs[[chr]]), females)
          
          probs[[chr]][females,,] =  clean_one_chr(probs[[chr]][females,,], thr)
          probs[[chr]][males,,]   =  clean_one_chr(probs[[chr]][males,,],   thr)
        
        } # else
      
      } # else
      
    } else {

      probs[[chr]] = clean_one_chr(probs[[chr]], thr)

    } # else

  } # for(chr)
  
  return(probs)

} # clean_genoprobs_by_col()


# Test code:
probs_file = '/projects/korstanje-lab/Pureplex/AnalyzedData/TumorStudy_combined/results/quilt/20250421_tumorstudy_combined/2000/geno_probs/complete_genoprobs.rds'
probs = readRDS(probs_file)

new_probs = clean_genoprobs_by_col(probs, 20)

# Check Chr each chromosome. This should be 1, 1.
for(i in seq_along(probs)) {

  print(paste('CHR:', i, ' : ', range(apply(probs[[i]], 3, rowSums))))

} # for(i)

for(i in seq_along(probs)) {

   cs = apply(probs[[i]][,c('AC', 'AD', 'BC', 'BD'),], 2:3, sum)

   print(range(apply(cs, 2, min)))

} # for(i)


