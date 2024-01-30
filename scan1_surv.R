################################################################################
# Survival QTL maping in the DO using qtl2.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-12-15
################################################################################

# Required libraries.
require(survival)
require(BiocParallel)

# Arguments:
# genoprobs: qtl2-style genoprobs object with founder allele probabilities.
#            List containing 20 elements, each of which is a 3 dimensional
#            numeric array with samples in rows, founders in columns and
#            markers in slices.
# pheno:     data.frame containing at least two columns:
#            'survival': Number of days that mouse lived.
#            'event':    Whether diabetes was observed.
#                        0 means that diabetes WAS NOT observed.
#                        1 means that diabetes WAS observed.
# addcovar:  Matrix containing additive covariates for the mapping model.
# cores:     Number of cores to use in attempting parallel computation. 
#            Not implemented yet.
scan1_surv = function(genoprobs, pheno, addcovar = NULL, 
                      intcovar = NULL, cores = 1) {
  
  register(BPPARAM = SnowParam(workers = cores, progressbar = TRUE))
  
  # Declare return value.
  retval = NULL

  pheno_surv        = Surv(time = pheno$survival, event = pheno$event)
  names(pheno_surv) = rownames(pheno)
  
  lod = vector("list", length(probs))
  
  # Null model log_likelihood.
  null_ll = coxph(pheno_surv ~ addcovar)$loglik[2]

  if(is.null(intcovar)) {
    
    ### Additive covariates ###
    
    # For each chromosome...
    lod = bplapply(seq_along(probs), scan1_surv_addcovar,
                   pr = probs, pheno = pheno_surv, addcovar = addcovar,
                   BPOPTIONS = bpoptions(packages = c('stats', 'survival')))
  
  } else {
    
    ### Interactive covariate ###
    
    # For each chromosome...
    lod = bplapply(seq_along(probs), scan1_surv_intcovar,
                   pr = probs, pheno = pheno_surv, addcovar = addcovar,
                   intcovar = intcovar, 
                   BPOPTIONS = bpoptions(packages = c('stats', 'survival')))

  } # else

  markers = unlist(sapply(probs, function(z) { dimnames(z)[[3]] }))
  
  retval = matrix(unlist(lod), ncol = 1, dimnames = list(markers, "survival"))
  # Convert log-likelihood to LOD.
  retval[,1] = (retval[,1] - null_ll) / log(10)
  class(retval) = c("scan1", "matrix")
  
  return(retval)
  
} # scan1_surv()


# Internal function used by scan1_surv() for addtive covariates.
scan1_surv_addcovar = function(chr, pr, pheno, addcovar) {

  pr = pr[[chr]][rownames(pheno),-1,, drop = FALSE]

  n_markers = dim(pr)[[3]]
  
  lod = data.frame(lod       = rep(0, n_markers),
                   row.names = dimnames(pr)[[3]])
  
  # Get log-likelihood for each marker.
  for(j in 1:n_markers) {
    
    lod[j,1] = coxph(pheno ~ addcovar + pr[,,j])$loglik[2]
    
  } # for(j)
  
  lod
  
} # scan1_surv_addcovar()


# Internal function used by scan1_surv() for addtive and interactive covariates.
scan1_surv_intcovar = function(chr, pr, pheno, addcovar, intcovar) {
  
  pr = pr[[chr]][rownames(pheno),-1,, drop = FALSE]
  
  n_markers = dim(pr)[[3]]
  
  lod = data.frame(lod       = rep(0, n_markers),
                   row.names = dimnames(pr)[[3]])
  
  # Get log-likelihood for each marker.
  for(j in 1:n_markers) {
    
    lod[j,1] = coxph(pheno ~ addcovar + pr[,,j] + intcovar:pr[,,j])$loglik[2]
    
  } # for(j)
  
  lod
  
} # scan1_surv_intcovar()


# DMG: NOT WORKING YET!!!
# Plot the allele effects at one marker.
# This works for something like a backcross or F2. Not sure about DO yet.
# pr: qtl2-style genoprobs object.
# pheno: data.frame containing a survival object. Must have sample IDs in
#       "id" column and survival in "surv" column. 
# addcovar: matrix of additive covariates.
# mkr: name of marker in pr at which to plot effects.
plot_surv_coef = function(pr, pheno, addcovar, mkr) {

  # Get genotype calls.
  gt = maxmarg(pr = pr)
  
  # Get the chromosome where the markers lies.
  dn  = lapply(pr, dimnames)
  dn  = lapply(dn, '[[', 3)
  chr = grep(mkr, dn)
  
  # Get genotypes at this marker.
  gt        = gt[[chr]][,mkr]
  names(gt) = rownames(pr[[chr]])
  
  stopifnot(all(names(gt) == pheno$id))
  
  surv = data.frame(id   = rownames(cross$covar), 
                    surv = Surv(time = cross$pheno[,'age'], event = cross$pheno[,'censored']),
                    treatment = cross$covar[,'treatment'],
                    gt  = gt)
  
  # Fit models for plotting.
  mod_coxph = coxph(surv ~ gt,   data = surv)
  p.value   = summary(mod_coxph)$waldtest['pvalue']
  mod       = survfit(surv ~ gt, data = surv)
  
  # Make the plot.
  plot(mod, conf.int = TRUE, col = 1:2, las = 1, lwd = 2, ...)
  legend('bottomleft', legend = c('DD', 'BD'), col = 1:2, lty = 1)
  text(x = 15, y = 0.1, labels = paste('Cox-PH p =', format(p.value, digits = 3)))

  
} # plot_surv_coef()


# DMG: NOT WORKING YET!!!
scan1_surv_coef = function(genoprobs, pheno, addcovar = NULL, 
                           intcovar = NULL, marker = NULL) {
  
  if(is.null(marker)) {
    
    stop('marker name is a required argument.')
    
  } # if(is.null(marker))
  
  # Declare return value.
  retval = NULL
  
  # Create survival object.
  pheno_surv        = Surv(time = pheno$survival, event = pheno$event)
  names(pheno_surv) = rownames(pheno)
  
  # Get probs at marker.
  dn  = lapply(genoprobs, dimnames)
  dn  = lapply(dn, '[', 3)
  chr = grep(marker, dn)
  idx = grep(marker, dimnames(genoprobs[[chr]])[[3]])
  pr  = genoprobs[[chr]][,,idx]

  if(is.null(intcovar)) {
    
    ### Additive covariates ###
    
    mod = coxph(pheno_surv ~ addcovar + pr)
    eff = coef(mod)[2:3]
    eff[2] = 0
    eff = eff - mean(eff)
    
  } else {
    
    ### Interactive covariate ###
    mod = coxph(pheno_surv ~ addcovar + pr + intcovar:pr)
    
  } # else
  

} # scan1_surv_coef()


### Test code.
#covar = covar[rownames(pheno),,drop = FALSE]

#qtl_surv = suppressWarnings(survscan(pheno_surv, probs, covar))

#qtl_surv = qtl_surv[rownames(qtl),,drop=FALSE]

#qtl = cbind(qtl, qtl_surv)



