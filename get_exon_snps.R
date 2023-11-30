################################################################################
# Given a single genomic region and a set of strain names, return SNPs which
# intersect with exons.
# The user must provide the full path to the Sanger VCF file.
# We assume Sanger VCF file format. Other VCFs may not work.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-11-29
################################################################################
options(stringsAsFactors = FALSE)

##### LIBRARIES #####

require(GenomicRanges)
require(VariantAnnotation)

# Return the available strain names in the VCF file.
# Arguments:
# vcf_file: character string containing the full path to the Sanger VCF file.
# Returns:
# character vector containing strain names in VCF file.
get_strain_names = function(vcf_file) {
  
  if(!file.exists(vcf_file)) {
    
    stop(paste('Cannot find VCF file:', vcf_file))
    
  } # if(!file.exists(vcf_file))
  
  hdr = scanVcfHeader(vcf_file)
  
  return(samples(hdr))
  
} # get_strain_names()


# Given a genomic region, strain names, and a VCF file path,
# return the polymorphic SNPs that intersect with exons.
# Arguments:
# region: GRanges object containing one genomic range.
# strains: character vector containing the strains to query. You may 
# include C57BL/6J. If C57BL/6J is not included, we will not use it
# to find polymorphic SNPs.
get_exon_snps = function(region = NULL, strains = NULL, vcf_file) {
  
  # Error handling.
  if(!file.exists(vcf_file)) {
    
    stop(paste('Cannot find VCF file:', vcf_file))
    
  } # if(!file.exists(vcf_file))
  
  if(is.null(region)) {
    
    stop('Region cannot be empty.')
    
  } # if(is.null(region))
  
  if(is.null(strains)) {
    
    stop('Strains cannot be empty.')
    
  } # if(is.null(strains))
  
  all_strains = get_strain_names(vcf_file)
  
  if(!all(strains %in% all_strains)) {
    
    wh = which(!strains %in% all_strains)
    
  } # if(all(!strains %in% all_strains))

  param_strains = strains[strains != 'C57BL_6J']
  param = ScanVcfParam(samples = param_strains, which = region)
  vcf   = readVcf(file   = vcf_file,
                  genome = '',
                  param = param)
  rm(param_strains)
  
  # Get SNPs from VCF.
  gt = geno(vcf)$GT
  
  # Add in C57BL/6J, if needed.
  if('C57BL_6J' %in% strains) {
    gt = cbind(C57BL_6J = '0/0', gt)
  } # if('C57BL_6J' %in% strains)
  
  # Subset to retain polymorphic SNPs.
  vcf = subset(vcf, rowMeans(gt == '0/0') < 1.0)
  vcf = genotypeCodesToNucleotides(vcf)
  vcf = expand(vcf)
  
  # Get SNPs from VCF again.
  gt = geno(vcf)$GT
  
  # Get SNP consequences.
  csq = info(vcf)$CSQ
  csq = lapply(csq, strsplit, split = '\\|')
  csq = lapply(csq, function(z) { matrix(unlist(z), ncol = length(z[[1]]), byrow = TRUE) })
  csq_type = lapply(csq, function(z) { unique(z[,2]) })
  
  # Filter to retain SNPs in coding exons.
  csq = lapply(csq, function(z) { z[z[,8] == 'protein_coding',] })
  vcf = subset(vcf, sapply(csq, length) > 0)

  # Get SNP consequences again.
  csq = info(vcf)$CSQ
  csq = lapply(csq, strsplit, split = '\\|')
  csq = lapply(csq, function(z) { matrix(unlist(z), ncol = length(z[[1]]), byrow = TRUE) })
  csq = lapply(csq, function(z) { z[z[,8] == 'protein_coding',,drop = FALSE] })
  
  csq_codes = 'intron|prime|missense|regulatory|splice|stop|synonymous'
  
  csq = lapply(csq, function(z) { z[grepl(csq_codes, z[,2]),, drop = FALSE] })
  
  # Filter VCF again.
  vcf = subset(vcf, sapply(csq, length) > 0)
  
  # The VCF is an "expandedVCF", which means that SNPs which are tri- or tetra-
  # allelic show up in multiple rows. Retain rows for which the observed
  # alleles match the alternate allele.
  gt  = geno(vcf)$GT
  
  # Add in C57BL/6J, if needed.
  if('C57BL_6J' %in% strains) {
    ref = as.character(fixed(vcf)$REF)
    gt = cbind(C57BL_6J = paste(ref, ref, sep = '/'), gt)
  } # if('C57BL_6J' %in% strains)
  
  ref = as.character(fixed(vcf)$REF)
  alt = as.character(fixed(vcf)$ALT)
  obs_alleles = apply(gt, 1, strsplit, split = '/')
  obs_alleles = lapply(lapply(obs_alleles, unlist), unique)
  vcf_alleles = cbind(as.character(ref), as.character(alt))
  
  keep = rep(TRUE, length(vcf))
  for(i in seq_along(obs_alleles)) {
    
    keep[i] = all(vcf_alleles[i,] %in% obs_alleles[[i]])

  } # for(i)
  
  # Subset the VCF to retain SNPs for which we observe the alleles in
  # the VCF.
  vcf = subset(vcf, keep)
  
  # Get the genotypes again.
  gt  = geno(vcf)$GT
  
  # Add in C57BL/6J, if needed.
  if('C57BL_6J' %in% strains) {

    ref = as.character(fixed(vcf)$REF)
    gt = cbind(C57BL_6J = paste(ref, ref, sep = '/'), gt)

  } # if('C57BL_6J' %in% strains)

  # Get the SNP postions, ref, & alt alleles.
  rr  = rowRanges(vcf)
  
  # Format the results into a table.
  retval = as.data.frame(rr)[,c('seqnames', 'start', 'REF', 'ALT')]
  retval = cbind(retval, gt)

  csq = info(vcf)$CSQ
  csq = lapply(csq, strsplit, split = '\\|')
  csq = lapply(csq, function(z) { matrix(unlist(z), ncol = length(z[[1]]), byrow = TRUE) })
  
  # Go through each consequence and retain values for the current ALT allele
  # and retain unique values.
  csq = lapply(csq, function(z) { z[grepl(csq_codes, z[,2]),, drop = FALSE] })
  for(i in seq_along(csq)) {
    
    csq[[i]] = csq[[i]][csq[[i]][,1] == retval$ALT[i],,drop = FALSE]
    csq[[i]] = data.frame(consequence = paste(unique(csq[[i]][,2]), collapse = ';'),
                          gene_id     = paste(unique(csq[[i]][,5]), collapse = ';'),
                          symbol      = paste(unique(csq[[i]][,4]), collapse = ';'))
    
  } # for(i)
  
  csq = do.call(rbind, csq)
  
  stopifnot(nrow(csq) == nrow(retval))
  
  retval = cbind(retval, csq)
  
  colnames(retval)[1:4] = c('chr', 'pos', 'ref', 'alt')
  
  return(retval)
  
} # get_exon_snps()


# Test code.
#vcf_file = '/projects/omics_share/mouse/GRCm39/genome/annotation/snps_indels/rel_2112_v8/mgp_REL2021_snps.vcf.gz'

#get_strain_names(vcf_file)

#snps = get_exon_snps(region   = GRanges('1', IRanges(3e6, 5e6)),
#                     strains  = c('C57BL_6J', 'DBA_2J', 'CAST_EiJ'),
#                     vcf_file = vcf_file)
#stopifnot(nrow(snps) != 11683)

#snps = get_exon_snps(region   = GRanges('1', IRanges(3e6, 5e6)),
#                     strains  = c('C57BL_6J', 'DBA_2J'),
#                     vcf_file = vcf_file)
#stopifnot(nrow(snps) != 5715)




