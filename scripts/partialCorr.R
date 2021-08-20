#####################
### partialCorr.R ###
#####################

# caculates partial correlation between all permutations of exons present in data
partialCorr <- function(data, exons, covariates){
  
  suppressPackageStartupMessages( library(data.table) )
  suppressPackageStartupMessages( library(ppcor) )
  
  # get all exon-exon pairs
  perm <- gtools::permutations(n = length(exons), r = 2, v = exons, repeats.allowed = T)
  
  # create empty data table for partial correlations
  pcorr <- data.table(matrix(ncol=8, nrow=0))
  colnames(pcorr) <- c("x", "y", "estimate", "p.value", "statistic", "n", "gp", "Method")
  
  # loop over exon-exon pairs
  for(i in 1:nrow(perm)){
    
    # if both exons are the same, simply save a correlation of 1
    if(perm[i,1] == perm[i,2]){
      pcorr <- rbind(use.names = FALSE, pcorr, 
                     data.table(perm[i,1], perm[i,2], 1.0, 0, NA, nrow(data), 2, "spearman"))
      
      # if not, calculate partial correlation and take into account the covariates
    } else { 
      pcorr <- rbind(use.names = FALSE, pcorr, 
                     cbind(perm[i,1], perm[i,2],
                           pcor.test(x = data[ ,perm[i,1], with=F], 
                                     y = data[ ,perm[i,2], with=F], 
                                     z = data[ ,covariates, with=F], method = "spearman")))
    }
  }
  
  # create a matrix with all partial correlations
  pcorr_r <- matrix(pcorr$estimate, 
                    nrow = length(exons), ncol = length(exons), 
                    dimnames = list(levels(pcorr$x), levels(pcorr$y)))
  
  # change the exon order
  pcorr_r <-pcorr_r[exons, exons]
  
  # calculate the FDR-corrected p value for the partial correlations and save them in another matrix
  pcorr_p <- matrix(p.adjust(pcorr$p.value, method = "fdr"), 
                    nrow = length(exons), ncol = length(exons), 
                    dimnames = list(levels(pcorr$x), levels(pcorr$y)))
  
  # change the exon order
  pcorr_p <-pcorr_p[exons, exons]
  
  return(list(r = pcorr_r, p = pcorr_p))
}