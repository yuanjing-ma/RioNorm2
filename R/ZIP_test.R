#' A differential abundance test function
#'
#' This function is for differential abundance test using zero-inflated Poisson model
#' @param nondisp_OTU non-overdispersed OTU table in form of OTUs x samples
#' @param class labels for each sample
#' @param log_sizefactor log of size factors
#' @keywords differential abundance test
#' @keywords non-overdispersion
#' @keywords Benjamini-Hochberg
#' @export
#' @examples
ZIP_test <- function(nondisp_OTU, class, log_sizefactor){
  library(pscl)
  library(MASS)
  # input:
  # table for non-overdispersed OTUs (OTUs x samples)
  # log_sizefactor: log of size factors used for normalization
  df_nondisp = data.frame(t(nondisp_OTU), class, log_sizefactor)
  ID_nondisp = rownames(nondisp_OTU)
  colnames(df_nondisp) = c(ID_nondisp, 'class', 'logsizefactor')

  pval_ZIP_foldchange = rep(0, length = length(ID_nondisp))
  pval_ZIP_zeroinf = rep(0, length = length(ID_nondisp))
  # if num of zero counts for OTU > 0, use ZIP
  # otherwise, use normal poisson regression
  for (i in 1:length(ID_nondisp)){
    if (sum(df_nondisp[,i]==0) >= 1) {
      tryCatch({
        m1 <- zeroinfl(df_nondisp[,i] ~ df_nondisp$class + df_nondisp$logsizefactor | df_nondisp$logsizefactor )
        pval_ZIP_foldchange[i] = summary(m1)$coefficients[[1]][2,4]
        pval_ZIP_zeroinf[i] = summary(m1)$coefficients[[2]][2,4]
      }, error = function(e){})
    }else{
      tryCatch({
        m1 <- glm(df_nondisp[,i] ~ df_nondisp$class + df_nondisp$logsizefactor, family = poisson)
        pval_ZIP_foldchange[i] = summary(m1)$coefficients[2,4]
        pval_ZIP_zeroinf[i] = 10
      }, error = function(e){})
    }
  }
  
  pval_ZIP_foldchange[is.na(pval_ZIP_foldchange)] = 1
  # Benjamini-Hochberg correction
  pval_ZIP_foldchange_adj = p.adjust(pval_ZIP_foldchange, method = 'BH', n = length(pval_ZIP_foldchange))

  return(list("padj" = pval_ZIP_foldchange_adj, "pvalue" = pval_ZIP_foldchange, "id" = ID_nondisp))
}
