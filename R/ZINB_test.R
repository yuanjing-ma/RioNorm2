#' A differential abundance test function
#'
#' This function is for differential abundance test using zero-inflated Negative Binomial model
#' @param disp_OTU overdispersed OTU table in form of OTUs x samples
#' @param class labels for each sample
#' @param log_sizefactor log of size factors
#' @keywords differential abundance test
#' @keywords over-dispersion
#' @keywords Benjamini-Hochberg
#' @export
#' @examples
ZINB_test <- function(disp_OTU, class, log_sizefactor){
  library(MASS)
  library(pscl)
  df_disp = data.frame(t(disp_OTU), class, log_sizefactor)
  ID_disp = rownames(disp_OTU)
  colnames(df_disp) = c(ID_disp, 'class', 'logsizefactor')

  pval_ZINB_foldchange = rep(10, length = length(ID_disp))
  pval_ZINB_zeroinf = rep(10, length = length(ID_disp))
  for (i in 1:length(ID_disp)){
    if (sum(df_disp[,i]==0) >= 1) {
      tryCatch({
        m1 <- zeroinfl(df_disp[,i] ~ df_disp$class + df_disp$logsizefactor | df_disp$logsizefactor, dist = "negbin")
        pval_ZINB_foldchange[i] = summary(m1)$coefficients[[1]][2,4]
        pval_ZINB_zeroinf[i] = summary(m1)$coefficients[[2]][2,4]
      }, error = function(e){})
    }else{
      tryCatch({
        m1 <- glm.nb(df_disp[,i] ~ df_disp$class + df_disp$logsizefactor)
        pval_ZINB_foldchange[i] = summary(m1)$coefficients[2,4]
        pval_ZINB_zeroinf[i] = 10
      }, error = function(e){})
    }
  }
  pval_ZINB_foldchange[is.na(pval_ZINB_foldchange)] = 1

  # Benjamini-Hochberg correction
  pval_ZINB_foldchange_adj = p.adjust(pval_ZINB_foldchange, method = "BH", n = length(pval_ZINB_foldchange))  # 8 pval<0.05

  return(list("padj" = pval_ZINB_foldchange_adj, "pvalue" = pval_ZINB_foldchange, "id" = ID_disp))
}
