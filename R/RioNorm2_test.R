#' Function for RioNorm2 differential abundance test
#'
#' This function allows you to apply RioNorm2 framework and detect differentially abundant OTUs.
#' @param OTU_table pool for OTUs which will be tested for differential abundance in form of OTUs x Samples.
#' @param class dummy variables for class information.
#' @param FDR level for FDR control
#' @keywords RioNorm2
#' @keywords Differential abundance analysis
#' @export
#' @examples
## load in required package

RioNorm2_test <- function(OTU_table, class, FDR = 1) {

  # Step 1: apply "hk_find" function in the package
  # to detect a group of relatively invariant OTUs (riOTUs)
  # use riOTUs to calculate size factors
  size_factor = hk_find(OTU_table, min_avg_counts = 5)$size_factor
  log_sizefactor = log(size_factor)

  # Step 2: apply "overdisp_scoretest" function in the package for over-dispersion test
  # output will be two lists of OTUs, one with dispersed OTUs, the other with non-overdispersed OTUs
  scoretest = overdisp_scoretest(OTU_table, class, log_sizefactor)
  ID_nondisp = scoretest$ID_nondisp
  ID_disp = scoretest$ID_disp

  # Step 3: apply "ZIP_test" function in the package to test differential abundance for non-overdispersed OTUs
  if(length(ID_nondisp)>0){
  nondisp_OTU = OTU_table[ID_nondisp,]
  nondisp_res = ZIP_test(nondisp_OTU, class, log_sizefactor)
  }else{
  nondisp_res = NULL}

  # Step 4: apply "ZINB_test" function in the package to test differential abundance for overdispersed OTUs
  if(length(ID_disp)>0){
  disp_OTU = OTU_table[ID_disp,]
  disp_res = ZINB_test(disp_OTU, class, log_sizefactor)
  }else{disp_res = NULL}

  # combine test results from ZIP and ZINB
  combined_res = apply(cbind(disp_res, nondisp_res),1,unlist)
  # if needed, rownames of combined_res can reflect overdispersion or not

  output = combined_res[as.numeric(combined_res[,1]) <= FDR,]
  return(output)
}


# output is differentially abundant OTUs
# first column is the adjusted p-value
# second column is the raw p-value
# third column is the OTU id

