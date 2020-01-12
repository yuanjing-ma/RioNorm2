#' A normalization function
#'
#' This function allows you to calculate normalized OTU table.
#' @param OTU_table pool for OTUs which will be tested for differential abundance in form of OTUs x Samples
#' @param size_factor size factors of all samples for normalization use
#' @keywords RiOTUs a group of relatively invariant OTUs
#' @keywords size factors
#' @keywords normalization
#' @export
#' @examples
hknorm <- function(OTU_table, size_factor){
  # input:
  # OTU_table: for later on different abundance testing
  # size_factor: normalization size factor for each sample
  # output:
  # normalized OTU_table
  norm_table = t(t(OTU_table)/size_factor) # OTUs x samples
  return(norm_table)
}
