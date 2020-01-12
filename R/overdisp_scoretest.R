#' A overdispersion test function
#'
#' This function allows you to test for over-dispersion in zero-inflated count regression model
#' @param OTU_table pool for OTUs which will be tested for differential abundance in form of OTUs x Samples
#' @param class labels for each sample
#' @param log_sizefactor log of size factors
#' @keywords overdispersion
#' @keywords score test
#' @export
#' @examples
overdisp_scoretest <- function(OTU_table, class, log_sizefactor){
  library(pscl)
  library(MASS)
  # input:
  # OTU_table: pool for OTUs which will be tested for differential abundance
  # OTU_table are sorted from most abundant to least abundant
  # outputs:
  # IDs(rownames) for non-dispersed OTUs and dispersed OTUs separatelyß
  samobs = apply(OTU_table, 1, function(x) sum(x != 0))
  sums = apply(OTU_table, 1, sum)
  otudf = data.frame(prev = samobs, sums = sums)
  otudf = otudf[order(-otudf$prev, -otudf$sums),]
  OTU_table = OTU_table[rownames(otudf),]

  nsamples = dim(OTU_table)[2]
  ntest_pool = dim(OTU_table)[1] # number of OTUs tested for DA
  test_OTU_ID = rownames(OTU_table)
  df = data.frame(t(OTU_table), class, log_sizefactor)
  colnames(df) =c(test_OTU_ID, 'class', 'logsizefactor')
  # store coefficients from logPoisson regression
  para_logPoi = matrix(0, nrow = ntest_pool, ncol = 3)
  para_zeroinf = matrix(0, nrow = ntest_pool, ncol = 2)
  num_nonzero_OTU = sum(apply(OTU_table, 1, function(x) {sum(x==0) == 0}))
  for (i in 1:ntest_pool){
    if (sum(df[,i]==0) >= 1) {
      tryCatch({
        m1 <- zeroinfl(df[,i] ~ df$class + df$logsizefactor | df$logsizefactor )
        para_logPoi[i,] = summary(m1)$coefficients[[1]][,1]
        para_zeroinf[i,] = summary(m1)$coefficients[[2]][,1]
      }, error = function(e){})
    }else{
      tryCatch({
        m1 <- glm(df[,i] ~ df$class + df$logsizefactor, family = poisson)
        para_logPoi[i,] = summary(m1)$coefficients[,1]
        para_zeroinf[i,] = c(0,0)
      }, error = function(e){})
    }
  }
  rownames(para_logPoi) = test_OTU_ID
  rownames(para_zeroinf) = test_OTU_ID

  # MLE of mu
  B = rbind(rep(1, nsamples), class, log_sizefactor)
  mu = para_logPoi %*% B
  mu = exp(mu)
  # MLE of pi
  A = rbind(rep(1, nsamples), log_sizefactor)
  expo = exp(para_zeroinf %*% A)
  pi = expo / (expo + 1)
  # expo can be Inf in some cases, cause pi be NA
  pi[is.na(pi)] = 1
  # adjust for OTUs without 0 count
  pi[1:num_nonzero_OTU,] = 0

  # Compute score function for each OTU
  index_0 = OTU_table
  index_0[index_0 == 0] <- 0.1
  index_0[index_0 >= 1] <- 0
  index_0 = index_0*10

  U = ((OTU_table-mu)^2-OTU_table) - index_0*mu^2*pi/(pi+(1-pi)*exp(-mu) + 1e-300)
  U = 0.5 * apply(U,1,sum)
  # compute inverse of Fisher Information matrix
  inv_I = mu^2*(2*(1-pi)-mu^2*pi*(1-pi/(pi+(1-pi)*exp(-mu) + 1e-300)))
  inv_I = 0.25*apply(inv_I,1,sum)

  # test statistic: sqrt of Score statistic
  test_stat = U*sqrt(inv_I)
  # calculate p-values
  overdisp_p = 1-pnorm(test_stat)

  # Extract OTUs with and without over-dispersion
  nondisp_OTU = OTU_table[overdisp_p > 0.1,]
  ID_nondisp = test_OTU_ID[overdisp_p > 0.1]
  disp_OTU = OTU_table[overdisp_p <= 0.1,]
  ID_disp = test_OTU_ID[overdisp_p <= 0.1]
  return(list(ID_nondisp = ID_nondisp, ID_disp = ID_disp))
}
