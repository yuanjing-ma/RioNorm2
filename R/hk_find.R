#' Function for finding a group of relatively invariant OTUs 
#'
#' This function allows you to find a group of relatively invariant OTUs and calculate size factors for normalization.
#' @param OTU_table pool for OTUs which will be tested for differential abundance in form of OTUs x Samples and rownames of OTU_table are OTU ids
#' @param min_avg_counts minimum average counts required for riOTUs
#' @keywords relatively invariant OTUs
#' @keywords size factors
#' @keywords normalization
#' @export
#' @examples
hk_find <- function(OTU_table, min_avg_counts = 5){
  # output:
  # A list of riOTUs' id and size factor used for normalization
  # result = hk_find(OTU_table, min_avg_counts)
  # riOTUs_id = result$riOTUs_ID
  # size_factor = result$size_factor
  
  nsamples = dim(OTU_table)[2]
  samobs = apply(OTU_table, 1, function(x) sum(x != 0))
  hk_pool = OTU_table[samobs >= nsamples * 0.8,]
  avg_count = apply(hk_pool, 1, mean)
  hk_pool = hk_pool[avg_count >= min_avg_counts,]
    
  nOTUs = dim(hk_pool)[1]
  OTU_ID = rownames(hk_pool)
  # create symmetric distance matrix between selected OTUs
  ratio_var = matrix(0, nrow = nOTUs, ncol = nOTUs)
  for (i in 1:nOTUs){
    for (j in 1:nOTUs){
      mul = (hk_pool[i,]*hk_pool[j,])
      ind = unname(unlist(mul)) != 0
      ratio_var[i,j] = var(unlist(log(hk_pool[i,][ind]/hk_pool[j,][ind])))
    }
  }
  ratio_var[lower.tri(ratio_var, diag = TRUE)] <- 0
  dist = ratio_var[ratio_var > 0]
  
  First = TRUE
  library(igraph)
  quans = seq(0.01,0.1,0.005)
  for (q in 1:length(quans)){
    print(q)
    nodes = data.frame(seq(1, nOTUs, 1), OTU_ID)
    colnames(nodes) = c("order","OTU_ID")
    # Build links dataset
    links = matrix(0, nrow = 1, ncol=3)
    dist = ratio_var[ratio_var>0]
    h = quantile(dist, probs = quans[q])
    order = seq(1,nOTUs,1)
    for (i in 1:dim(ratio_var)[1]){
      check = ratio_var[i,]
      ind = order[check < h & check > 0]
      if (length(ind)>0){
        dist = check[ind]
        rep = replicate(length(ind),i)
        record = cbind(rep,ind,dist)
        links = rbind(links,record)
      }
    }
    links = links[-1,]
    
    # plot the network
    net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
    
    # Find cliques
    list_cliques = cliques(net) 
    clique_size = sapply(list_cliques, length) 
    curr_largest_cliques = largest_cliques(net) 
    curr_length = length(curr_largest_cliques[[1]])
    curr_largest_cliques = matrix(unlist(curr_largest_cliques), ncol = curr_length, byrow = T)
    
    if (curr_length < 3) {
      next
    }
    
    if (First == TRUE) {
      prev_largest_cliques = curr_largest_cliques
      prev_length = dim(prev_largest_cliques)[2]
      First = FALSE
      next
    }
    
    if (curr_length == prev_length){
      riOTUs = prev_largest_cliques
      break  
    }
    
    print(prev_largest_cliques)
    print(curr_largest_cliques)
    updated_cliques = matrix(0, nrow = 1, ncol = curr_length)
    for (i in 1:dim(curr_largest_cliques)[1]) {
      for (j in 1:dim(prev_largest_cliques)[1]) {
        if (sum(prev_largest_cliques[j,] %in% curr_largest_cliques[i,]) == prev_length) {
          updated_cliques = rbind(updated_cliques, curr_largest_cliques[i,])
          break
        }
      }
    }
    
    if (dim(updated_cliques)[1] == 1) {
      riOTUs = prev_largest_cliques
      break
    } else if (q == length(quans) & dim(updated_cliques)[1] == 2) {
      riOTUs = matrix(updated_cliques[-1,], nrow = 1)
      break 
    } else if (dim(updated_cliques)[1] == 2) {
      prev_largest_cliques = matrix(updated_cliques[-1,], nrow = 1)
      prev_length = dim(prev_largest_cliques)[2]
    } else{
      prev_largest_cliques = updated_cliques[-1,]
      prev_length = dim(prev_largest_cliques)[2]
    }
  }
  riOTUs = riOTUs[order(riOTUs)]
  riOTUs_ID = rownames(hk_pool)[riOTUs]
  riOTUs_count = hk_pool[riOTUs,]
  size_factor = colSums(riOTUs_count)
  return(list(riOTUs_ID = riOTUs_ID, size_factor = size_factor))
}
