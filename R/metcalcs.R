metcalcs <- function(networks, indices, network_level = 'both'){
  #' @title calculate variance in network level metrics caused by MOTU clustering level
  #'
  #' Acts as a wrapper for bipartite's networklevel function, calculating the difference between your given network-level metrics for a series of networks and clustering thresholds
  #'
  #' @param networks A nested list of networks, with the first list level corresponding to MOTU clustering thresholds, the second level corresponding to individual networks
  #' @param indices A vector of indices to be calculated. See the function networklevel in package bipartite for details
  #' @param network_level The network level to analyse
  #' @return Produces a dataframe showing which metrics are robust in your dataset to clustering-level effects
  #' @seealso \code{\link{line_plot}} for visualisation of the resulting data
  #' @export
  #' @examples metcalcs(networks= batnets, indices = ind, network_level = 'higher')
  out <- list()
  for(i in 1:length(networks)){
    out[[names(networks)[i]]] <- lapply(networks[[i]], function(x) as.matrix(bipartite::networklevel(web = x, index = indices, level = network_level)))
  }


  m <- reshape2::melt(out)
  m <- m[,-2] #column two is just a tonne of 1s from the lapply etc, discard it
  colnames(m) <- c('metric', 'value', 'network', 'clustering')
  m$clustering <- as.numeric(m$clustering)
  m$metric <- gsub('_', ', ',m$metric)
  return(m)
}
