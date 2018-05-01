metcalcs <- function(networks, indices, network_level = 'both'){
  out <- list()
  for(i in 1:length(networks)){
    out[[names(networks)[i]]] <- lapply(networks[[i]], function(x) as.matrix(networklevel(web = x, index = indices, level = network_level)))
  }


  m <- reshape2::melt(out)
  m <- m[,-2] #column two is just a tonne of 1s from the lapply etc, discard it
  colnames(m) <- c('metric', 'value', 'network', 'clustering')
  m$clustering <- as.numeric(m$clustering)
  m$metric <- gsub('_', ', ',m$metric)
  return(m)
}
