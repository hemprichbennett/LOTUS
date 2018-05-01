metcalcs <- function(lis, indexes, netlevel = 'both'){
  out <- list()
  for(i in 1:length(lis)){
    out[[names(lis)[i]]] <- lapply(lis[[i]], function(x) as.matrix(networklevel(web = x, index = indexes, level = netlevel)))
  }


  m <- reshape2::melt(out)
  m <- m[,-2] #column two is just a tonne of 1s from the lapply etc, discard it
  colnames(m) <- c('metric', 'value', 'network', 'clustering')
  m$clustering <- as.numeric(m$clustering)
  m$metric <- gsub('_', ', ',m$metric)
  return(m)
}
