metcalcs <- function(networks, indices, network_level = 'both', list_format = 'clust_net'){
  #' calculate variance in network level metrics caused by MOTU clustering level
  #'
  #' Acts as a wrapper for bipartite's networklevel and computeModules functions, calculating the difference between your given network-level metrics for a series of networks and clustering thresholds
  #'
  #' @param networks A list of networks
  #' @param indices A vector of indices to be calculated. See the functions networklevel and computeModules in package bipartite for details
  #' @param network_level The network level to analyse
  #' @param list_format The input must be a nested list of networks, with the levels either being clustering, then network identity ('clust_net'), or network identity,
  #' then clustering ('net_clust').
  #' @return Produces a dataframe showing which metrics are robust in your dataset to clustering-level effects
  #' @seealso \code{\link{line_plot}} for visualisation of the resulting data
  #' @export
  #' @examples metcalcs(networks= batnets, indices = ind, network_level = 'higher')

  list_rev <-  function(ll) {#Function for standardising the input data
    nms <- unique(unlist(lapply(ll, function(X) names(X))))
    ll <- lapply(ll, function(X) setNames(X[nms], nms))
    ll <- apply(do.call(rbind, ll), 2, as.list)
    lapply(ll, function(X) X[!sapply(X, is.null)])
  }#Function sourced from https://stackoverflow.com/questions/15263146/revert-list-structure

  #Set things up
  if(!list_format %in% c('clust_net', 'net_clust', 'clust_only')){
    stop("list_format must either be \'clust_net\', \'net_clust\' or \'clust_only\'")
  }
  if(list_format=='net_clust'){ #standardise the format
    networks <- list_rev(networks)
  }
  if(list_format=='clust_only'){
    networks <- list(networks)
  }

  out <- list()
  modularity = F
  if('modularity' %in% indices){ #Modularity-only calculating session, as it uses a different command in bipartite to all other metrics
    indices <- indices[-which(indices=='modularity')]
    modularity = T
    if(list_format == 'clust_only'){
      for(i in 1:length(networks)){
        out[[i]] <- lapply(networks[[i]], function(x) slot(bipartite::computeModules(web = x), 'likelihood'))
      }
      mods <- reshape2::melt(out)
      mods <- cbind(mods, rep('modularity', nrow(mods)))
      mods <- mods[-3]# This value is meaningless if theres no network identity to work with
      colnames(mods) <- c('value', 'clustering', 'metric')
      mods <- mods[,c(3,1,2)] #Reordering the columns to keep it uniform with the networklevel outputs
      mods$clustering <- as.numeric(mods$clustering)
    }
    if(list_format != 'clust_only'){
      for(i in 1:length(networks)){
        out[[names(networks)[i]]] <- lapply(networks[[i]], function(x) slot(bipartite::computeModules(web = x), 'likelihood'))
      }
      mods <- reshape2::melt(out)
      mods <- cbind(mods, rep('modularity', nrow(mods)))
      colnames(mods) <- c('value', 'network', 'clustering', 'metric')
      mods <- mods[,c(4,1,2,3)] #Reordering the columns to keep it uniform with the networklevel outputs
      mods$clustering <- as.numeric(mods$clustering)
    }


  }


  if(length(indices>0)){ #The rest of the metrics
    if(list_format =='clust_only'){
        for(i in 1:length(networks)){
          out[[i]] <- lapply(networks[[i]], function(x) as.matrix(bipartite::networklevel(web = x, index = indices, level = network_level)))
        }
      m <- reshape2::melt(out)
      m <- m[-c(2,5)]# This value is meaningless if theres no network identity to work with
      colnames(m) <- c('metric',  'value', 'clustering')
      m$clustering <- as.numeric(m$clustering)

      }
    if(list_format != 'clust_only'){
      for(i in 1:length(networks)){
        out[[names(networks)[i]]] <- lapply(networks[[i]], function(x) as.matrix(bipartite::networklevel(web = x, index = indices, level = network_level)))
      }
      m <- reshape2::melt(out)
      m <- m[,-2] #column two is just a tonne of 1s from the lapply etc, discard it
      colnames(m) <- c('metric', 'value', 'network', 'clustering')
      m$clustering <- as.numeric(m$clustering)
      m$metric <- gsub('_', ', ',m$metric)
      }

    }


if(modularity ==T & length(indices>0)){
  m <- rbind(mods, m)
  return(m)
}else if(modularity ==T){
  return(mods)
}else if(length(indices>0)){
  return(m)
}



}
