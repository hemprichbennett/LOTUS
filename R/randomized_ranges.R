randomized_ranges <- function(networks, input_format = 'clust_net', indices, network_level = 'both', n_perm=1000, sums_to_preserve='both', summarise=T, quantiles_to_return=c(0.025, 0.975), out_format='data.frame', actual_vals=F){
  #' Identify meaningful trends emerging from MOTU clustering thresholds
  #'
  #'
  #' Acts as a wrapper for the \code{\link{metcalcs}} function, randomising a networks across a series of clustering levels
  #'
  #' @param networks A nested list of networks
  #' @param input_format The input must be a list of networks, either a single network identity generated at different clustering levels ('clust_only'),
  #' or a nested list of networks, with the levels either being clustering, then network identity ('clust_net'), or network identity,
  #' then clustering ('net_clust').
  #' @param indices A vector of indices to be calculated. See the functions networklevel and computeModules in package bipartite for details
  #' @param network_level The network level to analyse
  #' @param n_perm The number of permutations to run
  #' @param sums_to_preserve The sums to be preserved in randomisation. Possible values are 'none', 'rows', 'columns', or 'both'
  #' @param summarise Should the randomised values be summarised to the upper and lower quantiles, or should all randomised values be returned?
  #' @param quantiles_to_return The quantiles desired if summarising the output
  #' @param out_format The format for the data to be output in. Either a dataframe ('data.frame') or list ('list')
  #' @param actual_vals Should the actual values for each network be calculated?
  #' @return Produces either a dataframe or list for your desired metrics and randomisations
  #' @seealso \code{\link{metcalcs}}
  #' @export
  #' @examples

  list_rev <-  function(ll) {#Function for standardising the input data
    nms <- unique(unlist(lapply(ll, function(X) names(X))))
    ll <- lapply(ll, function(X) setNames(X[nms], nms))
    ll <- apply(do.call(rbind, ll), 2, as.list)
    lapply(ll, function(X) X[!sapply(X, is.null)])
  }#Function sourced from https://stackoverflow.com/questions/15263146/revert-list-structure

  #Set things up
  if(!out_format %in% c('list', 'data.frame')){
    stop('out_format is incorrect, can be either \'list\' or \'data.frame\'')
  }
  if(!is.na(quantiles_to_return) && summarise==F){
    stop('you have requested specific quantiles to be returned but then asked for the raw data. You can only have one of these. Perhaps return the raw data and then use the \'quantiles\' function on it yourself?')
  }
  if(!sums_to_preserve %in% c('none', 'rows', 'columns', 'both')){
    stop('incorrect value for sums_to_preserve, acceptable values are \'none\', \'rows\', \'columns\', \'both\'')
  }
  if(!input_format %in% c('clust_net', 'net_clust', 'clust_only')){
    stop("input_format must either be \'clust_net\', \'net_clust\' or \'clust_only\'")
  }
  if(input_format=='net_clust'){ #standardise the format for use in metcalcs
    networks <- list_rev(networks)
    input_format <- 'clust_net'
  }
  if(input_format=='clust_only'){
    networks <- list(networks)
  }


  if(out_format=='data.frame'){
    if(actual_vals==T){
      outmat <- matrix(nrow = 0, ncol = 4+length(quantiles_to_return))
    }else{
      outmat <- matrix(nrow = 0, ncol = 3+length(quantiles_to_return))
    }

  }
  if(out_format=='list'){
    penultimate_list <- list()
    out_list <- list()
  }

  if(actual_vals==T){
    actual <- metcalcs(networks=networks, indices= indices, network_level = network_level, list_format = input_format)
    actual$metric <- gsub('\\.', ' ', actual$metric)
    actual$metric <- gsub(' HL', '', actual$metric)
    actual$metric <- gsub(' LL', '', actual$metric)
    actual$metric <- gsub('ISA', 'interaction strength asymmetry', actual$metric)
  }

#####Do the network generation and index calculation ####
  for(index in 1:length(indices)){
    index_used <- indices[index]
    rand_list <- list()
    for(i in 1:length(networks)){
      nam <- names(networks)[i]
      #cat(nam, index_used, '\n')
      if(!is.na(network_level)){
        rand_list[[i]] <- lapply(networks[[i]], function(x)
          replicate(1000, bipartite::networklevel(permatfull(x, fixedmar=sums_to_preserve,mtype="count",times=1)$perm[[1]],
                                       index = index_used, level = network_level)))
      }else{
        rand_list[[i]] <- lapply(networks[[i]], function(x)
          replicate(1000, bipartite::networklevel(permatfull(x, fixedmar=sums_to_preserve,mtype="count",times=1)$perm[[1]],
                                       index = index_used)))
      }

      names(rand_list)[i] <- nam
    }

    if(out_format=='data.frame'){
      if(actual_vals==T){
        #Sort the matrix size
      mat <- matrix(nrow = 0, ncol = 2+length(quantiles_to_return))
      }else{
        mat <- matrix(nrow = 0, ncol = 1+length(quantiles_to_return))
      }
      for(a in 1:length(rand_list)){
        #This creates a matrix for all networks of clustering level a, index i
        m <- cbind(rep(names(rand_list)[a], length(rand_list[[a]])), t(sapply(rand_list[[a]], function(x) quantile(x, probs=quantiles_to_return))))
        if(actual_vals==T){# Here we need to add to m the appropriate real values
          m <- cbind(m, rep(NA, nrow(m)))
          mat <- rbind(mat, m)

        }else{
          mat <- rbind(mat, m)
        }


      }
      if(actual_vals==T){
        print(index_used)
        print(actual)

        mat[,ncol(mat)] <- actual[actual$metric==index_used,'value']
      }


      mat <- cbind(rownames(mat),  rep(index_used, nrow(mat)), mat)
      outmat <- rbind(outmat, mat)
      outmat <- as.data.frame(outmat)
      for(d in 4:ncol(outmat)){
        outmat[,d] <- as.numeric(as.character(outmat[,d]))
      }

    }
    if(out_format=='list'){
      penultimate_list[[index_used]] <- rand_list
    }
  }

  if(out_format=='data.frame'){
    rownames(outmat) <- seq(1,nrow(outmat)) #Otherwise it has annoying rownames
    if(actual_vals==T){
      #outmat <- cbind(outmat, actual$value)
      colnames(outmat)[seq(1,3)] <- c('network', 'metric', 'clustering')
      colnames(outmat)[ncol(outmat)] <- 'actual'
      return(outmat)
    }else{
      colnames(outmat)[seq(1,3)] <- c('network', 'metric', 'clustering')
      return(outmat)
    }

  }
  if(out_format=='list'){
    if(actual_vals==T){
      out_list$actual <- actual
      out_list$randomized <- penultimate_list
      return(out_list)
    }else{
      return(penultimate_list)
    }

  }
}
