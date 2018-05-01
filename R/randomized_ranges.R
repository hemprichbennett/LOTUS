randomized_ranges <- function(networks, indices, network_level = 'both', n_perm=1000, sums_to_preserve='both', out_format='data.frame', quantiles_to_return=c(0.025, 0.975), summarise=T){
  if(!out_format %in% c('list', 'data.frame')){
    stop('out_format is incorrect, can be either \'list\' or \'data.frame\'')
  }
  if(!is.na(quantiles_to_return) && summarise==F){
    stop('you have requested specific quantiles to be returned but then asked for the raw data. You can only have one of these. Perhaps return the raw data and then use the \'quantiles\' function on it yourself?')
  }
  if(!sums_to_preserve %in% c('none', 'rows', 'columns', 'both')){
    stop('incorrect value for sums_to_preserve, acceptable values are \'none\', \'rows\', \'columns\', \'both\'')
  }


  if(out_format=='data.frame'){
    outmat <- matrix(nrow = 0, ncol = 3+length(quantiles_to_return))
  }
  if(out_format=='list'){
    outlist <- list()
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
      mat <- matrix(nrow = 0, ncol = 1+length(quantiles_to_return))
      for(a in 1:length(rand_list)){

        m <- cbind(rep(names(rand_list)[a], length(rand_list[[a]])), t(sapply(rand_list[[a]], function(x) quantile(x, probs=quantiles_to_return))))
        mat <- rbind(mat, m)
        #print(a)
        #print(mat)
      }

      mat <- cbind(rownames(mat),  rep(index_used, nrow(mat)), mat)
      outmat <- rbind(outmat, mat)
      outmat <- as.data.frame(outmat)
      for(d in 4:ncol(outmat)){
        outmat[,d] <- as.numeric(as.character(outmat[,d]))
      }

    }
    if(out_format=='list'){
      outlist[[index_used]] <- rand_list
    }
  }

  if(out_format=='data.frame'){
    colnames(outmat)[seq(1,3)] <- c('network', 'metric', 'clustering')
    return(outmat)
  }
  if(out_format=='list'){
    return(outlist)
  }
}
