
#Make a nested list, where the main list has a value for each clustering level, and each item within it is a distinct network generated at that level #####

library(here)
library(bipartite)
library(ggplot2)

setwd('~/Dropbox/Education/PhD/Bat-diet/')
getwd()
r_network_gen <- function(input_network, collapse_species = T, desired_species = NULL, filter_species = F, include_malua=F, lulu= F){

  if(collapse_species==T && !is.null(desired_species)){
    break('Cannot have false for collapsing species AND have a species desired for selection')
  }
  if(interactive()==TRUE){
    library('bipartite')
    library('stringr')
    library('igraph')
    library('reshape')
  }else{
    require(methods)
    library(network, lib.loc = '/data/home/btw863/r_packages/')
    library(statnet.common, lib.loc = '/data/home/btw863/r_packages/')
    library(sna, lib.loc = '/data/home/btw863/r_packages/')
    library(igraph, lib.loc = '/data/home/btw863/r_packages/')
    library(permute, lib.loc = '/data/home/btw863/r_packages/')
    library(vegan, lib.loc = '/data/home/btw863/r_packages/')
    library(bipartite, lib.loc = '/data/home/btw863/r_packages/')
    library(stringr, lib.loc = '/data/home/btw863/r_packages/')
    library(reshape, lib.loc = '/data/home/btw863/r_packages/')
  }

  source('scripts/r/The.matrix.reloader.R')
  source('scripts/r/hernani_comparisons.R')



  #####Reading in data, formatting it for analysis ####

  all_interactions <- #read.table('data/processed_dna_data/25_april_strict_lengths/for_r/95/all_post_QC_otus.txt.table_binary.out', sep = '\t', header = F, stringsAsFactors = F)#
    all_interactions <- input_network

  desired_colnames <- c("Rhbo","Kein","Hice","Kepa","Rhse","Hiri","Rhtr", "Keha", "Hidi", "Hidy") #If we want to filter out species


  #This finds all samples with 'GC' in the name and gives them a useful name
  gc <- grep('GC',all_interactions[1,])
  for(i in 1:length(gc)){
    #print(str_split(all_interactions[1,gc[i]], pattern = '\\.')[[1]][7])
    temp <- str_split(all_interactions[1,gc[i]], pattern = '\\.')[[1]][7]
    all_interactions[1,gc[i]] <- str_split(temp, pattern='_')[[1]][1]
  }

  #This finds all samples without 'GC' in the name and gives them a useful name
  non_gc <- seq(1, ncol(all_interactions))[-gc]
  str_split(all_interactions[1,non_gc[2]], pattern = '\\.')[[1]][1] #Finish this
  for(i in 1:length(non_gc)){
    all_interactions[1,non_gc[i]] <- str_split(all_interactions[1,non_gc[i]], pattern = '\\.')[[1]][1]
  }

  badcols <- c('1774','4437', '2070', '2275', '4260', '4531', "1004", "1007", "1107", "1134", "1165", "1180",  "198",  "209",  "210","387",  "426",  "459",  "497",  "536",  "541",  "567",  "591",  "689","796",  "806",  "822",  "841",  "843",  "899",  "910",  "918",  "986","996", "3712", "4341", "4361",'1774','4437', '2070', '2275', '4260', '4531', '841', '843')#Sadly these columns match two different samples, so must be removed for now until checked against the field data

  all_interactions <- all_interactions[,-which(all_interactions[1,] %in% badcols)]

  field_data <- read.csv('data/Edited_all_files_with_lat_long_VKedits.csv')

  field_data$Site <- gsub('DVCA', 'DANUM', field_data$Site)
  field_data$Site <- gsub('MALUA', 'SBE', field_data$Site)
  field_data$Faeces_no1 <- gsub('T', '', field_data$Faeces_no1)
  field_data$Faeces_no2 <- gsub('T', '', field_data$Faeces_no2)

  field_data$SiteAndYear <- paste(field_data$Site, field_data$Year, sep = ', ')

  #badsites <- c('DV88, 2016', 'MALUA, 2016', 'DV89, 2016') #These ones won't make viable networks for various reasons
  if(include_malua==F){
    badsites <- c('DV88, 2016', 'SBE, 2016', 'DV89, 2016') #These ones won't make viable networks for various reasons
  }else{
    badsites <- c('DV88, 2016', 'DV89, 2016') #These ones won't make viable networks for various reasons
  }

  field_data <- field_data[-which(field_data$SiteAndYear %in% badsites),]

  if(collapse_species==T){

    nets <- lapply(unique(field_data$SiteAndYear), function(i) the.matrix.reloader(master.data = field_data, ID.column.1 = "Faeces_no1", ID.column.2 = 'Faeces_no2', species.column = "Species", split.by.column = "SiteAndYear", split.by.var = i, OTU.matrix = all_interactions))

    names(nets) <- unique(field_data$SiteAndYear)
    for(i in 1:length(nets)){
      print(names(nets)[i])
      print(colnames(nets[[i]]))
    }


    if(filter_species==T){
      for(i in 1: length(nets)){
        if(length(which(!colnames(nets[[i]]) %in% desired_colnames))>0){#If there are species that we don't want
          to_remove <- which(!colnames(nets[[i]]) %in% desired_colnames)
          nets[[i]] <- nets[[i]][,-to_remove]
        }
      }
    }



    return(nets)
  }
  else if(collapse_species==F){
    all_interactions <- rbind(all_interactions[1,], all_interactions)
    if(is.null(desired_species)){
      if(filter_species==F){
        for(i in 1:ncol(all_interactions)){
          if(as.character(all_interactions[1,i]) %in% field_data$Faeces_no1){
            row <- field_data[which(field_data$Faeces_no1==as.character(all_interactions[1,i])),]
            site <- row$SiteAndYear
            all_interactions[1,i] <- as.character(site)
          }else if(as.character(all_interactions[1,i]) %in% field_data$Faeces_no2){
            row <- field_data[which(field_data$Faeces_no2==as.character(all_interactions[1,i])),]
            site <- row$SiteAndYear
            all_interactions[1,i] <- as.character(site)
          }
        }
      }else if(filter_species==T){
        badcols <- c()
        for(i in 1:ncol(all_interactions)){
          if(as.character(all_interactions[1,i]) %in% field_data$Faeces_no1){
            row <- field_data[which(field_data$Faeces_no1==as.character(all_interactions[1,i])),]
            site <- row$SiteAndYear
            all_interactions[1,i] <- as.character(site)
            if(!row$Species %in% desired_colnames){
              badcols <- c(badcols,i)
            }
          }else if(as.character(all_interactions[1,i]) %in% field_data$Faeces_no2){
            row <- field_data[which(field_data$Faeces_no2==as.character(all_interactions[1,i])),]
            site <- row$SiteAndYear
            all_interactions[1,i] <- as.character(site)
            if(!row$Species %in% desired_colnames){
              badcols <- c(badcols,i)
            }
          }
        }
        all_interactions <- all_interactions[,-badcols]
      }


      # if(filter_species==T){
      #
      #   for(i in 1: length(nets)){
      #     if(length(which(!colnames(nets[[i]]) %in% desired_colnames))>0){#If there are species that we don't want
      #       to_remove <- which(!colnames(nets[[i]]) %in% desired_colnames)
      #       nets[[i]] <- nets[[i]][,-to_remove]
      #     }
      #   }
      # }
      return(all_interactions)
    }else{
      all_interactions_with_extra <- rbind(all_interactions[1,], all_interactions)
      all_interactions_with_extra <- the.matrix.reloader(master.data = field_data, ID.column.1 = "Faeces_no1", ID.column.2 = 'Faeces_no2',species.column = "SiteAndYear", split.by.column = "Species", split.by.var = desired_species, OTU.matrix = all_interactions_with_extra, collapse_top_species = F)
      return(all_interactions_with_extra)
    }


  }

}




filenames <- list.files(pattern = 'all_post_QC_otus.txt.table_binary.out', recursive = T)

filenames <- filenames[-grep('galaxy', filenames)]
filenames <- filenames[-grep('lulu', filenames)]

rawnets <- lapply(filenames, read.table,  sep = '\t', header = F, stringsAsFactors = F)
names(rawnets) <- gsub('.*\\/', '', dirname(filenames))
netlists <- lapply(rawnets, r_network_gen, collapse_species = T, filter_species = T)
names(netlists) <- names(rawnets)


#Specify which index(s) to use

ind <- c('connectance', 'web asymmetry')

calctime <- c()
for(i in 1:length(ind)){
  s <- Sys.time()
  networklevel(netlists$`91`$`SAFE, 2015`, index = ind[i])
  f <- Sys.time()
  calctime <- c(calctime, f-s)
}

timings <- data.frame(ind, calctime)

sums_to_preserve = 'both'

outmat <- matrix(nrow = 0, ncol = 4)
for(index in 1:length(ind)){
  index_used <- ind[index]
  rand_list <- list()
  for(i in 1:length(netlists)){
    nam <- names(netlists)[i]
    print(nam)
    rand_list[[i]] <- lapply(netlists[[i]], function(x)
      replicate(1000, networklevel(permatfull(x, fixedmar=sums_to_preserve,mtype="count",times=1)$perm[[1]],
                                   index = index_used)))
    names(rand_list)[i] <- nam
  }

  mat <- matrix(nrow = 0, ncol = 2)
  for(a in 1:length(rand_list)){
    print(names(rand_list)[a])
    mat <- rbind(mat, t(sapply(rand_list[[i]], function(x) quantile(x, probs=c(0.025, 0.975)))))
  }
  mat <- cbind(rownames(mat), rep(index_used, nrow(mat)), mat)
  outmat <- rbind(outmat, mat)
}





rand_list <- list()
for(i in 1:length(netlists)){
  nam <- names(netlists)[i]
  print(nam)
  rand_list[[i]] <- lapply(netlists[[i]], function(x)
    replicate(1000, networklevel(permatfull(x, fixedmar=sums_to_preserve,mtype="count",times=1)$perm[[1]],
                                 index = ind)))
  names(rand_list)[i] <- nam
}

a <- lapply(netlists[[1]], function(x)
  replicate(1000, networklevel(permatfull(x, fixedmar=sums_to_preserve,mtype="count",times=1)$perm[[1]],
                               index = ind)))
#SOMETHINGS NOT WORKED the values are the same for each clustering level

mat <- matrix(nrow = 0, ncol = 2)
for(i in 1:length(rand_list)){
  print(names(rand_list)[i])
  mat <- rbind(mat, t(sapply(rand_list[[i]], function(x) quantile(x, probs=c(0.025, 0.975)))))
}
mat <- cbind(rownames(mat), rep(ind[1], nrow(mat)), mat) ###CHANGE THE IND WHEN THIS LOOPS
#sapply(rand_list[[1]], function(x) quantile(x, probs=c(0.025, 0.975)))
