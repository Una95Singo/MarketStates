# ASCP state discovery
# Unarine Singo (R recode from “Potts Model Clustering” Author : Lionel Yelibi, 2018 - https://github.com/tehraio/potts-model-clustering/blob/master/aspc_v2a.py)
# 8 June 2020

# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/segmentation.R")
source("R/simulate.R")
source("R/stylisedFacts.R")
source("R/05_ICC.R")
source("R/NohAnzatsStateSimulation.R")

# libraries ---------------------
library('tidyverse')
library(readxl)
library(MASS)
library(fossil)
library(markovchain)
library(depmixS4)
library(mclust)
library(xts)
library(PerformanceAnalytics)
library(hash)
library(ramify)
library(stringr)


# cluster function: compute liklihood which occurs when objects are clustered
clus_lc = function(gij, gii, gjj, ns=2){
  if (ns == 1 ){
    return(0)
  }
 # intracluster correlation 
  cs = 2 * gij + gii + gjj
  if (cs<=ns){
    return(0)
  }
  return(0.5 * (log(ns/cs) + (ns - 1) * log(ns^2 - ns) / log(ns^2 - cs) ) )  
}


# ASCP - rquires a correlation matrix as input: which is converted to a dictionary for convenience
agglo_spc = function(G, cn = NA){
  #G = G_
  #cn = 10
  N = nrow(G)
  # correlation matrix converted to hash table - might change this later
  gdic = hash()
  for (i in 1:N) {
    gdic[[paste(i)]] = hash()
    for (j in 1:N){
      gdic[[paste(i)]][[paste(j)]] = G[i,j]
    }
  }
  
  # tracker hash table
  tracker = hash()
  for (i in 1:N){
    tracker[[paste(i)]] = paste(i)
  }
  
  # cluster size ns is stored in ns_ array
  ns_ = replicate(N, 1)
  
  #''' Create a list of object indices 'other_keys': at every iteration one object
  #   is clustered and removed from the list. It is also removed if no suitable
  #optimal destination is found.'''
  other_keys = paste(1:N)
  
  #''' the operation stocks once there is only one object left to cluster as we
  # need two objects at the very least.'''
  
  while (length(keys(tracker)) != cn){
    #  ''' a random initialization:
    #  pick a object 'node' at random to start clustering,
    #  this might have a consequence on the final result depending on the data.
    #  then loop through the other objects using 'nbor' and costs to store
    #  the likelihood resulting from clustering 'node' to the objects in
    #  'nbor'.
    #  indices: stores the indices which are combinations of node and others.
    #  costs: stores the cost which compute the difference between the likelihood
    #  of the resulting cluster and the sum of the two individual objects forming
    #  the result cluster.
    # '''
    #  ''' the routine uses other_keys and removes elements everytime they are clustered
    # or can't be clustered anymore. If a cluster number is not provided the routine
    #  stops there. If one is then it continues by looking at the elements in the
    #  optimal cluster solution (tracker) and continues merging until the preset
    #  number of clusters is met'''
     
    if (length(other_keys) > 1){
      node = sample(other_keys, 1, replace = T) # ask lionel if sampling with replacement was deliberate
    }else{
      node = sample(keys(tracker), 1, replace = T) # ask lionel if sampling with replacement was deliberate
    }
    nbor = keys(tracker)
    nbor = nbor[nbor!=node]
    costs = replicate(length(nbor), 0)
    indices = list()
    for (i in 1:length(nbor)){
      indices[[i]] = c(node, nbor[i])
    }
    node_lc = clus_lc(0, gdic[[node]][[node]], 0, ns = ns_[as.numeric(node)]) # check why gii and gjj is set to zero
    k = 1
    for( m in indices){
      i = m[1]
      j = m[2]
      costs[k] = clus_lc(gdic[[i]][[j]], gdic[[i]][[i]], gdic[[j]][[j]], ns=ns_[as.numeric(i)]+ns_[as.numeric(j)] ) - (node_lc + clus_lc(0, gdic[[j]][[j]], 0, ns = ns_[as.numeric(j)]))
      k = k + 1  
    }
    
    #''' find the optimal cost which will be the object clustered with node'''
    next_merge = argmax(mat(costs, nrow = 1))
    
    # stopping conditions
    if (costs[next_merge]<=0){
      if(length(other_keys)>1){
        #''' if no cost is positive then this node cannot be clustered further
        #          and must be removed from the list'''
        other_keys = other_keys[other_keys != node]
        next
      }else if (is.na(cn)){
        #''' if no cluster number is provided then the routine has completed
        #          and tracker is the final solution'''
        cn = length(tracker)
        next
      }
    }
    #''' on the other hand, the largest positive cost is the designated
    #    object 'label_b' clustered to node which here is stored as 'label_a'.
    #new clusters 'new_label' take values superior to N.
    #tracker, as previously explained, stores joined strings of the clusters
    #contents'''
    label_a = node
    label_b = indices[[next_merge]][2]
    new_label = paste(length(ns_) + 1)
      
    # '' removes merged elements and update others with the new cluster.
    #only do it when a positive cost is found.'''
    if(costs[next_merge]>0){
      other_keys = keys(tracker)
      other_keys = other_keys[other_keys != label_a]
      other_keys = other_keys[other_keys != label_b]
      other_keys = c(other_keys, new_label)
    }
    
    # ''' Once a cluster is formed, the correlation matrix gdic and tracker need to
    # be updated with the new cluster and the cluster size must be updated with ns_'''
    nbor = nbor[nbor != label_b]
    tracker[[new_label]] = paste(tracker[[label_a]],"_",tracker[[label_b]], sep = "")
    gdic[[new_label]] = hash()
    gdic[[new_label]][[new_label]] = 2*gdic[[label_a]][[label_b]] + gdic[[label_a]][[label_a]] + gdic[[label_b]][[label_b]]
    ns_ = c(ns_, ns_[as.numeric(label_a)] + ns_[as.numeric(label_b)])
    
    for (key in nbor){
      gdic[[new_label]][[key]] = gdic[[label_a]][[key]] + gdic[[label_b]][[key]]
      gdic[[key]][[new_label]] = gdic[[label_a]][[key]] + gdic[[label_b]][[key]]
    }
    
    del(label_a, tracker)
    del(label_b, tracker)
  }
  #''' create the final clustering array:
  #      tracker contains the cluster memberships but as a dictionary
  #we create a numpy array where stocked are labeled with the same number
  #if they belong to the same cluster, and 0 if unclustered'''
  solution = replicate(N, expr = 0)
  k = 1
  for (cluster in keys(tracker)){
    cluster_members = as.numeric(unlist(str_split(tracker[[cluster]], "_")))
    solution[cluster_members] = k
    k = k + 1
  }
  return(solution)
}



#test variables / Identify sectors
data = mat(rnorm(50*500),nrow = 500, ncol = 50)
G_ = cor(t(data))
dim(G_)
#G = matrix(runif(100, min = -1, max = 1), nrow = 100, ncol = 100)
cn = 5
test = agglo_spc(G = t(G_), cn = 5)
test


# identify states on real data
# identify estimate parameters from real data using a hidden markov model
#helper functions -----
naSums = function(x){sum(is.na(x))}
#set.seed(1)

# SnP 500 data clean -------------------------------------
allStocks = read_csv(file ="data/sandp500/all_stocks_5yr.csv")
allStocks = allStocks  %>% dplyr::select(date, close, Name)
flatStocks = allStocks %>% spread(key = Name, value = close, fill = NA ) # explode matrix
survivorStocks =  flatStocks %>% select_if(apply(flatStocks,naSums, MARGIN = 2) == 0) # removed stocks that did not trade in the whole period
# move this to a new document
#log returns
smaller = sample(2:400, size = 100)
GRet = survivorStocks[, c(smaller)]
GRet = diff(as.matrix(log(GRet)))

G_ = cor((Noh$stockMatrix))
dim(G_)
#G = matrix(runif(100, min = -1, max = 1), nrow = 100, ncol = 100)
cn = 2
test = agglo_spc(G = t(G_))
test
adj.rand.index(test,twoStateSequence)
