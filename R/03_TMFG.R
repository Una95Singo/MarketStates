#' TMFG
#'  Computes a Planar Maximally Filtered Graph (PMFG) starting from 
#'  a tetrahedron and inserting recursively vertices inside 
#'  existing triangles (T2 move) in order to approximate a
#'  maxiaml planar graph with the largest total weight - non negative weights
#' 
# Una Singo
# 10 July 2020
#'
#' @param W: dense pxp matrix with positive weights (covariances)
#' @return A: A sparse matrix, P, a filtered version of W fulfilling the planarity constraint
#' @return Cliques: 
#' 
#'% Reference: http://arxiv.org/pdf/1505.02445.pdf
#' Guido Previde Massara,  T. Di Matteo, and Tomaso Aste. "Network Filtering
#' for Big Data: Triangulated Maximally Filtered Graph." 
#'  arXiv preprint arXiv:1505.02445 (2015).

library(NetworkToolbox)
library(igraph)

TMFG.filter = function(W){
  N = nrow(W)
  if (N < 9) { print('W Matrix too small') }
  #if (any(W<0)) { print('W Matrix has negative elements')}
  
  A = matrix(0, nrow = N, ncol = N) # initialise adjacency matrix
  in_v = matrix(0, nrow = N, ncol = 1) # initialise list of inserted vertices
  tri = matrix(0, nrow = 2*N-4, ncol = 3) # initialise list of triangles
  separators = matrix(0, nrow = N-4, ncol = 3) # initialise list of 3-cliques (non face-triangles)
  
  ## find 3 vertices with largest strength
  s = rowSums(W*(W > mean(W)))
  j = sort(s, decreasing = T, index.return = T)$ix
  in_v[1:4,1] = j[1:4] 
  ou_v = setdiff(1:N, in_v) # list of vertices not inserted yet
  
  ## build the tetrahedron with largest strength
  tri[1, ] = in_v[c(1,2,3)]
  tri[2, ] = in_v[c(2,3,4)]
  tri[3, ] = in_v[c(1,2,4)]
  tri[4, ] = in_v[c(1,3,4)]
  A[in_v[1], in_v[2]] = 1 
  A[in_v[1], in_v[3]] = 1 
  A[in_v[1], in_v[4]] = 1 
  A[in_v[2], in_v[3]] = 1 
  A[in_v[2], in_v[4]] = 1 
  A[in_v[3], in_v[4]] = 1 
  
  ## build initial gain table
  gain = matrix(-Inf, N, 2*N-4)
  gain[ou_v, 1] = rowSums(W[ou_v, tri[1,]])
  gain[ou_v, 2] = rowSums(W[ou_v, tri[2,]])
  gain[ou_v, 3] = rowSums(W[ou_v, tri[3,]])
  gain[ou_v, 4] = rowSums(W[ou_v, tri[4,]])
  
  kk = 4 # number of triangles
  # helper variables
  gij = matrix(nrow=1,ncol=ncol(gain))
  v = matrix(nrow=1,ncol=ncol(gain))
  ve = array()
  tr = 0
  
  for (k in 5:N){
    ## find best vertex to add in a triangle
    if (length(ou_v) == 1){ # spevial case for the last vertex
      ve = ou_v
      v = 1
      w = 1
      tr = which.max(gain[ou_v,])
    }
    else{
      for(q in 1:ncol(gain))
      {
        gij[,q] =  max(gain[ou_v,q])
        v[,q] = which.max(gain[ou_v,q])
        tr = which.max(gij)
      }
      ve = ou_v[v[tr]]
      w  = v[tr] # use w instead of v
    }
    
    #update vertex lists / cache
    ou_v = ou_v[-w]
    in_v[k] = ve
    
    # update adjacency matrix
    A[ve, tri[tr,]] = 1
    
    # update 3-clique list
    separators[k-4, ] = tri[tr, ] 
    
    # update triangle list replacing 1 and adding 2 triangles 
    tri[kk+1, ] = c(tri[tr,c(1,3)],ve) # add
    tri[kk+2, ] = c(tri[tr,c(2,3)],ve) # add
    tri[tr, ]   = c(tri[tr,c(1,2)],ve) # replace
    
    # update gain table
    if( length(ou_v) == 1) {
      gain[ve, ] = 0 
      gain[ou_v, tr]  = sum(W[ou_v, tri[tr, ] ])
      gain[ou_v,kk+1] = sum(W[ou_v, tri[kk+1, ] ]) # sum(W(ou_v,tri(kk+1,:)),2);
      gain[ou_v,kk+2] = sum(W[ou_v, tri[kk+2, ] ]) #sum(W(ou_v,tri(kk+2,:)),2);
      
    }else{
      gain[ve, ] = 0 
      gain[ou_v, tr]  = rowSums(W[ou_v, tri[tr, ] ])
      gain[ou_v,kk+1] = rowSums(W[ou_v, tri[kk+1, ] ]) # sum(W(ou_v,tri(kk+1,:)),2);
      gain[ou_v,kk+2] = rowSums(W[ou_v, tri[kk+2, ] ]) #sum(W(ou_v,tri(kk+2,:)),2);
    }  
    # update number of triangles
    kk = kk+2 
    #% if mod(k,1000)==0,fprintf('TMFG - T2 only: %0.2f per-cent done\n',k/N*100);end
  }
  
  A = W * ( (A + t(A)) == 1)
  
  cliques = rbind(in_v[1:4],(cbind(separators,in_v[5:ncol(W)])))

  return(list('A' = A , 'cliques' = cliques, 'separators' = separators))
}





# test____
#n = 50
#W = matrix(runif(n*n, min = 0.01, max = 1), nrow = n)


#una = TMFG.filter(W)

#una$separators



