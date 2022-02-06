#' LoGo
#' computes sparse inverse covariance, J, from a clique tree made of cliques 
#' and separators
#' S is the complete covariance matrix
#' separators is the list of separators 
#' clique is the list of cliques
#' separators and cliques can be the outputs of TMFG.m function 
#' see: Wolfram Barfuss, Guido Previde Massara, Tiziana Di Matteo, 
#' and Tomaso Aste  "Parsimonious modeling with information filtering 
#' networks." Physical Review E 94.6 (2016): 062306.
#' TA May 2017
#' 
#' Una Singo 
#' # Una Singo
# 11 July 2020
#'
#' @param S dense pxp matrix with positive weights (covariances)
#' @param separators 
#' @param cliques
#' @return J: A sparse inverse covariance matrix


source('R/03_TMFG.R')
LoGo.solve = function(S, separators = NULL , cliques = NULL){
  
  if ( is.null(separators) || is.null(cliques) ){
    tmfg = TMFG.filter(S)
    separators = tmfg$separators
    cliques = tmfg$cliques
  }
  
  N = nrow(S)
  JLoGo = matrix(0, nrow = N, ncol = N)
  for(i in 1:nrow(cliques))
  {
    v = cliques[i,]
    JLoGo[v,v] = JLoGo[v,v]+solve(S[v,v])
  }
  
  for(i in 1:nrow(separators))
  {
    v = separators[i,]
    JLoGo[v,v] = JLoGo[v,v]-solve(S[v,v])
  }
  return(JLoGo)
}


