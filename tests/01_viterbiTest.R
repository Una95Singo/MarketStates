# Viterbi algorithm test script
# una singo, SNGUNA003
# 9 January 2020
# Test cases for the Viterbi algorithm

# source functions
source('R/01_viterbi.R')

# libraries
library('tidyverse')
library(readxl)
library(fossil)
library(xtable)


# Test set 1 ----------------------
# data import
InputData = read_excel("data/Test/viterbiDataset.xlsx",sheet = 'Test1')
DistanceMatrix = InputData[,2:5]

OptimalViterbi = viterbi(D = DistanceMatrix, K = 4, gamma = 0, debug = F )
viterbi(D = DistanceMatrix, K = 4, gamma = 0, debug = F )
InputData$`True path`
adj.rand.index(OptimalViterbi$Final_Path, InputData$`True path`)
#OptimalViterbi$Final_Path
#InputData$`True path`


# Test set 2 ----------------------
# data import

H.dist = list('A' = 0.2, 'C' = 0.3, 'G' = 0.3, 'T' = 0.2)
L.dist = list('A' = 0.3, 'C' = 0.2, 'G' = 0.2, 'T' = 0.3)

start = c(0.5, 0.5)

transition = matrix(c(0.5, 0.5, 0.4, 0.6), ncol = 2, nrow = 2,  byrow = T )

start = log2(start)
H.dist = lapply(H.dist, log2)
L.dist = lapply(L.dist, log2)
transition = log2(transition)

obs = c('G', 'G', 'C', 'A', 'C', 'T', 'G', 'A', 'A')




DistanceMatrix = matrix (NA, 2, 9)

DistanceMatrix[1,1] = -1+H.dist$G
DistanceMatrix[2,1] = -1+L.dist$G

temp = c(NA, NA)

for (i in 2:9){
  for (j in 1:2){
    if(j == 1){
      DistanceMatrix[j, i] = H.dist[[obs[i]]] + max(DistanceMatrix[1, i-1] + transition[1,1], DistanceMatrix[2, i-1] + transition[2,1])
    } else{
      DistanceMatrix[j, i] = L.dist[[obs[i]]] + max(DistanceMatrix[1, i-1] + transition[1,2], DistanceMatrix[2, i-1] + transition[2,2])
    }
  }
}

#xtable(DistanceMatrix)

OptimalViterbi = viterbi(D = t(DistanceMatrix), K = 2, gamma = 0, debug = F )

truePath = c(2,2,2,1,1,1,1,1) # taken from toy example

adj.rand.index(OptimalViterbi$Final_Path, truePath)


#Viterbi function operates as needed!

