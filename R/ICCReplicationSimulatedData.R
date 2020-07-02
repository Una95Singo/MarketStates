# ICC replication on simulated datasets
# Una Singo, SNGUNA003
# 5 June 2020


# import R functions ------------
source("R/plotter.R")
source("R/simulate.R")
source("R/segmentation.R")
source("R/simulate.R")
source("R/stylisedFacts.R")

# libraries ---------------------
library('tidyverse')
library(readxl)
library(MASS)
library(fossil)
library(markovchain)
library(igraph)

set.seed(1)

# using markovchain package, simulate 1000 datasets per state
# Simulate date--------------


# 2 state simulation
marketStates = c("1", "2")
byRow = TRUE
tMatrix = matrix(data = c(0.95, 0.05,
                           0.05, 0.95), byrow = byRow, nrow = 2, dimnames = list(marketStates, marketStates))
twoState =  new("markovchain", states = marketStates, byrow = byRow, transitionMatrix = tMatrix, name = "twoStateMarket")


show(twoState)
#mcDF =  as(twoState, 'igraph')
#mcDF
twoStateSequence = as(rmarkovchain(n = 2000, object = twoState, t0 = '1'), 'numeric')
plotter.states(twoStateSequence, title = 'simulated states', thickness = 4)

