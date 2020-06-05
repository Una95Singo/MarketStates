# Stylised facts
# R script takes a timeseries and produces financial timeseries stylised facts
#Una Singo, SNGUNA003
# 2 June 2020


# libraries ---------------------
library('tidyverse')
library(readxl)
library(MASS)
library(fossil)
library(timeSeries)
library(ape)
library(Rcpp)
library(umap)
library(emstreeR)
library(poweRlaw)
set.seed(1)


stylisedFacts.plot = function(tSeries, title = "Historgram of daily returns"){
  
  par(mfrow =  c(2,2))
  #par(mfrow =  c(1,1))
  hist(tSeries, freq = F, main = title, xlab = 'Daily return')
  qqnorm(tSeries)
  qqline(tSeries, col=2, lwd=2, lty=2 )
  
  lower = round(length(tSeries)*0.05)
  upper = round(length(tSeries)*0.95)
  
  lowerTail = abs(tSeries[1:lower])
  upperTail = abs(tSeries[upper:length(tSeries)])
  
 
  # Left tail plot
  m_pl = conpl$new(lowerTail)
  est = estimate_xmin(m_pl)
  m_pl$setXmin(est)
  plot(m_pl, main = 'Left tail power law CDF')
  lines(m_pl, col = 2)
  
  # Right tail plot
  m_pl = conpl$new(upperTail)
  est = estimate_xmin(m_pl)
  m_pl$setXmin(est)
  plot(m_pl, main = 'Right tail power law CDF')
  lines(m_pl, col = 2) 
}
