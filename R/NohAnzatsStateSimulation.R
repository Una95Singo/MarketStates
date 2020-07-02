# Noh Anzats State Simulation
# Unarine Singo
# 8 June 2020

NohAnzats.simulate = function (stateSequence ='', g = c(0.1, 0.1), clusters = 2, stocks = 100, gaussian = T){
  # step 1 Cluster numbers --------------------------------------------------------
  stateSequence = twoStateSequence
  spinLabels = stateSequence
  C = clusters  # two clusters of stocks
  s = max(table(twoStateSequence)) # clusters size is 50 each
  N = length(spinLabels) # number of days
  D = stocks  # number of stocks
  # step 2 simulate cluster effects --------------------------------------------------------
 
  eta = matrix( NA, C, D, byrow = T) 
  if (gaussian){
    for (i in 1:C){
      eta[i, ] = rnorm (D, 0, 1)
    }  
  } else{
    # cluster effects use t-distribution with 3 degrees of freedom
    for (i in 1:C){
      eta[i, ] = rt (D, df = 3)
    }
  }
  
  # step 3 simulate random effects ---------------------------------------------------------
  eps = matrix(NA, N, D) # specific random effects
  # expsilon
  for (n in 1:N) {
    eps[n, ] = rnorm(D, 0,1)
  }
  
  # test market return
  #par(mfrow =  c(1,1))
  #plot(cumsum(colMeans(t(eps))), type='l', main ='Cummulative plot of epsilon', ylab= 'price')
  # step 4 fix Intra-cluster binding strength --------------------------------------------------------
  gs = g # same cluster binding strength
  
  # step 5 compute daily returns --------------------------------------------------------
  xi = matrix(NA, N, D) # matrix to store daily returns

  # double for loop for simualation
  for (n in 1:N){
    for (d in 1:D){
      cluster = spinLabels[n]
      xi[n, d] = (sqrt(gs[cluster]) * eta[cluster, d] + eps[n, d]) / (sqrt(1+gs[cluster]))
    }
  }
  
  return(list( 'stockMatrix' = t(xi), 'stateSeq' = spinLabels, 'eta' = eta, 'eps' = eps, 'intraBinding' =gs))
}
