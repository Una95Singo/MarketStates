# R script to plot unqiue plots that pertain to the project 
# una singo, SNGUNA003
# 21 January 2020

# 
library(scales)



# takes a series of states and plots
plotter.states = function(state_series, title = 'no title', thickness = 1){
  image(matrix(state_series, ncol=1), main = title)
  
}
# plots probability tranisiton matrix using heat colours
# @param trans_matrix: transition probability matrix
# @param col: number of heat colours to use
plotter.transitions = function(trans_matrix, cols=5, title='no title') {
  image(trans_matrix, col = heat.colors(cols), main = title)
}

# plots series and underlying states
# @param S0: starting price. if not NA plots price path
# @param series: timeseries of returns
# @param states: timeseires of states 
plotter.series = function(series, states, S0 = NA, title =""){
  if (is.na(S0)) {
 
    Min = min(series)
    Max = max(series)
    Time = length(series)
    ## set up the plot region:
    plot(y = c(Min, Max),x = c(0, Time), type = "n",
         main = title, 
         xlab = 'Time', ylab ='Returns')
    
    K = length(unique(states))
    Colours = c("#4C00FFFF", "#00E5FFFF", 'red', 'green' )
    Colours = alpha(Colours, alpha = 0.3)
    for (i in 1:(K)){
      rect(which(states==i)-0.5, Min, which(states==i)+0.5, Max,  col=Colours[i], border =NA )
    } # plot rectangles by state 
    lines(series, col='Black', lwd=1.25)
  }
  else{
    series = S0*c(cumsum(series))+S0 # convert returns to prices
    Min = min(series)
    Max = max(series)
    Time = length(series)
    ## set up the plot region:
    plot(y = c(Min, Max),x = c(0, Time), type = "n",
         main = title, 
         xlab = 'Time', ylab ='Returns')
    K = length(unique(states))
    Colours = c("#4C00FFFF", "#00E5FFFF", 'red', 'green' )
    Colours = alpha(Colours, alpha = 0.3)
    for (i in 1:(K)){
      rect(which(states==i)-0.5, Min, which(states==i)+0.5, Max,  col=Colours[i], border =NA )
    } # plot rectangles by state 
    lines(series, col='Black', lwd=1.25)
  }
}























#'  ------------------------------------------------------------------------------------------
#' depricated functions  ---------------------------------------------------------------------
#' 
#' plotter.states = function(state_series, title = 'no title', thickness = 1){
#' plot(x = 1:length(state_series), y = replicate(length(state_series),1), col = state_series, type ='h', lwd = thickness, main =title, xlab='index', ylab = '', ylim = c(0,1))
#' }
#' 
#' 





