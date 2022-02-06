# russell
library(tidyverse)


data = read.csv("~/Downloads/RIY2_sheet5.csv", sep = ';' )

#data = data[1:5,1:5]
data
data[data == "#N/A"] = NA
(data[,-1])
t_data = apply(data[,-1], 2 , as.numeric)

n_rows = dim(t_data)[1]

survive = t_data[,colSums(is.na(t_data))<200]

gRet = diff(log(survive))

market = rowMeans(gRet, na.rm = T)
length(market)

market = market[!(is.na(market))]


plot(x = 1:length(market), y = market, main = "market return", type ='l' )

plot(x = 1:length(market), y = cumsum(market), main = "market return", type ='l' )

