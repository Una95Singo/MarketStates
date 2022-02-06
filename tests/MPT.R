# Libraries
library(openxlsx)
library(Rcpp)
library(timeSeries)

# paths 

fileName = "PT-TAA-JSE-Daily-1994-2017.xlsx"
filePath = "/Users/singo/Library/Mobile Documents/com~apple~CloudDocs/Advanced Portfolio Theory/R/Assignment 1/data"
ffilen = paste(filePath, fileName, sep="")


## 4. load the dataset by sheet
dfS <- list()
for (i in 1:4){
  dfS[[i]] <- read.xlsx(ffilen, sheet = i,detectDates = TRUE)
}
dim(dfS[[3]])
# if you need to convert column ???A??? to date then use
# df$A <- as.POSIXct(df$A,format="%H:%M:%S")



## 5. Keep only the specified list of Tickers and create timeSeries object
Entities = c('X1','STEFI','ALBI','J203','J500',sprintf("J5%d",seq(10,90, by = 10)))
Items = c('Date','TRI','Stefi')
# find Tickers in colnames, and
# TRI at the attribute type and reference and join
for (i in c(1,2,3,4)){
  # logical FALSE vector
  tI0 <- logical(length = length(colnames(dfS[[i]])))
  tI1 <- tI0
  # find the Entities in the data frame
  for (j in 1:length(Entities)){
    tI0 <- tI0 | grepl(Entities[j],colnames(dfS[[i]]))
  }
  # find the Items in the data frame
  for (k in 1:length(Items)){
    tI1 <- tI1 | grepl(Items[[k]],dfS[[i]][2,])
  }
  # combined the logical indices
  tI <- tI0 & tI1
  # remove the columns not required
  dfS[[i]] <- dfS[[i]][,tI]
  # remove the first two rows (as they are not dates)
  dfS[[i]] <- dfS[[i]][-c(1,2),]
  # rename the first column name to Dates
  names(dfS[[i]])[1] <- "Date"
  # clean up the remaining column names
  newColNames <- strsplit(colnames(dfS[[i]]),":")
  for (m in 2:length(newColNames)){
    names(dfS[[i]])[m] <- newColNames[[m]][1]
  }
}


## 6. Clean and convert into a single timeSeries object
# 6.1. Initialise the timeSeries object with first data frame
iN <- 1
tsTAA <- timeSeries(dfS[[iN]][,2:ncol(dfS[[iN]])],as.Date(dfS[[iN]][,1]))
print(dim(tsTAA)) # print dimensions to console
# correct the column names
# 6.3. Concatenate additional timeSeries columns on to the object
for (i in c(2,3,4)){
  # consider iterative merging using inherited zoo properties
  # the first column is the Date column the rest are features we do this
  # to ensure that the dates are correctly aligned when time series are
  merged
  tsTAA <- cbind(tsTAA,timeSeries(dfS[[i]][,2:ncol(dfS[[i]])],as.Date(dfS[[i
                                                                           ]][,1])))
  print(dim(tsTAA))
  print(colnames(tsTAA))
}
# 6.4 Set the units to TRI
setFinCenter(tsTAA) <- "Johannesburg"
# 6.5 Fix the colname errors introduce during the cbind
names(tsTAA)[grep("TS.1.1",names(tsTAA))] <- "ALBI"
names(tsTAA)[grep("TS.1.2",names(tsTAA))] <- "STEFI"
names(tsTAA)[grep("TS.1",names(tsTAA))] <- "ALSI"


## 7. Convert from Daily Sampled Data to Monthly Sampled Data
# 7.1. Decimate the daily data to monthly data
tsTAA <- daily2monthly(tsTAA)
# 7.2 Visualise the data on a single plot
# combine and visualise and prettify the y-axis
plot(tsTAA,plot.type = c("single"),
     format = "auto",
     at=pretty(tsTAA),
     ylab = "Returns",
     main = "TRI for sectors")