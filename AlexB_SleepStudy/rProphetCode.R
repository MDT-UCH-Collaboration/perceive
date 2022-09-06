library('prophet')
library(ggplot2)

#df2 <- read.csv('https://raw.githubusercontent.com/facebook/prophet/main/examples/example_yosemite_temps.csv')

fileLIST <- c('SPPD11_L_Prophet.csv','SPPD11_R_Prophet.csv','SPPD12_L_Prophet.csv',
              'SPPD12_R_Prophet.csv','SPPD13_L_Prophet.csv','SPPD13_R_Prophet.csv')
fileSAVE <- c('SPPD11_L_forcast.csv','SPPD11_R_forcast.csv','SPPD12_L_forcast.csv',
              'SPPD12_R_forcast.csv','SPPD13_L_forcast.csv','SPPD13_R_forcast.csv')

for (fi in 1:length(fileLIST)) {
  
  tmpfile <- fileLIST[fi]
  df <- read.csv(tmpfile)
  
  testCP <- as.POSIXct(df$ds[1],format = "%d-%b-%Y %H:%M:%S")
  
  if (is.na(testCP)){
    df$ds <- as.POSIXct(df$ds,format = "%m/%d/%Y %H:%M",tz="MST")
  }else{
    df$ds <- as.POSIXct(df$ds,format = "%d-%b-%Y %H:%M:%S",tz="MST")
  }
  
  m <- prophet(df, changepoint.prior.scale=0.01)
  future <- make_future_dataframe(m, periods = 300, freq = 60 * 60)
  fcst <- predict(m, future)
  fcstdf = as.data.frame(fcst)
  
  
  
  write.csv(fcstdf,file = fileSAVE[fi])
  
  
}



df <- read.csv('SPPD3_L_Prophet.csv')

df$ds <- as.POSIXct(df$ds,format = "%d-%b-%Y %H:%M:%S")

m <- prophet(df, changepoint.prior.scale=0.01)
future <- make_future_dataframe(m, periods = 300, freq = 60 * 60)
fcst <- predict(m, future)
plot(m, fcst)

ggplot(df,aes(x = ds, y = y)) +
  geom_point() + 
  geom_line(data = fcst, aes(x = ds, y = yhat), color = "red")

fcstdf = as.data.frame(fcst)
write.csv(fcstdf,file = 'testFCST.csv')
