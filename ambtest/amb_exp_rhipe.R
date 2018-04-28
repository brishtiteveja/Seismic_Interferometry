# con <- file("hadoop_output_log.txt")
# sink(con, append=TRUE)
# sink(con, append=TRUE, type="message")
# 
# # This will echo all input and not truncate 150+ character lines...
# source("amb_exp_rhipe.R", echo=TRUE, max.deparse.length=10000)
# 
# # Restore output to console
# sink() 
# sink(type="message")


library(RSEIS)
library(RPMG)
library(Rwave)
library(datadr)
library(signal)
library(pracma)
library(datadr)

# plot in a pdf file
pdf()
# set the margin of the plot
par(mar=c(1,1,1,1))

# get the file names from the directory in the server
fnames = list.files(path='./data/',
                    pattern=NULL, full.names=TRUE)

# get the number of the stations
num_station = length(fnames)

Amp <- list()

# create the user directory and move there
user <- 'eergun/'
user_dir <- paste('/user/', user, sep="")
hdfs.setwd(user_dir)

# Create the project directory 
dir <- 'interferometry_small'
project_dir <- paste(user_dir, dir, sep="")
if (!rhexists(project_dir))
  rhmkdir(project_dir)

# change to the project directory
hdfs.setwd(project_dir)

data_dir <- "data"
if (!rhexists(data_dir))
   rhmkdir(data_dir)

# change to the data directory
hdfs.setwd(data_dir)

# remove all data
rhdel("*")

# Read the SAC formatted seismic data for each station from the data directory
# and write them as key value pair file on hadoop
n = 2000
# Divide each station into 2 files because otherwise Hadoop 64 MB block size exceeds 
D <- NULL # or 16, needs to divide N-1


for (i in 1:num_station) {
  fn = fnames[i]
  st_num <- as.character(paste('100', i, sep="")) 
  station <- read1sac(fn , Iendian = 1 , HEADONLY=FALSE, BIGLONG=FALSE)
  
  N_orig <- length(station$amp)
  
  nn <- as.integer(N_orig/n)
  
  N <- nn * n - 1 
  N_del <- N_orig - N 
  
  station$station <- rep(st_num, N_orig)
  stationDF <- data.frame(station=station['station'], amp=station['amp'])
  stationDF <- stationDF[1:N, ]
  
  start_idx <- 1 
  
  k <- 10
  num_rows <- (N + 1) / k
  end_idx <- as.integer(num_rows)
 
  j <- 1 
  while(1) {
    if (end_idx > N) {
      end_idx = N
    }
    jN <- length(start_idx:end_idx)
    stationKV <- list(list(st_num, stationDF[start_idx:end_idx, ])) 
   
    stationDataFile <- paste(i, "_", j, "station", "_Data", sep="")
    j <- j + 1
    print(paste("Writing ", stationDataFile, " in HDFS."))
    if (!rhexists(stationDataFile))
      rhwrite(stationKV, file=stationDataFile)
   
    if (end_idx == N)
      break
    
    start_idx <- end_idx + 1 
    end_idx <- start_idx + num_rows - 1
  }
}  

hdfs_dir <- paste(user_dir, dir,"/", data_dir, sep="")
seismHDFSconn <- hdfsConn(hdfs_dir, autoYes = TRUE)

datakvDdf <- ddf(seismHDFSconn) 
datakvDdf  <- updateAttributes(datakvDdf)

byamp <- divide(datakvDdf, 
                by ="station",
                spill = n,
                #output = hdfsConn(dirs, autoYes=TRUE),
                update=TRUE)

dt = 0.004
fs = 1/dt
b1= butter(2, c(0.01/(fs/2), 3/(fs/2)))
signbit <- function(data){
  
  for (i in seq(1,length(data))){
    if (data[i] < 0) {
      data[i] = -1
    } else if (data[i] > 0) {
      data[i] = 1
    } else
      data[i] =0
  }  
  return(data)
}


time=(0:(n-1))*dt
l = n/2
proccc <- addTransform(byamp, function(v) {
  #print(length(v$amp))
  a = v$amp -  mean(v$amp)
  a = detrend(a)
  a = filtfilt(b1, a, type="pass")
  b = signbit(a)
  
  au_sta_22  = acf(b,lag.max = l-1, type = c("correlation"))
  # #print(length(au_sta_22$acf))
  vcrit = sqrt(2)*erfinv(0.99)
  lconf = -vcrit/sqrt(n);
  upconf = vcrit/sqrt(n);
  ind_22 = (au_sta_22$acf >=lconf & au_sta_22$acf <= upconf)
  au_sta_22$acf[ind_22=="TRUE"] = 0
  
  fit.loess22 <- loess(au_sta_22$acf ~ time[1:l], span=0.15, degree=2)
  predict.loess22 <- predict(fit.loess22, time[1:l], se=TRUE)
  a_22 <- ts(au_sta_22$acf, frequency = fs) # tell R the sampling frequency
  a_22_spec <- spec.pgram(a_22, demean = FALSE, detrend = TRUE,plot = TRUE)
  s_22 <- ts(predict.loess22$fit, frequency = fs) # tell R the sampling frequency
  s_22_spec <- spec.pgram(s_22, demean = FALSE, detrend = TRUE,plot = TRUE)
  # spectral whitening can be done dividing the power spectrum of autocorrelated data to smoothed data . add a little damping to the denominator
  wh_sta_22 = a_22_spec$spec / (s_22_spec$spec + 0.00001)
  wh_sta_22_time = abs(ifft((wh_sta_22)))
  b2= butter(2, c(10/(fs/2), 20/(fs/2)))
  result_station_22 <- filtfilt(b2, wh_sta_22_time, type="pass")
})
last = recombine(proccc, combRbind)

hdfs.setwd(project_dir)

# Now each station data will be saved as chunk
m <- n/4
for (i in 1:num_station) {
  fn = fnames[i]
  st_num <- as.character(paste('100', i, sep="")) 
  station = subset(last, last$station == st_num)
  
  NN <- length(station$val)
  stationDF <- data.frame(station = station$station, amp = station$val)
 
  output_dir <- paste('Station_', i, '_Output', sep="")
  if (rhexists(output_dir))
    rhdel(output_dir)  
  rhmkdir(output_dir)
  
  hdfs.setwd(output_dir)
  
  fD <- as.integer(NN / m)
  
  start_idx <- 1 
  end_idx <- as.integer(m)
  for (j in 1:fD) {
    if (j == fD)
      end_idx <- NN
    
    jN <- length(start_idx:end_idx)
    stationKV <- list(list(key=st_num, value=stationDF[start_idx:end_idx, ])) 
    
    stationDataFile <- paste(j, "_station_output_", "Data", sep="")
    print(paste("Writing ", stationDataFile, " output file in HDFS."))
    rhwrite(stationKV, file=stationDataFile)
    
    start_idx <- end_idx + 1 
    end_idx <- start_idx + as.integer(m) - 1
    
  }
  
  hdfs.setwd(project_dir)
}  

st_sum = list()

for (i in 1:num_station) { 
  print(paste("Processing station output", i))
  tryCatch({
    output_dir <- paste(user_dir, dir, "/", 'Station_', i, '_Output', sep="")
    #print(rhls(output_dir))
    outHDFSconn <- hdfsConn(output_dir, autoYes = TRUE)
    outputkvDdf <- ddf(outHDFSconn)
    outputkvDdf <- updateAttributes(outputkvDdf)
    
    byStation <- divide(outputkvDdf,
                        by = 'station',
                        spill = m
                       )
    byStationModified <- ddf(byStation[1:fD])
    sumReduce <- addTransform(byStationModified, function(v) {
          v$amp 
    })
    
    st_sum[[i]] = recombine(sumReduce, combMean)
  },
  error = function(err){
     print(err)
  })
}
hdfs.setwd(project_dir)

# # replace NA with 0
for (i in 1:num_station) {
  st_sum[[i]] <- st_sum[[i]] * fD
}

par(mar=c(4,4,4,4))
#par(mfrow=c(num_stations,1))
time = (0:(n/4 - 1)) * dt
for (i in 1:num_station) {
  t = paste("Summed Amplitude of noise from station", i)
  plot(rev(st_sum[[i]]), time, type='l', col=i, ylab='Time(s)', xlab='Amp', main=t, xlim=c(-max(st_sum[[i]]), max(st_sum[[i]])))
  abline(v=0, lty=2)
}
dev.off()