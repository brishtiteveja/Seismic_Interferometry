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
# set the margin of the plot
par(mar=c(1,1,1,1))

data_file_name = './data2/'
dir <- 'data2' # project directory

# get the file names from the directory in the server
fnames = list.files(path=data_file_name,
                    pattern=NULL, full.names=TRUE)

# get the number of the stations
num_station = length(fnames)

Amp <- list()

# create the user directory and move there
user <- 'eergun/'
user_dir <- paste('/user/', user, sep="")
hdfs.setwd(user_dir)

# Create the project directory 
project_dir <- paste(user_dir, dir, sep="")
rhdel(project_dir)
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
n = 8000
# Divide each station into 2 files because otherwise Hadoop 64 MB block size exceeds 
D <- NULL # or 16, needs to divide N-1


for (i in 1:num_station) {
  fn = fnames[i]
  st_num <- as.character(paste('100', i, sep="")) 
  station <- read1sac(fn , Iendian = 1 , HEADONLY=FALSE, BIGLONG=FALSE)
  
  N_orig <- length(station$amp)
  
  nn <- as.integer(N_orig/n)
  
  N <- N_orig - 2 #nn * n - 1 
  #N_del <- N_orig - N 
  
  station$station <- rep(st_num, N_orig)
  stationDF <- data.frame(station=station['station'], amp=station['amp'])
  stationDF <- stationDF[1:N, ]
  
  start_idx <- 1 
  
  k <- 10
  num_rows <- 2*n
  end_idx <- as.integer(num_rows)
 
  j <- 1 
  breakout = FALSE
  while(1) {
    print(start_idx)
    print(end_idx)
    print("")
    if (end_idx > N) {
      end_idx = N
      breakout = TRUE
    }
    jN <- length(start_idx:end_idx)
    stationKV <- list(list(st_num, stationDF[start_idx:end_idx, ])) 
   
    stationDataFile <- paste(i, "_", j, "station", "_Data", sep="")
    j <- j + 1
    print(paste("Writing ", stationDataFile, " in HDFS."))
    if (!rhexists(stationDataFile))
      rhwrite(stationKV, file=stationDataFile)
   
    if (breakout)
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
  v$amp[is.na(v$amp)] <- 0
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
  b2= butter(2, c(6/(fs/2), 12/(fs/2)))
  result_station_22 <- filtfilt(b2, wh_sta_22_time, type="pass")
})
last = recombine(proccc, combRbind)

hdfs.setwd(project_dir)

cor_threshhold <- 0.90
percent_subset_match <- 0.50

station = list()
station_m = list()
station_corm = list()
selected_subsets_per_station = list()
n_selected_subsets_per_station = c()

for (k in 1:num_station) {
  st_number = as.numeric(paste("100", k, sep=""))
  station[[k]] = subset(last, last$station == st_number)
  
  m = n/4
  v = station[[k]]$val
  station_m[[k]] = matrix(v, nrow=m, byrow=FALSE)
  station_m[[k]] = station_m[[k]][1:(m/2 - 1), ]
  station_corm[[k]] = cor(station_m[[k]])
  #print(station_corm[[i]])
  # rgb.palette <- colorRampPalette(c("blue", "yellow"), space = "rgb")
  # # levelplot(station_corm[[i]], 
  # #           main=paste("Station ", i, sep=""), 
  # #           xlab="", ylab="", 
  # #           col.regions=rgb.palette(120), 
  # #           cuts=100, at=seq(0,1,0.01))
  
  selected_stations = list()
  n_selected_stations = c()
  for(i in 1:nn) {
    selecteds = list()
    for(j in 1:nn) {
      if (i==j)
        next 
      if(station_corm[[k]][i,j] > cor_threshhold) 
        selecteds = c(selecteds, j)
    }
    selected_stations[[i]] = selecteds
    n_selected_stations = c(n_selected_stations, length(selecteds))
  }
  
  selected_subsets = rep(0, nn)
  ns <- length(n_selected_stations)
  for (i in 1:ns) {
    c = n_selected_stations[i]
    if (c > as.integer(percent_subset_match * nn)) {
      #print(c)
      selected_subsets[[i]]  = 1
    }
  }
  selected_subsets_per_station[[k]] = selected_subsets
  SS <- selected_subsets_per_station[[k]]
  nSS <- length(SS[SS == 1])
  n_selected_subsets_per_station = c(n_selected_subsets_per_station, nSS)
}

# 
# 
# for (i in 1:num_station) {
#   t = paste("Amplitude of noise from station", i)
#   plot(station[[i]]$val[1:l], type='l', col=i, ylab='Amp', main=t)
# }
# 


m = n/4
mm = (m/2)-1
st_sum = list()  
for(j in 1:num_station) {
  st = station[[j]]
  cur_st_sum = rep(0, m)
  nV <- length(st$val)
  for (i in 1:nV) {
    stv <- st$val[i]
    #print(i)
    subset_idx = as.integer((i-1)/m) + 1
    #print(subset_idx)
    if (selected_subsets_per_station[[j]][subset_idx] == 1) {
      idx = i %% m + 1
      cur_st_sum[[idx]] = cur_st_sum[[idx]] + stv
    }
  }
  st_sum[[j]] = cur_st_sum
}

pdf()
#par(mar=c(4,4,4,4))
#par(mfrow=c(num_stations,1))
time = (0:(n/4 - 1)) * dt
for (i in 1:num_station) {
  t = paste("Summed Amplitude of noise from station", i)
  plot(rev(st_sum[[i]][1:mm]), time, type='l', col=i, ylab='Time(s)', xlab='Amp', main=t, xlim=c(-max(st_sum[[i]]), max(st_sum[[i]])))
  abline(v=0, lty=2)
}
dev.off()

# # Now each station data will be saved as chunk
# m <- n/4
# for (i in 1:num_station) {
#   fn = fnames[i]
#   st_num <- as.character(paste('100', i, sep="")) 
#   station = subset(last, last$station == st_num)
#   
#   NN <- length(station$val)
#   stationDF <- data.frame(station = station$station, amp = station$val)
#  
#   output_dir <- paste('Station_', i, '_Output', sep="")
#   if (rhexists(output_dir))
#     rhdel(output_dir)  
#   rhmkdir(output_dir)
#   
#   hdfs.setwd(output_dir)
#   
#   fD <- as.integer(NN / m)
#   
#   start_idx <- 1 
#   end_idx <- as.integer(m)
#   for (j in 1:fD) {
#     if (j == fD)
#       end_idx <- NN
#     
#     jN <- length(start_idx:end_idx)
#     stationKV <- list(list(key=st_num, value=stationDF[start_idx:end_idx, ])) 
#     
#     stationDataFile <- paste(j, "_station_output_", "Data", sep="")
#     print(paste("Writing ", stationDataFile, " output file in HDFS."))
#     rhwrite(stationKV, file=stationDataFile)
#     
#     start_idx <- end_idx + 1 
#     end_idx <- start_idx + as.integer(m) - 1
#     
#   }
#   
#   hdfs.setwd(project_dir)
# }  
# 
# st_sum = list()
# 
# for (i in 1:num_station) { 
#   print(paste("Processing station output", i))
#   tryCatch({
#     output_dir <- paste(user_dir, dir, "/", 'Station_', i, '_Output', sep="")
#     #print(rhls(output_dir))
#     outHDFSconn <- hdfsConn(output_dir, autoYes = TRUE)
#     outputkvDdf <- ddf(outHDFSconn)
#     outputkvDdf <- updateAttributes(outputkvDdf)
#     
#     byStation <- divide(outputkvDdf,
#                         by = 'station',
#                         spill = m
#                        )
#     byStationModified <- ddf(byStation[1:fD])
#     sumReduce <- addTransform(byStationModified, function(v) {
#           v$amp 
#     })
#     
#     st_sum[[i]] = recombine(sumReduce, combMean)
#   },
#   error = function(err){
#      print(err)
#   })
# }
# hdfs.setwd(project_dir)
# 
# 
# for (i in 1:num_station) {
#   st_sum[[i]] <- st_sum[[i]] * fD
# }
# 
# pdf()
# par(mar=c(4,4,4,4))
# #par(mfrow=c(num_stations,1))
# time = (0:(n/4 - 1)) * dt
# for (i in 1:num_station) {
#   t = paste("Summed Amplitude of noise from station", i)
#   plot(rev(st_sum[[i]]), time, type='l', col=i, ylab='Time(s)', xlab='Amp', main=t, xlim=c(-max(st_sum[[i]]), max(st_sum[[i]])))
#   abline(v=0, lty=2)
# }
# dev.off()