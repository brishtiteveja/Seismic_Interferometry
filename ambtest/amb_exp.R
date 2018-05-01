library(RSEIS)
library(RPMG)
library(Rwave)
library(datadr)
library(signal)
library(pracma)

fnames = list.files(path='./data2',
                    pattern=NULL, full.names=TRUE)
num_station = length(fnames)

stations = list()
dt = 0.004
fs=1/dt;
n = 8000

for (i in 1:num_station) {
  fn = fnames[i]
  station =  read1sac(fn , Iendian = 1 , HEADONLY=FALSE, BIGLONG=FALSE)
  station = data.frame(station = station$HEAD$kstnm, amp = station$amp)
  N_orig <- length(station$amp)
  nn <- as.integer(N_orig/n)
  NN <- nn * n # Taking a multiple of n
  station = station[1:NN,]
  #station <- station[-nrow(station),]
  station <- station[-nrow(station),]
  #station$HEAD.time = (1:length(station$HEAD.npts)-1)*dt
  station <- station[,c('station', 'amp')]
  stations[[i]] <- station
} 


dataKV <- list()

for (i in 1:num_station) {
  st_name <- paste('station', i, sep="")
  st_id <- as.numeric(paste('100', i, sep=""))
  dataKV[[i]] <- kvPair(st_name, subset(stations[[i]], station == st_id))
}



datakvDdf <- ddf(dataKV)

byamp <- divide(datakvDdf, 
                by ="station",
                spill = n,
                #output = hdfsConn(dirs, autoYes=TRUE),
                update=TRUE)

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

cor_threshhold <- 0.60
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
  nn <- nrow(station_corm[[k]]) # number of subsets in each station
  #print(station_corm[[i]])
  # rgb.palette <- colorRampPalette(c("blue", "yellow"), space = "rgb")
  # # levelplot(station_corm[[i]], 
  # #           main=paste("Station ", i, sep=""), 
  # #           xlab="", ylab="", 
  # #           col.regions=rgb.palette(120), 
  # #           cuts=100, at=seq(0,1,0.01))
  
  selected_stations = list() # selected subsets for each station
  n_selected_stations = c()
  for(i in 1:nn) {
    selecteds = c()
    for(j in 1:nn) {
      if (i==j)
        next 
      if(abs(station_corm[[k]][i,j]) > cor_threshhold) 
        selecteds = c(selecteds, j)
    }
    selected_stations[[i]] = selecteds #selected for station i
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
  selected_subsets_per_station[[k]] <- selected_subsets
  SS <- selected_subsets_per_station[[k]]
  nSS <- length(SS[SS == 1])
  n_selected_subsets_per_station <- c(n_selected_subsets_per_station, nSS)
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
  plot(rev(st_sum[[i]]), time, type='l', col=i, ylab='Time(s)', xlab='Amp', main=t, xlim=c(-max(st_sum[[i]]), max(st_sum[[i]])))
  abline(v=0, lty=2)
}
dev.off()


########
#### Simple plotting with taking correlation of the subsets for each station
# station = list()
# 
# for (i in 1:num_station) {
#   st_number = as.numeric(paste("100", i, sep=""))
#   station[[i]] = subset(last, last$station == st_number)
# }
# 
# m = n/4
# st_sum = list()
# j = 1
# for(st in station) {
#   cur_st_sum = rep(0, m)
#   i = 1
#   for (stv in st$val) {
#     idx = i %% m + 1
#     cur_st_sum[[idx]] = cur_st_sum[[idx]] + stv
#     #print(idx)
#     #print("")
#     i = i+1
#   }
#   st_sum[[j]] = cur_st_sum
#   j = j+1
# }
# 
# pdf()
# #par(mar=c(4,4,4,4))
# #par(mfrow=c(num_stations,1))
# time = (0:(n/4 - 1)) * dt
# for (i in 1:num_station) {
#   t = paste("Summed Amplitude of noise from station", i)
#   plot(rev(st_sum[[i]]), time, type='l', col=i, ylab='Time(s)', xlab='Amp', main=t, xlim=c(-max(st_sum[[i]]), max(st_sum[[i]])))
#   abline(v=0, lty=2)
# }
# dev.off()