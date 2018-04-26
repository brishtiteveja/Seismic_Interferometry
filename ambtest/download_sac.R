#url https://service.iris.edu/irisws/timeseries/1/query?net=YW&sta=1001&cha=DPZ&start=2016-07-01T00:00:00&end=2016-07-15T00:00:00&output=miniseed&loc=--

base.url <- "http://service.iris.edu/fdsnws/dataselect/1/"

parameters <- list(
  name = c(
    "net",
    "sta",
    "loc",
    "cha",
    "starttime",
    "endtime"
  ),
  value = c(
    "YW",
    1001,
    '--',
    "DPZ",
    "2016-07-01T00:00:00",
    "2016-07-1T00:00:00"
  ) #End time
)

library(rFDSN)
library(RCurl)
download.file('http://service.iris.edu/fdsnws/dataselect/1/query?net=YW&sta=1001&loc=--&cha=DPZ&starttime=2016-07-01T00:00:00&endtime=2016-07-15T00:00:00', destfile='data.mseed', method='libcurl')

networks <- FDSNGetNetworks(base.url, parameters)

base.url <- "http://service.iris.edu/fdsnws/dataselect/1/"
parameters <- list(
  name = c(
    "sta",
    "cha",
    "start",
    "end"),
  value = c(
    "Y22D", #Station name
    "BDF", #SEED channel name
    "2014-08-09T00:00:00.000", #Midnight GMT
    "2014-08-09T23:59:59.999") #Midnight GMT
)
miniseed.file <- FDSNGetTimeSeries(base.url, parameters, save.file = "Y22D_BDF.mseed")


HADOOP_DIR <- '/user/azehady/seism/'

setwd('./data2')
fnames = list.files(path=getwd(),
                    pattern=NULL, full.names=FALSE )
for (f in fnames) {
  fnc <- paste(getwd(), f, sep="/")
  fn <- paste(HADOOP_DIR, f, sep="")
  rhput(fnc, fn)
}

setwd('..')

