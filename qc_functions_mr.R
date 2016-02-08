# Function that organizes the data by the date recorded in the file name. 
date.organize <- function(x) {
  splice.date <- str_split_fixed(x,'-',2)
  dts <- splice.date[,2]
  dts <- sub('.csv','',dts)
  date <- mdy(dts,tz='America/New_York')
  df <- (cbind(x,date))
  organized <- df[order(df[,2]),]
  ordered <- organized[,1]
  return(ordered)
}

# must be used on already ordered data. 
fix.savings <- function(x) {
  splice.date <- str_split_fixed(x,'-',2)
  dts <- splice.date[,2]
  dts <- sub('.csv','',dts)
  date <- mdy(dts,tz='America/New_York')
  d.times <- c('11/02/2014 02:00:00','03/15/2015 02:00:00')
  d.savings <- mdy_hms(d.times,tz='America/New_York')
  d.f <- data.frame(x,date,stringsAsFactors=F)
  d.f <- d.f[order(d.f$date),]
  df.fix <- d.f[c(1,which(d.f$date < d.savings[1] | d.f$date >  d.savings[2])+1),1]
  return(df.fix)
}



round_minute<-function(x,precision){
  m<-minute(x)+second(x)/60
  m.r<- round(m/precision)*precision
  minute(x)<-m.r
  second(x)<-0
  x
}

jumpr <- function(x,sites){
  icol<- paste(sites,'iter',sep='_')
  icm <- paste(sites,'_cm',sep='')
  j.dat <- as.list(sites)
  names(j.dat) <- sites
  hld.list <- j.dat
  beg <- function(y){ beg <- median(y[1:6],na.rm=T) }
  ends <- function(y){ end <- median(y[seq(length(y)-6,length(y),by=1)],na.rm=T)}
  first <- function(y){fst <- y[1]}
  last <- function(y){lst <- y[length(y)]}
  
  for(i in 1:length(sites)){
    dat <- x[,c('min10',icm[i],icol[i],'index')]
    #dat <- dat[!(is.na(dat[,icol[i]])),]
    dat$fact <- factor(dat[,icol[i]])
    b <- tapply(dat[,icm[i]],dat$fact,beg)
    e <- tapply(dat[,icm[i]],dat$fact,ends)
    firsts <- ymd_hms(tapply(as.character(dat$min10),dat$fact,first),tz=mytimezone)
    lasts <- ymd_hms(tapply(as.character(dat$min10),dat$fact,last),tz=mytimezone)
    firsts <- firsts[-1]
    lasts <- lasts[-1]
    times <- interval(firsts,lasts)
    
    b1 <- c(b,NA)
    e1 <- c(NA,e)
    diff <- (b1-e1)
    diff <- diff[-c(1,length(diff))]
    jump.data <- data.frame(jump=diff,interval=times,
                            min10start = firsts,
                            min10stop = lasts,
                            iteration = as.numeric(names(diff)))
    for(j in 1:nrow(jump.data)){
      naq <- dat[which(dat$min10 > jump.data[j-1,'min10stop']),c('min10',icm[i])]
      if(length(which(is.na(naq[1:60,icm[i]]))) > 59){
          jump.data[j,'jump'] <- 0.0000
        } 
      dat[which(dat$min10 < jump.data[j,'min10start']),icm[i]] <- 
      dat[which(dat$min10 < jump.data[j,'min10start']),icm[i]] + 
      jump.data[j,'jump']
    }
    names(dat) <- c('min10','cm','iter','index','fact')
    dat$site <- sites[i]
    jump.data$site <- sites[i]
    j.dat[[i]] <- jump.data
    hld.list[[i]] <- dat
  }
  return(list(j.dat,hld.list))
}

readr.c <- function(df){
  #Setup column names
  c.col.names1 <- c('Sample','Date.Time','Cond_low','Cond_high','Temp',
                    'Detached','Attached','Connected','Stopped','end',
                    'lubridate')
  c.col.names2 <- c('Sample','Date.Time','Cond_low','Temp','Detached',
                    'Attached','Connected','Stopped','end','lubridate')
  stor <- as.list(df)
  if(length(df)>0){
    for(w in 1:length(df)){
      #Read in data
      print(df[w])
      a <- read.csv(df[w],skip=2,header=F,stringsAsFactors=F)
      #Convert time column to posix for R to read time
      a$lubridate <- mdy_hms(a[,2])
      #Setup column names 
      if(length(colnames(a))>  10) colnames(a) <- c.col.names1 else { colnames(a) <- c.col.names2}
      if(length(colnames(a))==10) a$Cond_high <- 'NoData'
      if(df[w] %in% ord.it$og) a$iteration <- ord.it[which(ord.it$og %in% df[w]),'iter'] else {a$iteration <- NA}
      a$launch <- NA
      a$launch[seq(nrow(a)-10,nrow(a))] <- "Stopped"
      a$launch[1:10] <- "Started"
      stor[[w]] <- a
    }
  }
a <- do.call('rbind',stor)
return(a)
}


drift.c <- function(z){
  for(l in 1:length(unique(z$iteration))){
    if(l < length(unique(z$iteration))){
    dat1 <- tail(z[which(z$iteration == l),])
    dat2 <- head(z[which(z$iteration == l+1),])
    drift <- median(dat2$SC,na.rm=T) - median(dat1$SC,na.rm=T)
    drift.vect <- seq(0,drift,length.out=nrow(z[z$iteration==l,]))
    z[z$iteration ==l, 'SC'] <-  z[z$iteration ==l, 'SC'] + drift.vect
    }
  }
  return(z)
}






#Function that converts onset time stamp into a computer friendly version
chronify <- function(x) {
  # Substitute out 12: for 00: in order to convert to military time
  x <- sub(' 12:',' 00:',x)
  # Split out time stamp from onset output. 
  splice <- as.data.frame(str_split_fixed(x," ",3),stringsAsFactors=F)
  #convert PM to military time  
  indx <- which(splice[,3]== 'PM')
  pms <- splice[indx,]
  tm.split <- str_split_fixed(pms[,2],':',3)
  addtime <- as.numeric(tm.split[,1])+12
  pm.time <- paste(addtime,tm.split[,2],'00',sep=':')
  splice[indx,2] <- pm.time
  chron <- chron(splice[,1],splice[,2], format=c('m/d/y','h:m:s'))
  chron <- trunc(chron, times('00:10:00'))
  return(chron)
}