temp <- read.csv('~/Downloads/reclimatedata/tmean_yr_Prism_1900_full.csv')
precip <- read.csv('~/Downloads/reclimatedata/pr_monthly_Prism_1900_full.csv')
biomass_dat_est <- read.csv("/Users/paleolab/babySTEPPS/Data/biomass_prediction_v0.9-7_bam.csv")

head(temp)
head(biomass_dat_est)

temp <- as.data.frame(temp)
precip <- as.data.frame(precip)
biomass_dat_est <- as.data.frame(biomass_dat_est)

matches <- which(outer(temp$x,biomass_dat_est$x,'==')&
        outer(temp$y,biomass_dat_est$y,'=='),arr.ind = TRUE)

#temp[473,c('x','y')]
#biomass_dat_est[476,c('x','y')]

plot(temp[matches[,1],'Mean'],biomass_dat_est[matches[,2],'Total'],
     pch = 19, cex = .3,xlab = 'Mean Temperature',ylab = 'Total Biomass')

matches.precip <- which(outer(precip$x,biomass_dat_est$x,'==')&
                   outer(precip$y,biomass_dat_est$y,'=='),arr.ind = TRUE)
plot(precip[matches.precip[,1],'total'],biomass_dat_est[matches.precip[,2],'Total'],
     pch = 19, cex = .3,xlab = 'Total Precip',ylab = 'Total Biomass')


convert.lat.lon <- function(x,y){
  lat.long.reg <- cbind(as.numeric(as.character(x)),as.numeric(as.character(y)))
  lat.long.reg.df = data.frame(lat.long.reg)
  colnames(lat.long.reg.df) = c('x', 'y')
  
  coordinates(lat.long.reg.df) <- ~ x + y
  proj4string(lat.long.reg.df) <- CRS('+init=epsg:3175')#
  
  albers <- spTransform(lat.long.reg.df, CRS('+proj=longlat +ellps=WGS84'))
  albers <- as.matrix(data.frame(albers))
  return(albers)
}

breaks <-  seq(2,15,2)
colors <- rev(rainbow(length(breaks)-1,alpha = 1))
data_binned <-  cut(temp$Mean, c(breaks), include.lowest = FALSE, labels = FALSE)

bluefunc <- colorRampPalette(rev(c('red','orange','yellow','green','blue','purple')))
bluefuncs <- bluefunc(length(breaks))

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(albers[,1],albers[,2], pch=19,
       cex=.3,lwd=.2,col = bluefuncs[data_binned])
legend('topright',legend = round(breaks),col = bluefuncs, pch = rep(19,length(breaks)))
title('Average Temp in 1900')
plot(temp[matches[,1],'Mean'],biomass_dat_est[matches[,2],'Total'],
     pch = 19, cex = .3,xlab = 'Mean Temperature',ylab = 'Total Biomass')

breaks <-  seq(500,1141,50)
colors <- rev(rainbow(length(breaks)-1,alpha = 1))
data_binned <-  cut(precip$total, c(breaks), include.lowest = FALSE, labels = FALSE)

bluefunc <- colorRampPalette(rev(c('red','orange','yellow','green','blue','purple')),alpha=1)
bluefuncs <- bluefunc(length(breaks))

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(convert.lat.lon(precip$x,precip$y)[,1],
       convert.lat.lon(precip$x,precip$y)[,2], pch=19,
       cex=.5,lwd=.2,col = colors[data_binned])
legend('topright',legend = round(breaks),col = bluefuncs,
       pch = rep(19,length(breaks)))
title('Total Precip in 1900')

map('state', xlim=c(-98,-81), ylim=c(41.5,50))

max.total <- max(biomass_dat_est$Total,na.rm = TRUE)
alpha.use <- biomass_dat_est$Total[matches.precip[,2]]/max.total
alpha.use[alpha.use<.2]<-.2
lat.long.keep<-convert.lat.lon(precip$x,precip$y)
data_binned <-  cut(precip$total[matches.precip[,1]], c(breaks), 
                    include.lowest = FALSE, labels = FALSE)
for(i in 1:nrow(matches.precip)){
  #bluefunc <- colorRampPalette(rev(c('red','orange','yellow','green','blue','purple')),alpha = alpha.use)
  bluefuncs <- rev(rainbow(length(breaks)-1,start=0,end=1,alpha = alpha.use[i]))
  points(lat.long.keep[matches.precip[i,1],1],
         lat.long.keep[matches.precip[i,1],2], pch=19,
         cex=.5,lwd=.2,col = bluefuncs[data_binned[i]])
}

plot(precip[matches.precip[,1],'total'],biomass_dat_est[matches.precip[,2],'Total'],
     pch = 19, cex = .3,xlab = 'Total Precip',ylab = 'Total Biomass')

plot.stop <- cbind(precip[matches.precip[,1],'total'],biomass_dat_est[matches.precip[,2],c('Total','x','y')])
x <- convert.lat.lon(plot.stop$x,plot.stop$y)[,1]
y <- convert.lat.lon(plot.stop$x,plot.stop$y)[,2]
points(plot.stop[y>(45)&y<(46),1],plot.stop[y>(45)&y<(46),'Total'],pch = 19, cex = .5,col='red')

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(x[y>(45)&y<(46)],
       y[y>(45)&y<(46)],pch=19,
       cex=.2,col='gray')







