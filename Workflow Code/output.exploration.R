library(oce)

#####
##### Trace Plots #####
#####

for(t in 1:100){
plot(samplesList[[1]][,t],typ='l',main=t,ylim=c(0,150),col=rainbow(3,alpha = 1)[1])
points(samplesList[[2]][,t],typ='l',col=rainbow(3,alpha = 0.6)[2])
points(samplesList[[3]][,t],typ='l',col=rainbow(3,alpha = 0.6)[3])

}

for(i in 101){
  plot(samplesList[[1]][,i],typ='l',main='SIGMA',
       ylim=c(0,max(c(samplesList[[1]][,i],samplesList[[2]][,i],samplesList[[3]][,i]))),col=rainbow(3,alpha = 1)[1])
  points(samplesList[[2]][,i],typ='l',col=rainbow(3,alpha = 0.6)[2])
  points(samplesList[[3]][,i],typ='l',col=rainbow(3,alpha = 0.6)[3])
}

#####
##### Time Series #####
#####

plot_biomass_ts <- function(site_number, biomassCI){
  fig.mat <- matrix(1,27,1)
  fig.mat[1:6,]<-1
  fig.mat[7:27,]<-seq(2,22,1)
  #control.pts<-read.csv(text=getURL('https://raw.githubusercontent.com/PalEON-Project/stepps-baconizing/master/data/bacon_age_markers_v1.csv'))
  
  layout(fig.mat)
  par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),tcl=.5)
  
  breaks <-  c(seq(0,50,10),seq(75,200,25))
  colors <- rev(terrain.colors(length(breaks)))
  
  data_binned <-  cut(biomassCI[2,], c(breaks), include.lowest = FALSE, labels = FALSE)
  
  plot(seq(100,10000,100),biomassCI[2,],
       cex=.1,ylim=c(0,150),xlim=c(10000,-10),ylab="Biomass (Mg / Ha)",
       xlab="Years Before Present",main=NA, xaxt='n')
  
  title(x.meta[x.meta$site.id==site_number,'site.name'][1],outer=TRUE)
  axis(3,at=seq(0,10000,1000),labels=seq(0,10000,1000),padj=1)
  
  ciEnvelope(seq(100,10000,100),biomassCI[1,],biomassCI[3,],col="gray")
  
  points(seq(100,10000,100),biomassCI[2,],cex=.8,pch=16,col = colors[data_binned])
  
  keep.dataset.id <- unique(x.meta[x.meta$site.id==site_number,4])
  
  rug(x.meta[x.meta[,1]==site_number,]$age_bacon,lwd=2)
  #rug(control.pts[which(control.pts[,1]%in%keep.dataset.id),]$geo_age,lwd=3,col="red")
  
  ten.count.use = ten.count[which(x.meta[,1]==site_number),]
  prop.use <- prop.table(as.matrix(ten.count.use),margin=1)    
  
  for(p in rev(match(names(sort(colMeans(prop.use))),colnames(prop.use)))){
    prop.plot<- cbind(as.vector(x.meta[which(x.meta$site.id==site_number),]$age_bacon),as.matrix(prop.use[,p]))      	
    prop.plot<-prop.plot[order(prop.plot[,1]),]
    plot(x=prop.plot[,1],y=prop.plot[,2],type="l",xlim=c(10000,-10),
         ylim=c(0,max(prop.use[,p])),ylab=NA,yaxt='n', xaxt='n')
    #axis(2,at=signif(max(prop.use)/2),labels=signif(max(prop.use[,p])/2,digits=1))
    ciEnvelope(prop.plot[,1],rep(0,nrow(prop.use)),prop.plot[,2],col="darkblue")
    legend('topleft',colnames(prop.use)[p])
    #legend('topleft',paste(signif(mean(prop.use[,p]),digits=2)),bty="n")
  } 
  
}

plot_biomass_ts(site_number=unique(x.meta$site.id)[182], biomassCI=biomassCI[[182]])

##### Average time series 
pdf('increase.decrease.pdf')
breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))
plot(seq(100,10000,100),biomassCI[[6]][2,],cex=.1,ylim=c(0,150),
     xlim=c(10000,-10),ylab="Biomass (Mg / Ha)",xlab="Years BP",main='Biomass')
for(i in 6:182){
  if(length(biomassCI[[i]])>1){
    biomassCI.use <- biomassCI[[i]]
    #data_binned <-  cut(biomassCI.use[2,], c(breaks), include.lowest = FALSE, labels = FALSE)
    breaks <- c(-10000000,0,10000000)
    data_binned <-  cut(diff(rev(biomassCI.use[2,])), c(breaks), include.lowest = FALSE, labels = FALSE)
    points(seq(100,10000,100),biomassCI.use[2,],col=c(1,rev(data_binned)),pch=19)
  }
}
dev.off()

bio.quant <- apply(matrix(unlist(lapply(biomassCI,function(x) x[2,])),183,100,byrow = TRUE),2,quantile,c(0.025,0.5,0.975))
data_binned <-  cut(bio.quant[2,], c(breaks), include.lowest = FALSE, labels = FALSE)
diff.list <- list()
for(i in 1:182){
  if(length(biomassCI[[i]])>1 & i !=103 & i != 81){
    diff.list[[i]] <- diff(rev(biomassCI[[i]][2,]),lag=1)
  }else{
    diff.list[[i]] <- rep(0,99)
  }
}

diff.mat <- matrix(unlist(diff.list),182,99,byrow=TRUE)
diff.prop <- numeric(ncol(diff.mat))
for(i in 1:ncol(diff.mat)){
  diff.prop[i] <- length(which(diff.mat[,i]>0))/44
}
length(which(diff.mat[,1]!=0))

pdf('average.biomass.pdf')
#quartz()
zones=matrix(c(2,1), ncol=1, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
par(mar=c(5,4,1,1))
plot(seq(100,10000,100),bio.quant[2,],
     col=colors[data_binned],pch=19,xlim=c(10000,-10),
     ylim=c(0,150),xlab = 'Year Before Present',ylab='Biomass (Mg/ha)')
for(i in 1:182){  
  if(length(biomassCI[[i]])>1 & i !=103 & i != 81){
    lines(seq(100,10000,100),biomassCI[[i]][2,],col='gray')
  }
}
points(seq(100,10000,100),bio.quant[2,],
       col=colors[data_binned],pch=19)
par(mar=c(0,4,1,1))
barplot(height=diff.prop, ylim=c(0, 1), 
        space=0,xlim=c(0,100),ylab = 'Proportion Increasing',
        cex.lab=.6, cex.axis = .5)
abline(h=.5)
#title('Average Biomass Across Sites')
dev.off()

#####
##### Hemlock Decline #####
#####
pdf('increase.decrease.hemlock.pdf')
breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))
plot(seq(100,10000,100),biomassCI[[6]][2,],cex=.1,ylim=c(0,150),
     xlim=c(10000,-10),ylab="Biomass (Mg / Ha)",xlab="Years BP",main='Biomass',col='white')
for(i in 6:182){

  if(length(biomassCI[[i]])>1){
    
      site_number<-unique(x.meta$site.id)[i]
      ten.count.use = ten.count[which(x.meta[,1]==site_number),]
      prop.use <- prop.table(as.matrix(ten.count.use),margin=1) 
      prop.plot<- cbind(as.vector(x.meta[which(x.meta$site.id==site_number),]$age_bacon),as.matrix(prop.use[,11]))      	
      

  if(sum(prop.plot[prop.plot[,1]>4000&prop.plot[,1]<6000,2])>.04){
    biomassCI.use <- biomassCI[[i]]
    #data_binned <-  cut(biomassCI.use[2,], c(breaks), include.lowest = FALSE, labels = FALSE)
    breaks <- c(-10000000,0,10000000)
    data_binned <-  cut(diff(rev(biomassCI.use[2,])), c(breaks), include.lowest = FALSE, labels = FALSE)
    points(seq(100,10000,100),biomassCI.use[2,],col=c(1,rev(data_binned)),pch=19)
  }
  }
}
dev.off()




#####
##### biomassCI exploration #####
#####

##### biomass CI is ordered by unique(x.meta$site.id)

name.vec <- seq(100,10000,100)
only.means.all<-lapply(biomassCI,function(x){return(x[seq(2,299,3)])})
get.lat.long<-matrix(0,183,2)
for(i in 2:183){
  get.lat.long[i,1]<- unique(x.meta[x.meta[,1]==unique(x.meta[,1])[i],c(3)])
  get.lat.long[i,2]<- unique(x.meta[x.meta[,1]==unique(x.meta[,1])[i],c(2)])
}


names(only.means.all)<-unique(x.meta$site.id)[1:182]

breaks <-  c(seq(0,50,10),seq(75,200,25))
colors <- rev(terrain.colors(length(breaks)))

pdf(paste0('pred.points.map',Sys.Date(),'.pdf'))
for(r in seq(1,99,1)){
  
  only.means <- unlist(lapply(only.means.all,function(x){return(x[r])}))
  
  data_binned <-  cut(only.means, c(breaks), include.lowest = FALSE, labels = FALSE)
  
  long.keep <- list()
  lat.keep <- list()
  for(i in 1:length(only.means)){
    long.keep[[i]] <- x.meta[x.meta[,1]==names(only.means)[i],'long'][1]
    lat.keep[[i]] <- x.meta[x.meta[,1]==names(only.means)[i],'lat'][1]
  }
  
  
  map('state', xlim=c(-98,-81), ylim=c(41.5,50))
  points(unlist(long.keep),unlist(lat.keep), pch=21,
         cex=1.1, bg=colors[data_binned],lwd=.2)
  plotInset(-90,47,-82.5,50,
            expr={
              keep.col<-unique(data_binned)
              keep.col<-keep.col[!is.na(keep.col)]
              keep.col<-sort(keep.col)
              is.na.vec <- rep(NA,10)
              is.na.vec[keep.col]<-colors[keep.col]
              
              hist(data_binned,col=is.na.vec,xaxt="n",xlab=NA,
                   ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,
                   xlim=c(0,length(breaks)),ylim=c(0,20),breaks=seq(0,12,1))
              
              axis(side=1,breaks,at = seq(0,11,1),cex.axis = .5,las=2,line=0)
              mtext(side = 1, "Biomass (Mg/ha)", line = 1.5,cex=.5)
              mtext(side = 2, "Frequency", line = 1.7,cex=.5)
            })
  title(paste("Biomass @",name.vec[r]))
}
dev.off()

#####
##### Second Derivative #####
#####

calc.second.deriv <- function(biomassCI,h,second.deriv){
  for(i in 1:length(biomassCI)){
    if(length(biomassCI[[i]])>10){
      T <- dim(biomassCI[[i]][,1:11])[2] - h
      t <- h + 1 
      biomassCI[[i]]<-(biomassCI[[i]][,1:11]) #log or no log here?
      second.deriv[[i]]<-sum(((biomassCI[[i]][2,(t:T)+h]-2*biomassCI[[i]][2,(t:T)]+biomassCI[[i]][2,(t:T)-h])/((h*100)^2))^2)
      
    }else{
      second.deriv[[i]]<-NA
    }
  }
  return(second.deriv)
}

pdf(paste('second.deriv.map',Sys.Date(),'.pdf'))
for(i in c(2)){
  h=i
  second.deriv<-calc.second.deriv(biomassCI=biomassCI,h=h,second.deriv=list())
  
  #hist(unlist(second.deriv))
  names(second.deriv) <- unique(x.meta[,1])[1:182]
  second.deriv.unlist <- unlist(second.deriv)
  
  second.deriv.unlist <- na.omit(second.deriv.unlist)
  
  
  long.keep <- list()
  lat.keep <- list()
  settlebiomass <- list()
  site.name.list <- list()
  site.id.list <- list()
  for(i in 1:length(second.deriv.unlist)){
    long.keep[[i]] <- x.meta[x.meta[,1]==names(second.deriv.unlist)[i],'long'][1]
    lat.keep[[i]] <- x.meta[x.meta[,1]==names(second.deriv.unlist)[i],'lat'][1]
    settlebiomass[[i]] <- cast.x[cast.x$site.id == names(second.deriv.unlist)[i],ncol(cast.x)][1]
    site.name.list[[i]] <- x.meta[x.meta[,1]==names(second.deriv.unlist)[i],'site.name'][1]
    site.id.list[[i]] <- x.meta[x.meta[,1]==names(second.deriv.unlist)[i],'site.id'][1]
  }

  tot.change.mat <- data.frame(unlist(site.id.list),unlist(long.keep),unlist(lat.keep),unlist(settlebiomass),second.deriv.unlist)
  colnames(tot.change.mat)<-c('site.id','long','lat','settlebiomass','second_deriv1000')
  #write.csv(tot.change.mat,file='Second_Deriv_with_meta.csv')
  
  breaks <-  c(0,quantile(probs=c(.05,.25,.5,.75,.975),second.deriv.unlist,na.rm=TRUE),max(second.deriv.unlist)+1)
  
  colors <- colorRampPalette(c("white","yellow",'green','blue','purple'))(length(breaks)-1)
  
  data_binned <-  cut(second.deriv.unlist, c(breaks), include.lowest = FALSE, labels = FALSE)
  
  map('state', xlim=c(-98,-81), ylim=c(41.5,50),bg='grey')
  points(long.keep, lat.keep, pch=19,
         cex=1.3, col=colors[data_binned],lwd=.2)
  title(paste("Second Derivative Sum Point Estimates h = ",h*100))
  
  plotInset(-90,47,-82.5,50,
            expr={
              hist(data_binned,col=colors,xaxt="n",xlab=NA,
                   ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,breaks=0:length(unique(data_binned)))
              axis(side=1,signif(breaks,digits=2),at = 0:length(unique(data_binned)),cex.axis = .5,las=2,line=0)
              mtext(side = 1, "Squared Sum", line = 1.5,cex=.5)
              mtext(side = 2, "Frequency", line = 1.7,cex=.5)
            })
  
}
dev.off()

mat.bio.sec <- matrix(NA,ncol = length(second.deriv.unlist),
                      nrow = 5)

mat.bio.sec <- cbind(na.omit(second.deriv.unlist[match(names(second.deriv.unlist),
                                                       cast.x$site.id)]),
                     na.omit(cast.x[match(cast.x$site.id,
                                          as.numeric(names(second.deriv.unlist))),
                                    ncol(cast.x)]))

plot(mat.bio.sec[1,],mat.bio.sec[,2])

data.dir = c("/Users/paleolab/babySTEPPS/Data/")
biomass_dat_est <- read.csv(paste0(data.dir,"biomass_prediction_v0.9-7_bam.csv"))

lat.long.reg.df = data.frame(biomass_dat_est$x,biomass_dat_est$y)
colnames(lat.long.reg.df) = c('x', 'y')

coordinates(lat.long.reg.df) <- ~ x + y
proj4string(lat.long.reg.df) <-  CRS('+init=epsg:3175')

albers <- spTransform(lat.long.reg.df,CRS('+proj=longlat +ellps=WGS84'))
albers <- as.matrix(data.frame(albers))

sum.save <- numeric(182)
for(i in 1:182){
  sum.save[i]<-sum(biomassCI[[i]][2,])
}

diff.sum <- numeric(182)
for(i in 1:182){
  diff.sum[i]<-sum(diff(rev(biomassCI[[i]][2,])))
}

breaks <-  c(0,quantile(probs=c(.05,.25,.5,.75,.975),second.deriv.unlist,na.rm=TRUE),max(second.deriv.unlist)+1)
colors <- colorRampPalette(c("white","yellow",'green','blue','purple'))(length(breaks)-1)
data_binned <-  cut(second.deriv.unlist, c(breaks), include.lowest = FALSE, labels = FALSE)


quartz()
pdf('wiggle.map1.pdf')
map('state', xlim=c(-98,-81), ylim=c(41.5,50),col='darkgray',bg='lightgrey')
points(albers[which(biomass_dat_est$Hemlock>1),],pch=11,cex=.1,col='darkgrey')
points(albers[which(biomass_dat_est$Beech>1),],pch=19,cex=.26)
points(unlist(long.keep), unlist(lat.keep), pch=19,
       cex=1, col='darkgrey',lwd=.2)
#title(paste("Second Derivative Sum Point Estimates h = ",h*100))
title('Time Series Map')

for(i in 1:182){
  if(length(biomassCI[[i]])>1){
    long.use <- x.meta[x.meta[,1]==names(biomassCI)[i],'long'][1]
    lat.use<- x.meta[x.meta[,1]==names(biomassCI)[i],'lat'][1]
    color.use <- data_binned[which(names(second.deriv.unlist)==names(biomassCI)[i])]
    
    plotInset(long.use-3,lat.use-1,long.use,lat.use+1,
              expr={  
                plot(seq(100,10000,100),
                     biomassCI[[i]][2,],
                     typ='l',xaxt = 'n',
                     yaxt = 'n',xlab = NA,ylab=NA,xlim=c(10000,0),
                     bty = 'n',ylim=c(0,150),col=colors[color.use],lwd=3)
              })
  }
}
plotInset(-90,47,-82.5,50,
          expr={
            hist(data_binned,col=colors,xaxt="n",xlab=NA,
                 ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,breaks=0:length(unique(data_binned)))
            axis(side=1,signif(breaks,digits=2),
                 at = 0:length(unique(data_binned)),
                 cex.axis = .5,las=2,line=0)
            mtext(side = 1, "Squared Sum", line = 1.5,cex=.5)
            mtext(side = 2, "Frequency", line = 1.7,cex=.5)
          })
dev.off()

breaks <-  c(0,quantile(probs=c(.05,.25,.5,.75,.975),second.deriv.unlist,na.rm=TRUE),max(second.deriv.unlist)+1)
colors <- colorRampPalette(c("white","yellow",'green','blue','purple'))(length(breaks)-1)
data_binned <-  cut(second.deriv.unlist, c(breaks), include.lowest = FALSE, labels = FALSE)


zones=matrix(c(2,1), ncol=1, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
par(mar=c(5,4,1,1))
plot(seq(100,10000,100),bio.quant[2,],
     col='white',pch=19,xlim=c(10000,-10),
     ylim=c(0,150),xlab = 'Year Before Present',ylab='Biomass (Mg/ha)')
for(i in 1:182){  
  if(length(biomassCI[[i]])>1 & i !=103 & i != 81){
    color.use <- data_binned[which(names(second.deriv.unlist)==names(biomassCI)[i])]
    lines(seq(100,10000,100),biomassCI[[i]][2,],col=colors[color.use],lwd=2)
  }
}
par(mar=c(0,4,1,1))
barplot(height=diff.prop, ylim=c(0, 1), 
        space=0,xlim=c(0,100),ylab = 'Proportion Increasing',
        cex.lab=.6, cex.axis = .5)
abline(h=.5)


pdf('stability.group.ts.100.pdf')
par(mfrow=c(2,3))
for(count in 1:6){
  plot(seq(100,10000,100),seq(100,10000,100),
       cex=.1,ylim=c(0,150),xlim=c(10000,-10),ylab="Biomass (Mg / Ha)",
       xlab="Years Before Present",main=paste('Stability Bin',count))
  for(r in names(second.deriv.unlist[data_binned==count])){
    biomassCI.use <- biomassCI[[which(unique(x.meta$site.id)==as.numeric(r))]]
    print(x.meta[x.meta$site.id==r,'site.name'])[[1]]
    points(seq(100,10000,100),biomassCI.use[2,],
           main=(x.meta[x.meta$site.id==r,'site.name'])[[1]],
           ylim=c(0,150), pch = 21,bg=colors[count])
  }
}
dev.off()

#####
##### Principle Component Analysis #####
#####

prop.ten.count = prop.table(ten.count,margin = 1)

#whole sample #oak,pine,prairie
pca=princomp(prop.ten.count) 
summary(pca) 
biplot(pca)

#pca by stability group #number is time step
pdf('stability.group.pca.2000.pdf')
for(count in 1:6){
  pca.new <- numeric(ncol(prop.ten.count))
  for(b in as.numeric(names(second.deriv.unlist[data_binned==count]))){
    pca.use <- prop.ten.count[x.meta$site.id==b,]
    #print(pca.use)
    pca.new <- rbind(pca.use,pca.new)
  }
pca=princomp(pca.new) 
#summary(pca)
biplot(pca,main=paste('Stability Bin',count))
}
dev.off()

library(vegan)

save.cats <- list()
for(i in 1:nrow(prop.ten.count)){
  #save.cats[[i]] <- length(which(prop.ten.count[i,] > .1))
  save.cats[[i]] <- diversity(prop.ten.count[i,])
}

x.meta.new <- cbind(x.meta, unlist(save.cats))
x.meta.new <- cbind(x.meta.new, numeric(nrow(x.meta)))
colnames(x.meta.new) <- c(colnames(x.meta),c('divIndex','estBiomass'))

names(biomassCI) <- unique(x.meta$site.id)[1:182]

pdf('diversity.biomass.pdf')
plot(save.cats[which(x.meta$site.id==site.num)],
     biomassCI[[182]][2,x.meta[x.meta$site.id==site.num,'age_bacon']/100],
     pch = 19,ylab = 'Biomass',xlab = 'Diversity Index',
     ylim = c(0,150),xlim=c(min(unlist(save.cats)),max(unlist(save.cats))))

for(i in 1:182){
  if(length(biomassCI[[i]])>1){
    site.num <- names(biomassCI[i])
    points(save.cats[which(x.meta$site.id==site.num)],
         biomassCI[[i]][2,x.meta[x.meta$site.id==site.num,'age_bacon']/100],
         pch = 19)
    x.meta.new[which(x.meta$site.id==site.num),'estBiomass'] <- biomassCI[[i]][2,x.meta[x.meta$site.id==site.num,'age_bacon']/100]
    
 }
}
dev.off()


new.site.locs <- cbind(x.meta.new$lon,x.meta.new$lat)
centers_pol = data.frame(new.site.locs)
colnames(centers_pol) = c('x', 'y')

coordinates(centers_pol) <- ~ x + y
proj4string(centers_pol) <- CRS('+proj=longlat +ellps=WGS84')

centers_polA <- spTransform(centers_pol, CRS('+init=epsg:3175'))
centers_polA <- as.matrix(data.frame(centers_polA))

all.preds1 = cbind(x.meta.new,centers_polA)

all.preds2 <- all.preds1[all.preds1$estBiomass!=0,]

b <- gam(log(estBiomass) ~ te(x,y,age_bacon, d = c(2,1),bs = c("tp","cr"),k=30), data = as.data.frame(all.preds2))


#### Vegetation Proportions and Biomass
pdf('estimate.props.pdf')
par(mfrow=c(3,3))
for(i in 1:21){
  plot(x.meta.new$estBiomass,prop.ten.count[,i],
       main=colnames(ten.count)[i],
       ylab='Pollen Proportion',
       xlab='Mean Biomass Estimate (Mg/ha)',
       pch = 19,
       ylim=c(0,1),
       xlim=c(10,150),
       cex = .5)
}
dev.off()

hem.beech.bio <- x.meta.new$estBiomass[which(prop.ten.count[,11]>.05|prop.ten.count[,19]>.05)]

site.id.save <- unique(x.meta$site.id[which(prop.ten.count[,11]>.05|prop.ten.count[,19]>.05)])

hist(x.meta[which(prop.ten.count[,11]>.05|prop.ten.count[,19]>.05),
       c('age_bacon')]) #hemlock decline

site.age <- x.meta[which(prop.ten.count[,11]>.05|prop.ten.count[,19]>.05),
       c('site.id','age_bacon')]

for(s in site.id.save){
  site.age[site.age$site.id==s,]
  biomassCI[[which(match(unique(x.meta$site.id),s)!=0)]]
}

site.get <- which(match(unique(x.meta$site.id),site.id.save)!=0)

diff.sum <- numeric(182)
for(i in 1:182){
  diff.sum[i]<-sum(diff(rev(biomassCI[[i]][2,])))
}
sum(diff.sum[diff.sum!=0])
hist(diff.sum[diff.sum!=0])


map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(x.meta.new[which(prop.ten.count[,11]>.05|prop.ten.count[,19]>.05),'long'],
       x.meta.new[which(prop.ten.count[,11]>.05|prop.ten.count[,19]>.05),'lat'])

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(x.meta[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass>100),'long'],
       x.meta[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass>100),'lat'],
       pch=19,
       cex=1.3,lwd=.2)
hist(x.meta[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass>100),'age_bacon'])
pca.new1 <- prop.ten.count[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass>100),]
pca1=princomp(pca.new)
summary(pca1)
biplot(pca1,choices = 2:3)

map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(x.meta[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass<40),'long'],
       x.meta[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass<40),'lat'],
       pch=19,
       cex=1.3,lwd=.2)
hist(x.meta[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass<40),'age_bacon'])
pca.new <- prop.ten.count[which(x.meta.new$divIndex>1.75&x.meta.new$estBiomass<40),]
pca=princomp(pca.new)
summary(pca)
biplot(pca,choices = 2:3)

pdf('diversity.lobe.pies.pdf')
par(mfrow=c(1,2))
pie(colMeans(pca.new),main='low biomass')
pie(colMeans(pca.new1),main = 'high biomass')
dev.off()


list.biomass <- lapply(biomassCI,function(x) x[2,])

mat.biomass <- data.frame(c(list.biomass[[6]],
                        as.data.frame(x.meta[which(x.meta$site.id==names(list.biomass)[[6]])[1],])))

for(i in 7:182){
  if(length(list.biomass[[i]]) > 2){
    mat.biomass <- rbind(mat.biomass,
                         data.frame(c(list.biomass[[i]],
                           as.data.frame(x.meta[which(x.meta$site.id==names(list.biomass)[[i]])[1],]))))
  }
}

library(corrplot)
rownames(mat.biomass) <- mat.biomass$site.name
pdf('lake.ts.corrs.last1000.pdf')
corrplot(cor(t(mat.biomass[,1:11])),order='AOE')
dev.off()
cor.use <- cor(t(mat.biomass[,1:100]))

pdf('corr.map.holocene.pdf')
map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(x.meta.new[which(x.meta.new$estBiomass!=0),'long'],
       x.meta.new[which(x.meta.new$estBiomass!=0),'lat'],pch=19,col='gray')
for(i in c('Seidel','Rutz Lake','Kimble Pond','Pogonia Bog Pond', 'Wolsteld Lake','Sharkey Lake','Kirchner Marsh')){
  points(x.meta[which(x.meta$site.name==i),'long'],
         x.meta[which(x.meta$site.name==i),'lat'],pch=19)
}
title('Points Correlated Over Holocene')
dev.off()

pdf('corr.map.1000.pdf')
map('state', xlim=c(-98,-81), ylim=c(41.5,50))
points(x.meta.new[which(x.meta.new$estBiomass!=0),'long'],
       x.meta.new[which(x.meta.new$estBiomass!=0),'lat'],pch=19,col='gray')
for(i in c('Pogonia Bog Pond','Capitola Lake','East Soldier Lake',
           'Devils Lake','Lake Minnie','Jay Lake',
           'Camp 11 Lake','Penegor Lake','Kellys Hollow','Wood Lake')){
  points(x.meta[which(x.meta$site.name==i),'long'],
         x.meta[which(x.meta$site.name==i),'lat'],pch=19)
}

for(i in c('Rutz Lake','Lorraine Lake','Spirit Lake','Cub Lake',
           'Kellys Hollow','Wood Lake','Steel Lake')){
  points(x.meta[which(x.meta$site.name==i),'long'],
         x.meta[which(x.meta$site.name==i),'lat'],pch=19,col='black')
}
title('Points Correlated Over Last 1000')
dev.off()

plot.cor <- colSums(cor.use)

breaks <-  c(min(plot.cor)-.01,quantile(plot.cor,seq(.1,.9,.2)),max(plot.cor)+.01)
colors <- rainbow(length(breaks))#colorRampPalette(c("yellow",'green','blue','purple'))(length(breaks)-1)
data_binned <-  cut(plot.cor, c(breaks), include.lowest = FALSE, labels = FALSE)

zones=matrix(c(2,1), ncol=1, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
par(mar=c(5,4,1,1))
plot(seq(100,10000,100),bio.quant[2,],
     col='white',pch=19,xlim=c(10000,-10),
     ylim=c(0,150),xlab = 'Year Before Present',ylab='Biomass (Mg/ha)')
for(i in 1:length(plot.cor)){
    color.use <- data_binned[i]
    lines(seq(100,10000,100),mat.biomass[i,1:100],col=colors[color.use],lwd=3)
}
par(mar=c(0,4,1,1))
barplot(height=diff.prop, ylim=c(0, 1), 
        space=0,xlim=c(0,100),ylab = 'Proportion Increasing',
        cex.lab=.6, cex.axis = .5)
abline(h=.5)

pdf('lines.by.corrs.sums.pdf')
map('state', xlim=c(-98,-81), ylim=c(41.5,50),col='darkgray',bg='lightgrey')
points(albers[which(biomass_dat_est$Hemlock>1),],pch=11,cex=.1,col='darkgrey')
points(albers[which(biomass_dat_est$Beech>1),],pch=19,cex=.26)
points(unlist(long.keep), unlist(lat.keep), pch=19,
       cex=1, col='darkgrey',lwd=.2)
title('Time Series Map')
for(i in 1:length(plot.cor)){
  long.use <- mat.biomass$long[i]
  lat.use <- mat.biomass$lat[i]
  color.use <- data_binned[i]
  
  plotInset(long.use-3,lat.use-1,long.use,lat.use+1,
            expr={  
              plot(seq(100,10000,100),
                   mat.biomass[i,1:100],
                   typ='l',xaxt = 'n',
                   yaxt = 'n',xlab = NA,ylab=NA,xlim=c(10000,0),
                   bty = 'n',ylim=c(0,150),col=colors[color.use],lwd=3)
            })
}
plotInset(-90,47,-82.5,50,
          expr={
            hist(data_binned,col=colors,xaxt="n",xlab=NA,
                 ylab=NA,main=NA,cex.lab=.5,cex.axis=.5,breaks=0:length(unique(data_binned)))
            axis(side=1,signif(breaks,digits=2),
                 at = 0:length(unique(data_binned)),
                 cex.axis = .5,las=2,line=0)
            mtext(side = 1, "Squared Sum", line = 1.5,cex=.5)
            mtext(side = 2, "Frequency", line = 1.7,cex=.5)
          })
zones=matrix(1:6, ncol=1, byrow=TRUE)
layout(zones)
for(r in 1:max(data_binned)){
  par(mar=rep(.5,4))
  plot(seq(100,10000,100),bio.quant[2,],
       col='white',pch=19,xlim=c(10000,-10),
       ylim=c(0,150),xlab = 'Year Before Present',ylab='Biomass (Mg/ha)')
  for(i in 1:length(plot.cor)){
    if(data_binned[i]==r){
      color.use <- data_binned[i]
      lines(seq(100,10000,100),mat.biomass[i,1:100],
            col=colors[color.use],lwd=2)
    }
  }
}

dev.off()


breaks <-  c(0,quantile(probs=c(.05,.25,.5,.75,.975),second.deriv.unlist,na.rm=TRUE),max(second.deriv.unlist)+1)
colors <- colorRampPalette(c("white","yellow",'green','blue','purple'))(length(breaks)-1)
data_binned <-  cut(second.deriv.unlist, c(breaks), include.lowest = FALSE, labels = FALSE)

fig.mat <- matrix(1,16,1)
fig.mat[1:6,]<-1
fig.mat[7:16,]<-seq(2,11,1)
#control.pts<-read.csv(text=getURL('https://raw.githubusercontent.com/PalEON-Project/stepps-baconizing/master/data/bacon_age_markers_v1.csv'))

pdf('second.deriv.with.pollen.pdf')
layout(fig.mat)
par(mar = c(0, 4, 0, 3), oma = c(5, 0, 4, 0),tcl=.5,bg='thistle')

for(r in 1:max(data_binned)){
  plot(seq(100,10000,100),bio.quant[2,],
       col='thistle',pch=19,xlim=c(10000,-10),
       ylim=c(0,150),xlab = 'Year Before Present',ylab='Biomass (Mg/ha)')
  for(i in 1:length(plot.cor)){
    if(data_binned[i]==r){
      color.use <- data_binned[i]
      lines(seq(100,10000,100),mat.biomass[i,1:100],
            col=colors[color.use],lwd=2)
    }
  }
  
  ten.count.use1 = ten.count[match(x.meta[,1],mat.biomass[which(data_binned==r),'site.id']),]
  prop.use1 <- prop.table(as.matrix(ten.count.use),margin=1)  
  p.do <- c(11,19,match(names(sort(colMeans(prop.use1),decreasing = TRUE)[1:8]),colnames(prop.use1)))
  
  for(p in p.do){
    plot(x=1:10000,y=1:10000,type="l",xlim=c(10000,-10),
         ylim=c(-.01,max(prop.stop[,p])),ylab=NA,yaxt='n', xaxt='n',col='thistle')
  for(i in 1:length(plot.cor)){
    if(data_binned[i]==r){
      ten.count.use = ten.count[which(x.meta[,1]==mat.biomass$site.id[i]),]
      prop.use <- prop.table(as.matrix(ten.count.use),margin=1)    
      
        prop.plot<- cbind(as.vector(x.meta[which(x.meta[,1]==mat.biomass$site.id[i]),]$age_bacon),as.matrix(prop.use[,p]))      	
        prop.plot<-prop.plot[order(prop.plot[,1]),]
        #lines(x=prop.plot[,1],y=prop.plot[,2],type="l")
        #axis(2,at=signif(max(prop.use)/2),labels=signif(max(prop.use[,p])/2,digits=1))
        lines(prop.plot[,1],prop.plot[,2],col='darkblue',lwd=1)
        legend('topleft',colnames(prop.use)[p])
        #legend('topleft',paste(signif(mean(prop.use[,p]),digits=2)),bty="n")
      }
    }
  }

}
dev.off()


