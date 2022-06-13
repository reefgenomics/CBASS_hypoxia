library(plyr)
library(gplots)
library(lattice)
library(emmeans)
library(drc)

######### Dose-response curves for FvFm over temperature, determines ED50 #########

pamdata<-read.delim("./Input_files/Hypoxia_PAM_data_infile_Jun6.txt", header=T, sep="\t")   

pamdata$PAM<-pamdata$PAM
pamdata$Treatment<-as.factor(pamdata$Treatment)
pamdata$TreatmentTemp<-paste0(pamdata$Treatment,pamdata$Temperature)
pamdata$RepTemp<-paste0(pamdata$Temperature, '_', pamdata$Replicate)
pamdata$Replicate2<-paste0(pamdata$Treatment, '_', pamdata$Replicate)
pamdata$FacTemp<-as.factor(pamdata$Temperature)
aggregate(PAM~Treatment + Temperature, data=pamdata, summary )

#Checking sample sizes
aggregate(PAM~Treatment + Temperature, data=pamdata, length)
aggregate(PAM ~ Treatment + FacTemp, data=pamdata, length)
aggregate(PAM ~ Treatment + FacTemp + Replicate2, data=pamdata, length)

#DRC Curve Fitting
#getMeanFunctions()
DRCpam = drm(PAM ~ Temperature, data = pamdata, curveid = Treatment,
             fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpam)
compParm(DRCpam, 'ed50')
compParm(DRCpam, 'ed50', "-")
plot(DRCpam)
points(pamdata$Temperature, pamdata$PAM)
ED(DRCpam, c(50))

#fit to each treatment individually
#H0
DRCpamH0 = drm(PAM ~ Temperature, data = pamdata[pamdata$Treatment=="H0",], 
               fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamH0)
DRCpamH0$coefficients[3]
ED(DRCpamH0, c(50))[,1]

#replicate-specific curve fits
DRCpamH0rep = drm(PAM ~ Temperature, data = pamdata[pamdata$Treatment=="H0",], curveid=Replicate2,
                  fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamH0rep)

compParm(DRCpamH0rep, 'ed50')
compParm(DRCpamH0rep, 'ed50', "-")

DRCpamH0rep$coefficients[7:9]
ED(DRCpamH0rep, c(50))[,1]

#Hd
DRCpamHd = drm(PAM ~ Temperature, data = pamdata[pamdata$Treatment=="Hd",],
               fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamHd)
DRCpamHd$coefficients[3]
ED(DRCpamHd, c(50))

#replicate-specific curve fits
DRCpamHdrep = drm(PAM ~ Temperature, data = pamdata[pamdata$Treatment=="Hd",], curveid=Replicate2,
                  fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamHdrep)

compParm(DRCpamHdrep, 'ed50')
compParm(DRCpamHdrep, 'ed50', "-")

DRCpamHdrep$coefficients[7:9]
ED(DRCpamHdrep, c(50))[,1]

#Merging Coeffs
Coeffs<-c(DRCpamH0$coefficients[3],DRCpamHd$coefficients[3])
RepCoeffs<-data.frame("ED50"=c(DRCpamH0rep$coefficients[7:9],DRCpamHdrep$coefficients[7:9]),"Treatment"=c(rep("H0",3), rep("Hd",3)))

aggregate(ED50 ~ Treatment, data=RepCoeffs, mean)

#ED50 ttest
aggregate(ED50 ~ Treatment, data=RepCoeffs, FUN= function(x) shapiro.test(x)$p.value)
bartlett.test(ED50 ~ Treatment, data=RepCoeffs)

ED50.ttest<- t.test(ED50 ~ Treatment, data=RepCoeffs, alternative= "greater")
ED50.ttest

aggregate(ED50 ~ Treatment, data=RepCoeffs, summary)

RepCoeffs
SumaryStats<-data.frame("Treatment"=names(tapply(RepCoeffs$ED50,RepCoeffs$Treatment, mean)), "MeanED50"=tapply(RepCoeffs$ED50,RepCoeffs$Treatment, mean), "ED50StdDev"=tapply(RepCoeffs$ED50,RepCoeffs$Treatment, sd), "ED50StdErr"=tapply(RepCoeffs$ED50,RepCoeffs$Treatment, sd)/sqrt(tapply(RepCoeffs$ED50,RepCoeffs$Treatment, length)))
SumaryStats<-SumaryStats[order(SumaryStats$MeanED50),]
#write.table(data.frame("Replicate"=gsub("ed50:", "", row.names(RepCoeffs)), RepCoeffs), file="./ReplicateED50s.txt", quote=F, sep="\t", row.names=F)
#write.table(SumaryStats, file="./Heat+Hypoxia_ED50_SummaryStats.txt", quote=F, sep="\t", row.names=F)

#Plotting curves

Sytes<-c("H0","Hd")
Colorz<-c('#CCCCFF','#9933CC')
Syms<-c(19,17)
temp_x<- seq(29, 40, length = 100)
Denscity <- 45

pdf("./CBASS_Hd_H0_DRCCurves.pdf",10,7)

line_width=2 
offsets<-c(0.1875,0.0625,-0.0625,-0.1875)
i<-1 #H0
matplot(temp_x, predict(DRCpamH0, data.frame(Temperature = temp_x), interval= "confidence"),
        type="n",col=Colorz[i],lty=c(1,3,3),lwd=line_width,ylab="Fv/Fm",xlab="Temperature [°C]", xlim=c(29,40),ylim=c(0,0.65), cex.axis=2.5, cex.lab=2.5)
polygon(c(temp_x, rev(temp_x)),c(predict(DRCpamH0, data.frame(Temp = temp_x), interval="confidence")[,2],rev(predict(DRCpamH0, data.frame(Temp = temp_x), interval="confidence")[,3])), col=Colorz[i], density = Denscity)
matpoints(temp_x, predict(DRCpamH0, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col=Colorz[i],lty=c(1,3,3),lwd=line_width)
with(pamdata[pamdata$Treatment==Sytes[i],],matpoints(Temperature-offsets[i],PAM,pch=Syms[i], col=Colorz[i], cex=2.5))
i<-2 #Hd
polygon(c(temp_x, rev(temp_x)),c(predict(DRCpamHd, data.frame(Temp = temp_x), interval="confidence")[,2],rev(predict(DRCpamHd, data.frame(Temp = temp_x), interval="confidence")[,3])), col=Colorz[i], density = Denscity)
matpoints(temp_x, predict(DRCpamHd, data.frame(Temperature = temp_x), interval= "confidence"),
          type="l",col=Colorz[i],lty=c(1,3,3),lwd=line_width)
with(pamdata[pamdata$Treatment==Sytes[i],],matpoints(Temperature-offsets[i],PAM,pch=Syms[i], col=Colorz[i], cex=2.5))

legend("bottomleft",c("H0; ED50 = 36.12°C ","Hd; ED50 = 35.73°C"),pch=Syms, col=Colorz,pt.cex=2, bty="n",cex=2)
title(main="CBASS T1 H0 vs Hd")

abline(v= 36.12, col='#CCCCFF', lwd=3, lty=3)
abline(v= 35.73, col='#9933CC', lwd=3, lty=3)

dev.off()
