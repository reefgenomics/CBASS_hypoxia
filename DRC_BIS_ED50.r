library(plyr)
library(gplots)
library(lattice)
library(emmeans)
library(drc)

# upload data
bis<-read.delim("./Input_files/bleaching_scores_CBASS_Hd_T1_infile.txt", header=T, sep="\t") 

bis$Treatment<-as.factor(bis$Treatment)
bis$TreatmentTemp<-paste0(bis$Treatment,bis$Temperature)
bis$RepTemp<-paste0(bis$Temperature, '_', bis$Replicate)
bis$Replicate2<-paste0(bis$Treatment, '_', bis$Replicate)
bis$FacTemp<-as.factor(bis$Temperature)
aggregate(BIS~Treatment + Temperature, data=bis, summary )

# Checking sample sizes
aggregate(BIS~Treatment + Temperature, data=bis, length)
aggregate(BIS ~ Treatment + FacTemp, data=bis, length)
aggregate(BIS ~ Treatment + FacTemp + Replicate2, data=bis, length)


#### DRC Curve Fitting####
#getMeanFunctions()
DRCBIS = drm(BIS ~ Temperature, data = bis, curveid = Treatment,
             fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCBIS)
compParm(DRCBIS, 'ed50')
compParm(DRCBIS, 'ed50', "-")
plot(DRCBIS)
points(bis$Temperature, bis$BIS)
ED(DRCBIS, c(50))

# fit to each treatment individually
#### H0####
DRCBISH0 = drm(BIS ~ Temperature, data = bis[bis$Treatment=="H0",], #change to Hd/H0
               fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCBISH0)
DRCBISH0$coefficients[3]
ED(DRCBISH0, c(50))[,1]

# replicate-specific curve fits
DRCBISH0rep = drm(BIS ~ Temperature, data = bis[bis$Treatment=="H0",], curveid=Replicate2,
                  fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCBISH0rep)
DRCBISH0rep$coefficients[9:12] #changed row range as now 4 rather than 3 reps
ED(DRCBISH0rep, c(50))[,1]

#### Hd####
DRCBISHd = drm(BIS ~ Temperature, data = bis[bis$Treatment=="Hd",],
               fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCBISHd)
DRCBISHd$coefficients[3]
ED(DRCBISHd, c(50))

# replicate-specific curve fits
DRCBISHdrep = drm(BIS ~ Temperature, data = bis[bis$Treatment=="Hd",], curveid=Replicate2,
                  fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCBISHdrep)
DRCBISHdrep$coefficients[9:12]
ED(DRCBISHdrep, c(50))[,1]

#### Merging Coeffs ####
Coeffs<-c(DRCBISH0$coefficients[3],DRCBISHd$coefficients[3])
RepCoeffs<-data.frame("ED50"=c(DRCBISH0rep$coefficients[9:12],DRCBISHdrep$coefficients[9:12]),"Treatment"=c(rep("H0",4), rep("Hd",4))) #changed rep from 3 to 4 

aggregate(ED50 ~ Treatment, data=RepCoeffs, mean)

#### ED50 ttest
aggregate(ED50 ~ Treatment, data=RepCoeffs, FUN= function(x) shapiro.test(x)$p.value)
bartlett.test(ED50 ~ Treatment, data=RepCoeffs)

ED50.ttest<- t.test(ED50 ~ Treatment, data=RepCoeffs, alternative= "greater")
ED50.ttest

aggregate(ED50 ~ Treatment, data=RepCoeffs, summary)

RepCoeffs
SumaryStats<-data.frame("Treatment"=names(tapply(RepCoeffs$ED50,RepCoeffs$Treatment, mean)), "MeanED50"=tapply(RepCoeffs$ED50,RepCoeffs$Treatment, mean), "ED50StdDev"=tapply(RepCoeffs$ED50,RepCoeffs$Treatment, sd), "ED50StdErr"=tapply(RepCoeffs$ED50,RepCoeffs$Treatment, sd)/sqrt(tapply(RepCoeffs$ED50,RepCoeffs$Treatment, length)))
SumaryStats<-SumaryStats[order(SumaryStats$MeanED50),]
write.table(data.frame("Replicate"=gsub("ed50:", "", row.names(RepCoeffs)), RepCoeffs), file="./ReplicateED50s_BIS.txt", quote=F, sep="\t", row.names=F)
write.table(SumaryStats, file="./Heat+Deoxy_ED50_BIS_SummaryStats.txt", quote=F, sep="\t", row.names=F)

#### Plotting ####

Sytes<-c("H0","Hd")
Colorz<-c('#CCCCFF','#9933CC')
Syms<-c(19,17)
temp_x<- seq(29, 40, length = 100)
Density <- 45

pdf("./CBASS_Hd_H0_DRCCurves_BIS.pdf",9,6)
par(mar=c(5,6,4,1)+.1)

line_width=2 
offsets<-c(0.1875,0.0625,-0.0625,-0.1875)

i<-1 #H0
matplot(temp_x, predict(DRCBISH0, data.frame(Temperature = temp_x), interval= "confidence"),
        type="n",col=Colorz[i],lty=c(1,3,3),lwd=line_width,ylab="Bleaching index score",xlab="Temperature [°C]", xlim=c(29,40),ylim=c(0,7), cex.axis=2.5, cex.lab=2.5)
polygon(c(temp_x, rev(temp_x)),c(predict(DRCBISH0, data.frame(Temp = temp_x), interval="confidence")[,2],rev(predict(DRCBISH0, data.frame(Temp = temp_x), interval="confidence")[,3])), col=Colorz[i], density = Density)
matpoints(temp_x, predict(DRCBISH0, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col=Colorz[i],lty=c(1,3,3),lwd=line_width)
with(bis[bis$Treatment==Sytes[i],],matpoints(Temperature-offsets[i],BIS,pch=Syms[i], col=Colorz[i], cex=2.5))
i<-2 #Hd
polygon(c(temp_x, rev(temp_x)),c(predict(DRCBISHd, data.frame(Temp = temp_x), interval="confidence")[,2],rev(predict(DRCBISHd, data.frame(Temp = temp_x), interval="confidence")[,3])), col=Colorz[i], density = Density)
matpoints(temp_x, predict(DRCBISHd, data.frame(Temperature = temp_x), interval= "confidence"),
          type="l",col=Colorz[i],lty=c(1,3,3),lwd=line_width)
with(bis[bis$Treatment==Sytes[i],],matpoints(Temperature-offsets[i],BIS,pch=Syms[i], col=Colorz[i], cex=2.5))

legend("bottomleft",c("H0; BIS = 37.11°C ","Hd; BIS = 36.02°C"),pch=Syms, col=Colorz,pt.cex=2, bty="n",cex=2)
title(main="CBASS T1 H0 vs Hd")

abline(v= 37.11, col='#CCCCFF', lwd=3, lty=3)
abline(v= 36.02, col='#9933CC', lwd=3, lty=3)

dev.off()
