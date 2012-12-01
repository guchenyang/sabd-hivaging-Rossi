##% HIV-Aging: Motor Domain
##% Brown HIV Aging

## # Preprocessing #
d <- read.csv("data.rossi.041312.csv")
d$id <- factor(d$id)

d <- data.frame(d)
##Raw scores only from now on
##groupings for raw scores
grp.exe <- c(10:14)
grp.mem <- c(17:20)
grp.mot <- c(21:22)
grp.spd <- c(23:26)
##Remove 40minus subjects.
d <- d[d$age>=40,]    # dim=280*33

##Remove subjects missing either baseline or 12 month visit.
d1 <- d[d$visit==0, ]
d2 <- d[d$visit==12,]
dd <- merge(d1, d2, by=("id"))    # dim=85*65
dl <- rbind(d1, d2)


##calculate changes in the motor domain
aa <- data.frame(as.matrix(dd[,grp.mot+32]-dd[,grp.mot]))
names(aa) <- paste(names(aa), "change", sep="")
dd <- cbind(dd, aa)

##calculate relative changes in the ececutive domain
aa <- data.frame(as.matrix( (dd[,grp.mot+32]-dd[,grp.mot] )/(dd[,grp.mot]+1)   ))
names(aa) <- paste(names(aa), "rchange", sep="")
dd <- cbind(dd, aa)


##Separate into 3 age groups: 40-49, 50-59, and 60+, with group sizes
##40-49: 39; 50-59: 35; 60+: 10
##Use groups instead of numeric values of age
dd$age.grp <- NA
dd$age.grp[dd$age.x>=40] <- "40"
dd$age.grp[dd$age.x>=50] <- "50"
dd$age.grp[dd$age.x>=60] <- "60"
dd$age.grp <- factor(dd$age.grp)


##Separate into 2 age groups: 40-55, and 55+, with group sizes
##40-49: 39; 50-59: 35; 60+: 10
##Use groups instead of numeric values of age
dd$age.grp55 <- NA
dd$age.grp55[dd$age.x>=40] <- "40"
dd$age.grp55[dd$age.x>=55] <- "55"
dd$age.grp55 <- factor(dd$age.grp55)


##Univaraite approaches to yield easier interpretation

##Add clinical information: duration, CD4, CD4 Nadir, HCV, HIV/RAV4, CART
dc <- read.csv("AgeEffectsMASTER031212.csv", sep=";")
sel <- c(1,2,13:22)
dc <- data.frame(dc[,sel])
dc1 <- dc[dc$Visit==0,]
dc2 <- dc[dc$Visit==12,]
dcd <- merge(dc1, dc2, by=("ID"))
names(dcd)[1] <- "id"
################ merge all together, non-hiv=2, hiv>100 = 1, hiv<100 = 0
dall <- merge(dd, dcd, by="id")
dall$CD4current100.x <-  1*(dall$CD4current.x>=100)+2*(dall$hiv.x==0)
dall$CD4nadir100.x <-  1*(dall$CD4nadir.x>=100)+2*(dall$hiv.x==0)
dall$CD4nadir100.y <-  1*(dall$CD4nadir.y>=100)+2*(dall$hiv.x==0)

##Dictomize age by 55 age cut-off
dall$age.grp2 <- 0 
dall$age.grp2[dall$age.x>=55] <- 1

##Tranform age  by logarithm
dall$age.log <- log(dall$age.x)

##Tranform age  by square
dall$age.sq <- (dall$age.x)^2

################################################################
##Long format 
names(dc)[1:2] <- c("id", "visit")
dl$visit <- factor(dl$visit)
dc$visit <- factor(dc$visit)
dc$id <- factor(dc$id)
dlall <- merge(dl, dc, by=c("id", "visit"))
##Remove missing visits 
a <- table(dlall$id, dlall$visit)
selids <- row.names(a)[apply(a, 1, prod) == 1]
sel <- dlall$id %in% selids
dlall <- dlall[sel,]
dlall$id <- factor(dlall$id)

dlall$CD4nadir100 <-  1*(dlall$CD4nadir>=100)+2*(dlall$hiv==0)


## source CorrectFunctions
source("/Users/Chenyang/Desktop/CorrectingFunctions.R")
##Long format: correct for age, sex, education for non-HIV, 
##Correction for all clinical variables and demographical variables, at baseline and 12month
aa <- data.frame(as.matrix(dlall[,c(grp.mot)]))
for (jj in 1:ncol(aa)) {
  ## visit 0
  vsel <-  dlall$visit == 0 
  sel <- dlall$hiv==0  & vsel
  ##correcting for education only, and apply the same correction to the second visit.  
  tmp <- correctEffect(aa[sel, jj], data.frame(education=dlall[sel,]$education), data.frame(education=13))
  aa[sel, jj] <- tmp$ce
  
  ##correct for hiv by the same formula
  sel <- dlall$hiv==1 & vsel
  rr <- applyCorrect(tmp$fit, aa[sel, jj], data.frame(education=dlall[sel,]$education), data.frame(education=13))
  aa[sel, jj] <- rr
  
  ## correct for HIV status, set to HIVduration=0, detectable=0, cart=0, HCV=0, CD4nadir100=1 (>100), like controls would have
  ## fit <- lm(aa[sel,jj]~HIVduration+HIVRNA_Detectable+cART+HCVcurrent+CD4nadir100, data=dlall[sel,])
  
  ##repeat above for second visit
  vsel <-  dlall$visit == 12
  sel <- dlall$hiv==0  & vsel
  rr <- applyCorrect(tmp$fit, aa[sel, jj],  data.frame(education=dlall[sel,]$education), data.frame(education=13))
  aa[sel,jj] <- rr

  ##correct for hiv by the same formula
  sel <- dlall$hiv==1 & vsel 
  rr <- applyCorrect(tmp$fit, aa[sel, jj], data.frame(education=dlall[sel,]$education), data.frame(education=13))
  aa[sel,jj] <- rr
}

names(aa) <- paste(names(aa), "C", sep="")
dlall <- cbind(dlall, aa)

dlall$age2 <- (dlall$age)^2 

uid <- unique(dlall$id)
dlall$age.grp <-  NA
for (jjj in uid) {
  sel <- dlall$id == jjj
  age <- min(dlall$age[sel])
  if (age <50) {
    dlall$age.grp[sel] <- "40"
  } else if (age>=50 && age<60 ) {
    dlall$age.grp[sel] <- "50"
  } else {
    dlall$age.grp[sel] <- "60"
  }
}
dlall$age.grp <- factor(dlall$age.grp)





##Select color for plots
library(RColorBrewer)
source("/Users/Chenyang/Desktop/MyPlot.R")
mycols <-  brewer.pal(11, "RdBu")
col.n <- mycols[9:11]
col.p <- rev(mycols[1:3])
col.all <- as.vector(rbind(col.n, col.p))
mytheme()


##Correction for all clinical variables and demographical variables, at baseline and 12month
aa <- data.frame(as.matrix(dall[,c(grp.mot, grp.mot+32)]))
ashift <- (ncol(aa)/2)
for (jj in 1:(ncol(aa)/2) ) {
  sel <- dall$hiv.x==0
  ##correcting for education only 
  tmp <- correctEffect(aa[sel, jj], data.frame(education=dall[sel,]$education.x), data.frame(education=13))
  aa[sel, jj] <- tmp$ce
  
  ##correct this for the second visit
  aa[sel, jj+ashift] <- applyCorrect(tmp$fit, aa[sel, jj+ashift], data.frame(education=dall[sel,]$education.y), data.frame(education=13))


  ##correct this for hiv positives as well
  sel <- dall$hiv.x==1
  aa[sel, jj] <- applyCorrect(tmp$fit, aa[sel, jj], data.frame(education=dall[sel,]$education.x), data.frame(education=13))

  aa[sel, jj+ashift] <- applyCorrect(tmp$fit, aa[sel, jj+ashift], data.frame(education=dall[sel,]$education.y), data.frame(education=13))
}

names(aa) <- paste(names(aa), "C", sep="")
dall <- cbind(dall, aa)

##Add age^2 terms
dall <- data.frame(dall, age2 = (dall$age.x)^2)

# dim(dall)
#[1]  84 104

# dall[1,100:103]
#  np.motor.groovedpeg_d.xC np.motor.groovedpeg_nd.xC np.motor.groovedpeg_d.yC np.motor.groovedpeg_nd.yC
#                       88                        91                       75                        74


################################################################
##  #Baseline Analysis
################################################################
## ## MANOVA ##
##Baseline:
##Without Age.sq
fit <- manova(as.matrix(dall[,100:101])~hiv.x*age.x, data=dall)
summary(fit, test="Wilks")
#            Df   Wilks approx F num Df den Df  Pr(>F)  
#hiv.x        1 0.99709   0.1137      2     78 0.89264  
#age.x        1 0.90567   4.0621      2     78 0.02098 *
#hiv.x:age.x  1 0.99414   0.2298      2     78 0.79523  
#Residuals   79  


## ## HIV Positive vs Negative, Univariate ##
##Use regression to compare each raw T score. Not considering CD4current because not significant, consistent with past analyses.  

## ### np.exec.actionflu_total ###
##3 age groups, before AIC:
fit <- lm(np.motor.groovedpeg_d.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  41.3293    23.2226   1.780    0.079 .
#hiv.x       -15.3954    31.7698  -0.485    0.629  
#age.x         0.6653     0.4306   1.545    0.126  
#hiv.x:age.x   0.3944     0.6124   0.644    0.521  

postscript(file="MotBaselineGroovedpegD.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.motor.groovedpeg_d.xC~age.grp*hiv.x, data=dall, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="BD")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()

fit <- lm(np.motor.groovedpeg_d.xC~hiv.x*age.grp, data=dall)
summary(fit)
##                Estimate Std. Error t value Pr(>|t|)    
##(Intercept)       71.406      6.435  11.097   <2e-16 ***
##hiv.x              2.840      7.462   0.381    0.705    
##age.grp50          3.812      8.713   0.438    0.663    
##age.grp60         15.607     10.028   1.556    0.124    
##hiv.x:age.grp50    2.657     10.441   0.254    0.800    
##hiv.x:age.grp60    9.386     15.901   0.590    0.557  

## Uncorrected
fit <- lm(np.motor.groovedpeg_d.x~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  43.7041    23.2703   1.878   0.0641 .
##hiv.x       -17.1210    31.8351  -0.538   0.5922  
##age.x         0.6105     0.4315   1.415   0.1610  
##hiv.x:age.x   0.4385     0.6137   0.714   0.4770  


## ### np.exec.cowat_fas ###
##3 age groups, before AIC:
fit <- lm(np.motor.groovedpeg_nd.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  56.8233    25.4224   2.235   0.0282 *
#hiv.x       -16.7046    34.7792  -0.480   0.6323  
#age.x         0.5073     0.4714   1.076   0.2851  
#hiv.x:age.x   0.3673     0.6704   0.548   0.5854   

postscript(file="MotBaselineGroovedpegND.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.motor.groovedpeg_nd.xC~age.grp*hiv.x, data=dall, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="BS")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()

fit <- lm(np.motor.groovedpeg_nd.x~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  58.6611    25.4256   2.307   0.0237 *
##hiv.x       -18.0400    34.7837  -0.519   0.6055  
##age.x         0.4649     0.4714   0.986   0.3270  
##hiv.x:age.x   0.4014     0.6705   0.599   0.5511  


################################################################
## # Change 12-0
################################################################

## ## MANOVA ##
##change:
fit <- manova(as.matrix(dall[,102:103]-dall[,100:101])~hiv.x*age.x, data=dall)
summary(fit, test="Wilks")
#            Df   Wilks approx F num Df den Df Pr(>F)
#hiv.x        1 0.99311  0.27058      2     78 0.7637
#age.x        1 0.97184  1.13001      2     78 0.3283
#hiv.x:age.x  1 0.99546  0.17786      2     78 0.8374


## ## HIV Positive vs Negative ##
##Use baseline variables to predict 12-0 changes

## ### np.motor.groovedpeg_d ###
##3 age groups, before AIC:

# remove outliers
#dif <- dall$np.motor.groovedpeg_d.yC - dall$np.motor.groovedpeg_d.xC
#outs <- which(dif < -40)
#dall <- dall[-outs,]

fit <- lm(np.motor.groovedpeg_d.yC-np.motor.groovedpeg_d.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  16.3467    18.3194   0.892    0.375
#hiv.x       -13.6410    25.0619  -0.544    0.588
#age.x        -0.3802     0.3397  -1.119    0.266
#hiv.x:age.x   0.2756     0.4831   0.570    0.570


fit <- lm(np.motor.groovedpeg_d.y-np.motor.groovedpeg_d.x~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)
##(Intercept)  16.3467    18.3194   0.892    0.375
##hiv.x       -13.6410    25.0619  -0.544    0.588
##age.x        -0.3802     0.3397  -1.119    0.266
##hiv.x:age.x   0.2756     0.4831   0.570    0.570


## ### np.motor.groovedpeg_nd ###
##3 age groups, before AIC:
fit <- lm(np.motor.groovedpeg_nd.yC-np.motor.groovedpeg_nd.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  -9.1837    18.2715  -0.503    0.617
#hiv.x        -2.6920    24.9964  -0.108    0.915
#age.x         0.1252     0.3388   0.370    0.713
#hiv.x:age.x   0.0350     0.4819   0.073    0.942

fit <- lm(np.motor.groovedpeg_nd.y-np.motor.groovedpeg_nd.x~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)
##(Intercept)  -9.1837    18.2715  -0.503    0.617
##hiv.x        -2.6920    24.9964  -0.108    0.915
##age.x         0.1252     0.3388   0.370    0.713
##hiv.x:age.x   0.0350     0.4819   0.073    0.942

fname <- "MotChangeGroovedpegNDContCorrectEdu.eps"
postscript(file=fname, width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
yc <-  dall$np.motor.groovedpeg_nd.yC-dall$np.motor.groovedpeg_nd.xC
nonNA <- !is.na(yc)
sel <- dall$hiv.x == 1 & nonNA
plot(dall$age.x[sel], yc[sel], col=mycols[1], xlim=range(dall$age.x), ylim=range(yc[nonNA]), xlab="Age", ylab="GroovedpegND Change", pch=2)
sel <- dall$hiv.x == 0 & nonNA
points(dall$age.x[sel], yc[sel], col=mycols[2])
dnew <- data.frame(age.x=c(40, 75), hiv.x=1)
y1 <- predict(fit, dnew)
lines(c(40, 75), y1, col=mycols[1], lwd=3)
dnew <- data.frame(age.x=c(40, 75), hiv.x=0)
y1 <- predict(fit, dnew)
lines(c(40, 75), y1, col=mycols[2], lwd=3, lty=2)
legend("topright", c("HIV", "Non-HIV"), lty=c(1,2), pch=c(2,1), lwd=3, col=mycols, merge=F, seg.len=2)
dev.off(); system(paste("evince ", fname))


##3 age groups: hiv.x*age.grp
fit <- lm(np.motor.groovedpeg_d.yC-np.motor.groovedpeg_d.xC~hiv.x*age.grp, data=dall)
summary(fit)
#(Intercept)      -2.8000     5.0347  -0.556    0.580
#hiv.x             0.2483     5.8386   0.043    0.966
#age.grp50         1.8000     6.8170   0.264    0.792
#age.grp60        -7.6286     7.8460  -0.972    0.334
#hiv.x:age.grp50  -1.2028     8.1691  -0.147    0.883
#hiv.x:age.grp60   4.8470    12.4416   0.390    0.698

fit <- lm(np.motor.groovedpeg_nd.yC-np.motor.groovedpeg_nd.xC~hiv.x*age.grp, data=dall)
summary(fit)
#(Intercept)      -2.6000     4.9817  -0.522    0.603
#hiv.x            -3.1241     5.7772  -0.541    0.590
#age.grp50         0.5167     6.7453   0.077    0.939
#age.grp60        -0.5429     7.7635  -0.070    0.944
#hiv.x:age.grp50   4.5711     8.0832   0.566    0.573
#hiv.x:age.grp60  -5.0663    12.3108  -0.412    0.682


##2 age groups: hiv.x*age.grp55
fit <- lm(np.motor.groovedpeg_d.yC-np.motor.groovedpeg_d.xC~hiv.x*age.grp55, data=dall)
summary(fit)
##                  Estimate Std. Error t value Pr(>|t|)  
##(Intercept)         -1.300      3.489  -0.373   0.7104  
##hiv.x               -2.550      4.273  -0.597   0.5523  
##age.grp5555         -8.367      6.262  -1.336   0.1854  
##hiv.x:age.grp5555   13.717      7.917   1.732   0.0871 .


fit <- lm(np.motor.groovedpeg_nd.yC-np.motor.groovedpeg_nd.xC~hiv.x*age.grp55, data=dall)
summary(fit)
##                  Estimate Std. Error t value Pr(>|t|)
##(Intercept)         -2.200      3.425  -0.642    0.522
##hiv.x               -4.400      4.194  -1.049    0.297
##age.grp5555         -1.022      6.148  -0.166    0.868
##hiv.x:age.grp5555   11.194      7.773   1.440    0.154