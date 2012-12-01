##% HIV-Aging: Speed Domain
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
aa <- data.frame(as.matrix(dd[,grp.spd+32]-dd[,grp.spd]))
names(aa) <- paste(names(aa), "change", sep="")
dd <- cbind(dd, aa)

##calculate relative changes in the ececutive domain
aa <- data.frame(as.matrix( (dd[,grp.spd+32]-dd[,grp.spd] )/(dd[,grp.spd]+1)   ))
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
aa <- data.frame(as.matrix(dlall[,c(grp.spd)]))
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
aa <- data.frame(as.matrix(dall[,c(grp.spd, grp.spd+32)]))
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
#[1]  84 112

# dall[1,103:111]
#  age.sq np.speed.stroop_c.xC np.speed.trail_a_time.xC np.speed.wais_digsym.xC np.speed.wais_symsrch.xC np.speed.stroop_c.yC
#   2500                   NA                 30.62764                41.66504                 4.102439             53.60014
#  np.speed.trail_a_time.yC np.speed.wais_digsym.yC np.speed.wais_symsrch.yC
#                 41.62764                37.66504                 21.10244


################################################################
##  #Baseline Analysis
################################################################
## ## MANOVA ##
##Baseline:
##Without Age.sq
fit <- manova(as.matrix(dall[,104:107])~hiv.x*age.x, data=dall)
summary(fit, test="Wilks")
#            Df   Wilks approx F num Df den Df  Pr(>F)  
#hiv.x        1 0.94021   0.8108      4     51 0.52412  
#age.x        1 0.81183   2.9552      4     51 0.02849 *
#hiv.x:age.x  1 0.92030   1.1041      4     51 0.36475  
#Residuals   54     



## ## HIV Positive vs Negative, Univariate ##
##Use regression to compare each raw T score. Not considering CD4current because not significant, consistent with past analyses. 

## ### np.speed.stroop_c ###
##3 age groups, before AIC:
fit <- lm(np.speed.stroop_c.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 83.15154   15.85644   5.244 2.68e-06 ***
#hiv.x       -3.51389   23.61772  -0.149    0.882    
#age.x       -0.30645    0.28447  -1.077    0.286    
#hiv.x:age.x -0.04643    0.45553  -0.102    0.919    

postscript(file="SpdcBaselineSpeedStroopC.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.speed.stroop_c.xC~age.grp*hiv.x, data=dall, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="BD")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()

## Uncorrected
fit <- lm(np.speed.stroop_c.x~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)    
##(Intercept)  78.4378    16.6343   4.715 1.74e-05 ***
##hiv.x         1.7158    24.7764   0.069    0.945    
##age.x        -0.1659     0.2984  -0.556    0.581    
##hiv.x:age.x  -0.2043     0.4779  -0.427    0.671  


## ### np.speed.trail_a_time ###
##3 age groups, before AIC:
fit <- lm(np.speed.trail_a_time.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  20.0962    10.6737   1.883   0.0634 .
#hiv.x       -18.4402    14.6028  -1.263   0.2103  
#age.x         0.1761     0.1979   0.890   0.3763  
#hiv.x:age.x   0.4157     0.2815   1.477   0.1436  

## Uncorrected
fit <- lm(np.speed.trail_a_time.x~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  24.2615    10.6888   2.270   0.0259 *
##hiv.x       -21.4684    14.6235  -1.468   0.1460  
##age.x         0.0802     0.1982   0.405   0.6868  
##hiv.x:age.x   0.4928     0.2819   1.748   0.0843 .

fit <- lm(np.speed.trail_a_time.x~hiv.x*age.grp, data=dall)
summary(fit)
##                Estimate Std. Error t value Pr(>|t|)    
##(Intercept)      28.5000     2.9763   9.576 8.29e-15 ***
##hiv.x             1.3621     3.4515   0.395   0.6942    
##age.grp50        -0.6538     3.9588  -0.165   0.8692    
##age.grp60         1.3571     4.6382   0.293   0.7706    
##hiv.x:age.grp50   1.7463     4.7701   0.366   0.7153    
##hiv.x:age.grp60  12.7808     7.3549   1.738   0.0862 .  

fit <- lm(np.speed.trail_a_time.x~hiv.x*age.grp55, data=dall)
summary(fit)
##                  Estimate Std. Error t value Pr(>|t|)    
##(Intercept)       28.52381    1.99384  14.306   <2e-16 ***
##hiv.x              0.30119    2.46221   0.122   0.9029    
##age.grp5555        0.03175    3.64024   0.009   0.9931    
##hiv.x:age.grp5555  8.71468    4.61535   1.888   0.0626 .  


## ### np.speed.wais_digsym ###
##3 age groups, before AIC:
fit <- lm(np.speed.wais_digsym.xC~hiv.x*age.x, data=dall)
summary(fit)
#(Intercept)  91.6204    15.0992   6.068 4.09e-08 ***
#hiv.x        -2.9564    20.6574  -0.143   0.8866    
#age.x        -0.5555     0.2799  -1.984   0.0506 .  
#hiv.x:age.x   0.0505     0.3982   0.127   0.8994    

fit <- lm(np.speed.wais_digsym.xC~hiv.x*age.grp, data=dall)
summary(fit)
##                Estimate Std. Error t value Pr(>|t|)    
##(Intercept)       69.268      4.149  16.693   <2e-16 ***
##hiv.x             -3.497      4.812  -0.727   0.4696    
##age.grp50         -8.812      5.519  -1.597   0.1144    
##age.grp60        -14.653      6.466  -2.266   0.0262 *  
##hiv.x:age.grp50    5.010      6.650   0.753   0.4535    
##hiv.x:age.grp60    5.659     10.254   0.552   0.5826    

fit <- lm(np.speed.wais_digsym.x~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)    
##(Intercept)  76.1247    16.3433   4.658 1.25e-05 ***
##hiv.x         8.3094    22.3594   0.372    0.711    
##age.x        -0.1988     0.3030  -0.656    0.514    
##hiv.x:age.x  -0.2363     0.4310  -0.548    0.585    


## ### np.speed.wais_symsrch ###
##3 age groups, before AIC:
fit <- lm(np.speed.wais_symsrch.xC~hiv.x*age.x, data=dall)
summary(fit)
#(Intercept)  34.0699     7.4291   4.586 1.65e-05 ***
#hiv.x        13.1526    10.1639   1.294    0.199    
#age.x        -0.1142     0.1377  -0.829    0.410    
#hiv.x:age.x  -0.2881     0.1959  -1.470    0.145 

## Uncorrected
fit <- lm(np.speed.wais_symsrch.x~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)    
##(Intercept)  28.1133     7.7371   3.634 0.000492 ***
##hiv.x        17.4831    10.5853   1.652 0.102526    
##age.x         0.0229     0.1434   0.160 0.873544    
##hiv.x:age.x  -0.3983     0.2041  -1.952 0.054430 . 

fit <- lm(np.speed.wais_symsrch.x~hiv.x*age.grp, data=dall)
summary(fit)
##                Estimate Std. Error t value Pr(>|t|)    
##(Intercept)      28.0000     2.1000  13.334   <2e-16 ***
##hiv.x             0.4483     2.4352   0.184    0.854    
##age.grp50         2.9231     2.7932   1.046    0.299    
##age.grp60         0.2857     3.2725   0.087    0.931    
##hiv.x:age.grp50  -4.8259     3.3656  -1.434    0.156    
##hiv.x:age.grp60 -11.4007     5.1894  -2.197    0.031 *  


################################################################
## # Change 12-0
################################################################

## ## MANOVA ##
##change:
fit <- manova(as.matrix(dall[,108:111]-dall[,104:107])~hiv.x*age.x, data=dall)
summary(fit, test="Wilks")
#            Df   Wilks approx F num Df den Df  Pr(>F)  
#hiv.x        1 0.80311  3.12569      4     51 0.02246 *
#age.x        1 0.92756  0.99579      4     51 0.41842  
#hiv.x:age.x  1 0.86749  1.94757      4     51 0.11671  


## ## HIV Positive vs Negative ##
##Use baseline variables to predict 12-0 changes

## ### speed.stroop_c ###
##3 age groups, before AIC:
fit <- lm(np.speed.stroop_c.yC-np.speed.stroop_c.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)  
#(Intercept) 15.54215    9.27143   1.676   0.0995 .
#hiv.x       -4.49986   13.80953  -0.326   0.7458  
#age.x       -0.28194    0.16633  -1.695   0.0958 .
#hiv.x:age.x  0.08613    0.26635   0.323   0.7477 

fname <- "SpdChangeStroopCContCorrectEdu.eps"
postscript(file=fname, width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
yc <-  dall$np.speed.stroop_c.yC-dall$np.speed.stroop_c.xC
nonNA <- !is.na(yc)
sel <- dall$hiv.x == 1 & nonNA
plot(dall$age.x[sel], yc[sel], col=mycols[1], xlim=range(dall$age.x), ylim=range(yc[nonNA]), xlab="Age", ylab="StroopC Change", pch=2)
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

## hiv.x*age.grp
fit <- lm(np.speed.stroop_c.yC-np.speed.stroop_c.xC~hiv.x*age.grp, data=dall)
summary(fit)
#                Estimate Std. Error t value Pr(>|t|)  
#(Intercept)        6.200      2.895   2.142   0.0369 *
#hiv.x             -3.200      3.253  -0.984   0.3298  
#age.grp50         -6.783      3.445  -1.969   0.0543 .
#age.grp60         -9.629      3.790  -2.541   0.0141 *
#hiv.x:age.grp50    3.855      4.131   0.933   0.3551  
#hiv.x:age.grp60   -1.371      7.646  -0.179   0.8583  

## Uncorrected 
fit <- lm(np.speed.stroop_c.y-np.speed.stroop_c.x~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)  
##(Intercept) 15.54215    9.27143   1.676   0.0995 .
##hiv.x       -4.49986   13.80953  -0.326   0.7458  
##age.x       -0.28194    0.16633  -1.695   0.0958 .
##hiv.x:age.x  0.08613    0.26635   0.323   0.7477   

fname <- "SpdChangeStroopCContUncorrectEdu.eps"
postscript(file=fname, width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
yc <-  dall$np.speed.stroop_c.y-dall$np.speed.stroop_c.x
nonNA <- !is.na(yc)
sel <- dall$hiv.x == 1 & nonNA
plot(dall$age.x[sel], yc[sel], col=mycols[1], xlim=range(dall$age.x), ylim=range(yc[nonNA]), xlab="Age", ylab="StroopC Change", pch=2)
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

fit <- lm(np.speed.stroop_c.y-np.speed.stroop_c.x~hiv.x*age.grp, data=dall)
summary(fit)
##                Estimate Std. Error t value Pr(>|t|)  
##(Intercept)        6.200      2.895   2.142   0.0369 *
##hiv.x             -3.200      3.253  -0.984   0.3298  
##age.grp50         -6.783      3.445  -1.969   0.0543 .
##age.grp60         -9.629      3.790  -2.541   0.0141 *
##hiv.x:age.grp50    3.855      4.131   0.933   0.3551  
##hiv.x:age.grp60   -1.371      7.646  -0.179   0.8583  


## ### np.speed.trail_a_time ###
##3 age groups, before AIC:
fit <- lm(np.speed.trail_a_time.yC-np.speed.trail_a_time.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  16.2470    13.5025   1.203    0.232
#hiv.x       -20.8354    18.4729  -1.128    0.263
#age.x        -0.2725     0.2503  -1.088    0.280
#hiv.x:age.x   0.4021     0.3561   1.129    0.262 

fit <- lm(np.speed.trail_a_time.yC-np.speed.trail_a_time.xC~hiv.x*age.grp, data=dall)
summary(fit)
#                Estimate Std. Error t value Pr(>|t|)
#(Intercept)        4.300      3.625   1.186    0.239
#hiv.x             -4.955      4.203  -1.179    0.242
#age.grp50         -2.992      4.821  -0.621    0.537
#age.grp60         -5.443      5.649  -0.964    0.338
#hiv.x:age.grp50    9.466      5.809   1.629    0.107
#hiv.x:age.grp60    2.431      8.957   0.271    0.787

## Uncorrected
fit <- lm(np.speed.trail_a_time.y-np.speed.trail_a_time.x~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)
##(Intercept)  16.2470    13.5025   1.203    0.232
##hiv.x       -20.8354    18.4729  -1.128    0.263
##age.x        -0.2725     0.2503  -1.088    0.280
##hiv.x:age.x   0.4021     0.3561   1.129    0.262


## ### np.speed.wais_digsym ###
##3 age groups, before AIC:
fit <- lm(np.speed.wais_digsym.yC-np.speed.wais_digsym.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  -0.2563     9.1981  -0.028    0.978
#hiv.x        11.0075    12.5840   0.875    0.384
#age.x         0.1012     0.1705   0.593    0.555
#hiv.x:age.x  -0.3564     0.2426  -1.469    0.146

fit <- lm(np.speed.wais_digsym.yC-np.speed.wais_digsym.xC~hiv.x*age.grp, data=dall)
summary(fit)
#                Estimate Std. Error t value Pr(>|t|)  
#(Intercept)       5.2000     2.5178   2.065   0.0422 *
#hiv.x            -5.3724     2.9198  -1.840   0.0696 .
#age.grp50        -0.5077     3.3490  -0.152   0.8799  
#age.grp60         0.6571     3.9238   0.167   0.8674  
#hiv.x:age.grp50  -2.6381     4.0353  -0.654   0.5152  
#hiv.x:age.grp60  -7.8181     6.2220  -1.257   0.2127 

fit <- lm(np.speed.wais_digsym.y-np.speed.wais_digsym.x~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)
##(Intercept)  -0.2563     9.1981  -0.028    0.978
##hiv.x        11.0075    12.5840   0.875    0.384
##age.x         0.1012     0.1705   0.593    0.555
##hiv.x:age.x  -0.3564     0.2426  -1.469    0.146


## ### np.speed.wais_symsrch ###
##3 age groups, before AIC:
fit <- lm(np.speed.wais_symsrch.yC-np.speed.wais_symsrch.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  1.11715    6.35444   0.176    0.861
#hiv.x        7.25553    8.69358   0.835    0.406
#age.x        0.02095    0.11781   0.178    0.859
#hiv.x:age.x -0.16348    0.16758  -0.976    0.332

fit <- lm(np.speed.wais_symsrch.yC-np.speed.wais_symsrch.xC~hiv.x*age.grp, data=dall)
summary(fit)
#                Estimate Std. Error t value Pr(>|t|)
#(Intercept)       1.8000     1.7499   1.029    0.307
#hiv.x             0.2000     2.0293   0.099    0.922
#age.grp50         0.8154     2.3276   0.350    0.727
#age.grp60         0.3429     2.7271   0.126    0.900
#hiv.x:age.grp50  -2.5427     2.8046  -0.907    0.367
#hiv.x:age.grp60   0.3238     4.3244   0.075    0.941


fit <- lm(np.speed.wais_symsrch.y-np.speed.wais_symsrch.x~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)
##(Intercept)  1.11715    6.35444   0.176    0.861
##hiv.x        7.25553    8.69358   0.835    0.406
##age.x        0.02095    0.11781   0.178    0.859
##hiv.x:age.x -0.16348    0.16758  -0.976    0.332