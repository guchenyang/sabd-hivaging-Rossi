##% HIV-Aging: Memory Domain
##% Brown HIV Aging

########### multiple (6+) # will not show
########### one leading space will lead to formatted code block (adding empty lines before and after), except the space before #, +
########### lines starting without ## will be replaced with empty lines
################################################################
## # Version Information #
##+ Version 10: Add analysis only using subjects with 1st visit
##+ Version 9: uncorrected for other clinical factors
##+ Version 7: Remove correcting for age 
##+ Version 6: correcting to have age=50, other covariates=mean, new linear mixed effect age model, retain age for two visits 
##+ Version 5: Continuous age and correction
##+ Version 4: abandoned
##+ Version 3: new way for corrrection
##+ Version 2: add Manova after clinical variables, plots, 100 threhold for nadir
##+ Version 1: Adapted from ExecutiveDomain3
################################################################


################################################################
## # Result Summary #
##+ interaction significant at baseline: bvmt_delay and bvmt_sum
##+ interaction significant at change: hvlt_delay, especially using the variance stablization transformed response, see SQRT change
################################################################

## # Preprocessing #
d <- read.csv("data.rossi.041312.csv")
d$id <- factor(d$id)
##Replace age var with initial age 
## for (ji in unique(d$id)) {
##   sel <- d$id==ji
##   d[sel,]$age <- min(d[sel,]$age)
## }

d <- data.frame(d)
##Raw scores only from now on
##groupings for raw scores
grp.exe <- c(10:14)
grp.mem <- c(17:20)
grp.mot <- c(21:22)
grp.spd <- c(23:26)
##Remove 40minus subjects.
d <- d[d$age>=40,]

##Remove subjects missing either baseline or 12 month visit.
d1 <- d[d$visit==0, ]                  # dim(d1)=132*33
d2 <- d[d$visit==12,]                  # dim(d2)=92*33
dd <- merge(d1, d2, by=("id"))         # dim(dd)=85*65
dl <- rbind(d1, d2)

##calculate changes in the Exe domain
aa <- data.frame(as.matrix(dd[,grp.exe+32]-dd[,grp.exe]))
names(aa) <- paste(names(aa), "change", sep="")
dd <- cbind(dd, aa)
##calculate changes in the memory domain
aa <- data.frame(as.matrix(dd[,grp.mem+32]-dd[,grp.mem]))
names(aa) <- paste(names(aa), "change", sep="")
dd <- cbind(dd, aa)

##calculate relative changes in the memory domain
aa <- data.frame(as.matrix( (dd[,grp.mem+32]-dd[,grp.mem] )/(dd[,grp.mem]+1)   ))
names(aa) <- paste(names(aa), "rchange", sep="")
dd <- cbind(dd, aa)

##calculate sqrt changes in the memory domain
aa <- data.frame(as.matrix( ( sqrt(dd[,grp.mem+32])-sqrt(dd[,grp.mem]) )  ))
names(aa) <- paste(names(aa), "sqchange", sep="")
dd <- cbind(dd, aa)

##calculate sqrt changes in the memory domain
aa <- data.frame(as.matrix( ( log(1+dd[,grp.mem+32])-log(1+dd[,grp.mem]) )  ))
names(aa) <- paste(names(aa), "logchange", sep="")
dd <- cbind(dd, aa)                                                                   # dim(dd)=85*86

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
dd$age.grp55 <- factor(dd$age.grp55)    #dim(dd)=85*88




##Univaraite approaches to yield easier interpretation

##Add clinical information: duration, CD4, CD4 Nadir, HCV, HIV/RAV4, CART
dc <- read.csv("AgeEffectsMASTER031212.csv", sep=";")      # dim(dc)=425*411
sel <- c(1,2,13:22)
dc <- data.frame(dc[,sel])                                 # dim(dc)=425*12
# dc[1,]
#  ID Visit HIVduration HIVRNA_Detectable HIVRNA_value CD4current CD4nadir HAART ARV cART HCVlifetime HCVcurrent
#  1     0          22                 0            0        394      201     1   1    1           1          1

dc1 <- dc[dc$Visit==0,]
dc2 <- dc[dc$Visit==12,]
dcd <- merge(dc1, dc2, by=("ID"))
names(dcd)[1] <- "id"    # dim(dcd)=144*23
################ merge all together, non-hiv=2, hiv>100 = 1, hiv<100 = 0
dall <- merge(dd, dcd, by="id")  # dim(dcd)=84*110
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

# dim(dall)=84*116
################################################################
##Long format 
names(dc)[1:2] <- c("id", "visit")
dl$visit <- factor(dl$visit)   #dl <- rbind(d1, d2)
dc$visit <- factor(dc$visit)   #dc <- data.frame(dc[,sel])
dc$id <- factor(dc$id)
dlall <- merge(dl, dc, by=c("id", "visit")) #dim(dlall)=222*43
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
aa <- data.frame(as.matrix(dlall[,c(grp.mem)]))
for (jj in 1:ncol(aa)) {
  ## visit 0
  vsel <-  dlall$visit == 0 
  sel <- dlall$hiv==0  & vsel
  ##Version 9: correcting for education only, and apply the same correction to the second visit.  
  ##Version 7: not correcting for age 
  ##Version 5: correct for continous age, age=50, education=13, sex=1, or next year 51 for second visit
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
#mytheme()


##Correction for all clinical variables and demographical variables, at baseline and 12month
aa <- data.frame(as.matrix(dall[,c(grp.mem, grp.mem+32)]))
ashift <- (ncol(aa)/2)
for (jj in 1:(ncol(aa)/2) ) {
  sel <- dall$hiv.x==0
  ##Version 9: correcting for education only 
  ##Version 5: correct for continous age, age=50, education=13, sex=1, or next year 51 for change
  tmp <- correctEffect(aa[sel, jj], data.frame(education=dall[sel,]$education.x), data.frame(education=13))
  aa[sel, jj] <- tmp$ce
  
  ##correct this for the second visit
  aa[sel, jj+ashift] <- applyCorrect(tmp$fit, aa[sel, jj+ashift], data.frame(education=dall[sel,]$education.y), data.frame(education=13))

  ##correct this for hiv positives as well
  #sel <- dall$hiv.x==0, here should be 1
  sel <- dall$hiv.x==1
  aa[sel, jj] <- applyCorrect(tmp$fit, aa[sel, jj], data.frame(education=dall[sel,]$education.x), data.frame(education=13))

  aa[sel, jj+ashift] <- applyCorrect(tmp$fit, aa[sel, jj+ashift], data.frame(education=dall[sel,]$education.y), data.frame(education=13))
}

names(aa) <- paste(names(aa), "C", sep="")
dall <- cbind(dall, aa)

##Add age^2 terms
dall <- data.frame(dall, age2 = (dall$age.x)^2)


################################################################
##  #Baseline Analysis
################################################################
## ## MANOVA ##
##Baseline:
##Without Age.sq
fit <- manova(as.matrix(dall[,117:120])~hiv.x*age.x, data=dall)
summary(fit, test="Wilks")
##            Df   Wilks approx F num Df den Df Pr(>F)
##hiv.x        1 0.98794  0.23202      4     76 0.9196
##age.x        1 0.99450  0.10502      4     76 0.9804
##hiv.x:age.x  1 0.97458  0.49557      4     76 0.7390



## ## HIV Positive vs Negative, Univariate ##
##Use regression to compare each raw T score. Not considering CD4current because not significant, consistent with past analyses.  Also not selected by AIC.  Consider different ways to dictomizing/transforming age, including log, two groups, and square.  Square seems to work the best.

##HCV current/lifetime don't change the answer much.  Stick to current in the following analysis.

## ### np.learnmem.bvmt_delay ###
##3 age groups, before AIC:
fit <- lm(np.learnmem.bvmt_delay.xC~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  7.18417    3.25774   2.205   0.0303 *
##hiv.x        3.81574    4.45695   0.856   0.3945  
##age.x        0.01413    0.06040   0.234   0.8157  
##hiv.x:age.x -0.06633    0.08591  -0.772   0.4423   

postscript(file="MemBaselineBvmtDelay.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.learnmem.bvmt_delay.xC~age.grp*hiv.x, data=dall, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="BD")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()

## ### np.learnmem.bvmt_sum ###
##3 age groups, before AIC:
fit <- lm(np.learnmem.bvmt_sum.xC~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)  
##(Intercept) 16.47227    7.31564   2.252   0.0271 *
##hiv.x       13.28566   10.00861   1.327   0.1881  
##age.x        0.07099    0.13563   0.523   0.6022  
##hiv.x:age.x -0.25501    0.19293  -1.322   0.1900  

postscript(file="MemBaselineBvmtSum.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.learnmem.bvmt_sum.xC~age.grp*hiv.x, data=dall, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="BS")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()



## ### np.learnmem.hvlt_delay ###
##3 age groups, before AIC:
fit <- lm(np.learnmem.hvlt_delay.xC~hiv.x*age.x , data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  7.09250    2.75078   2.578   0.0118 *
##hiv.x        2.31028    3.77322   0.612   0.5421  
##age.x        0.01679    0.05100   0.329   0.7428  
##hiv.x:age.x -0.04660    0.07283  -0.640   0.5241    



## ### np.learnmem.hvlt_sum ###
##3 age groups, before AIC:
fit <- lm(np.learnmem.hvlt_sum.xC~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)    
##(Intercept) 22.46230    5.24647   4.281 5.12e-05 ***
##hiv.x        3.58627    7.17776   0.500    0.619    
##age.x        0.01531    0.09727   0.157    0.875    
##hiv.x:age.x -0.06517    0.13836  -0.471    0.639  



################################################################
## # Change 12-0
################################################################

## ## MANOVA ##
##change:
fit <- manova(as.matrix(dall[,121:124]-dall[,117:120])~hiv.x*age.x, data=dall)
summary(fit, test="Wilks")
##            Df   Wilks approx F num Df den Df Pr(>F)
##hiv.x        1 0.95486  0.89823      4     76 0.4694
##age.x        1 0.92156  1.61729      4     76 0.1786
##hiv.x:age.x  1 0.94568  1.09140      4     76 0.3669
##Residuals   79    


## ## HIV Positive vs Negative ##
##Use baseline variables to predict 12-0 changes

## ### np.learnmem.bvmt_delay ###
##3 age groups, before AIC:
fit <- lm(np.learnmem.bvmt_delay.yC-np.learnmem.bvmt_delay.xC~hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)  2.27287    2.73385   0.831    0.408
## hiv.x        3.62412    3.74021   0.969    0.335
## age.x       -0.02953    0.05069  -0.583    0.562
## hiv.x:age.x -0.08500    0.07210  -1.179    0.242


## ### np.learnmem.bvmt_sum ###
##3 age groups, before AIC:
fit <- lm(np.learnmem.bvmt_sum.yC-np.learnmem.bvmt_sum.xC~hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)   8.4627     6.6411   1.274    0.206
## hiv.x         1.7591     9.0858   0.194    0.847
## age.x        -0.1332     0.1231  -1.082    0.283
## hiv.x:age.x  -0.0640     0.1751  -0.365    0.716

## ### np.learnmem.hvlt_delay ###
##continuous age, uncorrected
fit <- lm(np.learnmem.hvlt_delay.y-np.learnmem.hvlt_delay.x~hiv.x*age.x, data=dall)
summary(fit)
##              Estimate Std. Error t value Pr(>|t|)  
## (Intercept)  0.100440   2.650568   0.038   0.9699  
## hiv.x        5.806967   3.635758   1.597   0.1142  
## age.x        0.009378   0.049142   0.191   0.8491  
## hiv.x:age.x -0.121971   0.070175  -1.738   0.0861 .
fname <- "MemChangeHvltDelayContUncorrected.eps"
postscript(file=fname, width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
yc <-  dall$np.learnmem.hvlt_delay.y-dall$np.learnmem.hvlt_delay.x
nonNA <- !is.na(yc)
sel <- dall$hiv.x == 1 & nonNA
plot(dall$age.x[sel], yc[sel], col=mycols[1], xlim=range(dall$age.x), ylim=range(yc[nonNA]), xlab="Age", ylab="HVLT Delay Change", pch=2)
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



fit <- lm(np.learnmem.hvlt_delay.yC-np.learnmem.hvlt_delay.xC~hiv.x*age.x, data=dall)
summary(fit)
##              Estimate Std. Error t value Pr(>|t|)  
## (Intercept)  0.100440   2.650568   0.038   0.9699  
## hiv.x        5.806967   3.635758   1.597   0.1142  
## age.x        0.009378   0.049142   0.191   0.8491  
## hiv.x:age.x -0.121971   0.070175  -1.738   0.0861 .
fname <- "MemChangeHvltDelayContCorrectEdu.eps"
postscript(file=fname, width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
yc <-  dall$np.learnmem.hvlt_delay.yC-dall$np.learnmem.hvlt_delay.xC
nonNA <- !is.na(yc)
sel <- dall$hiv.x == 1 & nonNA
plot(dall$age.x[sel], yc[sel], col=mycols[1], xlim=range(dall$age.x), ylim=range(yc[nonNA]), xlab="Age", ylab="HVLT Delay Change", pch=2)
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



fit <- lm(np.learnmem.hvlt_delay.yC-np.learnmem.hvlt_delay.xC~hiv.x*age.grp, data=dall)
summary(fit)
##                 Estimate Std. Error t value Pr(>|t|)  
## (Intercept)       0.5000     0.7157   0.699   0.4869  
## hiv.x             0.6034     0.8299   0.727   0.4694  
## age.grp50         0.1154     0.9519   0.121   0.9038  
## age.grp60         0.2143     1.1153   0.192   0.8481  
## hiv.x:age.grp50  -1.5522     1.1518  -1.348   0.1817  
## hiv.x:age.grp60  -3.3177     1.7685  -1.876   0.0644 .
fname <- "MemChangeHvltDelayTrendErrorBarCorrectEdu.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
par(mfrow = (c(1,2)))
##non hiv plot
plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delayC)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##non-hiv
  yh <- dlall$np.learnmem.hvlt_delayC
  age <- dlall$age.grp
  visit <- dlall$visit
  hiv <- dlall$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==0], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==0 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==0]),  std.err(yh[visit == 12 & age == gc & hiv==0 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.n[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.n[gg])
}
legend("topleft", c("Non-HIV 40", "Non-HIV 50", "Non-HIV 60"), col = col.n, lty=1, merge=F, lwd=2.5)

plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delayC)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##hiv
  yh <- dlall$np.learnmem.hvlt_delayC
  age <- dlall$age.grp
  visit <- dlall$visit
  hiv <- dlall$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==1], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==1 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==1]),  std.err(yh[visit == 12 & age == gc & hiv==1 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.p[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.p[gg])
}
legend("topleft", c("HIV 40", "HIV 50", "HIV 60"), col = col.p, lty=1, merge=F, lwd=2.5)
dev.off(); system(paste("evince ", fname))



fit <- lm(np.learnmem.hvlt_delay.y-np.learnmem.hvlt_delay.x~hiv.x*age.grp, data=dall)
summary(fit)
##                 Estimate Std. Error t value Pr(>|t|)  
## (Intercept)       0.5000     0.7157   0.699   0.4869  
## hiv.x             0.6034     0.8299   0.727   0.4694  
## age.grp50         0.1154     0.9519   0.121   0.9038  
## age.grp60         0.2143     1.1153   0.192   0.8481  
## hiv.x:age.grp50  -1.5522     1.1518  -1.348   0.1817  
## hiv.x:age.grp60  -3.3177     1.7685  -1.876   0.0644 .
fname <- "MemChangeHvltDelayTrendErrorBarUncorrected.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
par(mfrow = (c(1,2)))
##non hiv plot
plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delay)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##non-hiv
  yh <- dlall$np.learnmem.hvlt_delay
  age <- dlall$age.grp
  visit <- dlall$visit
  hiv <- dlall$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==0], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==0 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==0]),  std.err(yh[visit == 12 & age == gc & hiv==0 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.n[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.n[gg])
}
legend("topleft", c("Non-HIV 40", "Non-HIV 50", "Non-HIV 60"), col = col.n, lty=1, merge=F, lwd=2.5)

plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delay)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##hiv
  yh <- dlall$np.learnmem.hvlt_delay
  age <- dlall$age.grp
  visit <- dlall$visit
  hiv <- dlall$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==1], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==1 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==1]),  std.err(yh[visit == 12 & age == gc & hiv==1 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.p[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.p[gg])
}
legend("topleft", c("HIV 40", "HIV 50", "HIV 60"), col = col.p, lty=1, merge=F, lwd=2.5)
dev.off(); system(paste("evince ", fname))


## ## Remove too old controls ## 
## remove >65 controls, actually increase the performance of the 60s group, not used
dallnew <-  dall[dall$age.x<=65, ]
fit <- lm(np.learnmem.hvlt_delay.yC-np.learnmem.hvlt_delay.xC~hiv.x*age.grp, data=dallnew)
summary(fit)

selid <- dall[dall$age.x<=65, ]$id
dlallnew <- dlall[!dlall$id %in% c(140, 152, 179, 197), ]

fname <- "MemChangeHvltDelayTrendErrorBarCorrectEdu40_65.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
par(mfrow = (c(1,2)))
##non hiv plot
plot(c(0, 12), range(na.omit(dlallnew$np.learnmem.hvlt_delayC)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##non-hiv
  yh <- dlallnew$np.learnmem.hvlt_delayC
  age <- dlallnew$age.grp
  visit <- dlallnew$visit
  hiv <- dlallnew$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==0], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==0 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==0]),  std.err(yh[visit == 12 & age == gc & hiv==0 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.n[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.n[gg])
}
legend("topleft", c("Non-HIV 40", "Non-HIV 50", "Non-HIV 60"), col = col.n, lty=1, merge=F, lwd=2.5)

plot(c(0, 12), range(na.omit(dlallnew$np.learnmem.hvlt_delayC)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##hiv
  yh <- dlallnew$np.learnmem.hvlt_delayC
  age <- dlallnew$age.grp
  visit <- dlallnew$visit
  hiv <- dlallnew$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==1], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==1 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==1]),  std.err(yh[visit == 12 & age == gc & hiv==1 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.p[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.p[gg])
}
legend("topleft", c("HIV 40", "HIV 50", "HIV 60"), col = col.p, lty=1, merge=F, lwd=2.5)
dev.off(); system(paste("evince ", fname))





##Two groups 
fit <- lm(np.learnmem.hvlt_delay.y-np.learnmem.hvlt_delay.x~hiv.x*age.grp55, data=dall)
summary(fit)
##                   Estimate Std. Error t value Pr(>|t|)  
## (Intercept)         0.4762     0.4811   0.990   0.3253  
## hiv.x               0.4488     0.5941   0.755   0.4522  
## age.grp5555         0.4127     0.8784   0.470   0.6398  
## hiv.x:age.grp5555  -2.7223     1.1256  -2.419   0.0179 *




postscript(file="MemChangeHvltDelay.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.learnmem.hvlt_delay.yC-np.learnmem.hvlt_delay.xC~age.grp*hiv.x, data=dall, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="HD Change")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()





fname <- "MemChangeHvltDelayTrend.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
pd <- na.omit(data.frame(dall$np.learnmem.hvlt_delay.xC, dall$np.learnmem.hvlt_delay.yC, dall$hiv.x, dall$age.grp))
mytheme()
par(oma=c(3,5,3,3), mfrow=c(2,3))
par(las=1, oma=c(9,5,1,5), pty="m")
par(mar=c(0,0,0,0),mgp=c(2,0.25,0),tcl=-0.15)
par(ps=12)
##Non-HIV
sel <- pd[,3]==0
yrange <- range(pd[sel,1:2])
par(mar=c(0,0,0,0))
sel <- pd[,3] ==  0 & pd[,4] == "40"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.n[1])
points(cbind(1, pd[sel,2]), type="p", col=col.n[1])
segments(rep(0, sum(sel)),  pd[sel,1], rep(1, sum(sel)),  pd[sel,2], col=col.n[1]!80, lty=2)

box("plot", col="black")
axis(2)
mtext("Non-HIV", side=2, line=2, outer=T, adj=0.75, cex=par("cex.lab"), las=0)
par(mar=c(0,0,0,0))
sel <- pd[,3] ==  0 & pd[,4] == "50"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.n[2])
points(cbind(1, pd[sel,2]), type="p", col=col.n[2])
box("plot", col="black")
par(mar=c(0,0,0,0))
sel <- pd[,3] ==  0 & pd[,4] == "60"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.n[3])
points(cbind(1, pd[sel,2]), type="p", col=col.n[3])
box("plot", col="black")
axis(4)
##HIV
sel <- pd[,3]==1
yrange <- range(pd[sel,1:2])
sel <- pd[,3] ==  1 & pd[,4] == "40"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.p[1])
points(cbind(1, pd[sel,2]), type="p", col=col.p[1])
box("plot", col="black")
axis(2); axis(1, at=c(0, 1))
sel <- pd[,3] ==  1 & pd[,4] == "50"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.p[2])
points(cbind(1, pd[sel,2]), type="p", col=col.p[2])
box("plot", col="black")
axis(1, at=c(0, 1))
sel <- pd[,3] ==  1 & pd[,4] == "60"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.p[3])
points(cbind(1, pd[sel,2]), type="p", col=col.p[3])
box("plot", col="black")
axis(4); axis(1, at=c(0, 1))
mtext("HIV", side=2, line=2, outer=T, adj=0.25, cex=par("cex.lab"), las=0)
mtext("40", side=1, line=3, outer=T, adj=0.16, cex=par("cex.lab"), las=0)
mtext("50", side=1, line=3, outer=T, adj=0.5, cex=par("cex.lab"), las=0)
mtext("60", side=1, line=3, outer=T, adj=0.84, cex=par("cex.lab"), las=0)
dev.off(); system(paste("evince ", fname))




fname <- "MemChangeHvltDelayTrend.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
pd <- na.omit(data.frame(dall$np.learnmem.hvlt_delay.xC, dall$np.learnmem.hvlt_delay.yC, dall$hiv.x, dall$age.grp))
mytheme()
par(oma=c(3,5,3,3), mfrow=c(2,3))
par(las=1, oma=c(9,5,1,5), pty="m")
par(mar=c(0,0,0,0),mgp=c(2,0.25,0),tcl=-0.15)
par(ps=12)
##Non-HIV
sel <- pd[,3]==0
yrange <- range(pd[sel,1:2])
par(mar=c(0,0,0,0))
sel <- pd[,3] ==  0 & pd[,4] == "40"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.n[1])
points(cbind(1, pd[sel,2]), type="p", col=col.n[1])
segments(rep(0, sum(sel)),  pd[sel,1], rep(1, sum(sel)),  pd[sel,2], col=col.n[1]!80, lty=2)

box("plot", col="black")
axis(2)
mtext("Non-HIV", side=2, line=2, outer=T, adj=0.75, cex=par("cex.lab"), las=0)
par(mar=c(0,0,0,0))
sel <- pd[,3] ==  0 & pd[,4] == "50"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.n[2])
points(cbind(1, pd[sel,2]), type="p", col=col.n[2])
box("plot", col="black")
par(mar=c(0,0,0,0))
sel <- pd[,3] ==  0 & pd[,4] == "60"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.n[3])
points(cbind(1, pd[sel,2]), type="p", col=col.n[3])
box("plot", col="black")
axis(4)
##HIV
sel <- pd[,3]==1
yrange <- range(pd[sel,1:2])
sel <- pd[,3] ==  1 & pd[,4] == "40"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.p[1])
points(cbind(1, pd[sel,2]), type="p", col=col.p[1])
box("plot", col="black")
axis(2); axis(1, at=c(0, 1))
sel <- pd[,3] ==  1 & pd[,4] == "50"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.p[2])
points(cbind(1, pd[sel,2]), type="p", col=col.p[2])
box("plot", col="black")
axis(1, at=c(0, 1))
sel <- pd[,3] ==  1 & pd[,4] == "60"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.p[3])
points(cbind(1, pd[sel,2]), type="p", col=col.p[3])
box("plot", col="black")
axis(4); axis(1, at=c(0, 1))
mtext("HIV", side=2, line=2, outer=T, adj=0.25, cex=par("cex.lab"), las=0)
mtext("40", side=1, line=3, outer=T, adj=0.16, cex=par("cex.lab"), las=0)
mtext("50", side=1, line=3, outer=T, adj=0.5, cex=par("cex.lab"), las=0)
mtext("60", side=1, line=3, outer=T, adj=0.84, cex=par("cex.lab"), las=0)
dev.off(); system(paste("evince ", fname))



## ### np.learnmem.hvlt_sum ###
##3 age groups, before AIC:
fit <- lm(np.learnmem.hvlt_sum.yC-np.learnmem.hvlt_sum.xC~hiv.x*age.x , data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)
## (Intercept)  3.26818    4.66880   0.700    0.486
## hiv.x        2.82656    6.38744   0.443    0.659
## age.x       -0.03369    0.08656  -0.389    0.698
## hiv.x:age.x -0.05570    0.12313  -0.452    0.652



################################################################
## # Mixed effects #
################################################################
library(lme4)
## ## Overall repeated meausres  ##
## ### np.learnmem.bvmt_delayC ###
summary(lmer(np.learnmem.bvmt_delayC~(1|id) + hiv*age + hiv*visit, data=dlall))
##               Estimate Std. Error t value
## (Intercept)  7.5011783  2.7105690   2.767
## hiv          9.9387720  3.6855934   2.697
## age          0.0002164  0.0501302   0.004
## visit12      1.2122161  0.3917784   3.094
## hiv:age     -0.1248662  0.0708543  -1.762
## hiv:visit12 -4.3725678  0.4907209  -8.910


## ### np.learnmem.bvmt_sumC ###
summary(lmer(np.learnmem.bvmt_sumC~(1|id) + hiv*age+hiv*visit, data=dlall))
##             Estimate Std. Error t value
## (Intercept) 20.83270    6.72326   3.099
## hiv         18.53673    9.14278   2.027
## age         -0.04122    0.12429  -0.332
## visit12      1.89995    1.04108   1.825
## hiv:age     -0.26574    0.17571  -1.512
## hiv:visit12 -2.73220    1.30329  -2.096


## ### np.learnmem.hvlt_delayC ###
summary(lmer(np.learnmem.hvlt_delayC~(1|id)+hiv*visit+hiv*age, data=dlall))
##             Estimate Std. Error t value
## (Intercept)  6.75846    2.51937   2.683
## hiv          5.32337    3.42769   1.553
## visit12      0.82634    0.41162   2.008
## age          0.01337    0.04656   0.287
## hiv:visit12  2.06612    0.51669   3.999
## hiv:age     -0.10334    0.06589  -1.568


## ### np.learnmem.hvlt_sumC ###
summary(lmer(np.learnmem.hvlt_sumC~(1|id) + hiv*age + hiv*visit, data=dlall))
##             Estimate Std. Error t value
## (Intercept) 24.07022    4.67522   5.148
## hiv          4.34707    6.35801   0.684
## age         -0.03778    0.08641  -0.437
## visit12      1.51408    0.75088   2.016
## hiv:age     -0.07257    0.12217  -0.594
## hiv:visit12  1.13213    0.93976   1.205


## ## change ##
## ### np.learnmem.bvmt_delayC ###
summary(lmer(log(np.learnmem.bvmt_delayC)~(1|id) +  hiv*age*visit, data=dlall))
##                  Estimate Std. Error t value
## (Intercept)      1.555065   0.536254   2.900
## hiv              1.380255   0.733286   1.882
## age              0.006339   0.009942   0.638
## visit12          0.532607   0.514643   1.035
## hiv:age         -0.017329   0.014134  -1.226
## hiv:visit12      0.202978   0.701556   0.289
## age:visit12     -0.005591   0.009447  -0.592
## hiv:age:visit12 -0.017232   0.013372  -1.289


## ### np.learnmem.bvmt_sumC ###
summary(lmer(np.learnmem.bvmt_sumC~(1|id) + hiv*age*visit, data=dlall))
##                 Estimate Std. Error t value
## (Intercept)     17.76669    7.44931   2.385
## hiv             19.78168   10.18467   1.942
## age              0.01634    0.13811   0.118
## visit12          8.22017    6.69768   1.227
## hiv:age         -0.28643    0.19631  -1.459
## hiv:visit12     -5.47834    9.13011  -0.600
## age:visit12     -0.11745    0.12294  -0.955
## hiv:age:visit12  0.04593    0.17402   0.264

## ### np.learnmem.hvlt_delayC ###
##Very similar results to the regression estimates, 
summary(lmer(np.learnmem.hvlt_delayC~(1|id) + hiv*age*visit, data=dlall))
##                  Estimate Std. Error t value
## (Intercept)      7.373245   2.815450   2.619
## hiv              2.476428   3.856633   0.642
## age              0.001828   0.052199   0.035
## visit12         -0.445647   2.611894  -0.171
## hiv:age         -0.046456   0.074403  -0.624
## hiv:visit12      7.676671   3.567421   2.152
## age:visit12      0.023637   0.047944   0.493
## hiv:age:visit12 -0.110618   0.068068  -1.625

## ### np.learnmem.hvlt_sumC ###
summary(lmer(np.learnmem.hvlt_sumC~(1|id) + hiv*age*visit, data=dlall))
##                 Estimate Std. Error t value
## (Intercept)     23.23462    5.22389   4.448
## hiv              3.86930    7.14271   0.542
## age             -0.02209    0.09685  -0.228
## visit12          3.24632    4.85614   0.668
## hiv:age         -0.06166    0.13768  -0.448
## hiv:visit12      1.98726    6.61981   0.300
## age:visit12     -0.03219    0.08914  -0.361
## hiv:age:visit12 -0.01958    0.12618  -0.155

## ## age group ##
summary(lmer(np.learnmem.hvlt_delayC~(1|id) + hiv*age.grp*visit, data=dlall))
##                       Estimate Std. Error t value
## (Intercept)            7.48265    0.78107   9.580
## hiv                    0.06903    0.90578   0.076
## age.grp50             -0.16629    1.03892  -0.160
## age.grp60              0.25730    1.21721   0.211
## visit12                0.54714    0.69808   0.784
## hiv:age.grp50          0.42177    1.25518   0.336
## hiv:age.grp60         -0.36992    1.93016  -0.192
## hiv:visit12            2.82331    0.80955   3.488
## age.grp50:visit12      0.38663    0.92854   0.416
## age.grp60:visit12      0.53968    1.08789   0.496
## hiv:age.grp50:visit12 -1.48161    1.12258  -1.320
## hiv:age.grp60:visit12 -3.26789    1.72509  -1.894

mytheme()
fit <- lmer(np.learnmem.hvlt_delayC~(1|id) + hiv*age.grp*visit, data=dlall)
visitx <- as.numeric(as.character(dlall$visit))
age.grpx <- as.numeric(dlall$age.grp)
uid.n <- dall$id[dall$hiv.x == 0 ]
uid.p <-  dall$id[dall$hiv.x == 1 ]
bb <- as.numeric(fit@fixef)
getModelVal <- function(bb, hiv, age.grp, visit) {
  if (age.grp=="40") {
    ad <- c(0, 0)
  } else if (age.grp=="50") {
    ad <- c(1, 0)
  } else {
    ad <- c(0, 1)
  }
  dc <- c(1, hiv, ad, visit, hiv*ad, hiv*visit, ad*visit, hiv*visit*ad)
  sum(bb*dc)
}


fname <- "MemChangeHvltDelayTrendGroupOverlay.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
par(mfrow = (c(1,2)))
## non hiv plot
plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delayC)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (jj in uid.n) {
  sel <- dlall$id == jj
  age.x <- min(age.grpx[sel])
  hiv <- 0
  ## month 0
  yy <- dlall$np.learnmem.hvlt_delayC[sel]
  yv <- visitx[sel]
  ## ccols <- ifelse(rep(hiv, 3), col.p, col.n)
  print(age.x)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.n[age.x], lwd=0.5)
  points(yv, yy, col=col.n[age.x])
}

for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ## non-hiv
  yy <-  c( getModelVal(bb, 0, gc, 0), getModelVal(bb, 0, gc, 1))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.n[gg], lwd=10)
}
legend("topleft", c("Non-HIV 40", "Non-HIV 50", "Non-HIV 60"), col = col.n, lty=1, pch=1, merge=F, lwd=2.5)

## hiv plot
plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delayC)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (jj in uid.p) {
  sel <- dlall$id == jj
  age.x <- min(age.grpx[sel])
  hiv <- 0
  ## month 0
  yy <- dlall$np.learnmem.hvlt_delayC[sel]
  yv <- visitx[sel]
  ## ccols <- ifelse(rep(hiv, 3), col.p, col.n)
  print(age.x)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.p[age.x], lwd=0.5)
  points(yv, yy, col=col.p[age.x], pch=2)
}

for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ## hiv
  yy <-  c( getModelVal(bb, 1, gc, 0), getModelVal(bb, 1, gc, 1))
  yv <-  c(0, 12)
  segments(yv[1], yy[1],  yv[2], yy[2], col=col.p[gg], lwd=10)
}
legend("topleft", c("HIV 40", "HIV 50", "HIV 60"), col = col.p, lty=1, pch=2, merge=F, lwd=2.5)
dev.off(); system(paste("evince ", fname))

################################################################
## # Simple plot # 
################################################################


##Uncorrected, 3 age groups 
fname <- "MemChangeHvltDelayTrendErrorBarUncorrected.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
par(mfrow = (c(1,2)))
error.bars <- function (x, upper, lower, width = 0.02, ...) 
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}

std.err <- function(x) {
  x <- na.omit(x)
  sd(x)/sqrt(length(x))
}
##non hiv plot
plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delay)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##non-hiv
  yh <- dlall$np.learnmem.hvlt_delay
  age <- dlall$age.grp
  visit <- dlall$visit
  hiv <- dlall$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==0], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==0 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==0]),  std.err(yh[visit == 12 & age == gc & hiv==0 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.n[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.n[gg])
}
legend("topleft", c("Non-HIV 40", "Non-HIV 50", "Non-HIV 60"), col = col.n, lty=1, merge=F, lwd=2.5)

plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delay)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##hiv
  yh <- dlall$np.learnmem.hvlt_delayC
  age <- dlall$age.grp
  visit <- dlall$visit
  hiv <- dlall$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==1], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==1 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==1]),  std.err(yh[visit == 12 & age == gc & hiv==1 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.p[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.p[gg])
}
legend("topleft", c("HIV 40", "HIV 50", "HIV 60"), col = col.p, lty=1, merge=F, lwd=2.5)

dev.off(); system(paste("evince ", fname))

##                   Estimate Std. Error t value Pr(>|t|)  
## (Intercept)         0.4762     0.4811   0.990   0.3253  
## hiv.x               0.4488     0.5941   0.755   0.4522  
## age.grp5555         0.4127     0.8784   0.470   0.6398  
## hiv.x:age.grp5555  -2.7223     1.1256  -2.419   0.0179 *


fname <- "MemChangeHvltDelayCont.eps"
postscript(file=fname, width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
yc <-  dall$np.learnmem.hvlt_delay.yC-dall$np.learnmem.hvlt_delay.xC
nonNA <- !is.na(yc)
sel <- dall$hiv.x == 1 & nonNA
plot(dall$age.x[sel], yc[sel], col=mycols[1], xlim=range(dall$age.x), ylim=range(yc[nonNA]), xlab="Age", ylab="HVLT Delay Change", pch=2)
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


fname <- "MemChangeHvltDelayContUncorrected.eps"
postscript(file=fname, width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
yc <-  dall$np.learnmem.hvlt_delay.y-dall$np.learnmem.hvlt_delay.x
nonNA <- !is.na(yc)
sel <- dall$hiv.x == 1 & nonNA
plot(dall$age.x[sel], yc[sel], col=mycols[1], xlim=range(dall$age.x), ylim=range(yc[nonNA]), xlab="Age", ylab="HVLT Delay Change", pch=2)
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



postscript(file="MemChangeHvltDelay.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.learnmem.hvlt_delay.yC-np.learnmem.hvlt_delay.xC~age.grp*hiv.x, data=dall, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="HD Change")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()





fname <- "MemChangeHvltDelayTrend.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
pd <- na.omit(data.frame(dall$np.learnmem.hvlt_delay.xC, dall$np.learnmem.hvlt_delay.yC, dall$hiv.x, dall$age.grp))
mytheme()
par(oma=c(3,5,3,3), mfrow=c(2,3))
par(las=1, oma=c(9,5,1,5), pty="m")
par(mar=c(0,0,0,0),mgp=c(2,0.25,0),tcl=-0.15)
par(ps=12)
##Non-HIV
sel <- pd[,3]==0
yrange <- range(pd[sel,1:2])
par(mar=c(0,0,0,0))
sel <- pd[,3] ==  0 & pd[,4] == "40"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.n[1])
points(cbind(1, pd[sel,2]), type="p", col=col.n[1])
segments(rep(0, sum(sel)),  pd[sel,1], rep(1, sum(sel)),  pd[sel,2], col=col.n[1]!80, lty=2)

box("plot", col="black")
axis(2)
mtext("Non-HIV", side=2, line=2, outer=T, adj=0.75, cex=par("cex.lab"), las=0)
par(mar=c(0,0,0,0))
sel <- pd[,3] ==  0 & pd[,4] == "50"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.n[2])
points(cbind(1, pd[sel,2]), type="p", col=col.n[2])
box("plot", col="black")
par(mar=c(0,0,0,0))
sel <- pd[,3] ==  0 & pd[,4] == "60"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.n[3])
points(cbind(1, pd[sel,2]), type="p", col=col.n[3])
box("plot", col="black")
axis(4)
##HIV
sel <- pd[,3]==1
yrange <- range(pd[sel,1:2])
sel <- pd[,3] ==  1 & pd[,4] == "40"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.p[1])
points(cbind(1, pd[sel,2]), type="p", col=col.p[1])
box("plot", col="black")
axis(2); axis(1, at=c(0, 1))
sel <- pd[,3] ==  1 & pd[,4] == "50"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.p[2])
points(cbind(1, pd[sel,2]), type="p", col=col.p[2])
box("plot", col="black")
axis(1, at=c(0, 1))
sel <- pd[,3] ==  1 & pd[,4] == "60"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.p[3])
points(cbind(1, pd[sel,2]), type="p", col=col.p[3])
box("plot", col="black")
axis(4); axis(1, at=c(0, 1))
mtext("HIV", side=2, line=2, outer=T, adj=0.25, cex=par("cex.lab"), las=0)
mtext("40", side=1, line=3, outer=T, adj=0.16, cex=par("cex.lab"), las=0)
mtext("50", side=1, line=3, outer=T, adj=0.5, cex=par("cex.lab"), las=0)
mtext("60", side=1, line=3, outer=T, adj=0.84, cex=par("cex.lab"), las=0)
dev.off(); system(paste("evince ", fname))




fname <- "MemChangeHvltDelayTrend.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
pd <- na.omit(data.frame(dall$np.learnmem.hvlt_delay.xC, dall$np.learnmem.hvlt_delay.yC, dall$hiv.x, dall$age.grp))
mytheme()
par(oma=c(3,5,3,3), mfrow=c(2,3))
par(las=1, oma=c(9,5,1,5), pty="m")
par(mar=c(0,0,0,0),mgp=c(2,0.25,0),tcl=-0.15)
par(ps=12)
##Non-HIV
sel <- pd[,3]==0
yrange <- range(pd[sel,1:2])
par(mar=c(0,0,0,0))
sel <- pd[,3] ==  0 & pd[,4] == "40"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.n[1])
points(cbind(1, pd[sel,2]), type="p", col=col.n[1])
segments(rep(0, sum(sel)),  pd[sel,1], rep(1, sum(sel)),  pd[sel,2], col=col.n[1]!80, lty=2)

box("plot", col="black")
axis(2)
mtext("Non-HIV", side=2, line=2, outer=T, adj=0.75, cex=par("cex.lab"), las=0)
par(mar=c(0,0,0,0))
sel <- pd[,3] ==  0 & pd[,4] == "50"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.n[2])
points(cbind(1, pd[sel,2]), type="p", col=col.n[2])
box("plot", col="black")
par(mar=c(0,0,0,0))
sel <- pd[,3] ==  0 & pd[,4] == "60"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.n[3])
points(cbind(1, pd[sel,2]), type="p", col=col.n[3])
box("plot", col="black")
axis(4)
##HIV
sel <- pd[,3]==1
yrange <- range(pd[sel,1:2])
sel <- pd[,3] ==  1 & pd[,4] == "40"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.p[1])
points(cbind(1, pd[sel,2]), type="p", col=col.p[1])
box("plot", col="black")
axis(2); axis(1, at=c(0, 1))
sel <- pd[,3] ==  1 & pd[,4] == "50"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.p[2])
points(cbind(1, pd[sel,2]), type="p", col=col.p[2])
box("plot", col="black")
axis(1, at=c(0, 1))
sel <- pd[,3] ==  1 & pd[,4] == "60"
plot(c(0,1), yrange, type="n", ylab="", axes=F)
points(cbind(0, pd[sel,1]), type="p", col=col.p[3])
points(cbind(1, pd[sel,2]), type="p", col=col.p[3])
box("plot", col="black")
axis(4); axis(1, at=c(0, 1))
mtext("HIV", side=2, line=2, outer=T, adj=0.25, cex=par("cex.lab"), las=0)
mtext("40", side=1, line=3, outer=T, adj=0.16, cex=par("cex.lab"), las=0)
mtext("50", side=1, line=3, outer=T, adj=0.5, cex=par("cex.lab"), las=0)
mtext("60", side=1, line=3, outer=T, adj=0.84, cex=par("cex.lab"), las=0)
dev.off(); system(paste("evince ", fname))



## ### np.learnmem.hvlt_sum ###
##3 age groups, before AIC:
fit <- lm(np.learnmem.hvlt_sum.yC-np.learnmem.hvlt_sum.xC~hiv.x*age.x , data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)
## (Intercept)  3.26818    4.66880   0.700    0.486
## hiv.x        2.82656    6.38744   0.443    0.659
## age.x       -0.03369    0.08656  -0.389    0.698
## hiv.x:age.x -0.05570    0.12313  -0.452    0.652



################################################################
## # Mixed effects #
################################################################
library(lme4)
## ## Overall repeated meausres  ##
## ### np.learnmem.bvmt_delayC ###
summary(lmer(np.learnmem.bvmt_delayC~(1|id) + hiv*age + hiv*visit, data=dlall))
##               Estimate Std. Error t value
## (Intercept)  7.5011783  2.7105690   2.767
## hiv          9.9387720  3.6855934   2.697
## age          0.0002164  0.0501302   0.004
## visit12      1.2122161  0.3917784   3.094
## hiv:age     -0.1248662  0.0708543  -1.762
## hiv:visit12 -4.3725678  0.4907209  -8.910


## ### np.learnmem.bvmt_sumC ###
summary(lmer(np.learnmem.bvmt_sumC~(1|id) + hiv*age+hiv*visit, data=dlall))
##             Estimate Std. Error t value
## (Intercept) 20.83270    6.72326   3.099
## hiv         18.53673    9.14278   2.027
## age         -0.04122    0.12429  -0.332
## visit12      1.89995    1.04108   1.825
## hiv:age     -0.26574    0.17571  -1.512
## hiv:visit12 -2.73220    1.30329  -2.096


## ### np.learnmem.hvlt_delayC ###
summary(lmer(np.learnmem.hvlt_delayC~(1|id)+hiv*visit+hiv*age, data=dlall))
##             Estimate Std. Error t value
## (Intercept)  6.75846    2.51937   2.683
## hiv          5.32337    3.42769   1.553
## visit12      0.82634    0.41162   2.008
## age          0.01337    0.04656   0.287
## hiv:visit12  2.06612    0.51669   3.999
## hiv:age     -0.10334    0.06589  -1.568


## ### np.learnmem.hvlt_sumC ###
summary(lmer(np.learnmem.hvlt_sumC~(1|id) + hiv*age + hiv*visit, data=dlall))
##             Estimate Std. Error t value
## (Intercept) 24.07022    4.67522   5.148
## hiv          4.34707    6.35801   0.684
## age         -0.03778    0.08641  -0.437
## visit12      1.51408    0.75088   2.016
## hiv:age     -0.07257    0.12217  -0.594
## hiv:visit12  1.13213    0.93976   1.205


## ## change ##
## ### np.learnmem.bvmt_delayC ###
summary(lmer(log(np.learnmem.bvmt_delayC)~(1|id) +  hiv*age*visit, data=dlall))
##                  Estimate Std. Error t value
## (Intercept)      1.555065   0.536254   2.900
## hiv              1.380255   0.733286   1.882
## age              0.006339   0.009942   0.638
## visit12          0.532607   0.514643   1.035
## hiv:age         -0.017329   0.014134  -1.226
## hiv:visit12      0.202978   0.701556   0.289
## age:visit12     -0.005591   0.009447  -0.592
## hiv:age:visit12 -0.017232   0.013372  -1.289


## ### np.learnmem.bvmt_sumC ###
summary(lmer(np.learnmem.bvmt_sumC~(1|id) + hiv*age*visit, data=dlall))
##                 Estimate Std. Error t value
## (Intercept)     17.76669    7.44931   2.385
## hiv             19.78168   10.18467   1.942
## age              0.01634    0.13811   0.118
## visit12          8.22017    6.69768   1.227
## hiv:age         -0.28643    0.19631  -1.459
## hiv:visit12     -5.47834    9.13011  -0.600
## age:visit12     -0.11745    0.12294  -0.955
## hiv:age:visit12  0.04593    0.17402   0.264

## ### np.learnmem.hvlt_delayC ###
##Very similar results to the regression estimates, 
summary(lmer(np.learnmem.hvlt_delayC~(1|id) + hiv*age*visit, data=dlall))
##                  Estimate Std. Error t value
## (Intercept)      7.373245   2.815450   2.619
## hiv              2.476428   3.856633   0.642
## age              0.001828   0.052199   0.035
## visit12         -0.445647   2.611894  -0.171
## hiv:age         -0.046456   0.074403  -0.624
## hiv:visit12      7.676671   3.567421   2.152
## age:visit12      0.023637   0.047944   0.493
## hiv:age:visit12 -0.110618   0.068068  -1.625

## ### np.learnmem.hvlt_sumC ###
summary(lmer(np.learnmem.hvlt_sumC~(1|id) + hiv*age*visit, data=dlall))
##                 Estimate Std. Error t value
## (Intercept)     23.23462    5.22389   4.448
## hiv              3.86930    7.14271   0.542
## age             -0.02209    0.09685  -0.228
## visit12          3.24632    4.85614   0.668
## hiv:age         -0.06166    0.13768  -0.448
## hiv:visit12      1.98726    6.61981   0.300
## age:visit12     -0.03219    0.08914  -0.361
## hiv:age:visit12 -0.01958    0.12618  -0.155

## ## age group ##
summary(lmer(np.learnmem.hvlt_delayC~(1|id) + hiv*age.grp*visit, data=dlall))
##                       Estimate Std. Error t value
## (Intercept)            7.48265    0.78107   9.580
## hiv                    0.06903    0.90578   0.076
## age.grp50             -0.16629    1.03892  -0.160
## age.grp60              0.25730    1.21721   0.211
## visit12                0.54714    0.69808   0.784
## hiv:age.grp50          0.42177    1.25518   0.336
## hiv:age.grp60         -0.36992    1.93016  -0.192
## hiv:visit12            2.82331    0.80955   3.488
## age.grp50:visit12      0.38663    0.92854   0.416
## age.grp60:visit12      0.53968    1.08789   0.496
## hiv:age.grp50:visit12 -1.48161    1.12258  -1.320
## hiv:age.grp60:visit12 -3.26789    1.72509  -1.894

mytheme()
fit <- lmer(np.learnmem.hvlt_delayC~(1|id) + hiv*age.grp*visit, data=dlall)
visitx <- as.numeric(as.character(dlall$visit))
age.grpx <- as.numeric(dlall$age.grp)
uid.n <- dall$id[dall$hiv.x == 0 ]
uid.p <-  dall$id[dall$hiv.x == 1 ]
bb <- as.numeric(fit@fixef)
getModelVal <- function(bb, hiv, age.grp, visit) {
  if (age.grp=="40") {
    ad <- c(0, 0)
  } else if (age.grp=="50") {
    ad <- c(1, 0)
  } else {
    ad <- c(0, 1)
  }
  dc <- c(1, hiv, ad, visit, hiv*ad, hiv*visit, ad*visit, hiv*visit*ad)
  sum(bb*dc)
}


fname <- "MemChangeHvltDelayTrendGroupOverlay.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
par(mfrow = (c(1,2)))
## non hiv plot
plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delayC)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (jj in uid.n) {
  sel <- dlall$id == jj
  age.x <- min(age.grpx[sel])
  hiv <- 0
  ## month 0
  yy <- dlall$np.learnmem.hvlt_delayC[sel]
  yv <- visitx[sel]
  ## ccols <- ifelse(rep(hiv, 3), col.p, col.n)
  print(age.x)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.n[age.x], lwd=0.5)
  points(yv, yy, col=col.n[age.x])
}

for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ## non-hiv
  yy <-  c( getModelVal(bb, 0, gc, 0), getModelVal(bb, 0, gc, 1))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.n[gg], lwd=10)
}
legend("topleft", c("Non-HIV 40", "Non-HIV 50", "Non-HIV 60"), col = col.n, lty=1, pch=1, merge=F, lwd=2.5)

## hiv plot
plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delayC)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (jj in uid.p) {
  sel <- dlall$id == jj
  age.x <- min(age.grpx[sel])
  hiv <- 0
  ## month 0
  yy <- dlall$np.learnmem.hvlt_delayC[sel]
  yv <- visitx[sel]
  ## ccols <- ifelse(rep(hiv, 3), col.p, col.n)
  print(age.x)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.p[age.x], lwd=0.5)
  points(yv, yy, col=col.p[age.x], pch=2)
}

for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ## hiv
  yy <-  c( getModelVal(bb, 1, gc, 0), getModelVal(bb, 1, gc, 1))
  yv <-  c(0, 12)
  segments(yv[1], yy[1],  yv[2], yy[2], col=col.p[gg], lwd=10)
}
legend("topleft", c("HIV 40", "HIV 50", "HIV 60"), col = col.p, lty=1, pch=2, merge=F, lwd=2.5)
dev.off(); system(paste("evince ", fname))

################################################################
## # Simple plot # 
################################################################
fname <- "MemChangeHvltDelayTrendErrorBar.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
par(mfrow = (c(1,2)))
error.bars <- function (x, upper, lower, width = 0.02, ...) 
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}

std.err <- function(x) {
  x <- na.omit(x)
  sd(x)/sqrt(length(x))
}
##non hiv plot
plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delayC)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##non-hiv
  yh <- dlall$np.learnmem.hvlt_delayC
  age <- dlall$age.grp
  visit <- dlall$visit
  hiv <- dlall$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==0], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==0 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==0]),  std.err(yh[visit == 12 & age == gc & hiv==0 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.n[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.n[gg])
}
legend("topleft", c("Non-HIV 40", "Non-HIV 50", "Non-HIV 60"), col = col.n, lty=1, merge=F, lwd=2.5)

plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delayC)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##hiv
  yh <- dlall$np.learnmem.hvlt_delayC
  age <- dlall$age.grp
  visit <- dlall$visit
  hiv <- dlall$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==1], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==1 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==1]),  std.err(yh[visit == 12 & age == gc & hiv==1 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.p[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.p[gg])
}
legend("topleft", c("HIV 40", "HIV 50", "HIV 60"), col = col.p, lty=1, merge=F, lwd=2.5)

dev.off(); system(paste("evince ", fname))


##Uncorrected, 3 age groups 
fname <- "MemChangeHvltDelayTrendErrorBarUncorrected.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
par(mfrow = (c(1,2)))
error.bars <- function (x, upper, lower, width = 0.02, ...) 
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}

std.err <- function(x) {
  x <- na.omit(x)
  sd(x)/sqrt(length(x))
}
##non hiv plot
plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delay)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##non-hiv
  yh <- dlall$np.learnmem.hvlt_delay
  age <- dlall$age.grp
  visit <- dlall$visit
  hiv <- dlall$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==0], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==0 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==0]),  std.err(yh[visit == 12 & age == gc & hiv==0 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.n[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.n[gg])
}
legend("topleft", c("Non-HIV 40", "Non-HIV 50", "Non-HIV 60" ), col = col.n, lty=1, merge=F, lwd=2.5)

plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delay)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##hiv
  yh <- dlall$np.learnmem.hvlt_delay
  age <- dlall$age.grp
  visit <- dlall$visit
  hiv <- dlall$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==1], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==1 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==1]),  std.err(yh[visit == 12 & age == gc & hiv==1 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.p[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.p[gg])
}
legend("topleft", c("HIV 40", "HIV 50", "HIV 60"), col = col.p, lty=1, merge=F, lwd=2.5)

dev.off(); system(paste("evince ", fname))



##Uncorrected, 2 age groups 
fname <- "MemChangeHvltDelayTrendErrorBarUncorrected2Groups.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
par(mfrow = (c(1,2)))
error.bars <- function (x, upper, lower, width = 0.02, ...) 
{
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
}

std.err <- function(x) {
  x <- na.omit(x)
  sd(x)/sqrt(length(x))
}
##non hiv plot
plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delay)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:2)  {
  gc <-   c("40", "55")[gg]
  ##non-hiv
  yh <- dlall$np.learnmem.hvlt_delay
  age <- dlall$age.grp55
  visit <- dlall$visit
  hiv <- dlall$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==0], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==0 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==0]),  std.err(yh[visit == 12 & age == gc & hiv==0 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.n[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.n[gg])
}
legend("topleft", c("Non-HIV 40", "Non-HIV 55"), col = col.n, lty=1, merge=F, lwd=2.5)

plot(c(0, 12), range(na.omit(dlall$np.learnmem.hvlt_delay)), type="n", ylab="HVLT Delay", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "55")[gg]
  ##hiv
  yh <- dlall$np.learnmem.hvlt_delay
  age <- dlall$age.grp
  visit <- dlall$visit
  hiv <- dlall$hiv 
  yy <-  c(mean(yh[visit == 0 & age == gc & hiv ==1], na.rm=T), mean(yh[visit == 12 & age == gc & hiv==1 ], na.rm=T))
  yse <-  c(std.err(yh[visit == 0 & age == gc& hiv==1]),  std.err(yh[visit == 12 & age == gc & hiv==1 ]))
  yv <-  c(0, 12)
  segments(yv[1], yy[1], yv[2], yy[2], col=col.p[gg], lwd=5)
  error.bars(yv, yy+yse, yy-yse, lwd=3, col=col.p[gg])
}
legend("topleft", c("HIV 40", "HIV 50", "HIV 60"), col = col.p, lty=1, merge=F, lwd=2.5)

dev.off(); system(paste("evince ", fname))











##regression lines
fname <-  "MemChangeHvltDelayRegErrorBar.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
fit <- lm(np.learnmem.hvlt_delay.yC-np.learnmem.hvlt_delay.xC~hiv.x*age.x, data=dall)
summary(fit)

newdata <- data.frame(hiv.x=rep(c(0, 1), each=100), age.x=c(seq(40, 65, length.out=100), seq(40, 65, length.out=100)))
dfit <- cbind(data.frame(predict(fit, newdata=newdata, interval="confidence")), newdata)

plot(range(dfit$age.x), range(c(dfit$lwr, dfit$upr)), type="n", xlab="Age", ylab="HVLT Delay Change")
##hiv
sel <- dfit$hiv.x==1
lines(dfit$age.x[sel], dfit$fit[sel], col=mycols[1], lwd=5)
lines(dfit$age.x[sel], dfit$upr[sel], col=mycols[1], lty=1)
lines(dfit$age.x[sel], dfit$lwr[sel], col=mycols[1], lty=1)

##non-hiv
sel <- dfit$hiv.x==0
lines(dfit$age.x[sel], dfit$fit[sel], col=mycols[2], lwd=5)
lines(dfit$age.x[sel], dfit$upr[sel], col=mycols[2], lty=1)
lines(dfit$age.x[sel], dfit$lwr[sel], col=mycols[2], lty=1)
legend("topright", c("Non-HIV", "HIV"), col=mycols[c(2,1)], lty=1, lwd=3)
dev.off(); system(paste("evince ", fname))







################################################################
## # Age trend plot # 
################################################################
## ### learnmem.hvlt_delay ###
fit <- lm(np.learnmem.hvlt_delay.yC-np.learnmem.hvlt_delay.xC~hiv.x*age.x, data=dall)
summary(fit)

y <- dddd$np.learnmem.hvlt_delay.yC-dddd$np.learnmem.hvlt_delay.xC

dddd <- dall[-26, ]
sel <- dddd$hiv.x==1
yhat <- fitted(fit)
age <- dddd$age.x



fname <- "MemChangeHvltDelayPoints.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
plot(age[sel], y[sel], col=mycols[1], ylim=range(y), xlim=range(age), xlab="Age",ylab="Yearly Decline")
points(age[!sel], y[!sel], col=mycols[2], pch=2)

dnew <- data.frame(age.x=c(40, 75), hiv.x=1)
y1 <- predict(fit, dnew)
lines(c(40, 75), y1, col=mycols[1])
dnew <- data.frame(age.x=c(40, 75), hiv.x=0)
y1 <- predict(fit, dnew)
lines(c(40, 75), y1, col=mycols[2])
dev.off(); system(paste("evince ", fname))




fname <- "MemChangeHvltDelayHypothetical.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
##HIV group
b <- (-3.56301 + 7.14195)
a <- (0.06689-0.14001)
A0 <- 40
c <- 8 - 0.5*a*A0^2 - b*A0
xx <- 40:75
yy <-  0.5*a*xx^2 + b*xx+c
plot(xx, yy, type="l", cols=mycols[1])
b <- (-3.56301)
a <- (0.06689)
A0 <- 40
c <- 8 - 0.5*a*A0^2 - b*A0
xx <- 40:75
yy <-  0.5*a*xx^2 + b*xx+c
lines(xx, yy, cols=mycols[2])


plot(age[sel], y[sel], col=mycols[1], ylim=range(y), xlim=range(age), xlab="Age",ylab="Yearly Decline")
points(age[!sel], y[!sel], col=mycols[2], pch=2)

dnew <- data.frame(age.x=c(40, 75), hiv.x=1)
y1 <- predict(fit, dnew)
lines(c(40, 75), y1, col=mycols[1])
dnew <- data.frame(age.x=c(40, 75), hiv.x=0)
y1 <- predict(fit, dnew)
lines(c(40, 75), y1, col=mycols[2])
dev.off(); system(paste("evince ", fname))


################################################################
## # Focus on 1st visit #
################################################################
d1 <- d[d$visit==0, ]                  

##Separate into 3 age groups: 40-49, 50-59, and 60+, with group sizes
##40-49: 39; 50-59: 35; 60+: 10
##Use groups instead of numeric values of age
d1$age.grp <- NA
d1$age.grp[dd$age>=40] <- "40"
d1$age.grp[dd$age>=50] <- "50"
d1$age.grp[dd$age>=60] <- "60"
d1$age.grp <- factor(d1$age.grp)


##Separate into 2 age groups: 40-55, and 55+, with group sizes
##40-49: 39; 50-59: 35; 60+: 10
##Use groups instead of numeric values of age
d1$age.grp55 <- NA
d1$age.grp55[d1$age>=40] <- "40"
d1$age.grp55[d1$age>=55] <- "55"
d1$age.grp55 <- factor(d1$age.grp55)       


##Univaraite approaches to yield easier interpretation

##Add clinical information: duration, CD4, CD4 Nadir, HCV, HIV/RAV4, CART
dc <- read.csv("AgeEffectsMASTER031212.csv", sep=";")      
sel <- c(1,2,13:22)
dc <- data.frame(dc[,sel])                                 
# dc[1,]
#  ID Visit HIVduration HIVRNA_Detectable HIVRNA_value CD4current CD4nadir HAART ARV cART HCVlifetime HCVcurrent
#  1     0          22                 0            0        394      201     1   1    1           1          1

dc1 <- dc[dc$Visit==0,]         
names(dc1)[1] <- "id"
################ merge all together, non-hiv=2, hiv>100 = 1, hiv<100 = 0
dalll <- merge(d1, dc1, by=("id"))    
dalll$CD4current100 <-  1*(dalll$CD4current >= 100) + 2*(dalll$hiv == 0)
dalll$CD4nadir100 <-  1*(dalll$CD4nadir >= 100) + 2*(dalll$hiv == 0)
  

##Dictomize age by 55 age cut-off
dalll$age.grp2 <- 0 
dalll$age.grp2[dalll$age>=55] <- 1

##Tranform age  by logarithm
dalll$age.log <- log(dalll$age)

##Tranform age  by square
dalll$age.sq <- (dalll$age)^2


## source CorrectFunctions
source("/Users/Chenyang/Desktop/CorrectingFunctions.R")
##Long format: correct for age, sex, education for non-HIV, 
##Correction for all clinical variables and demographical variables, at baseline
aa <- data.frame(as.matrix(dalll[,c(grp.mem)]))
for (jj in 1:ncol(aa)) {
  
  sel <- dalll$hiv==0  
  ## correcting for education only, and apply the same correction to the second visit.  
  tmp <- correctEffect(aa[sel, jj], data.frame(education=dalll[sel,]$education), data.frame(education=13))
  aa[sel, jj] <- tmp$ce
  
  ##correct for hiv by the same formula
  sel <- dalll$hiv==1 
  rr <- applyCorrect(tmp$fit, aa[sel, jj], data.frame(education=dalll[sel,]$education), data.frame(education=13))
  aa[sel, jj] <- rr
  
}

names(aa) <- paste(names(aa), "C", sep="")
dalll <- cbind(dalll, aa)

dalll$age2 <- (dalll$age)^2 

uid <- dalll$id
dalll$age.grp <-  NA
for (jjj in uid) {
  sel <- dalll$id == jjj
  age <- min(dalll$age[sel])
  if (age <50) {
    dalll$age.grp[sel] <- "40"
  } else if (age>=50 && age<60 ) {
    dalll$age.grp[sel] <- "50"
  } else {
    dalll$age.grp[sel] <- "60"
  }
}
dalll$age.grp <- factor(dalll$age.grp)    # dim(dalll)=131*56



##Select color for plots
library(RColorBrewer)
source("/Users/Chenyang/Desktop/MyPlot.R")
mycols <-  brewer.pal(11, "RdBu")
col.n <- mycols[9:11]
col.p <- rev(mycols[1:3])
col.all <- as.vector(rbind(col.n, col.p))


# dalll[1,52:55]
#  np.learnmem.bvmt_delayC np.learnmem.bvmt_sumC np.learnmem.hvlt_delayC np.learnmem.hvlt_sumC
#                5.690232              16.05586                10.63333              32.20038


################################################################
##  #Baseline Analysis
################################################################
## ## MANOVA ##
##Baseline:
##Without Age.sq
fit <- manova(as.matrix(dalll[,52:55])~hiv*age, data=dalll)
summary(fit, test="Wilks")
##           Df   Wilks approx F num Df den Df Pr(>F)
##hiv         1 0.97381  0.82693      4    123 0.5105
##age         1 0.96581  1.08843      4    123 0.3653
##hiv:age     1 0.98371  0.50922      4    123 0.7290



## ## HIV Positive vs Negative, Univariate ##

## ### np.learnmem.bvmt_delay ###
##3 age groups, before AIC:
fit <- lm(np.learnmem.bvmt_delayC~hiv*age, data=dalll)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)    
##(Intercept)  9.82725    2.61302   3.761 0.000257 ***
##hiv          0.58626    3.73525   0.157 0.875530    
##age         -0.03948    0.04948  -0.798 0.426426    
##hiv:age     -0.01166    0.07327  -0.159 0.873797     

postscript(file="MemBaselineBvmtDelay.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.learnmem.bvmt_delayC~age.grp*hiv, data=dalll, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="BD")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()

fit <- lm(np.learnmem.bvmt_delay~hiv*age, data=dalll)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)   
##(Intercept)  7.708254   2.783610   2.769  0.00646 **
##hiv          1.860910   3.979101   0.468  0.64082   
##age          0.006447   0.052715   0.122  0.90286   
##hiv:age     -0.043263   0.078050  -0.554  0.58034  

## ### np.learnmem.bvmt_sum ###
##3 age groups, before AIC:
fit <- lm(np.learnmem.bvmt_sumC~hiv*age, data=dalll)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)    
##(Intercept) 23.19419    5.68600   4.079 7.93e-05 ***
##hiv          6.47925    8.12799   0.797    0.427    
##age         -0.07011    0.10768  -0.651    0.516    
##hiv:age     -0.13259    0.15943  -0.832    0.407   


fit <- lm(np.learnmem.bvmt_sum~hiv*age, data=dalll)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)   
##(Intercept) 16.73573    6.22087   2.690   0.0081 **
##hiv         10.36424    8.89258   1.165   0.2460   
##age          0.06988    0.11781   0.593   0.5541   
##hiv:age     -0.22891    0.17443  -1.312   0.1918  

## ### np.learnmem.hvlt_delay ###
##3 age groups, before AIC:
fit <- lm(np.learnmem.hvlt_delayC~hiv*age, data=dalll)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)    
##(Intercept)  7.393721   2.121281   3.485 0.000677 ***
##hiv         -0.685933   3.038812  -0.226 0.821782    
##age          0.008036   0.040172   0.200 0.841775    
##hiv:age      0.011257   0.059661   0.189 0.850646  

fit <- lm(np.learnmem.hvlt_delay~hiv*age, data=dalll)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  4.88546    2.28306   2.140   0.0343 *
##hiv          0.94598    3.27057   0.289   0.7729  
##age          0.06240    0.04324   1.443   0.1514  
##hiv:age     -0.02904    0.06421  -0.452   0.6519  


## ### np.learnmem.hvlt_sum ###
##3 age groups, before AIC:
fit <- lm(np.learnmem.hvlt_sumC~hiv*age, data=dalll)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)    
##(Intercept) 21.332636   3.821835   5.582 1.38e-07 ***
##hiv          0.722418   5.463217   0.132    0.895    
##age          0.021202   0.072377   0.293    0.770    
##hiv:age      0.004844   0.107161   0.045    0.964 

fit <- lm(np.learnmem.hvlt_sumC~hiv*age.grp, data=dalll)
summary(fit)
##              Estimate Std. Error t value Pr(>|t|)    
##(Intercept)   20.95964    0.98810  21.212   <2e-16 ***
##hiv            1.90674    1.17608   1.621   0.1075    
##age.grp50      3.53520    1.45774   2.425   0.0167 *  
##age.grp60      0.86373    1.77370   0.487   0.6271    
##hiv:age.grp50 -2.44408    1.76735  -1.383   0.1692    
##hiv:age.grp60  0.02037    2.90423   0.007   0.9944    

fit <- lm(np.learnmem.hvlt_sum~hiv*age, data=dalll)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)    
##(Intercept) 15.86277    4.17948   3.795 0.000227 ***
##hiv          4.01273    5.97447   0.672 0.503029    
##age          0.13976    0.07915   1.766 0.079834 .  
##hiv:age     -0.07673    0.11719  -0.655 0.513793   

fit <- lm(np.learnmem.hvlt_sum~hiv*age.grp, data=dalll)
summary(fit)
##              Estimate Std. Error t value Pr(>|t|)    
##(Intercept)     20.200      1.057  19.114  < 2e-16 ***
##hiv              2.050      1.258   1.630 0.105664    
##age.grp50        5.800      1.559   3.720 0.000299 ***
##age.grp60        4.022      1.897   2.120 0.035956 *  
##hiv:age.grp50   -4.262      1.890  -2.255 0.025883 *  
##hiv:age.grp60   -1.522      3.106  -0.490 0.624942   




