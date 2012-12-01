##% HIV-Aging: Executive Domain
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
dd <- merge(d1, d2, by=("id"))
dl <- rbind(d1, d2)

## first remove the variable np.exec.stroop-interf
## since this variable contains a lot of missing values
## grp.exe <- c(10,11,13,14)

##calculate changes in the ececutive domain
aa <- data.frame(as.matrix(dd[,grp.exe+32]-dd[,grp.exe]))
names(aa) <- paste(names(aa), "change", sep="")
dd <- cbind(dd, aa)

##calculate relative changes in the ececutive domain
aa <- data.frame(as.matrix( (dd[,grp.exe+32]-dd[,grp.exe] )/(dd[,grp.exe]+1)   ))
names(aa) <- paste(names(aa), "rchange", sep="")
dd <- cbind(dd, aa)

##calculate sqrt changes in the ececutive domain
##aa <- data.frame(as.matrix( ( sqrt(dd[,grp.exe+32])-sqrt(dd[,grp.exe]) )  ))   ## negative values in dd[,grp.exe], NaNs produced
##names(aa) <- paste(names(aa), "sqchange", sep="")
##dd <- cbind(dd, aa)

##calculate sqrt changes in the ececutive domain
##aa <- data.frame(as.matrix( ( log(1+dd[,grp.exe+32])-log(1+dd[,grp.exe]) )  )) ## NaNs produced
##names(aa) <- paste(names(aa), "logchange", sep="")
##dd <- cbind(dd, aa)


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
## New CorrectFunctions
source("/Users/Chenyang/Desktop/CorrectingFunctions.R")
##Long format: correct for age, sex, education for non-HIV, 
##Correction for all clinical variables and demographical variables, at baseline and 12month
aa <- data.frame(as.matrix(dlall[,c(grp.exe)]))
#aa <- data.frame(as.matrix(dlall[,c(10,11,13,14)]))
#data_in <- data.frame(education=dlall[sel,]$education, HIVduration=dlall[sel,]$HIVduration, HIVRNA_Detectable=dlall[sel,]$HIVRNA_Detectable,
#           cART=dlall[sel,]$cART, HCVcurrent=dlall[sel,]$HCVcurrent, CD4nadir100=dlall[sel,]$CD4nadir100)
#data_n <- data.frame(education=13, HIVduration=17, HIVRNA_Detectable=0, cART=1, HCVcurrent=0, CD4nadir100=0)
for (jj in 1:ncol(aa)) {
  ## visit 0
  vsel <-  dlall$visit == 0 
  sel <- dlall$hiv==0  & vsel
  ##version1: correcting for education only, and apply the same correction to the second visit.
  ##version2: correcting for HIVduration=17, HIVRNA_Detectable=0, cART=1, HCVcurrent=0, CD4nadir100=0  
  tmp <- correctEffect(aa[sel, jj], data.frame(education=dlall[sel,]$education), data.frame(education=13))
  aa[sel, jj] <- tmp$ce
  #tmp <- correctEffect(aa[sel, jj], data.frame(education=dlall[sel,]$education), data.frame(education=13))
  #aa[sel, jj] <- tmp$ce
 
  ##correct for hiv by the same formula
  sel <- dlall$hiv==1 & vsel
  
#  rr <- applyCorrect(tmp$fit, aa[sel, jj], data.frame(education=dlall[sel,]$education), data.frame(education=13))
#  aa[sel, jj] <- rr
#  tmph <- correctEffect(aa[sel, jj], 
#                     data.frame(HIVduration=dlall[sel,]$HIVduration, 
#                                HIVRNA_Detectable=dlall[sel,]$HIVRNA_Detectable, cART=dlall[sel,]$cART, 
#                                HCVcurrent=dlall[sel,]$HCVcurrent, CD4nadir100=dlall[sel,]$CD4nadir100), 
#                     data.frame(HIVduration=17, HIVRNA_Detectable=0, cART=1, HCVcurrent=0, CD4nadir100=0))

  tmph <- correctEffect(aa[sel, jj], 
                     data.frame(education=dlall[sel,]$education, HIVduration=dlall[sel,]$HIVduration, 
                                HIVRNA_Detectable=dlall[sel,]$HIVRNA_Detectable, cART=dlall[sel,]$cART, 
                                HCVcurrent=dlall[sel,]$HCVcurrent, CD4nadir100=dlall[sel,]$CD4nadir100), 
                     data.frame(education=13, HIVduration=17, HIVRNA_Detectable=0, cART=1, HCVcurrent=0, CD4nadir100=0))
  aa[sel, jj] <- tmph$ce

   
  ##repeat above for second visit
  ## visit 12
  vsel <-  dlall$visit == 12
  sel <- dlall$hiv==0  & vsel
  ##version1: correcting for education only
  ##version2: correcting for HIVduration=17, HIVRNA_Detectable=0, cART=1, HCVcurrent=0, CD4nadir100=0  
  rr <- applyCorrect(tmp$fit, aa[sel, jj], data.frame(education=dlall[sel,]$education), data.frame(education=13))
  aa[sel, jj] <- rr
  
  ##correct for hiv by the same formula
  sel <- dlall$hiv==1 & vsel
  
#  rr <- applyCorrect(tmp$fit, aa[sel, jj], data.frame(education=dlall[sel,]$education), data.frame(education=13))
#  aa[sel, jj] <- rr
#  rrh <- applyCorrect(tmph$fit, aa[sel, jj], 
#                     data.frame(HIVduration=dlall[sel,]$HIVduration, 
#                                HIVRNA_Detectable=dlall[sel,]$HIVRNA_Detectable, cART=dlall[sel,]$cART, 
#                                HCVcurrent=dlall[sel,]$HCVcurrent, CD4nadir100=dlall[sel,]$CD4nadir100), 
#                     data.frame(HIVduration=17, HIVRNA_Detectable=0, cART=1, HCVcurrent=0, CD4nadir100=0))
  
  rrh <- applyCorrect(tmph$fit, aa[sel, jj], 
                     data.frame(education=dlall[sel,]$education, HIVduration=dlall[sel,]$HIVduration, 
                                HIVRNA_Detectable=dlall[sel,]$HIVRNA_Detectable, cART=dlall[sel,]$cART, 
                                HCVcurrent=dlall[sel,]$HCVcurrent, CD4nadir100=dlall[sel,]$CD4nadir100), 
                     data.frame(education=13, HIVduration=17, HIVRNA_Detectable=0, cART=1, HCVcurrent=0, CD4nadir100=0))
  aa[sel, jj] <- rrh 
  #rr <- applyCorrect(tmp$fit, aa[sel, jj], data.frame(education=dlall[sel,]$education), data.frame(education=13))
  #aa[sel, jj] <- rr
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
aa <- data.frame(as.matrix(dall[,c(grp.exe, grp.exe+32)]))
ashift <- (ncol(aa)/2)
for (jj in 1:(ncol(aa)/2) ) {
  sel <- dall$hiv.x==0
  ##correcting for education only 
  tmp <- correctEffect(aa[sel, jj], data.frame(education=dall[sel,]$education.x), data.frame(education=13))
  aa[sel, jj] <- tmp$ce
  
  ##correct this for the second visit
  aa[sel, jj+ashift] <- applyCorrect(tmp$fit, aa[sel, jj+ashift], data.frame(education=dall[sel,]$education.y), data.frame(education=13))

  ##correct this for hiv positives as well
  #sel <- dall$hiv.x==0    This should be sel <- dall$hiv.x==1
  sel <- dall$hiv.x==1
  
#  rr<- applyCorrect(tmp$fit, aa[sel, jj], data.frame(education=dall[sel,]$education.x), data.frame(education=13))
#  aa[sel, jj] <- rr
#  tmph <- correctEffect(aa[sel, jj], 
#                        data.frame(HIVduration=dall[sel,]$HIVduration.x, 
#                                   HIVRNA_Detectable=dall[sel,]$HIVRNA_Detectable.x,cART=dall[sel,]$cART.x, 
#                                   HCVcurrent=dall[sel,]$HCVcurrent.x, CD4nadir100=dall[sel,]$CD4nadir100.x), 
#                        data.frame(HIVduration=17, HIVRNA_Detectable=0, cART=1, HCVcurrent=0, CD4nadir100=0))
  tmph <- correctEffect(aa[sel, jj], 
                        data.frame(education=dall[sel,]$education.x, HIVduration=dall[sel,]$HIVduration.x, 
                                   HIVRNA_Detectable=dall[sel,]$HIVRNA_Detectable.x,cART=dall[sel,]$cART.x, 
                                   HCVcurrent=dall[sel,]$HCVcurrent.x, CD4nadir100=dall[sel,]$CD4nadir100.x), 
                        data.frame(education=13, HIVduration=17, HIVRNA_Detectable=0, cART=1, HCVcurrent=0, CD4nadir100=0))
  aa[sel, jj] <- tmph$ce

#  aa[sel, jj+ashift] <- applyCorrect(tmp$fit, aa[sel, jj+ashift], data.frame(education=dall[sel,]$education.y), data.frame(education=13))
#  aa[sel, jj+ashift] <- applyCorrect(tmph$fit, aa[sel, jj+ashift], 
#                        data.frame(HIVduration=dall[sel,]$HIVduration.y, 
#                                   HIVRNA_Detectable=dall[sel,]$HIVRNA_Detectable.y, cART=dall[sel,]$cART.y, 
#                                   HCVcurrent=dall[sel,]$HCVcurrent.y, CD4nadir100=dall[sel,]$CD4nadir100.y), 
#                        data.frame(HIVduration=17, HIVRNA_Detectable=0, cART=1, HCVcurrent=0, CD4nadir100=0))

  aa[sel, jj+ashift] <- applyCorrect(tmph$fit, aa[sel, jj+ashift], 
                        data.frame(education=dall[sel,]$education.y, HIVduration=dall[sel,]$HIVduration.y, 
                                   HIVRNA_Detectable=dall[sel,]$HIVRNA_Detectable.y, cART=dall[sel,]$cART.y, 
                                   HCVcurrent=dall[sel,]$HCVcurrent.y, CD4nadir100=dall[sel,]$CD4nadir100.y), 
                        data.frame(education=13, HIVduration=17, HIVRNA_Detectable=0, cART=1, HCVcurrent=0, CD4nadir100=0))
}

## Only correct for Education
#for (jj in 1:(ncol(aa)/2) ) {
#  sel <- dall$hiv.x==0
#  ##correcting for education only 
# 
#  tmp <- correctEffect(aa[sel, jj], data.frame(HIVduration=dlall[sel,]$HIVduration.x), data.frame(HIVduration=17))
#  aa[sel, jj] <- tmp$ce
#  
#  ##correct this for the second visit
#  aa[sel, jj+ashift] <- applyCorrect(tmp$fit, aa[sel, jj+ashift], data.frame(HIVduration=dall[sel,]$HIVduration.y), data.frame(HIVduration=17))
#
#  ##correct this for hiv positives as well
#  sel <- dall$hiv.x==0
#  aa[sel, jj] <- applyCorrect(tmp$fit, aa[sel, jj], data.frame(HIVduration=dall[sel,]$HIVduration.x), data.frame(HIVduration=17))
#
#  aa[sel, jj+ashift] <- applyCorrect(tmp$fit, aa[sel, jj+ashift], data.frame(HIVduration=dall[sel,]$HIVduration.y), data.frame(HIVduration=17))
#}



names(aa) <- paste(names(aa), "C", sep="")
dall <- cbind(dall, aa)

##Add age^2 terms
dall <- data.frame(dall, age2 = (dall$age.x)^2)     #dim(dall)=84*116

## dim(dall)   # after rmoving grp.exe[13], the dim is 84*112
## dall[1,104:111]
##  np.exec.actionflu_total.xC np.exec.cowat_fas.xC np.exec.trail_b_time.xC np.exec.wais_lnseq.xC np.exec.actionflu_total.yC
##                         20                   38                      90                    11                         15
##  np.exec.cowat_fas.yC np.exec.trail_b_time.yC np.exec.wais_lnseq.yC
##                   46                      99                     8

## dall[1,106:115]
## np.exec.actionflu_total.xC np.exec.cowat_fas.xC       np.exec.stroop_interf.xC  np.exec.trail_b_time.xC
##                         20                   38                             NA                       90
## np.exec.wais_lnseq.xC      np.exec.actionflu_total.yC np.exec.cowat_fas.yC      np.exec.stroop_interf.yC
##                    11                         15                        46                     -7.833333
## np.exec.trail_b_time.yC    np.exec.wais_lnseq.yC
##                      99                        8

################################################################
##  #Baseline Analysis
################################################################
## ## MANOVA ##
##Baseline:
##Without Age.sq
fit <- manova(as.matrix(dall[,106:110])~hiv.x*age.x, data=dall)
summary(fit, test="Wilks")
#            Df   Wilks approx F num Df den Df Pr(>F)
#hiv.x        1 0.87512  1.54114      5     54 0.1926
#age.x        1 0.93448  0.75728      5     54 0.5845
#hiv.x:age.x  1 0.86436  1.69486      5     54 0.1516

## ## HIV Positive vs Negative, Univariate ##
##Use regression to compare each raw T score. Not considering CD4current because not significant, consistent with past analyses.  Also not selected by AIC.  Consider different ways to dictomizing/transforming age, including log, two groups, and square.  Square seems to work the best.

##HCV current/lifetime don't change the answer much.  Stick to current in the following analysis.

## ### np.exec.actionflu_total ###
##3 age groups, before AIC:
fit <- lm(np.exec.actionflu_total.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 17.38683    5.55334   3.131  0.00243 **
#hiv.x        4.15111    7.59759   0.546  0.58633   
#age.x       -0.04856    0.10296  -0.472  0.63844   
#hiv.x:age.x -0.02185    0.14645  -0.149  0.88177  

postscript(file="ExecBaselineActionfluTotal.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.exec.actionflu_total.xC~age.grp*hiv.x, data=dall, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="BD")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()

fit <- lm(np.exec.actionflu_total.x ~ hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)
##(Intercept)   9.4231     6.5278   1.444    0.153
##hiv.x         6.9831     8.9308   0.782    0.437
##age.x         0.1347     0.1210   1.113    0.269
##hiv.x:age.x  -0.1381     0.1722  -0.802    0.425


## ### np.exec.cowat_fas ###
##3 age groups, before AIC:
fit <- lm(np.exec.cowat_fas.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  46.0149    11.3878   4.041 0.000123 ***
#hiv.x        -4.6600    15.6052  -0.299 0.766017    
#age.x        -0.1331     0.2111  -0.631 0.530162    
#hiv.x:age.x   0.1752     0.3011   0.582 0.562247  

postscript(file="ExeBaselineCowatFas.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.exec.cowat_fas.xC~age.grp*hiv.x, data=dall, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="BS")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()

fit <- lm(np.exec.cowat_fas.x~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  36.2383    12.8148   2.828  0.00594 **
#hiv.x       -11.7394    17.5607  -0.669  0.50576   
#age.x         0.0919     0.2376   0.387  0.69995   
#hiv.x:age.x   0.2110     0.3388   0.623  0.53520  


## ### np.exec.stroop_interf ###
##3 age groups, before AIC:
fit <- lm(np.exec.stroop_interf.xC~hiv.x*age.x , data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)   8.0733    16.2195   0.498   0.6205  
#hiv.x       -34.5733    23.2378  -1.488   0.1422  
#age.x        -0.1693     0.2932  -0.577   0.5659  
#hiv.x:age.x   0.7998     0.4445   1.799   0.0772 . 

postscript(file="ExeBaselineStroopInterf.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.exec.stroop_interf.xC~age.grp*hiv.x, data=dall, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="BS")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()

fit <- lm(np.exec.stroop_interf.xC~hiv.x*age.grp , data=dall)
summary(fit)
##                Estimate Std. Error t value Pr(>|t|)  
##(Intercept)       0.5085     5.1097   0.100   0.9211  
##hiv.x             3.0710     5.8260   0.527   0.6002  
##age.grp50        -1.1375     6.2581  -0.182   0.8564  
##age.grp60        -4.0968     6.9634  -0.588   0.5587  
##hiv.x:age.grp50   1.5956     7.5790   0.211   0.8340  
##hiv.x:age.grp60  20.7925    11.6039   1.792   0.0786 .

fit <- lm(np.exec.stroop_interf.x~hiv.x*age.x , data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)  
##(Intercept)   8.8696    16.9574   0.523   0.6029  
##hiv.x       -38.0481    24.2950  -1.566   0.1228  
##age.x        -0.1911     0.3065  -0.623   0.5355  
##hiv.x:age.x   0.7879     0.4648   1.695   0.0954 .


## ### np.exec.trail_b_time ###
##3 age groups, before AIC:
fit <- lm(np.exec.trail_b_time.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)  
#(Intercept)  77.07221   31.32588   2.460    0.016 *
#hiv.x       -52.73987   42.85728  -1.231    0.222  
#age.x         0.01342    0.58079   0.023    0.982  
#hiv.x:age.x   0.82098    0.82614   0.994    0.323     

postscript(file="ExeBaselineTrailBTime.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.exec.trail_b_time.xC~age.grp*hiv.x, data=dall, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="BS")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()

fit <- lm(np.exec.trail_b_time.x~hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  90.5824    35.3655   2.561   0.0123 *
##hiv.x       -21.7906    48.3840  -0.450   0.6537  
##age.x        -0.2975     0.6557  -0.454   0.6512  
##hiv.x:age.x   0.6010     0.9327   0.644   0.5212  


## ### np.exec.wais_lnseq ###
##3 age groups, before AIC:
fit <- lm(np.exec.wais_lnseq.xC~hiv.x*age.x, data=dall)
summary(fit)
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  9.968041   2.676644   3.724 0.000363 ***
#hiv.x        1.764386   3.661946   0.482 0.631252    
#age.x        0.001668   0.049626   0.034 0.973264    
#hiv.x:age.x -0.031309   0.070589  -0.444 0.658578  

postscript(file="ExeBaselineWaisLnseq.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.exec.wais_lnseq.xC~age.grp*hiv.x, data=dall, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="BS")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()

fit <- lm(np.exec.wais_lnseq.x~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)   
#(Intercept)  9.48245    2.95910   3.205  0.00194 **
#hiv.x        0.86983    4.04838   0.215  0.83043   
#age.x        0.01285    0.05486   0.234  0.81548   
#hiv.x:age.x -0.03048    0.07804  -0.391  0.69718   

################################################################
## # Change 12-0
################################################################

## ## MANOVA ##
##change:
fit <- manova(as.matrix(dall[,111:115]-dall[,106:110])~hiv.x*age.x, data=dall)
summary(fit, test="Wilks")
#            Df   Wilks approx F num Df den Df  Pr(>F)  
#hiv.x        1 0.80047   2.6921      5     54 0.03038 *
#age.x        1 0.84365   2.0015      5     54 0.09311 .
#hiv.x:age.x  1 0.87771   1.5047      5     54 0.20373  



## ## HIV Positive vs Negative ##
##Use baseline variables to predict 12-0 changes

## ### np.exec.actionflu_total ###
##3 age groups, before AIC:
fit <- lm(np.exec.actionflu_total.yC-np.exec.actionflu_total.xC~hiv.x*age.x, data=dall)
summary(fit)
#             Estimate Std. Error t value Pr(>|t|)
#(Intercept) -1.760307   5.783243  -0.304    0.762
#hiv.x       -2.166403   7.912119  -0.274    0.785
#age.x        0.051195   0.107223   0.477    0.634
#hiv.x:age.x  0.003821   0.152517   0.025    0.980


## ### np.exec.cowat_fas ###
##3 age groups, before AIC:
fit <- lm(np.exec.cowat_fas.yC-np.exec.cowat_fas.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  -9.2469     9.5873  -0.964    0.338
#hiv.x         7.8163    13.1380   0.595    0.554
#age.x         0.1605     0.1778   0.903    0.369
#hiv.x:age.x  -0.1587     0.2535  -0.626    0.533


## ### np.exec.stroop_interf ###
## Uncorrected
fit <- lm(np.exec.stroop_interf.y-np.exec.stroop_interf.x~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)  
#(Intercept) -16.8746    15.1553  -1.113   0.2701  
#hiv.x        48.8261    21.7131   2.249   0.0283 *
#age.x         0.2865     0.2740   1.046   0.3000  
#hiv.x:age.x  -0.9796     0.4154  -2.358   0.0217 *

fname <- "ExeChangeStroopInterfContUncorrected.eps"
postscript(file=fname, width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
yc <-  dall$np.exec.stroop_interf.y-dall$np.exec.stroop_interf.x
nonNA <- !is.na(yc)
sel <- dall$hiv.x == 1 & nonNA
plot(dall$age.x[sel], yc[sel], col=mycols[1], xlim=range(dall$age.x), ylim=range(yc[nonNA]), xlab="Age", ylab="Stroop Interf Change", pch=2)
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


##remove outliers
#dif <- dall$np.exec.stroop_interf.yC - dall$np.exec.stroop_interf.xC
#outl <- which(dif > 20)
##outs <- which(dif < -20)
#out <- c(outl, outs)
#dall <- dall[-out,]

## Corrected
fit <- lm(np.exec.stroop_interf.yC - np.exec.stroop_interf.xC ~ hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)  
#(Intercept) -16.8746    15.4803  -1.090   0.2802  
#hiv.x        46.4306    22.1787   2.093   0.0407 *
#age.x         0.2865     0.2798   1.024   0.3102  
#hiv.x:age.x  -0.9235     0.4243  -2.176   0.0336 *

#dall <- dall[-out,]
fname <- "ExeChangeStroopInterfContCorrectEdu.eps"
postscript(file=fname, width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
yc <-  dall$np.exec.stroop_interf.yC-dall$np.exec.stroop_interf.xC
nonNA <- !is.na(yc)
sel <- dall$hiv.x == 1 & nonNA
plot(dall$age.x[sel], yc[sel], col=mycols[1], xlim=range(dall$age.x), ylim=range(yc[nonNA]), xlab="Age", ylab="Stroop Interf Change", pch=2)
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



## ### np.exec.trail_b_time ###
##3 age groups, before AIC:
fit <- lm(np.exec.trail_b_time.yC-np.exec.trail_b_time.xC~hiv.x*age.x, data=dall)
summary(fit)
#             Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.14854   30.53239   0.005    0.996
#hiv.x       -47.06743   41.77170  -1.127    0.263
#age.x        -0.06286    0.56608  -0.111    0.912
#hiv.x:age.x   1.13157    0.80521   1.405    0.164


## ### np.exec.wais_lnseq ###
##3 age groups, before AIC:
fit <- lm(np.exec.wais_lnseq.yC-np.exec.wais_lnseq.xC~hiv.x*age.x, data=dall)
summary(fit)
#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)  1.38464    2.35313   0.588    0.558
#hiv.x       -0.30669    3.21934  -0.095    0.924
#age.x       -0.01786    0.04363  -0.409    0.683
#hiv.x:age.x -0.01402    0.06206  -0.226    0.822






##3 age groups, before AIC:
## hiv.x*age.grp
## ### np.exec.actionflu_total ###
fit <- lm(np.exec.actionflu_total.yC-np.exec.actionflu_total.xC~hiv.x*age.grp, data=dall)
summary(fit)
#                Estimate Std. Error t value Pr(>|t|)
#(Intercept)       -0.200      1.523  -0.131    0.896
#hiv.x             -2.012      1.766  -1.139    0.258
#age.grp50          1.661      2.026   0.820    0.415
#age.grp60          1.914      2.374   0.806    0.422
#hiv.x:age.grp50    1.240      2.441   0.508    0.613
#hiv.x:age.grp60   -5.153      3.764  -1.369    0.175


## ### np.exec.cowat_fas ###
fit <- lm(np.exec.cowat_fas.yC-np.exec.cowat_fas.xC~hiv.x*age.grp, data=dall)
summary(fit)
#                Estimate Std. Error t value Pr(>|t|)
#(Intercept)       -3.500      2.617  -1.337    0.185
#hiv.x              2.510      3.035   0.827    0.411
#age.grp50          4.731      3.481   1.359    0.178
#age.grp60          3.214      4.078   0.788    0.433
#hiv.x:age.grp50   -5.876      4.212  -1.395    0.167
#hiv.x:age.grp60   -1.441      6.467  -0.223    0.824


## ### np.exec.stroop_interf ###
fit <- lm(np.exec.stroop_interf.yC-np.exec.stroop_interf.xC~hiv.x*age.grp, data=dall)
summary(fit)
#                Estimate Std. Error t value Pr(>|t|)  
#(Intercept)       -6.687      4.824  -1.386   0.1712  
#hiv.x              6.896      5.500   1.254   0.2151  
#age.grp50          6.962      5.908   1.178   0.2436  
#age.grp60          7.626      6.574   1.160   0.2509  
#hiv.x:age.grp50   -9.807      7.155  -1.371   0.1760  
#hiv.x:age.grp60  -25.878     10.954  -2.362   0.0217 *

#dlall <- dlall[-out,]
fname <- "MemChangeStroopinterfTrendErrorBarCorrectEdu.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
par(mfrow = (c(1,2)))
##non hiv plot
plot(c(0, 12), range(na.omit(dlall$np.exec.stroop_interfC)), type="n", ylab="Stroop Interf", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##non-hiv
  yh <- dlall$np.exec.stroop_interfC
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

plot(c(0, 12), range(na.omit(dlall$np.exec.stroop_interfC)), type="n", ylab="Stroop Interf", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##hiv
  yh <- dlall$np.exec.stroop_interfC
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


fit <- lm(np.exec.stroop_interf.y-np.exec.stroop_interf.x~hiv.x*age.grp, data=dall)
summary(fit)
#                Estimate Std. Error t value Pr(>|t|)  
#(Intercept)       -6.687      4.732  -1.413   0.1632  
#hiv.x              7.117      5.396   1.319   0.1925  
#age.grp50          6.962      5.796   1.201   0.2347  
#age.grp60          7.626      6.449   1.183   0.2420  
#hiv.x:age.grp50  -11.264      7.019  -1.605   0.1142  
#hiv.x:age.grp60  -26.099     10.746  -2.429   0.0184 *

fname <- "MemChangeStroopinterfTrendErrorBarUncorrected.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
par(mfrow = (c(1,2)))
##non hiv plot
plot(c(0, 12), range(na.omit(dlall$np.exec.stroop_interfC)), type="n", ylab="Stroop Interf", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##non-hiv
  yh <- dlall$np.exec.stroop_interf
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

plot(c(0, 12), range(na.omit(dlall$np.exec.stroop_interfC)), type="n", ylab="Stroop Interf", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##hiv
  yh <- dlall$np.exec.stroop_interf
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


## ### np.exec.trail_b_time ###
fit <- lm(np.exec.trail_b_time.yC-np.exec.trail_b_time.xC~hiv.x*age.grp, data=dall)
summary(fit)
#                Estimate Std. Error t value Pr(>|t|)  
#(Intercept)       -1.400      8.140  -0.172   0.8639  
#hiv.x             -1.851      9.439  -0.196   0.8450  
#age.grp50         -6.985     10.827  -0.645   0.5207  
#age.grp60          5.257     12.685   0.414   0.6797  
#hiv.x:age.grp50   26.956     13.046   2.066   0.0421 *
#hiv.x:age.grp60   12.344     20.115   0.614   0.5412   


fname <- "MemChangeTrailBTimeCorrectEdu.eps"
postscript(file=fname, width=1.5*mywid, height=mywid)
mytheme()
par(mfrow = (c(1,2)))
##non hiv plot
plot(c(0, 12), range(na.omit(dlall$np.exec.trail_b_timeC)), type="n", ylab="Trail B Time", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##non-hiv
  yh <- dlall$np.exec.trail_b_timeC
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

plot(c(0, 12), range(na.omit(dlall$np.exec.trail_b_timeC)), type="n", ylab="Trail B Time", xlab="Follow-up Month" )
for (gg in 1:3)  {
  gc <-   c("40", "50", "60")[gg]
  ##hiv
  yh <- dlall$np.exec.trail_b_timeC
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


## ### np.exec.wais_lnseq ###
fit <- lm(np.exec.wais_lnseq.yC-np.exec.wais_lnseq.xC~hiv.x*age.grp, data=dall)
summary(fit)
#                Estimate Std. Error t value Pr(>|t|)
#(Intercept)       0.2000     0.6391   0.313    0.755
#hiv.x            -0.4477     0.7411  -0.604    0.547
#age.grp50         0.4154     0.8501   0.489    0.626
#age.grp60         0.2286     0.9959   0.230    0.819
#hiv.x:age.grp50  -1.1533     1.0242  -1.126    0.264
#hiv.x:age.grp60   0.7052     1.5793   0.447    0.656




## 2 age group
## ### np.exec.stroop_interf ###
## Uncorrected
fit <- lm(np.exec.stroop_interf.y-np.exec.stroop_interf.x~hiv.x*age.grp55, data=dall)
summary(fit)
##                  Estimate Std. Error t value Pr(>|t|)  
##(Intercept)         -2.264      2.865  -0.790   0.4327  
##hiv.x                2.549      3.592   0.710   0.4808  
##age.grp5555          2.929      4.775   0.613   0.5421  
##hiv.x:age.grp5555  -13.608      6.488  -2.098   0.0403 *

## Corrected
fit <- lm(np.exec.stroop_interf.yC - np.exec.stroop_interf.xC ~ hiv.x*age.grp55, data=dall)
summary(fit)
##                  Estimate Std. Error t value Pr(>|t|)  
##(Intercept)         -2.264      2.909  -0.778   0.4396  
##hiv.x                2.927      3.646   0.803   0.4253  
##age.grp5555          2.929      4.848   0.604   0.5481  
##hiv.x:age.grp5555  -13.594      6.586  -2.064   0.0435 *


