##% HIV-Aging: Additional Meansure
##% Brown HIV Aging


## # Preprocessing #
d <- read.csv("data.rossi.041312.csv")
d$id <- factor(d$id)
d <- data.frame(d)

## include exec Measure
library(xlsx)
da <- read.xlsx("data.rossi.stroop_CW.xls",1)
da[,3] <- as.numeric(as.character(da[,3]))
da$id <- factor(da$id)
da <- data.frame(da)

d <- merge(d, da, by=c("id","visit"))

## additional measures
## np.domain.ds.attwmexec np.domain.ds.learn   np.domain.ds.mem  np.domain.ds.motor
## np.domain.ds.speed     np.domain.ds.verbal  np.global.gds     np.exec.stroop_cw
grp.add <- c(27:34)

##Remove 40minus subjects.
d <- d[d$age>=40,]

##Remove subjects missing either baseline or 12 month visit.
d1 <- d[d$visit==0, ]                  # dim(d1)=132*34
d2 <- d[d$visit==12,]                  # dim(d2)=92*34
dd <- merge(d1, d2, by=("id"))         # dim(dd)=85*67
dl <- rbind(d1, d2)                    # dim(dl)=224*34 


##calculate changes in the additional measure
aa <- data.frame(as.matrix(dd[,grp.add+33]-dd[,grp.add]))
names(aa) <- paste(names(aa), "change", sep="")
dd <- cbind(dd, aa)

##calculate relative changes in the additional measure
aa <- data.frame(as.matrix( (dd[,grp.add+33]-dd[,grp.add] )/(dd[,grp.add]+1)   ))
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
dall <- merge(dd, dcd, by="id")   # dim(dall)=84*107
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
aa <- data.frame(as.matrix(dlall[,c(grp.add)]))
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
aa <- data.frame(as.matrix(dall[,c(grp.add, grp.add+33)]))
ashift <- (ncol(aa)/2)
for (jj in 1:(ncol(aa)/2) ) {
  sel <- dall$hiv.x==0
  ## correcting for education only 
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

##dim(dall)=84*130
##dall[1,114:129]
## np.domain.ds.attwmexec.xC  np.domain.ds.learn.xC     np.domain.ds.mem.xC     np.domain.ds.motor.xC     
## np.domain.ds.speed.xC      np.domain.ds.verbal.xC    np.global.gds.xC        np.exec.stroop_cw.xC      
## np.domain.ds.attwmexec.yC  np.domain.ds.learn.yC     np.domain.ds.mem.yC     np.domain.ds.motor.yC     
## np.domain.ds.speed.yC      np.domain.ds.verbal.yC    np.global.gds.yC        np.exec.stroop_cw.yC

################################################################
##  #Baseline Analysis
################################################################
## ## MANOVA ##
##Baseline:
fit <- manova(as.matrix(dall[,114:121])~hiv.x*age.x, data=dall)
summary(fit, test="Wilks")
#            Df   Wilks approx F num Df den Df Pr(>F)
#hiv.x        1 0.92305  0.85742      7     72 0.5442
#age.x        1 0.91153  0.99833      7     72 0.4397
#hiv.x:age.x  1 0.94136  0.64072      7     72 0.7208



## ## HIV Positive vs Negative, Univariate ##

## ### np.domain.ds.attwmexec ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.attwmexec.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  0.723224   0.415431   1.741   0.0855 .
##hiv.x       -0.369738   0.568356  -0.651   0.5172  
##age.x       -0.010574   0.007702  -1.373   0.1736  
##hiv.x:age.x  0.007535   0.010956   0.688   0.4936   

## Uncorrected
fit <- lm(np.domain.ds.attwmexec.x ~ hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  0.792394   0.418558   1.893    0.062 .
##hiv.x       -0.420026   0.572634  -0.733    0.465  
##age.x       -0.012166   0.007760  -1.568    0.121  
##hiv.x:age.x  0.008815   0.011038   0.799    0.427 

##Dictomize np.domain.ds.attwmexec.x by 0 cut-off
dall$np.domain.ds.attwmexec.xd <- 0 
dall$np.domain.ds.attwmexec.xd[dall$np.domain.ds.attwmexec.x>0] <- 1

## logistic regression
fit <- glm(np.domain.ds.attwmexec.xd ~ hiv.x*age.x, family = binomial, data= dall)
summary(fit)
##            Estimate Std. Error z value Pr(>|z|)
##(Intercept)  4.88926    3.51279   1.392    0.164
##hiv.x       -4.51454    4.32045  -1.045    0.296
##age.x       -0.11475    0.06993  -1.641    0.101
##hiv.x:age.x  0.08574    0.08661   0.990    0.322


## ### np.domain.ds.learn ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.learn.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)   
##(Intercept)  3.74533    1.26942   2.950  0.00416 **
##hiv.x       -2.47753    1.73671  -1.427  0.15760   
##age.x       -0.04524    0.02354  -1.922  0.05816 . 
##hiv.x:age.x  0.04051    0.03348   1.210  0.22976 

postscript(file="MemBaselineLearn.eps", width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
boxplot(np.domain.ds.learn.xC ~ age.grp*hiv.x, data=dall, col=c(col.n, col.p), names=rep(c("40", "50", "60"), 2), boxwex=0.6, at=c(1, 1.8, 2.6, 4.4, 5.2, 6), ylab="BD")
mtext("Non-HIV", side=1, line=2.5, outer=F, adj=0.2, cex=par("cex.lab"))
mtext("HIV", side=1, line=2.5, outer=F, adj=0.8, cex=par("cex.lab"))
mtext("Age", side=1, line=0.3, outer=F, adj=0, cex=par("cex"))
dev.off()

## 3 age groups, hiv.x*age.grp
fit <- lm(np.domain.ds.learn.xC ~ hiv.x*age.grp, data=dall)
summary(fit)
#                Estimate Std. Error t value Pr(>|t|)    
#(Intercept)       1.9782     0.3437   5.756 1.62e-07 ***
#hiv.x            -0.8553     0.3985  -2.146   0.0350 *  
#age.grp50        -0.7964     0.4571  -1.742   0.0854 .  
#age.grp60        -1.2743     0.5356  -2.379   0.0198 *  
#hiv.x:age.grp50   0.6237     0.5508   1.132   0.2610    
#hiv.x:age.grp60   0.9524     0.8493   1.121   0.2655    

## Uncorrected
fit <- lm(np.domain.ds.learn.x ~ hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)    
##(Intercept)  5.08340    1.38709   3.665 0.000444 ***
##hiv.x       -3.45034    1.89770  -1.818 0.072781 .  
##age.x       -0.07603    0.02572  -2.957 0.004088 ** 
##hiv.x:age.x  0.06528    0.03658   1.784 0.078135 .  

fit <- lm(np.domain.ds.learn.x ~ hiv.x*age.grp, data=dall)
summary(fit)
##                Estimate Std. Error t value Pr(>|t|)    
##(Intercept)       2.2000     0.3644   6.038 4.98e-08 ***
##hiv.x            -0.9241     0.4226  -2.187  0.03174 *  
##age.grp50        -1.6231     0.4847  -3.349  0.00125 ** 
##age.grp60        -1.9857     0.5678  -3.497  0.00078 ***
##hiv.x:age.grp50   1.2790     0.5840   2.190  0.03150 *  
##hiv.x:age.grp60   1.3765     0.9004   1.529  0.13038    


## ### np.domain.ds.mem ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.mem.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  2.69743    1.36890   1.971   0.0522 .
##hiv.x       -1.28581    1.87281  -0.687   0.4943  
##age.x       -0.02683    0.02538  -1.057   0.2936  
##hiv.x:age.x  0.02006    0.03610   0.556   0.5799   

fit <- lm(np.domain.ds.mem.x ~ hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)   
##(Intercept)  4.10564    1.48188   2.771  0.00696 **
##hiv.x       -2.30961    2.02738  -1.139  0.25801   
##age.x       -0.05924    0.02747  -2.156  0.03407 * 
##hiv.x:age.x  0.04612    0.03908   1.180  0.24140   

fit <- lm(np.domain.ds.mem.x ~ hiv.x*age.grp, data=dall)
summary(fit)
##                Estimate Std. Error t value Pr(>|t|)    
##(Intercept)       1.8000     0.3958   4.548 1.96e-05 ***
##hiv.x            -0.4207     0.4590  -0.917  0.36219    
##age.grp50        -1.0692     0.5264  -2.031  0.04566 *  
##age.grp60        -1.6571     0.6168  -2.687  0.00881 ** 
##hiv.x:age.grp50   0.5536     0.6343   0.873  0.38552    
##hiv.x:age.grp60   1.2778     0.9781   1.307  0.19522   


## ### np.domain.ds.motor ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.motor.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##              Estimate Std. Error t value Pr(>|t|)
##(Intercept)  0.4401086  1.0197890   0.432    0.667
##hiv.x       -0.2882270  1.3951271  -0.207    0.837
##age.x        0.0001742  0.0189080   0.009    0.993
##hiv.x:age.x  0.0058979  0.0268935   0.219    0.827 

fit <- lm(np.domain.ds.motor.x ~ hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept)  0.294334   1.025101   0.287    0.775
##hiv.x       -0.182304   1.402394  -0.130    0.897
##age.x        0.003539   0.019007   0.186    0.853
##hiv.x:age.x  0.003191   0.027034   0.118    0.906


## ### np.domain.ds.speed ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.speed.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept)  0.377655   0.532512   0.709    0.480
##hiv.x       -0.585448   0.728536  -0.804    0.424
##age.x       -0.004122   0.009873  -0.417    0.677
##hiv.x:age.x  0.013861   0.014044   0.987    0.327

fit <- lm(np.domain.ds.speed.x ~ hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept)  0.586028   0.527817   1.110    0.270
##hiv.x       -0.736939   0.722112  -1.021    0.311
##age.x       -0.008918   0.009786  -0.911    0.365
##hiv.x:age.x  0.017718   0.013920   1.273    0.207



## ### np.domain.ds.verbal ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.verbal.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept)  0.518677   0.484721   1.070    0.288
##hiv.x        0.119041   0.664237   0.179    0.858
##age.x       -0.004243   0.008987  -0.472    0.638
##hiv.x:age.x -0.005072   0.012817  -0.396    0.693

fit <- lm(np.domain.ds.verbal.x ~ hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept)  0.707518   0.499910   1.415    0.161
##hiv.x       -0.021390   0.685051  -0.031    0.975
##age.x       -0.008589   0.009269  -0.927    0.357
##hiv.x:age.x -0.001501   0.013218  -0.114    0.910


## ### np.global.gds ###
##3 age groups, before AIC:
fit <- lm(np.global.gds.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  1.318702   0.515598   2.558   0.0124 *
##hiv.x       -0.748882   0.705396  -1.062   0.2916  
##age.x       -0.014265   0.009559  -1.492   0.1396  
##hiv.x:age.x  0.012837   0.013598   0.944   0.3480 

fit <- lm(np.global.gds.x ~ hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)   
##(Intercept)  1.79739    0.54647   3.289   0.0015 **
##hiv.x       -1.09690    0.74764  -1.467   0.1463   
##age.x       -0.02528    0.01013  -2.495   0.0146 * 
##hiv.x:age.x  0.02170    0.01441   1.505   0.1362   

fit <- lm(np.global.gds.x ~ hiv.x*age.grp, data=dall)
summary(fit)
##                Estimate Std. Error t value Pr(>|t|)    
##(Intercept)       0.7920     0.1459   5.430 6.18e-07 ***
##hiv.x            -0.1730     0.1691  -1.023  0.30945    
##age.grp50        -0.4497     0.1940  -2.318  0.02307 *  
##age.grp60        -0.6277     0.2273  -2.762  0.00717 ** 
##hiv.x:age.grp50   0.2294     0.2338   0.981  0.32952    
##hiv.x:age.grp60   0.5221     0.3604   1.449  0.15148   


## ### np.exec.stroop_cw ###
##3 age groups, before AIC:
fit <- lm(np.exec.stroop_cw.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)    
##(Intercept)  45.9330    12.0522   3.811 0.000324 ***
##hiv.x       -12.1680    16.9781  -0.717 0.476301    
##age.x        -0.1902     0.2179  -0.873 0.386117    
##hiv.x:age.x   0.1880     0.3228   0.583 0.562307  

fit <- lm(np.exec.stroop_cw.xC ~ hiv.x*age.grp, data=dall)
summary(fit)
##                Estimate Std. Error t value Pr(>|t|)    
##(Intercept)      36.4433     3.7637   9.683 8.48e-14 ***
##hiv.x            -3.2816     4.2913  -0.765    0.447    
##age.grp50         0.5763     4.6096   0.125    0.901    
##age.grp60        -4.2366     5.1291  -0.826    0.412    
##hiv.x:age.grp50   0.7384     5.4973   0.134    0.894    
##hiv.x:age.grp60   2.3302     8.5472   0.273    0.786    


fit <- lm(np.exec.stroop_cw.x ~ hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)   
##(Intercept) 42.15413   12.72464   3.313  0.00156 **
##hiv.x       -9.58898   17.92532  -0.535  0.59464   
##age.x       -0.08694    0.23002  -0.378  0.70676   
##hiv.x:age.x  0.10520    0.34076   0.309  0.75860 

################################################################
## # Change 12-0
################################################################

## ## MANOVA ##
##change:
fit <- manova(as.matrix(dall[,117:123]-dall[,110:116])~hiv.x*age.x, data=dall)
summary(fit, test="Wilks")
##            Df   Wilks approx F num Df den Df  Pr(>F)  
##hiv.x        1 0.94783  0.56613      7     72 0.78097  
##age.x        1 0.84115  1.94238      7     72 0.07518 .
##hiv.x:age.x  1 0.92361  0.85068      7     72 0.54951  
##Residuals   78                                         


## ## HIV Positive vs Negative ##
##Use baseline variables to predict 12-0 changes

## ### np.domain.ds.attwmexec ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.attwmexec.yC - np.domain.ds.attwmexec.xC ~ hiv.x*age.x, data=dall)
summary(fit)
#             Estimate Std. Error t value Pr(>|t|)
#(Intercept) -0.128326   0.391570  -0.328    0.744
#hiv.x        0.119304   0.535711   0.223    0.824
#age.x        0.002196   0.007260   0.303    0.763
#hiv.x:age.x -0.001834   0.010327  -0.178    0.860 


## ### np.domain.ds.learn ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.learn.yC-np.domain.ds.learn.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##              Estimate Std. Error t value Pr(>|t|)  
##(Intercept) -2.3539062  1.1468903  -2.052   0.0434 *
##hiv.x        0.4855828  1.5690734   0.309   0.7578  
##age.x        0.0363687  0.0212638   1.710   0.0911 .
##hiv.x:age.x  0.0001478  0.0302461   0.005   0.9961  

fname <- "AddChangeLearnContcorrected.eps"
postscript(file=fname, width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
yc <-  dall$np.domain.ds.learn.yC-dall$np.domain.ds.learn.xC
nonNA <- !is.na(yc)
sel <- dall$hiv.x == 1 & nonNA
plot(dall$age.x[sel], yc[sel], col=mycols[1], xlim=range(dall$age.x), ylim=range(yc[nonNA]), xlab="Age", ylab="Learn Change", pch=2)
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


## 3 age group, hiv.x*age.grp
fit <- lm(np.domain.ds.learn.yC-np.domain.ds.learn.xC ~ hiv.x*age.grp, data=dall)
summary(fit)
##                Estimate Std. Error t value Pr(>|t|)   
##(Intercept)      -0.9500     0.3127  -3.038  0.00324 **
##hiv.x             0.6741     0.3626   1.859  0.06679 . 
##age.grp50         0.7192     0.4159   1.729  0.08772 . 
##age.grp60         0.9500     0.4873   1.950  0.05483 . 
##hiv.x:age.grp50  -0.2616     0.5011  -0.522  0.60322   
##hiv.x:age.grp60  -0.5075     0.7727  -0.657  0.51328   



## Uncorrected
fit <- lm(np.domain.ds.learn.y-np.domain.ds.learn.x ~ hiv.x*age.x, data=dall)
summary(fit)
#              Estimate Std. Error t value Pr(>|t|)  
#(Intercept) -2.3539062  1.1468903  -2.052   0.0434 *
#hiv.x        0.4855828  1.5690734   0.309   0.7578  
#age.x        0.0363687  0.0212638   1.710   0.0911 .
#hiv.x:age.x  0.0001478  0.0302461   0.005   0.9961  

fname <- "AddChangeLearnContUncorrected.eps"
postscript(file=fname, width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
yc <-  dall$np.domain.ds.learn.y-dall$np.domain.ds.learn.x
nonNA <- !is.na(yc)
sel <- dall$hiv.x == 1 & nonNA
plot(dall$age.x[sel], yc[sel], col=mycols[1], xlim=range(dall$age.x), ylim=range(yc[nonNA]), xlab="Age", ylab="Learn Change", pch=2)
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

## Uncorrected, 3 age group
fit <- lm(np.domain.ds.learn.y-np.domain.ds.learn.x ~ hiv.x*age.grp, data=dall)
summary(fit)
##                Estimate Std. Error t value Pr(>|t|)   
##(Intercept)      -0.9500     0.3127  -3.038  0.00324 **
##hiv.x             0.6741     0.3626   1.859  0.06679 . 
##age.grp50         0.7192     0.4159   1.729  0.08772 . 
##age.grp60         0.9500     0.4873   1.950  0.05483 . 
##hiv.x:age.grp50  -0.2616     0.5011  -0.522  0.60322   
##hiv.x:age.grp60  -0.5075     0.7727  -0.657  0.51328   




## ### np.domain.ds.mem ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.mem.yC-np.domain.ds.mem.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)  
##(Intercept) -1.47507    1.16238  -1.269   0.2081  
##hiv.x       -2.57289    1.59027  -1.618   0.1096  
##age.x        0.02143    0.02155   0.995   0.3229  
##hiv.x:age.x  0.05734    0.03065   1.870   0.0651 .

fname <- "AddChangeMemContCorrected.eps"
postscript(file=fname, width=mywid, height=mywid)
mytheme()
par(oma=c(1,0,0,0))
yc <-  dall$np.domain.ds.mem.yC-dall$np.domain.ds.mem.xC
nonNA <- !is.na(yc)
sel <- dall$hiv.x == 1 & nonNA
plot(dall$age.x[sel], yc[sel], col=mycols[1], xlim=range(dall$age.x), ylim=range(yc[nonNA]), xlab="Age", ylab="Mem Change", pch=2)
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


fit <- lm(np.domain.ds.mem.yC-np.domain.ds.mem.xC ~ hiv.x*age.grp, data=dall)
summary(fit)
##                Estimate Std. Error t value Pr(>|t|)  
##(Intercept)     -0.55000    0.31448  -1.749   0.0842 .
##hiv.x           -0.08793    0.36469  -0.241   0.8101  
##age.grp50        0.20385    0.41829   0.487   0.6274  
##age.grp60        0.55000    0.49008   1.122   0.2652  
##hiv.x:age.grp50  0.75227    0.50401   1.493   0.1396  
##hiv.x:age.grp60  1.08793    0.77713   1.400   0.1655 

fit <- lm(np.domain.ds.mem.yC-np.domain.ds.mem.xC ~ hiv.x*age.grp55, data=dall)
summary(fit)
##                  Estimate Std. Error t value Pr(>|t|)  
##(Intercept)       -0.47619    0.21247  -2.241   0.0278 *
##hiv.x             -0.01131    0.26238  -0.043   0.9657  
##age.grp5555        0.47619    0.38791   1.228   0.2232  
##hiv.x:age.grp5555  0.79702    0.49182   1.621   0.1090 


## ### np.domain.ds.motor ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.motor.yC-np.domain.ds.motor.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept)  0.260801   0.785531   0.332    0.741
##hiv.x       -0.008481   1.074650  -0.008    0.994
##age.x       -0.009432   0.014565  -0.648    0.519
##hiv.x:age.x  0.001699   0.020716   0.082    0.935


## ### np.domain.ds.speed ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.speed.yC-np.domain.ds.speed.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept)  0.100928   0.563983   0.179    0.858
##hiv.x       -0.107411   0.771591  -0.139    0.890
##age.x       -0.002108   0.010456  -0.202    0.841
##hiv.x:age.x  0.001868   0.014874   0.126    0.900

## ### np.domain.ds.verbal ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.verbal.yC-np.domain.ds.verbal.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept) -0.527841   0.561575  -0.940    0.350
##hiv.x        0.781034   0.769554   1.015    0.313
##age.x        0.008032   0.010412   0.771    0.443
##hiv.x:age.x -0.012402   0.014849  -0.835    0.406

## ### np.global.gds ###
##3 age groups, before AIC:
fit <- lm(np.global.gds.yC-np.global.gds.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept) -0.631421   0.438391  -1.440    0.154
##hiv.x       -0.194338   0.599768  -0.324    0.747
##age.x        0.008756   0.008128   1.077    0.285
##hiv.x:age.x  0.006865   0.011561   0.594    0.554


## ### np.exec.stroop_cw ###
##3 age groups, before AIC:
fit <- lm(np.exec.stroop_cw.yC-np.exec.stroop_cw.xC ~ hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)
##(Intercept)  1.51600    8.79625   0.172    0.864
##hiv.x       14.17457   12.39136   1.144    0.257
##age.x       -0.02187    0.15901  -0.138    0.891
##hiv.x:age.x -0.26981    0.23556  -1.145    0.257

fit <- lm(np.exec.stroop_cw.y-np.exec.stroop_cw.x ~ hiv.x*age.x, data=dall)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)
##(Intercept)  1.51600    8.79625   0.172    0.864
##hiv.x       14.17457   12.39136   1.144    0.257
##age.x       -0.02187    0.15901  -0.138    0.891
##hiv.x:age.x -0.26981    0.23556  -1.145    0.257







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
aa <- data.frame(as.matrix(dalll[,c(grp.add)]))
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
dalll$age.grp <- factor(dalll$age.grp)   


################################################################
##  #Baseline Analysis
################################################################
## ## MANOVA ##
##Baseline:
##Without Age.sq
fit <- manova(as.matrix(dalll[,53:60])~hiv*age, data=dalll)
summary(fit, test="Wilks")

## ## HIV Positive vs Negative, Univariate ##

## ### np.domain.ds.attwmexec ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.attwmexecC ~ hiv*age, data=dalll)
summary(fit)
##              Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  0.5207222  0.3109987   1.674   0.0965 .
##hiv          0.0206206  0.4445649   0.046   0.9631  
##age         -0.0062585  0.0058896  -1.063   0.2900 
##hiv:age     -0.0007752  0.0087201  -0.089   0.9293    

## Uncorrected
fit <- lm(np.domain.ds.attwmexec ~ hiv*age, data=dalll)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  0.688818   0.314408   2.191   0.0303 *
##hiv         -0.080495   0.449438  -0.179   0.8581  
##age         -0.009902   0.005954  -1.663   0.0988 .
##hiv:age      0.001732   0.008816   0.196   0.8446  



## ### np.domain.ds.learn ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.learnC ~ hiv*age, data=dalll)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)   
##(Intercept)  3.13294    1.04800   2.989  0.00336 **
##hiv         -1.57644    1.49809  -1.052  0.29466   
##age         -0.03066    0.01985  -1.545  0.12485   
##hiv:age      0.02350    0.02938   0.800  0.42544  

## 3 age groups, hiv.x*age.grp
fit <- lm(np.domain.ds.learnC ~ hiv*age.grp, data=dalll)
summary(fit)
##              Estimate Std. Error t value Pr(>|t|)    
##(Intercept)     1.9667     0.2756   7.137 6.87e-11 ***
##hiv            -0.6878     0.3280  -2.097   0.0380 *  
##age.grp50      -0.6715     0.4066  -1.652   0.1011    
##age.grp60      -0.9199     0.4947  -1.860   0.0653 .  
##hiv:age.grp50   0.5431     0.4929   1.102   0.2726    
##hiv:age.grp60   0.4151     0.8100   0.513   0.6092    

## Uncorrected
fit <- lm(np.domain.ds.learn ~ hiv*age, data=dalll)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)    
##(Intercept)  4.63302    1.14490   4.047 8.97e-05 ***
##hiv         -2.47880    1.63660  -1.515  0.13236    
##age         -0.06318    0.02168  -2.914  0.00422 ** 
##hiv:age      0.04587    0.03210   1.429  0.15551  

fit <- lm(np.domain.ds.learn ~ hiv*age.grp, data=dalll)
summary(fit)
##              Estimate Std. Error t value Pr(>|t|)    
##(Intercept)     2.1750     0.2956   7.359 2.16e-11 ***
##hiv            -0.7271     0.3518  -2.067  0.04081 *  
##age.grp50      -1.2926     0.4360  -2.965  0.00363 ** 
##age.grp60      -1.7861     0.5305  -3.367  0.00101 ** 
##hiv:age.grp50   1.0417     0.5286   1.971  0.05099 .  
##hiv:age.grp60   0.8382     0.8687   0.965  0.33646    



## ### np.domain.ds.mem ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.memC ~ hiv*age, data=dalll)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept)  1.521279   1.083825   1.404    0.163
##hiv          0.040084   1.549300   0.026    0.979
##age         -0.004456   0.020525  -0.217    0.828
##hiv:age     -0.001029   0.030389  -0.034    0.973 

fit <- lm(np.domain.ds.mem ~ hiv*age, data=dalll)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  2.77303    1.16849   2.373   0.0191 *
##hiv         -0.71289    1.67032  -0.427   0.6702  
##age         -0.03159    0.02213  -1.427   0.1559  
##hiv:age      0.01764    0.03276   0.538   0.5913 

fit <- lm(np.domain.ds.mem ~ hiv*age.grp, data=dalll)
summary(fit)
##              Estimate Std. Error t value Pr(>|t|)    
##(Intercept)     1.6750     0.3040   5.509 1.97e-07 ***
##hiv            -0.1438     0.3619  -0.397   0.6919    
##age.grp50      -0.8515     0.4486  -1.898   0.0600 .  
##age.grp60      -1.1750     0.5458  -2.153   0.0332 *  
##hiv:age.grp50   0.5475     0.5438   1.007   0.3160    
##hiv:age.grp60   0.3938     0.8936   0.441   0.6603   


## ### np.domain.ds.motor ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.motorC ~ hiv*age, data=dalll)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept)  0.295660   0.889696   0.332    0.740
##hiv          0.569242   1.271786   0.448    0.655
##age          0.006846   0.016857   0.406    0.685
##hiv:age     -0.012325   0.024952  -0.494    0.622

fit <- lm(np.domain.ds.motor ~ hiv*age, data=dalll)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept)  0.262445   0.890136   0.295    0.769
##hiv          0.589224   1.272416   0.463    0.644
##age          0.007566   0.016866   0.449    0.654
##hiv:age     -0.012821   0.024964  -0.514    0.608


## ### np.domain.ds.speed ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.speedC ~ hiv*age, data=dalll)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)
##(Intercept)  0.438333   0.388846   1.127    0.262
##hiv         -0.255421   0.555846  -0.460    0.647
##age         -0.005109   0.007364  -0.694    0.489
##hiv:age      0.006835   0.010903   0.627    0.532

fit <- lm(np.domain.ds.speed.x ~ hiv.x*age.x, data=dall)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  0.657810   0.389415   1.689   0.0936 .
##hiv         -0.387444   0.556660  -0.696   0.4877  
##age         -0.009866   0.007375  -1.338   0.1833  
##hiv:age      0.010108   0.010919   0.926   0.3563  



## ### np.domain.ds.verbal ###
##3 age groups, before AIC:
fit <- lm(np.domain.ds.verbalC ~ hiv*age, data=dalll)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  0.630363   0.339219   1.858   0.0655 .
##hiv          0.011423   0.485581   0.024   0.9813  
##age         -0.007064   0.006424  -1.100   0.2736  
##hiv:age     -0.002422   0.009532  -0.254   0.7998  

fit <- lm(np.domain.ds.verbal.x ~ hiv.x*age.x, data=dall)
summary(fit)
##              Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  0.7327231  0.3428429   2.137   0.0345 *
##hiv         -0.0513381  0.4907683  -0.105   0.9169  
##age         -0.0092831  0.0064927  -1.430   0.1553  
##hiv:age     -0.0008668  0.0096334  -0.090   0.9285  


## ### np.global.gds ###
##3 age groups, before AIC:
fit <- lm(np.global.gdsC ~ hiv*age, data=dalll)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)  
##(Intercept)  1.035095   0.417280   2.481   0.0144 *
##hiv         -0.188251   0.596492  -0.316   0.7528  
##age         -0.007887   0.007902  -0.998   0.3202  
##hiv:age      0.002349   0.011700   0.201   0.8412 

fit <- lm(np.global.gds ~ hiv*age, data=dalll)
summary(fit)
##             Estimate Std. Error t value Pr(>|t|)    
##(Intercept)  1.537797   0.443825   3.465 0.000724 ***
##hiv         -0.490644   0.634436  -0.773 0.440750    
##age         -0.018783   0.008405  -2.235 0.027185 *  
##hiv:age      0.009847   0.012444   0.791 0.430270    

fit <- lm(np.global.gds ~ hiv*age.grp, data=dalll)
summary(fit)
##              Estimate Std. Error t value Pr(>|t|)    
##(Intercept)     0.8245     0.1144   7.204 4.84e-11 ***
##hiv            -0.1151     0.1362  -0.845  0.39966    
##age.grp50      -0.4239     0.1688  -2.511  0.01333 *  
##age.grp60      -0.5456     0.2054  -2.656  0.00894 ** 
##hiv:age.grp50   0.2027     0.2047   0.990  0.32395    
##hiv:age.grp60   0.2412     0.3364   0.717  0.47463   


## ### np.exec.stroop_cw ###
##3 age groups, before AIC:
fit <- lm(np.exec.stroop_cwC ~ hiv*age, data=dalll)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)    
##(Intercept)  43.1086     8.2887   5.201 1.14e-06 ***
##hiv         -10.8189    12.0232  -0.900    0.370    
##age          -0.1471     0.1519  -0.968    0.335    
##hiv:age       0.1535     0.2322   0.661    0.510   

fit <- lm(np.exec.stroop_cw ~ hiv*age, data=dalll)
summary(fit)
##            Estimate Std. Error t value Pr(>|t|)    
##(Intercept) 38.56872    8.71019   4.428 2.54e-05 ***
##hiv         -9.32503   12.63456  -0.738    0.462    
##age         -0.04107    0.15966  -0.257    0.798    
##hiv:age      0.10352    0.24401   0.424    0.672    









   