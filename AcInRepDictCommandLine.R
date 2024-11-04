library(ppcor)
library(MASS)
library(berryFunctions)
library(ggplot2)
library(nlme )  # library nlme includes the ctrl <- lmeControl(opt='optim') option, 
# which deals with some convergence issues, as per
# https://stats.stackexchange.com/questions/40647/lme-error-iteration-limit-reached
# Hence use function lme with option control=ctrl, where 
# ctrl <- lmeControl(opt='optim') .


try( source("/media/sf_mmpsy/Dropbox/FIL_aux/R_scripts/nspn_utils.R"))
try( source( "C:/Users/mmpsy/Dropbox/FIL_aux/R_scripts/nspn_utils.R"))


### ########  putting together the key data object, abFit ######### ####
##            ------>  code block moved to nr. end of file ...      ####

## Exploring Liam idea - pts. w. diff.erent aInitEv at bsl -------------
bslKeyTaskHd <- c("UserID","pH0.A", "pS0.A", "dInitEv.A", "aInitEv.A", "alphaPrec.A", 
                  "mem.A", "wH.A" , "wS.A" , "w0.A" ,     "othLR.A"  ,  "BIC.A",           # first 12 are ID and model measures
                  "pred1.A", "HI1.A",  "SI1.A",   "predAv.A",   "HIAv.A",  "SIAv.A");      # 13-18 are descriptives.
fuKeyTaskHd  <-  c("UserID","pH0.B", "pS0.B", "dInitEv.B", "aInitEv.B", "alphaPrec.B", 
                   "mem.B", "wH.B" , "wS.B" , "w0.B" ,     "othLR.B"  , "BIC.B",
                   "pred1.B","HI1.B",  "SI1.B",   "predAv.B",   "HIAv.B",  "SIAv.B" );

V <- 'pHO';
#  yw is the wide format data, could refresh from allD$wiDat[[1]]
VA <- paste(V,'A',sep='.'); VB <- paste(V,'B',sep='.'); 
medVA  <- median(na.omit(yw[,VA]))
keyHi <- vecTRUE(yw[,VA] >= medVA); 
keyLo <- vecTRUE(yw[,VA] < medVA); 
drugHi <- vecTRUE((yw[,VA] >= medVA)  & (yw[,'placebo0SSRI1.A'] == 1)) 
drugLo <- vecTRUE((yw[,VA]  < medVA)  & (yw[,'placebo0SSRI1.A'] == 1)) 
placLo <- vecTRUE((yw[,VA]  < medVA)  & (yw[,'placebo0SSRI1.A'] == 0)) 
placHi <- vecTRUE((yw[,VA] >= medVA)  & (yw[,'placebo0SSRI1.A'] == 0)) 

boxplot(yw[drugLo,VA] ,
        yw[drugLo,VB] ,
        yw[drugHi,VA] ,
        yw[drugHi,VB] ,
        yw[placLo,VA] ,
        yw[placLo,VB] ,
        yw[placHi,VA] ,
        yw[placHi,VB] ,
        names= c('drgLo,bsl','drgLo,f/u', 'drgHi.,bsl','drgHi.,f/u',
                 'plcLo,bsl','plcLo,f/u', 'plcHi.,bsl','plcHi.,f/u'),
        col = c('green1','green1','green4','green4','cyan1','cyan1','cyan4','cyan4'),
        main = paste(V,'follow-up - ',V,'baseline'),
        notch=T)

boxplot(yw[drugLo | drugHi,VA],
        yw[placLo | placHi,VA],
        yw[drugHi,VB]-yw[drugHi,VA] ,
        yw[placHi,VB]-yw[placHi,VA] ,
        yw[drugLo,VB]-yw[drugLo,VA] ,
        yw[placLo,VB]-yw[placLo,VA] ,
        names= c("drgBsl",'plcBsl','drgHi,fu-bsl','drgLo,fu-bsl','plcHi,fu-bsl','plcLo,fu-bsl'),
        col=c('green3','cyan3','green1','green4','cyan1','cyan4'),
        main = paste(V,': baseline and change, SSRI vs. placebo'),
        notch=T); abline(v=2.5)

## Simple factor analysis and scores over psychometrics ----------------
library(psych)
library(missMDA) # This loads also required package FactoMineR
library(nFactors)
library(lavaan)
library(sem)
library(corrplot)
library(caret)
library(car)
library(GPArotation)
library(fungible)

# key (not all) self-report and dispositional headings
bslKeyPsHd <- c("UserID","PHQ.A","GAD.A","MASQ.A", "AMI.A", "PSS.A", "RSE.A");
fuKeyPsyHd  <- c("UserID","PHQ.B","GAD.B","MASQ.B", "AMI.B", "PSS.B", "RSE.B");
keyDispHd   <- c("UserID","HPS","subjSES","ethnCode","sexF1M2","Weight.A"); 

bslKeyPsy <- allD$psy[,bslKeyPsHd];   
bslKeyPsy <- data.frame(1:dim(bslKeyPsy)[1],  allD$psy[,bslKeyPsHd] )
colnames(bslKeyPsy)[1] <- 'keyPsyID';
# There is no point interpolating missing values :
ps1 <- na.omit(bslKeyPsy)
y <- ps1[,3:8]
# See CFA-decAc-UCHANGE.R for more on the below ...
ev <- eigen(cor(y)) # get eigenvalues
ap <- parallel(subject=nrow(y),var=ncol(y),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS, main='Scree Test - baseline state-like psy measures'); 

## Run the factor analysis on no. of factors suggested by the scree plot
factN <- nS[[1]]$noc; # see nS[[1]]
nocfit<-fa(r=cor(y), nfactors=factN,fm="minres"); nocfit$loadings
# Minimal EFA with just one factor:
fit1<-fa(r=cor(y), nfactors=1,fm="minres"); fit1$loadings

# Tidy up the scores, first baseline:
faScor1 <- faScores(y,Loadings = nocfit$loadings, Phi=nocfit$Phi);
bslFASc <- data.frame(ps1[,'keyPsyID'],faScor1$fscores);
colnames(bslFASc) <- c('keyPsyID', 'anxDep.A','masqAmi.A'); 
print(head(bslFASc));

# Now follow-up scores based on the same loadings etc:
fuKeyPsy <- allD$psy[,fuKeyPsyHd];   
fuKeyPsy <- data.frame(1:dim(fuKeyPsy)[1],  allD$psy[,fuKeyPsyHd] )
colnames(fuKeyPsy)[1] <- 'keyPsyID';
# There is no point interpolating missing values :
ps2 <- na.omit(fuKeyPsy);
faScor2 <- faScores(ps2[,3:8],Loadings = nocfit$loadings, Phi=nocfit$Phi);
fuFASc <- data.frame(ps2[,'keyPsyID'],faScor2$fscores);
colnames(fuFASc) <- c('keyPsyID', 'anxDep.B','masqAmi.B'); 
print(head(fuFASc)); 

keyPsy <- allD$psy[,c(bslKeyPsHd,fuKeyPsyHd[2:length(fuKeyPsyHd)], keyDispHd[2:length(keyDispHd)])]
keyPsy <- data.frame(1:dim(keyPsy)[1],keyPsy);
colnames(keyPsy)[1] <- 'keyPsyID';
test <- merge(keyPsy,bslFASc,by.x='keyPsyID',by.y='keyPsyID',all=TRUE);
test <- merge(test,fuFASc,by.x='keyPsyID',by.y='keyPsyID',all=TRUE);
psyD <- test;  remove(test); 

allD$psyFA$psyD <- psyD;  allD$psyFA$nocfit <- nocfit; 
allD$psyFA$faScore.A <- faScor1; allD$psyFA$faScore.B <- faScor2


## FA over params                                        ----------------
bslKeyTaskHd <- c("UserID","pH0.A", "pS0.A", "dInitEv.A", "aInitEv.A", "alphaPrec.A", 
                  "mem.A", "wH.A" , "wS.A" , "w0.A" ,     "othLR.A"  ,  "BIC.A",           # first 12 are ID and model measures
                  "pred1.A", "HI1.A",  "SI1.A",   "predAv.A",   "HIAv.A",  "SIAv.A");      # 13-18 are descriptives.
fuKeyTaskHd  <-  c("UserID","pH0.B", "pS0.B", "dInitEv.B", "aInitEv.B", "alphaPrec.B", 
                   "mem.B", "wH.B" , "wS.B" , "w0.B" ,     "othLR.B"  , "BIC.B",
                   "pred1.B","HI1.B",  "SI1.B",   "predAv.B",   "HIAv.B",  "SIAv.B" );

bslKeyMod <- allD$wiDat[[1]][,bslKeyTaskHd[1:11]];     # ID and model params
bslKeyMod <- data.frame(1:dim(bslKeyMod)[1],  allD$wiDat[[1]][,bslKeyTaskHd[1:11]] )
colnames(bslKeyMod)[1] <- 'keyPsyID';
# There is no point interpolating missing values :
md1 <- na.omit(bslKeyMod)
y <- md1[,3:dim(md1)[2]]
# See CFA-decAc-UCHANGE.R for more on the below ...
ev <- eigen(cor(y)) # get eigenvalues
ap <- parallel(subject=nrow(y),var=ncol(y),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS, main='Scree Test - baseline task parameters'); 

## Run the factor analysis on no. of factors suggested by the f/u scree plot
mFactN <- 2   # or nS[[1]]$noc; # see nS[[1]]
efaM<-fa(r=cor(y), nfactors=mFactN,fm="minres"); efaM$loadings
# Minimal EFA with just one factor:
fitM1<-fa(r=cor(y), nfactors=1,fm="minres"); fitM1$loadings

# Tidy up the scores, first baseline:
faMScor1 <- faScores(y,Loadings = efaM$loadings, Phi=efaM$Phi);
bslMFASc <- data.frame(md1[,'keyPsyID'],faMScor1$fscores);
colnames(bslMFASc) <- c('keyPsyID', 'wHwS','pS0pH0'); 
print(head(bslMFASc));

# Now follow-up scores based on own or baseline loadings etc:
fuKeyMod <- allD$wiDat[[1]][,fuKeyTaskHd[1:11]];     # ID and model params
fuKeyMod <- data.frame(1:dim(fuKeyMod)[1],  allD$wiDat[[1]][,fuKeyTaskHd[1:11]] )
colnames(fuKeyMod)[1] <- 'keyPsyID';
# There is no point interpolating missing values :
md2 <- na.omit(fuKeyMod);

y <- md2[,3:dim(md2)[2]]
# See CFA-decAc-UCHANGE.R for more on the below ...
ev <- eigen(cor(y)) # get eigenvalues
ap <- parallel(subject=nrow(y),var=ncol(y),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS, main='Scree Test - follow-up task parameters'); 

## Run the factor analysis on no. of factors suggested by the f/u scree plot
mFactN2 <- 2; # nS[[1]]$noc; # see nS[[1]]
efaM2<-fa(r=cor(y), nfactors=mFactN2,fm="minres"); efaM2$loadings
# Minimal EFA with just one factor:
fitM1fu<-fa(r=cor(y), nfactors=1,fm="minres"); fitM1fu$loadings



faMScor2b <- faScores(md2[,3:dim(md2)[2]],Loadings = efaM$loadings, Phi=efaM$Phi);
fuMFAScb <- data.frame(md2[,'keyPsyID'],faMScor2b$fscores);
colnames(fuMFAScb) <- c('keyPsyID', 'wHwS','pS0pH0');  # c('keyPsyID', 'OthPol.B','pSH0.B','BelUpd.B','Mem.B'); 
print(head(fuMFAScb)); 

keyPsy <- allD$psy[,c(bslKeyPsHd,fuKeyPsyHd[2:length(fuKeyPsyHd)], keyDispHd[2:length(keyDispHd)])]
keyPsy <- data.frame(1:dim(keyPsy)[1],keyPsy);
colnames(keyPsy)[1] <- 'keyPsyID';
test <- merge(keyPsy,bslFASc,by.x='keyPsyID',by.y='keyPsyID',all=TRUE);
test <- merge(test,fuFASc,by.x='keyPsyID',by.y='keyPsyID',all=TRUE);
psyD <- test;  remove(test); 

allD$psyFA$psyD <- psyD;  allD$psyFA$nocfit <- nocfit; 
allD$psyFA$faScore.A <- faScor1; allD$psyFA$faScore.B <- faScor2


## Statistics and plots for Results ------------------------------------
# plot margins for if the get screwed up:
par(mar=c(5,5,4,2))  #  down, left, up, right (ffs ...)
bakPar <- par()

# Citalopram and attributions
taskPlacebs <- vecTRUE(allD$clean[,"placebo0SSRI1","qA"] == 0)
taskSSRIs <- vecTRUE(allD$clean[,"placebo0SSRI1","qA"] == 1)
mAPredAv <- mean(na.omit(allD$clean[,"predAv","qA"]));
sdAPredAv <-  sd(na.omit(allD$clean[,"predAv","qA"]));
mApH0 <-  mean(na.omit(allD$clean[,"pH0","qA"]));
sdApH0 <- sd(na.omit(allD$clean[,"pH0","qA"]));
boxplot((allD$clean[,"predAv","qA"] - mAPredAv)/sdAPredAv, 
        (allD$clean[taskSSRIs,"predAv","qA"] - mAPredAv)/sdAPredAv, 
        (allD$clean[taskPlacebs,"predAv","qB"] - mAPredAv)/sdAPredAv, 
        (allD$clean[taskSSRIs,"predAv","qB"] - mAPredAv)/sdAPredAv, 
        
        # (allD$clean[,"pH0","qA"] -  mApH0)/sdApH0,
        # (allD$clean[taskSSRIs,"pH0","qA"] -  mApH0)/sdApH0,
        # (allD$clean[taskPlacebs,"pH0","qB"] -  mApH0)/sdApH0,
        # (allD$clean[taskSSRIs,"pH0","qB"] -  mApH0)/sdApH0,
        # 
        # (allD$clean[,"pS0","qA"] -  mApH0)/sdApH0,
        # (allD$clean[taskSSRIs,"pS0","qA"] -  mApH0)/sdApH0,
        # (allD$clean[taskPlacebs,"pS0","qB"] -  mApH0)/sdApH0,
        # (allD$clean[taskSSRIs,"pS0","qB"] -  mApH0)/sdApH0,
        
        col=c('gold4','gold3','gold2','yellow'), # ,'cyan4','cyan4','cyan3','green','purple4','purple4','purple2','green3'),
        ylab = "Predicted average fairness (scaled to baseline)",
        names=rep(c('plac BSL','Cital BSL','plac FU','Cital FU'),1),
        # xlab = "    av. prediction#, #                                       pHI0                                                pSI0",
        main = 'Key measures of attributed motivation',
        notch = TRUE)
# abline(v=3.5); abline(v=6.5); 

# Exploring psychometrics
explorePsych <- 1;
#
if (explorePsych != 0){
  par(mar=c(4.1,8,2.5,1))
  DV <- c('AMI','PHQ','RSE','MASQ','GAD','PSS');
  IV <- c('sexF1M2','subjSES','HPS');
  ses <- allD$psy[,IV[2]];
  sex <- 1*(allD$psy[,"Gender.A"] == 'Female');
  hps <- allD$psy[,IV[3]];
  ssri <- allD$psy[,"placebo0SSRI1"];
  olsrs <- list();    # to hold a bunch of ols regressions
                      #  FollowUp ~ Baseline + ssri + ses + sex + hps
  placRo <- vecTRUE(allD$psy[,"placebo0SSRI1"]==0);
  ssriRo <- vecTRUE(allD$psy[,"placebo0SSRI1"]==1);
  dPlac <- matrix(NA,sum(placRo),length(DV));
  dSSRI <- matrix(NA,sum(ssriRo),length(DV));
  colnames(dPlac) <- paste('d',DV,'.plac',sep='');
  colnames(dSSRI) <- paste('d',DV,'.Cital',sep='');
  for (DVi in 1:length(DV)){ 
    v <- DV[DVi]; vA <- paste(v,'A',sep='.'); vB <- paste(v,'B',sep='.');
    sdA <- sd(na.omit(allD$psy[,vA]));
    dPlac[,DVi] <- (allD$psy[placRo,vB] - allD$psy[placRo,vA])/sdA; 
    dSSRI[,DVi] <- (allD$psy[ssriRo,vB] - allD$psy[ssriRo,vA])/sdA; 
    
    Baseline <- allD$psy[,vA];  FollowUp <-  allD$psy[,vB];
    olsrs[[DVi]] <- lm(FollowUp ~ Baseline + ssri + ses + sex + hps ); 
    olsrs[[DVi]]$name <- DV[DVi];
    ssriDependency <- (summary(olsrs[[DVi]])$coeff['ssri',])
    ti <- paste(paste(c('beta:','; SE_b:',';\n t_val:','; p:'),
                round(ssriDependency,4),sep=''), collapse='');
    ti <- paste(v,'~ drug:',ti,coll='')
    
    plot(Baseline,FollowUp,pch=21,bg='blue',col='cyan',main=ti); lines(Baseline[ssriRo],FollowUp[ssriRo],t='p',pch=21,bg='green',col='green4');  abline(0,1,lwd=3,col='gray30');
    
  }
  
  nam <- paste(repAdjVec(DV,2),
               rep(c('Cital.','Placebo'),length(DV)))
  boxplot( dSSRI[,1], dPlac[,1],
           dSSRI[,2], dPlac[,2],
           dSSRI[,3], dPlac[,3],
           dSSRI[,4], dPlac[,4],
           dSSRI[,5], dPlac[,5],
           dSSRI[,6], dPlac[,6],
           #  main=
          col = rep(c('green3','cyan4'),6),
          ylim = c(-2,2),
          horizontal = TRUE,
          notch=TRUE,
          xaxt='n',yaxt='n')
  axis(1); 
  axis(2,at=1:length(nam),labels=nam,las=2)
  title(main = 'Follow-up - Baseline change in psychometric scores', xlab='Change, in units of baseline SD')
  abline(v=0,col='gray30',lwd=3);
 
}

### Income Distribution
plotIncomeDistr <- 0; 
if (plotIncomeDistr){
  bakPar <- par()
  par(mar=c(3.1,8,2.5,1))
  incHist <- horizHist(SSRISelfRep1[,"income"], 
                       (breaks = 1:9-0.5),
                       xaxt='n',yaxt='n',
                       col=c(rep('gold3',7),'gray30'),
                       xlim=c(0,25),
                       main='Income Distribution'); 
  axis(1); axis(2,at=((1:8)-0.5),labels=incomeCats,las=2)
}

### Ethnicity
plotEthnDistr <- 1; 
if (plotEthnDistr){
  bakPar <- par()
  par(mar=c(3.1,13,2.5,1))  #margins: down, left, up, right (ffs ...)
  
  ethnHist <- horizHist(SSRISelfRep1[,"ethnCode"], 
                       breaks = (1:19)-0.5,
                       xaxt='n',yaxt='n',
                       col=c(rep('lavenderblush',3),
                             rep('gray85',2),
                             rep('wheat1',2),
                             rep('brown',4),
                             rep('wheat3',2),
                             rep('gray10',3),
                             rep('brown1',1),
                             rep('wheat2',1)
                       ),
                       xlim=c(0,20),
                       main='Ethnicity Distribution'); 
  axis(1); 
  axis(2,at=((1:18)-0.5),labels=ethnCat,las=2 )
  par <- backPar
}

### Subjective Socioecon
plotSubjSES <- 0; 
if (plotSubjSES){
  bakPar <- par()
  par(mar=c(6,5,3,5))   # margins for labels etc.
  ssesHist <- horizHist(SSRISelfRep1[,"subjSES"], 
                        breaks = (1:11)-0.5,
                        xaxt='n',yaxt='n',
                        col='orchid4',
                        main='Subjective Socio-economic status',
                        ylim=c(0.99,10.01),
                        xlab = 'participant count',
                        ylab='McArthur  ladder  rung'); 
  axis(1); axis(2,at=((1:11)-0.5),labels=1:11,las=2)
}
## Correlations between fitted parameters ----------------------------

### The loop below plots correlations between fitted parameters (here 
### using column names for model q
d1 <- 'qB'; 
print(round(spearAB[[d1]]$P,4)); 
hd <- dimnames(allD$clean)[[2]][c(8:11,13:18,3)]; 
for (j in 1:10) { 
  for (k in (j+1):11) { 
    v1 <- hd[k] ; v2 <- hd[j] ;  
    xyd <- allD$clean[,c(v2,v1,'placebo0SSRI1'),d1]; 
    robr <- rlm(xyd[,v2] ~ xyd[,v1] + xyd[,'placebo0SSRI1']);  
    surr <- summary(robr);
    vBeta <- surr$coeff[2,1]; drugBeta <- surr$coeff[3,1];
    vTVal <- surr$coeff[2,3]; drugTVal <- surr$coeff[3,3];
    mainTi <- paste(d1,' robust regr.:  ',paste(v2,'~',v1,'+','drug  Groups: red->SSRI, blue->placebo'),
                    '\n',v1,': beta=',format(vBeta,digits=4),
                    ', t=',format(vTVal,digits=4),
                    '     drug: beta=',format(drugBeta,digits=4),
                    ', t=',format(drugTVal,digits=4),
                    sep='')
    plot(xyd[,v1],xyd[,v2], main = mainTi,
         pch=21,bg='navy',col='cyan',xlab=v1,ylab=v2); 
    lines(xyd[vecTRUE(xyd[,3]==1),v1],xyd[vecTRUE(xyd[,3]==1),v2], 
          t='p',pch=21,bg='red',col='cyan');  abline(0,1,col='gray60',lwd=2); 
    cat(paste('\n',v2,'~',v1,'+','drug\n')); 
    b <- robr$coefficients; 
    lines(xyd[,v1],xyd[,v1]*b[2] + b[1],col='gray30',lwd=2); 
    print(summary(robr)$coeff[2:3,c(1,3)]) 
  } 
}
## plots of regrn / stability bet. par.  ----------------------------

print(dimnames(allD$clean)[[2]],qu=F);
v <- 'LL'; d1='xA'; d2='xB'; 
cat('\nNow to examine stability and change of:'); print('  v       d1    d2',qu=F);print(paste(v,'for',d1,'vs',d2),qu=F)

# Use data version cleaned by hand.
xygCor <- allD$clean[,c(v,'wave','placebo0SSRI1'),d1]
xygCor[,2] <- allD$clean[,v,d2]
colnames(xygCor)[c(1,2)] <- paste(v,1:2,sep='')
# All data, incl. outliers:
xygUnc <- allD$abWdrug[,c(v,'wave','placebo0SSRI1'),d1]
xygUnc[,2] <- allD$abWdrug[,v,d2]
colnames(xygUnc)[c(1,2)] <- paste(v,1:2,sep='')

xc <- xygCor[,1]; yc <- xygCor[,2]; gc <- xygCor[,3]; 
xu <- xygUnc[,1]; yu <- xygUnc[,2]; gu <- xygUnc[,3]; 

# Do the basic ols regressions. The robust regression is applied to the 
# data including outliers etc, but the ordinary regr to the cleaned up
# data, to see if the latter in some way cheats in terms of significance
# Extract # the overall p-value and adjusted proportion of variance explained: 

colsr <- lm(yc~xc+gc);    # Ordinary regression, cleaned data
robr  <- rlm(yu~xu+gu);   # robust regression, all data.

f <- summary(colsr)$fstatistic;
pcolsr <- round(pf(f[1],f[2],f[3],lower.tail=F),6); 
adjrsq <- round(summary(colsr)$adj.r.squared,3) ; 

xy <- as.data.frame(na.omit(xygUnc[,1:2]));   x <- xy[,1]; y <- xy[,2];
mainTi <- paste(v,'from models',d1,'(x) vs.', d2, '(y)');
subTi  <- paste('RSE incl. drug: ols:',round(summary(colsr)$sigma,3),
                ';  robust:',round(summary(robr)$sigma,3),
                '\np val. ols:',round(pcolsr,5),
                '; adj.R.sq:', adjrsq);
capt <- 'Blue: robust regression; grey: identity';
ggplot( xy , aes(x = x, y = y)) + 
  geom_point() +
  stat_smooth(method = "rlm")+
  labs(x = 'Baseline', y='Follow-up') +
  geom_abline(size=1.5, slope=1, intercept=0,  col='gray60',
              mapping = NULL, data = NULL,  na.rm = FALSE, 
              show.legend = NA  )+
  labs(title = mainTi,subtitle=subTi, caption = capt );

txtFileName <- paste(v,'_',d2,'_vs_',d1,'+drug.txt',sep='')
sink(txtFileName); cat(paste(v,d2,'vs',d1,'+ group;\nBelow, g for group (0=placebo); c for outliers excluded, u for all included.\n')); print(summary(colsr)); print(summary(robr)); sink();
cat('Robust Regression incl. outliers:\n'); print(summary(robr)); cat('---------------------\n'); cat('Ordinary regression excl. outliers:\n'); print(summary(colsr))


# Back up : 
save.image("/media/sf_mmpsy/Dropbox/BASOR/AcInSOR/AcInRepeatedDictator/analyses/.RData")

# Illustrate placebo vs. SSRI:
v1plac <- allD$clean[vecTRUE(allD$clean[,'placebo0SSRI1',d1] == 0),v,d1]; 
v2plac <- allD$clean[vecTRUE(allD$clean[,'placebo0SSRI1',d2] == 0),v,d2]; 
dvPlac <- v2plac - v1plac; 
hist(dvPlac,15,col='cyan',main=paste('change in',v,'under placebo')); 
dMedPlac <- median(na.omit(dvPlac)); 
abline(v=dMedPlac,lwd=3,col='red3'); 
v1ssri <- allD$clean[vecTRUE(allD$clean[,'placebo0SSRI1',d1] == 1),v,d1]; 
v2ssri <- allD$clean[vecTRUE(allD$clean[,'placebo0SSRI1',d2] == 1),v,d2]; 
dvSSRI <- v2ssri - v1ssri; 
hist(dvSSRI,15,col='pink3',main=paste('change in',v,'under SSRI')); 
dMedSSRI <- median(na.omit(dvSSRI)); 
abline(v=dMedSSRI,lwd=3,col='red3'); 

### BIC stuff for some key model checks and comparisons ------------------
Da <- 'kB'; Db <- 'qB';
v <- 'BIC';
D <- allD$abWdrug[,v,c(Da,Db)];
pwilc <- wilcox.test(D[,1]-D[,2])
mediDif <- median(na.omit(D[,1] - D[,2]))
pv <-pwilc$p.value
mainTi <- paste('Median dif. bsl- f/u:',round(mediDif,3),'; p val:',round(pv,6) )
plot(D[,1],D[,2],
     pch=21,bg='blue3',col='white',lwd=2,
     xlab=paste(v,Da),
     ylab=paste(v,Db),
     main = mainTi); 
abline(6,1,col='gray20'); abline(-6,1,col='grey20'); abline(0,1,col='red',lw=2)
hist(D[,1] - D[,2],30,col='cyan',main=paste('Difference in',v,'(',Da,'-',Db,'). Blue line=median.\nMedian=',round(mediDif,3),';  Wilcox p=',format(pv,digits=4)), xlab=paste('delta',v)); abline(v=mediDif,col='navy',lwd=5)


### ######  Basic POC work 
pocWork <- 0
if (pocWork != 0) {
v <- 'retAv';
yw <- yl[,c('intID','wave','othEthn','ethnOrd',v)]; 
yw <- {yw %>% pivot_wider(names_from=c('othEthn','wave','ethnOrd'), # columns whose combo of entries uniquely specifies measurement 
                          values_from = v)}; # the measurement in question.
yw$whitish0 <- (yw$whitish_0_1+yw$whitish_0_2)/2;  yw$whitish1 <- (yw$whitish_1_1+yw$whitish_1_2)/2
yw$whitish <- ( yw$whitish0 + yw$whitish1) / 2
for (ro in 1:dim(yw)[1]){
  if (is.na(yw[ro,'whitish0'])) {yw[ro,'whitish'] <- yw[ro,'whitish1'] }
  if (is.na(yw[ro,'whitish1'])) {yw[ro,'whitish'] <- yw[ro,'whitish0'] }   
}
yw$poc0 <- (yw$poc_0_1+yw$poc_0_2)/2;  yw$poc1 <- (yw$poc_1_1+yw$poc_1_2)/2
yw$poc <- ( yw$poc0 + yw$poc1) / 2
for (ro in 1:dim(yw)[1]){
  if (is.na(yw[ro,'poc0'])) {yw[ro,'poc'] <- yw[ro,'poc1'] }
  if (is.na(yw[ro,'poc1'])) {yw[ro,'poc'] <- yw[ro,'poc0'] }   
}
# Simple estimates of within-pt standard deviation at baseline
# and follow-up:
yw$sd0 <- 1*NA;      yw$sd1 <- 1*NA;
for (ro in 1:dim(yw)[1]){
  x <- as.numeric(yw[ro,c('whitish_0_1','whitish_0_2','poc_0_1','poc_0_2')])
  yw[ro,'sd0'] <- sqrt(var(x));
  x <- as.numeric(yw[ro,c('whitish_1_1','whitish_1_2','poc_1_1','poc_1_2')])
  yw[ro,'sd1'] <- sqrt(var(x));
}
yw$whiPocAvg0 <-(yw$whitish0 + yw$poc0)/2;   yw$whiPocAvg1 <-(yw$whitish1 + yw$poc1)/2;
yw <- as.data.frame(yw);
# View(yw)
stab <- data.frame(na.omit(yw[,c('whiPocAvg0','whiPocAvg1')]));  co <- pcor(stab); 
# stab <- data.frame(na.omit(yw[,c('whitish0','whitish1')]));  co <- pcor(stab); 
print(paste('Stability basics for ',v),quote=F); print('cor:      p:',quote=F); 
print(round(c(co$est[1,2],co$p.val[1,2]),4))
plot(stab[,1],stab[,2],xlab=paste(v,'bsl'),ylab=paste(v,'f/u')); abline(0,1)

boxplot(yw$whitish0,yw$poc0,yw$whitish1,yw$poc1,yw$whitish,yw$poc, 
        (yw$whitish0+yw$poc0)/2, (yw$whitish1+yw$poc1)/2, notch=T, 
        col=c('pink3','brown3','pink4','brown4','pink','brown','gray90','gray60'),
        names =c('whi1','poc1','whi2','poc2','whi','poc','wav1','wav2'),  main=v)

} # end if to do basic POC work

### #####  Demo of categorical learning resembling Alex Pike 
#          exponentially decreasing learning rate
catlrn <- 0            # lines so that rest of file
if (catLrn != 0) {     # can be run easily.

categoLrn <- function(bIn, ph=0.9, tau=2) {
  totT <- length(bIn)
  pl <- 1-ph    # ph is the prob. of +ve feedback if 
                # the best option is chosen.
  Ph <- rep(NA,totT+1);    Pl <- Ph;     
  res <- matrix(NA,totT,7);
  colnames(res) <- c('r','PhPri','dQ','act','PEL','PEH','lr');
  # initial values:
  Ph[1]  <- 0.5
  Rcor <- 1
  Rinc <- 0
  Q <- c(0.5,0.5); names(Q) <- c('lo', 'hi');
  
  if (sum(abs(bIn * (1-bIn))) > 1e-6){
    stop('bIn to have only 1s and 0s please');
  }
  
  for (t in 2:(totT+1)){
      # Action values before choice or outcome:
      Q['hi'] <-  Ph[t-1]   *Rcor + (1-Ph[t-1])*Rinc
      Q['lo'] <- (1-Ph[t-1])*Rcor +    Ph[t-1] *Rinc
      dQ <- diff(Q)
      # Resulting policy :
      piH <- 1/(1+exp(-dQ/tau))
      # Action chosen:
      act <- sample( c(1,2), 1, repl=T, c(1-piH,piH) )
      res[t-1,'act'] <- act
      
      # Action-value Prediction Errors - experienced and counterfactual:
      if (bIn[t-1] > 0.5){ # observation was H
           # experienced  (or counter-factual) PE for L
           PEL <- Rinc - Q['lo'] 
           # experienced  (or counter-factual) PE for H
           PEH <- Rcor - Q['hi'] 
      } else {             # observation was L
           PEL <- Rcor - Q['lo']
           PEH <- Rinc - Q['hi']
      }
      
      # Belief update
      if (bIn[t-1] > 0.5) {
        Ph[t] <- ph*Ph[t-1] / (ph*Ph[t-1] + pl*(1-Ph[t-1]))
       } else {
        Ph[t] <- pl*Ph[t-1] / (pl*Ph[t-1] + ph*(1-Ph[t-1]))
       }
 
      # Apparent learning rate is the one which would give 
      # the correct *next* difference in action values 
      # Q(hi,t)-Q(lo,t) = Q(hi,t-1) + lr*PEH - (Q(lo,t-1) + lr*PEL) therefore:
      lr <- ( Ph[t]*Rcor+(1-Ph[t])*Rinc - (1-Ph[t])*Rcor+Ph[t]*Rinc - dQ) /
            ( PEH - PEL)
      
      res[t-1,'r'] <- bIn[t-1]
      res[t-1,'PhPri'] <- Ph[t-1]
      res[t-1,'dQ'] <- dQ
      res[t-1,'act'] <- act
      res[t-1,'PEL'] <- PEL
      res[t-1,'PEH'] <- PEH
      res[t-1,'lr']  <- lr
      
    }
  
    return(res); 
}

} # end block of auxiliary if statement to skip block.

###  ########  putting together the key data object, abFit ######### ###
makeAllD <- 0  # auxiliary to skip block, which has a lot of 
# by-hand work. 
if (makeAllD != 0) {
  
  corrA09j <- importCSV('corrFit_a09j.csv')
  corrA09k <- importCSV('corrFit_a09k.csv')
  corrA09l <- importCSV('corrFit_a09l.csv')
  corrA09m <- importCSV('corrFit_a09m.csv')
  corrA09o <- importCSV('corrFit_a09o.csv')
  corrA09p <- importCSV('corrFit_a09p.csv')
  corrA09q <- importCSV('corrFit_a09q.csv')
  corrA09r <- importCSV('corrFit_a09r.csv')
  corrA09s <- importCSV('corrFit_a09s.csv')
  corrA09u <- importCSV('corrFit_a09u.csv')
  corrA09v <- importCSV('corrFit_a09v.csv')
  corrA10x <- importCSV('corrFit_a10x.csv')
  
  corrB09j <- importCSV('corrFit_b09j.csv')
  corrB09k <- importCSV('corrFit_b09k.csv')
  corrB09l <- importCSV('corrFit_b09l.csv')
  corrB09m <- importCSV('corrFit_b09m.csv')
  corrB09o <- importCSV('corrFit_b09o.csv')
  corrB09p <- importCSV('corrFit_b09p.csv')
  corrB09q <- importCSV('corrFit_b09q.csv')
  corrB09r <- importCSV('corrFit_b09r.csv')
  corrB09s <- importCSV('corrFit_b09s.csv')
  corrB09u <- importCSV('corrFit_b09u.csv')
  corrB09v <- importCSV('corrFit_b09v.csv')
  
  print(dimnames(abFit)[[3]]); abFit[,,'jA'] <- as.matrix(corrA09j)
  
  abFit[,,'oA'] <- as.matrix(corrA09o)
  abFit[,,'pA'] <- as.matrix(corrA09p)
  abFit[,,'qA'] <- as.matrix(corrA09q)
  abFit[,,'rA'] <- as.matrix(corrA09r)
  abFit[,,'sA'] <- as.matrix(corrA09s)
  abFit[,,'uA'] <- as.matrix(corrA09u)
  abFit[,,'vA'] <- as.matrix(corrA09v)
  abFit[,,'xA'] <- as.matrix(corrA10x)
  abFit[,,'jB'] <- as.matrix(corrB09j)
  abFit[,1:29,'kB'] <- as.matrix(corrB09k);
  abFit[,,'lB'] <- as.matrix(corrB09l)
  abFit[,,'mB'] <- as.matrix(corrB09m)
  abFit[,,'oB'] <- as.matrix(corrB09o)
  abFit[,1:29,'qB'] <- as.matrix(corrB09q);
  abFit[,,'rB'] <- as.matrix(corrB09r)
  abFit[,,'sB'] <- as.matrix(corrB09s)
  abFit[,,'vB'] <- as.matrix(corrB09v)
  abFit[,,'uB'] <- as.matrix(corrB09u)
  
  # Add columns with drug info etc. to abFit:
  test <- abFit
  test <- abind(abFit,array(NA,replace(dim(abFit),2,3)),along=2)
  colnames(test)[30:32] <- colnames(abIDmapDrug)[2:4]
  for (k in 1:24){ test[,30:32,k] <- as.matrix(abIDmapDrug[,2:4]) }
  for (k in 1:24){ test[,1,k] <- 1:74 }
  names(dimnames(test)) <- list('MMptID','var','model')
  abFit <- test
  allD <- list()
  allD[[1]] <- abFit   # To store cleaned-by-hand variables from
  allD[[2]] <- abFit
  names(allD) <- c('clean','abWdrug')
  
  ### baseline - follow-up correlations and data hygiene #######
  library(ggplot2)
  library(ppcor)
  print(dimnames(abFit)[[2]]); 
  # $MMptID
  #  "pt1"  "pt2" ... "pt74"
  # $var
  #  "MMptN"         "wave"          "LL"            "othEthn"       "othAge"          
  #  "LP"            "F"             "pH0"           "pS0"           
  #  "dInitEv"       "aInitEv"       "initEvRat"    
  #  "alphaPrec"     "mem"           "wH"            
  #   "wS"            "w0"            "othLR"        
  #  "pocB"          "dInitS"        "retAv"         "pred1"         "HI1"           "SI1"          
  #  "predAv"        "HIAv"          "SIAv"          "AIC"           "BIC"           "abPtNmap"     
  #  "placebo0SSRI1" "AnaisMegPtN"  
  # $model
  #  "jA" "kA" "lA" "mA" "oA" "pA" "qA" "rA" "sA" "uA" "vA" "xA" "jB" "kB" "lB" "mB" "oB" "pB" "qB" "rB" "sB" "uB"
  #  "vB" "xB"
  
  v <- 'othLR'; d1='xA'; d2='xB'; 
  
  xy<-allD$abWdrug[,v,c(d1,d2)]; 
  x <- xy[,1]; y<- xy[,2]; 
  
  hist(x,15,main=d1,xlab=v,col='gray90'); hist(y,15,main=d2,xlab=v,col='gray70'); 
  plot(x,y,main=paste(v,';',d1,'vs.',d2),xlab=paste(d1,v),ylab=paste(d2,v),pch=21,col='white',bg='navy'); abline(0,1); 
  olsr <- lm(y ~x); # plot(y,rstandard(olsr),ylab='std resid.',xlab=paste(d2,v))
  robr <- rlm(y~x); summary(robr);
  
  print(summary(olsr));
  plot(olsr);
  
  ### find and exclude outliers -------------------------------------------
  # Cleaned model q by 8 Feb 23 for:  pH0, pS0, dInitEv, aInitEv, 
  #  (initEvRat is fixed in q), alphaPrec, mem, wH, wS, w0, othLR.
  allD$clean[,v,c(d1,d2)] <- xy; # we may want to exclude no outliers by hand.
  # xy['pt31',1] <- NA;  # aInitEv outlier excluded by hand.
  # xy[ vecTRUE(xy[,2] < -9) , 2 ] <- NA; # exclude pH0 outlier by hand
  # Exclude outliers based on histograms and quality control plots above :
  #badRow <- (rownames(xy) %in% c('pt19','pt38','pt67')) * 1:dim(xy)[1]; 
  #badRow <- badRow[naRow > 0]; xy[badRow,]
  # badRow <- (rownames(xy) %in% c('pt67')) * 1:dim(xy)[1]; 
  # badRow <- c('pt41','pt46') # for w0
  # xy[badRow,]
  # xy[badRow,] <- NA
  # badRow <- vecTRUE(xy[,1] < -1) * 1:dim(xy)[1]; badRow <- badRow[badRow > 0.5]
  # xy[badRow,] <- NA
  allD$clean[,v,c(d1,d2)] <- xy; 
  
  # Back up : 
  save.image("/media/sf_mmpsy/Dropbox/BASOR/AcInSOR/AcInRepeatedDictator/analyses/.RData")
  save(file = 'allD.RData',allD);
  
} # End block putting allD$abWdrug 
# and allD$clean together.

### Posterior correlation matrices -------------------------------------------
bakPar <- par()   # settings for plotting
par(mar=c(4,5,4,2))  #  down, left, up, right (ffs ...)
load('allD.RData')

bslKeyTaskHd <- c("UserID","pH0.A", "pS0.A", "dInitEv.A", "aInitEv.A", "alphaPrec.A", 
                  "mem.A", "wH.A" , "wS.A" , "w0.A" ,     "othLR.A"  ,  "BIC.A",           # first 12 are ID and model measures
                  "pred1.A", "HI1.A",  "SI1.A",   "predAv.A",   "HIAv.A",  "SIAv.A");      # 13-18 are descriptives.
fuKeyTaskHd  <-  c("UserID","pH0.B", "pS0.B", "dInitEv.B", "aInitEv.B", "alphaPrec.B", 
                   "mem.B", "wH.B" , "wS.B" , "w0.B" ,     "othLR.B"  , "BIC.B",
                   "pred1.B","HI1.B",  "SI1.B",   "predAv.B",   "HIAv.B",  "SIAv.B" );

bslPa4cor <-  c("pH0.A", "pS0.A", "dInitEv.A", "aInitEv.A", "alphaPrec.A", 
                  "mem.A", "wH.A" , "wS.A" , "w0.A" ,     "othLR.A") 
fuPa4cor <-  c("pH0.B", "pS0.B", "dInitEv.B", "aInitEv.B", "alphaPrec.B", 
                "mem.B", "wH.B" , "wS.B" , "w0.B" ,     "othLR.B") 
bslPs4cor <- c("PHQ.A",  "GAD.A",  "MASQ.A",  "AMI.A",  "PSS.A",   "RSE.A", "HPS") 
fuPs4cor <- c("PHQ.B",  "GAD.B",  "MASQ.B",  "AMI.B",  "PSS.B",   "RSE.B","HPS")

dPs <- as.matrix(na.omit(allD$wiDat[['yw']][,c(bslPs4cor,"placebo0SSRI1.A")]));
dPs <- dPs[vecTRUE(dPs[,"placebo0SSRI1.A"]==1),1:(dim(dPs)[2]-1)];
psCor <- cor(dPs)
Bonf <-  (dim(psCor)[1] * (dim(psCor)[1]-1) / 2);
sigPsCor <- cor.mtest(dPs, conf.level=0.95);
corrplot(psCor,
         type='upper',       
         method = 'ellipse',
               p.mat=sigPsCor$p, insig='p-value',sig.level = 0.01,
               order='alphabet',addrect=2
        )

# if all together, uncomment:
dPa <- as.matrix(na.omit(allD$wiDat[['yw']][,c(fuPa4cor)]));
# if only Citalo group, uncomment:
# dPa <- as.matrix(na.omit(allD$wiDat[['yw']][,c(fuPa4cor,"placebo0SSRI1.A")]));
# dPa <- dPs[vecTRUE(dPa[,"placebo0SSRI1.A"]==1),1:(dim(dPa)[2]-1)];
paCor <- cor(dPa)
paBonf <-  (dim(paCor)[1] * (dim(paCor)[1]-1) / 2);
sigLev = 0.05 / paBonf; 
sigPaCor <- cor.mtest(dPa, conf.level=(1-sigLev));
corrplot(paCor,
         type='upper',       
         method = 'ellipse',
         p.mat=sigPaCor$p, insig='p-value',sig.level = sigLev,
         # order='AOE',addrect=2
         order='alphabet'
)

###  Mixed effects analyses ==========================================
y <- allD$clean[,,'qA']
l1 <- lme(SIAv ~ othEthn+retAv, random = ~1|ptN, data=y); summary(l1)

###  Gender analyses ==========================================
# Convenience copy of data, to be deleted:
yw <- allD$wiDat$yw  # Rem this is most data in wide form
# Merge with version with cleanest Female/Male data and do quick check:
y <- merge(allD$psy,yw,by.y='keyPsyID',by.x='psyID', all=TRUE); dim(y)

v1A <- 'HIAv.A'; v1B <- 'HIAv.B'; vlab <- 'Average Attributions'
v2A <- 'SIAv.A'; v2B <- 'SIAv.B';
boxplot(y[vecTRUE(y[,"Gender.A"]=='Female'),v1A],
        y[vecTRUE(y[,"Gender.A"]=='Male'),v1A],
        y[vecTRUE(y[,"Gender.A"]=='Female'),v1B],
        y[vecTRUE(y[,"Gender.A"]=='Male'),v1B],
        y[vecTRUE(y[,"Gender.A"]=='Female'),v2A],
        y[vecTRUE(y[,"Gender.A"]=='Male'),v2A],
        y[vecTRUE(y[,"Gender.A"]=='Female'),v2B],
        y[vecTRUE(y[,"Gender.A"]=='Male'),v2B],
        main = vLab,
        col = c('green','green4','gold','gold4','cyan','cyan4','brown','brown4'),
        names =c('HI F bsl','HI M bsl',
                 'HI F f/u','HI M f/ul',
                 'SI F bsl','SI M bsl',
                 'sI F f/u','SI M f/u'), 
        notch=TRUE)
abline(v=2.5); abline(v=4.5); abline(v=6.5); 

v1A <- 'pH0.A'; v1B <- 'pH0.B'; vlab <- 'Prior-mean parameters'
v2A <- 'pS0.A'; v2B <- 'pS0.B';
boxplot(y[vecTRUE(y[,"Gender.A"]=='Female'),v1A],
        y[vecTRUE(y[,"Gender.A"]=='Male'),v1A],
        y[vecTRUE(y[,"Gender.A"]=='Female'),v1B],
        y[vecTRUE(y[,"Gender.A"]=='Male'),v1B],
        y[vecTRUE(y[,"Gender.A"]=='Female'),v2A],
        y[vecTRUE(y[,"Gender.A"]=='Male'),v2A],
        y[vecTRUE(y[,"Gender.A"]=='Female'),v2B],
        y[vecTRUE(y[,"Gender.A"]=='Male'),v2B],
        main = vLab,
        col = c('green','green4','gold','gold4','cyan','cyan4','brown','brown4'),
        names =c('pH0 F bsl','pH0 M bsl',
                 'pH0 F f/u','pH0 M f/ul',
                 'pS0 F bsl','pS0 M bsl',
                 'pS0 F f/u','pS0 M f/u'), 
        notch=TRUE)
abline(v=2.5); abline(v=4.5); abline(v=6.5); 




remove(y)  # delete convenience copy
# tidy up =============================================================
par(backPar)
###                               eof                                ###
### ---------------------------------------------------------------- ###


