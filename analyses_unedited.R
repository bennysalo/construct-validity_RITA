rm(list = ls())

# Read data from SPSS-file
library(foreign, pos=4)
FinPrisDataM <- read.spss("C:/Users/benny_000/Dropbox/BennypåSarasToshiba/RITA-CFA/males_all_random01_1496.sav", 
                    use.value.labels=TRUE, max.value.labels=Inf, to.data.frame=TRUE)


names(FinPrisDataM)

#subset of data interesting right now
subsetData<-FinPrisDataM[c(1,7,8,32:83,129,132:134)]
names(subsetData)

# set NA on martial and parenting as 0 = no problem
subsetData$marital<-ifelse(is.na(subsetData$marital),0,subsetData$marital) 
subsetData$parenting<-ifelse(is.na(subsetData$parenting),0,subsetData$parenting)
summary(subsetData)

# Create collapsed variables, the need will be proven later
Housing12 = (FinPrisDataM$Housing1 + FinPrisDataM$Housing2)

Drug2345678 = FinPrisDataM$Drug2 + FinPrisDataM$Drug3 + FinPrisDataM$Drug4 + FinPrisDataM$Drug5 + 
  FinPrisDataM$Drug6 + FinPrisDataM$Drug7 + FinPrisDataM$Drug8

# Recode into scale 0,1,2 with the 0 as 0, up to a mean of 1 = 1, above = 2

Housing12 = ifelse(Housing12 == 2,1,Housing12)
Housing12 = ifelse(Housing12 >= 3,2,Housing12)
table(Housing12)

Drug2345678 = ifelse(Drug2345678 == 2,1,Drug2345678)
Drug2345678 = ifelse(Drug2345678 == 3,1,Drug2345678)
Drug2345678 = ifelse(Drug2345678 == 4,1,Drug2345678)
Drug2345678 = ifelse(Drug2345678 == 5,1,Drug2345678)
Drug2345678 = ifelse(Drug2345678 == 6,1,Drug2345678)
Drug2345678 = ifelse(Drug2345678 == 7,1,Drug2345678)
Drug2345678 = ifelse(Drug2345678 >= 8,2,Drug2345678)
table(Drug2345678)

#add to subset
subsetData<-data.frame(subsetData, Housing12, Drug2345678)

#create a subsetversion that remains numeric
subsetDataNum<-subsetData
#define all variables as ordered
names(subsetData)
subsetData[c(1:61)] <- lapply(subsetData[c(1:61)], ordered)
summary(subsetData)

# 'cvgroup' is a randomized group variable - 748 in each group
# create subgroups with FRNAF-variables
summary(subsetData$cvgroup)
FRNAF0<-subsetData[subsetData$cvgroup == "original",4:55]
FRNAF1<-subsetData[subsetData$cvgroup == "crossvalidated",4:55]
FRNAF<-subsetData[4:55]
summary(FRNAF0)
summary(FRNAF1)

#load neccessary packages
library(lavaan)
library(psych)


#use the polycoric correlation function in lavaan to calculate polychoric correlations
RNcor0<-lavCor(FRNAF0)

#Check square multiple correlation
smc(RNcor0)

#extreme multicollinerarity found
#find correlations over 0.8 (all other set to 0)
RNcor0_0.8<-ifelse(RNcor0<0.8,0,RNcor0)
RNcor0_0.8


#find same correlations in other sample
RNcor1<-lavCor(FRNAF1)
RNcor1_0.8<-ifelse(RNcor1<0.8,0,RNcor1)
RNcor1_0.8

#subset only drug variables
names(FRNAF0)
DRUG0<-FRNAF0[33:41]
DRUG1<-FRNAF1[33:41]

# polychoric correlations in drug variables
Dcor0<-lavCor(DRUG0)
Dcor1<-lavCor(DRUG1)

Dcor0
Dcor1

#subset Drug variables 2 to 8
D280<-DRUG0[2:8]
D281<-DRUG1[2:8]

D28cor0<-lavCor(D280)
D28cor1<-lavCor(D281)

alpha(D28cor0)
alpha(D28cor1)

# Replace relevant variables with Housing12 and Drug 2345678 
names(subsetData)
FRNAF0<-subsetData[subsetData$cvgroup == "original",c(60,7:36,61,44:55)]
FRNAF1<-subsetData[subsetData$cvgroup == "crossvalidated",c(60,7:36,61,44:55)]

#check smc
RNcor0<-lavCor(FRNAF0)
smc(RNcor0)

RNcor1<-lavCor(FRNAF1)
smc(RNcor1)

#Exlude Drug9 (44)
names(subsetData)
FRNAF0<-subsetData[subsetData$cvgroup == "original",c(60,7:36,61,45:55)]
FRNAF1<-subsetData[subsetData$cvgroup == "crossvalidated",c(60,7:36,61,45:55)]
summary(FRNAF0)
#check smc
RNcor0<-lavCor(FRNAF0)
sort(smc(RNcor0))

RNcor1<-lavCor(FRNAF1)
sort(smc(RNcor1))

#extreme multicollinearity still in crossvalidated group
#Exlude ALKO3 (30)
names(subsetData)
FRNAF0<-subsetData[subsetData$cvgroup == "original",c(60,7:29,31:36,61,45:55)]
FRNAF1<-subsetData[subsetData$cvgroup == "crossvalidated",c(60,7:29,31:36,61,45:55)]
FRNAFA<-subsetData[c(60,7:29,31:36,61,45:55)]
#check smc and KMO
RNcor0<-lavCor(FRNAF0)
sort(smc(RNcor0)) #OK
sort(KMO(RNcor0)$MSAi) #items with low MSA

RNcor1<-lavCor(FRNAF1)
sort(smc(RNcor1)) #some items with smc over .9
sort(KMO(RNcor1)$MSAi)#items with low MSA

RNcorA<-lavCor(FRNAFA)
sort(smc(RNcorA)) #OK
sort(KMO(RNcorA)$MSAi) #Peruskoulu has the lowest MSA (.507)


#Exclude Peruskoulu on basis of low KMO (11) in second sample
names(subsetData)
#define the set of variables to include for quicker coding later
set<-c(60,7:10,12:29,31:36,61,45:55)
FRNAF0<-subsetData[subsetData$cvgroup == "original",set]
FRNAF1<-subsetData[subsetData$cvgroup == "crossvalidated",set]
FRNAFA<-subsetData[set]
#check smc and KMO
RNcor0<-lavCor(FRNAF0)
sort(smc(RNcor0)) #OK
sort(KMO(RNcor0)$MSAi) #items with low MSA, marital .419, domvioPerp, domvioVict .472

RNcor1<-lavCor(FRNAF1)
sort(smc(RNcor1)) #OK
sort(KMO(RNcor1)$MSAi) #OK, marital has MSA .634, domvioVict .654 

RNcorA<-lavCor(FRNAFA)
sort(smc(RNcorA)) #OK
sort(KMO(RNcorA)$MSAi) #OK but marital lowest .634, domvioVict .654 

#Exclude marital on basis of low KMO (20)
names(subsetData)
set<-c(60,7:10,12:19,21:29,31:36,61,45:55)
FRNAF0<-subsetData[subsetData$cvgroup == "original",set]
FRNAF1<-subsetData[subsetData$cvgroup == "crossvalidated",set]
FRNAFA<-subsetData[set]
#check smc and KMO
RNcor0<-lavCor(FRNAF0)
sort(smc(RNcor0)) #OK
sort(KMO(RNcor0)$MSAi) #low MSA parenting .395, domvioPerp, domvioVict .511

RNcor1<-lavCor(FRNAF1)
sort(smc(RNcor1)) #OK
sort(KMO(RNcor1)$MSAi)  #low MSA domvioVict .427

RNcorA<-lavCor(FRNAFA)
sort(smc(RNcorA)) #OK
sort(KMO(RNcorA)$MSAi)

#Exclude domvioVict on basis of low KMO (22)
names(subsetData)
set<-c(60,7:10,12:19,21,23:29,31:36,61,45:55)
FRNAF0<-subsetData[subsetData$cvgroup == "original",set]
FRNAF1<-subsetData[subsetData$cvgroup == "crossvalidated",set]
FRNAFA<-subsetData[set]
#check smc and KMO
RNcor0<-lavCor(FRNAF0)
sort(smc(RNcor0))
sort(KMO(RNcor0)$MSAi)

RNcor1<-lavCor(FRNAF1)
sort(smc(RNcor1))
sort(KMO(RNcor1)$MSAi)

RNcorA<-lavCor(FRNAFA)
sort(smc(RNcorA))
sort(KMO(RNcorA)$MSAi)

#Exclude parenting on basis of low KMO (23)
names(subsetData)
set<-c(60,7:10,12:19,21,24:29,31:36,61,45:55)
FRNAF0<-subsetData[subsetData$cvgroup == "original",set]
FRNAF1<-subsetData[subsetData$cvgroup == "crossvalidated",set]
FRNAFA<-subsetData[set]
#check smc and KMO
RNcor0<-lavCor(FRNAF0)
sort(smc(RNcor0))
sort(KMO(RNcor0)$MSAi)

RNcor1<-lavCor(FRNAF1)
sort(smc(RNcor1))
sort(KMO(RNcor1)$MSAi)

RNcorA<-lavCor(FRNAFA)
sort(smc(RNcorA))
sort(KMO(RNcorA)$MSAi)

#Exclude LUKIproblem on basis of low KMO (12)
names(subsetData)
set<-c(60,7:10,13:19,21,24:29,31:36,61,45:55)
FRNAF0<-subsetData[subsetData$cvgroup == "original",set]
FRNAF1<-subsetData[subsetData$cvgroup == "crossvalidated",set]
FRNAFA<-subsetData[set]
#check smc and KMO
RNcor0<-lavCor(FRNAF0)
sort(smc(RNcor0))
sort(KMO(RNcor0)$MSAi)

RNcor1<-lavCor(FRNAF1)
sort(smc(RNcor1))
sort(KMO(RNcor1)$MSAi)

RNcorA<-lavCor(FRNAFA)
sort(smc(RNcorA))
sort(KMO(RNcorA)$MSAi)

#Exclude submissive on basis of low KMO (27)
names(subsetData)
set<-c(60,7:10,13:19,21,24:26,28,29,31:36,61,45:55)
FRNAF0<-subsetData[subsetData$cvgroup == "original",set]
FRNAF1<-subsetData[subsetData$cvgroup == "crossvalidated",set]
FRNAFA<-subsetData[set]
#check smc and KMO
RNcor0<-lavCor(FRNAF0)
sort(smc(RNcor0))
sort(KMO(RNcor0)$MSAi)

RNcor1<-lavCor(FRNAF1)
sort(smc(RNcor1))
sort(KMO(RNcor1)$MSAi)

RNcorA<-lavCor(FRNAFA)
sort(smc(RNcorA))
sort(KMO(RNcorA)$MSAi)

# This subset fills criteria. Smc under .9 and KMO over .6 in both samples

# create subset in numeric version
set<-c(60,7:10,13:19,21,24:26,28,29,31:36,61,45:55)
FRNAF0Num<-subsetDataNum[subsetData$cvgroup == "original",set]
FRNAF1Num<-subsetDataNum[subsetData$cvgroup == "crossvalidated",set]
FRNAFANum<-subsetDataNum[set]

#Change a 2 to a 1 in Tarkkis to get around fa.parallel.poly bug
FRNAF0Num$Tarkkis[6]<-1
FRNAF1Num$Tarkkis[6]<-1
#FRNAFANum$Tarkkis[6]<-1

# Do parallel analysis
parAn<-fa.parallel.poly(FRNAF0Num, fm = "ml", fa = "fa", n.iter = 100)
parAn$fa.sim$mean+qnorm(0.95)*parAn$fa.sim$sd
parAn$fa.values

# compare to a permutation test
parres<-matrix(nrow = 100, ncol=36)
for (i in 1:100)
{
  parAn<-fa.parallel.poly(FRNAF0Num, fm = "ml", fa = "fa", n.iter = 1, global = FALSE)
  parres[i,]<-parAn$fa.sim$mean
}
parAn$fa.values[6]
sort(parres[,6], decreasing = T)
parAn$fa.values[7]
sort(parres[,7], decreasing = T)
parAn$fa.values[8]
sort(parres[,8], decreasing = T)

?sort
#check parallel analysis in other sample
parAn1<-fa.parallel.poly(FRNAF1Num, fm = "ML", fa = "fa", n.iter = 100, sim =T, global = FALSE)

parres1<-matrix(nrow = 100, ncol=36)
for (i in 1:100)
{
  parAn<-fa.parallel.poly(FRNAF1Num, fm = "ml", fa = "fa", n.iter = 1, global = FALSE)
  parres[i,]<-parAn$fa.sim$mean
}
parAn1$fa.values[6]
sort(parres[,6], decreasing = T)
parAn1$fa.values[7]
sort(parres[,7], decreasing = T)
parAn1$fa.values[8]
sort(parres[,8], decreasing = T)

#parAnA<-fa.parallel.poly(FRNAF0NumA, fm = "ML", fa = "fa", n.iter = 100, sim =T, global = FALSE)

# Do factor analysis in CFA framework in lavaan

library(semTools)
EFA06<-efaUnrotate(FRNAF0, nf = 6)
EFA06RotQ <- oblqRotate(EFA06, method="quartimin")
summary(EFA06RotQ, suppress = 0)
summary(EFA06, fit = T)

EFA05<-efaUnrotate(FRNAF0, nf = 5)
EFA05RotQ <- oblqRotate(EFA05, method="quartimin")


#Exclude familyties on basis of no loadings over .3
names(subsetData)
set<-c(60,7:10,13:18,21,24:26,28,29,31:36,61,45:55)
FRNAF0m<-subsetData[subsetData$cvgroup == "original",set]

EFA06m<-efaUnrotate(FRNAF0m, nf = 6)
EFA06RotQm <- oblqRotate(EFA06m, method="quartimin")
fitmeasures(EFA06m)

EFA06m<-efaUnrotate(FRNAF0m, nf = 6)

EFA08RotQ <- oblqRotate(EFA08, method="quartimin")
fitmeasures(EFA07m)

#check if any variables can be added
names(subsetData)
set<-c(60,7:10,13:18,21,24:26,28,29,31:36,61,44,45:55)
FRNAF0<-subsetData[subsetData$cvgroup == "original",set]
FRNAF1<-subsetData[subsetData$cvgroup == "crossvalidated",set]
FRNAFA<-subsetData[set]
#check smc and KMO
RNcor0<-lavCor(FRNAF0)
sort(smc(RNcor0))
sort(KMO(RNcor0)$MSAi)

RNcor1<-lavCor(FRNAF1)
sort(smc(RNcor1))
sort(KMO(RNcor1)$MSAi)



# Define CFA model based on EFA06rot
ModEFA06<-'
f1 =~ Finance2 + Housing3 + Finance3 + Housing 12 + Finance1 + jobseek
f2 =~ ALKO5 + ALKO7 + ALKO1 + ALKO6 + ALKO8 + ALKO2 + ALKO4
f3 =~ ALKO8 + motivation + insight + procrime + others + attitudestaff + manipulative + coping + socialskills + attitudesuperv + hostile + instrumental + jobattitude
f4 =~ procrime + Drug2345678 + Drug1 + criminalpeers + ALKO2 + riskseeking
f5 =~ others + attitudestaff + impulsive + ALKO4 + hostile + instrumental + domvioPerp
f6 =~ skillsneeds + employment + jobattitude + eduattitude + Tarkkis + jobseek'


CFA06pars<-cfa(ModEFA06, data = FRNAF0, std.lv = T, estimator = "WLSMV")
summary(CFA06pars, fit = T)
mi<-modindices(CFA06pars)
head(mi[order(modindices(CFA06pars)$mi.scaled, decreasing = T),],50)

#exclude coping
ModEFA06minC <-' 
f1 =~ Finance2 + Housing3 + Finance3 + Housing 12 + Finance1 + jobseek
f2 =~ ALKO5 + ALKO7 + ALKO1 + ALKO6 + ALKO8 + ALKO2 + ALKO4
f3 =~ ALKO8 + motivation + insight + procrime + others + attitudestaff + manipulative +  socialskills + attitudesuperv + hostile + instrumental + jobattitude
f4 =~ procrime + Drug2345678 + Drug1 + criminalpeers + ALKO2 + riskseeking
f5 =~ others + attitudestaff + impulsive + ALKO4 + hostile + instrumental + domvioPerp
f6 =~ skillsneeds + employment + jobattitude + eduattitude + Tarkkis + jobseek'

CFA06parsMinC<-cfa(ModEFA06minC, data = FRNAF0, std.lv = T, estimator = "WLSMV")
summary(CFA06parsMinC, fit = T)
mi<-modindices(CFA06parsMinC)
head(mi[order(modindices(CFA06parsMinC)$mi.scaled, decreasing = T),],50)


#Test in crossvalidating sample - this is Table 2 and 3 
CFA16parsMinC<-cfa(ModEFA06minC, data = FRNAF1, std.lv = T, estimator = "WLSMV")
summary(CFA16parsMinC, fit = T, standardized = T)

# Check reliabily coefficient omega. omega3 is reported.
require(semTools)
reliability(CFA16parsMinC)


