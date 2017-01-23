## data preparation for analysis.
## run in this directory
## ------------libraries--------
library(data.table)
library(reshape2)

## -------functions--------
oddit <- function(x) x/(1-x)
ioddit <- function(x) x/(1+x)
logit <- function(x) log(oddit(x))
ilogit <- function(x) ioddit(exp(x))

## --------work----------
B <- read.csv("country_incidence_disaggregated_age_sex_num.csv") #WHO 2016 incidence

B$mid <- B$inc.num.014
B$var <- ((B$inc.num.014.hi-B$inc.num.014.lo)/(2*1.96))^2


## sample from the estimates assuming gamma (mean=k scale, var = k scale^2)
B$S <- B$var/B$mid
B$k <- B$mid/B$S

## inspect if reasonable representation of burden
who <- 200+17                                #choose test country
png(paste0('../test/test_TB_',B$country[who],'.png'))
hist(rgamma(1e4,shape=B$k[who],scale=B$S[who]),main=B$country[who])
abline(v=B$mid[who],col=2);
abline(v=B$inc.num.014.hi[who],lty=2,col=2);
abline(v=B$inc.num.014.lo[who],lty=2,col=2);
dev.off()


## merge with notification data
N <- read.csv("TB_notifications_2016-10-18.csv") #current WHO notifications
nmz <- grep('14',names(N),value=TRUE)
nmz <- c(nmz,grep('04',names(N),value=TRUE))
N <- N[,c('country','iso3','year','g_whoregion',nmz)] #restsrict to relevant
N <- as.data.table(N)                   #convert to data.table
N <- N[year==2015,]                     #restrict current year

## all kids
N[,notifs:=newrel_m014 + newrel_f014 ]   #0-14 notifications, new system
addon <- N$newrel_sexunk014
addon[is.na(addon)] <- 0                #most are NA
N$notifs <- N$notifs + addon            #add on those that aren't

## young
N[,notifsY:=newrel_m04 + newrel_f04]
addon <- N$newrel_sexunk04
addon[is.na(addon)] <- 0                #most are NA
N$notifsY <- N$notifsY + addon            #add on those that aren't

## Old
N[,notifsO:=newrel_m514 + newrel_f514]
addon <- N$newrel_sexunk514
addon[is.na(addon)] <- 0                #most are NA
N$notifsO <- N$notifsO + addon            #add on those that aren't

## inspect
N[,list(notifs,notifsY,notifsO)]
N[,sum(is.na(notifs))]
N[,sum(is.na(notifsY))]
N[,sum(is.na(notifsO))]
keepBB <- c('iso3','year',
            "inc.num","inc.num.lo","inc.num.hi",
            "inc.num.014","inc.num.014.lo","inc.num.014.hi",
            "e.pop.num","e.pop.014",
            "mid","var","S",'k')
BB <- merge(B[,keepBB],N[,list(notifs,notifsY,notifsO,g_whoregion,iso3)],by='iso3') #merge with notifications
BB$id <- 1:nrow(BB)                     #for later


## ---- HIV data
load('HIVcleaning/paedHIV2016.Rdata')
H <- H1615c[,list(HIVmid=mid,HIVlo=lo,HIVhi=hi,iso3)]
BBA <- merge(BB,H,by='iso3',all.x=TRUE)
BB <- BBA

## ---- ART data
load('HIVcleaning/paedART2016.Rdata')
names(A)[2] <- 'ART'
BBA <- merge(BB,A,by='iso3',all.x=TRUE)
BB <- BBA

## tidying and saving
BB <- BB[,c('id','iso3','g_whoregion','e.pop.num','e.pop.014',
            'notifs','notifsY','notifsO',
            'inc.num','inc.num.lo','inc.num.hi',
            'inc.num.014','inc.num.014.lo','inc.num.014.hi',
            'S','k',
            'ART','HIVmid','HIVlo','HIVhi')]

## convert ART/HIV NAs to 0 (countries with no HIV data ~ no HIV)
BB$ART[is.na(BB$ART)] <- 0
BB$HIVmid[is.na(BB$HIVmid)] <- 0
BB$HIVlo[is.na(BB$HIVlo)] <- 0
BB$HIVhi[is.na(BB$HIVhi)] <- 0
## proportion ART/HIV
## make ART uncertainty derive from HIV denominator uncertainty!! NB
BB$artmid <- BB$ART/100
ses <- (BB$HIVhi-BB$HIVlo) / (BB$HIVmid+1e-6) 
BB$artlo <- ilogit( logit(BB$artmid) * (1+ses/2))
BB$arthi <- ilogit( logit(BB$artmid) * (1-ses/2))


## HIV
BB$hivmid <- BB$HIVmid/BB$e.pop.014
BB$hivlo <- BB$HIVlo/BB$e.pop.014
BB$hivhi <- BB$HIVhi/BB$e.pop.014


## distributions
BB$hS <- ((BB$hivhi-BB$hivlo)/3.92)^2/BB$hivmid
BB$hk <- BB$hivmid/BB$hS
bad <- is.nan(BB$hk) | is.nan(BB$hS) 
BB$hk[bad] <- BB$hS[bad] <- NA
BB$aS <- ((BB$arthi-BB$artlo)/3.92)^2/BB$artmid
BB$ak <- BB$artmid/BB$aS
bad <- is.nan(BB$ak) | is.nan(BB$aS) 
BB$ak[bad] <- BB$aS[bad] <- NA

## better country names
load('whokey.Rdata')
BB <- merge(BB,WHOkey[,c(1,3)],by='iso3')

## inspect if reasonable representation of burden
who <- 200+1*17                                #choose test country
png(paste0('../test/test_HIV_',BB$country[who],'.png'))
hist(rgamma(1e4,shape=BB$hk[who],scale=BB$hS[who]),main=BB$country[who])
abline(v=BB$hivmid[who],col=2);abline(v=(BB$hivlo)[who],lty=2);abline(v=(BB$hivhi)[who],lty=2);
dev.off()

png(paste0('../test/test_ART_',BB$country[who],'.png'))
hist(rgamma(1e4,shape=BB$ak[who],scale=BB$aS[who]),main=BB$country[who])
abline(v=BB$artmid[who],col=2);abline(v=(BB$artlo)[who],lty=2);abline(v=(BB$arthi)[who],lty=2);
dev.off()

## age splits from Dodd/Seddon model
load('AF.Rdata')
AF <- as.data.frame(AF)
AF <- AF[,c(1,4,5)]
names(AF)[2:3] <- c('aa','ab')
BA <- merge(BB,AF,by='iso3',all.x = TRUE)
BA$aa[is.na(BA$aa)] <- mean(BA$aa,na.rm=TRUE)
BA$ab[is.na(BA$ab)] <- mean(BA$ab,na.rm=TRUE)
BB <- BA


## save/load here!
save(BB,file='BB.Rdata')
