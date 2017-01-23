## ########
## paed HIV data cleaning
## ########
library(data.table)
## to be run in this directory

## ---------- new version 2/12/2016 ------
H <- read.csv('KHcleaner.csv',header=FALSE) #from running cleanHIV.sh
names(H) <- c('country',paste(c('mid','lo','hi'),rep(1990:2015,each=3),sep='_'))
H <- reshape2::melt(H,id='country')
tmp <- strsplit(as.character(H$variable),'_')
H$type <- unlist(lapply(tmp,function(x)x[1]))
H$year <- as.numeric(unlist(lapply(tmp,function(x)x[2])))
H <- reshape2::dcast(H, country + year ~ type, value=value)

##deal with <
lts1 <- grep('<',H$mid)
lts2 <- grep('<',H$lo)
lts3 <- grep('<',H$hi)
H$mid[lts1] <- substr(H$mid[lts1],start=2,stop=nchar(H$mid[lts1]))
H$lo[lts2] <- substr(H$lo[lts2],start=2,stop=nchar(H$lo[lts2]))
H$hi[lts3] <- substr(H$hi[lts3],start=2,stop=nchar(H$hi[lts3]))

## numeric
for(i in 3:5)
    H[,i] <- as.numeric(H[,i])

## his with < are same (no action)
## los with < are 0
H$lo[lts2] <- 0
## mids with < are halfway
H$mid[lts1] <- .5*(H$lo[lts1] + H$hi[lts1])
## sames are 50% ramge
sames <- which(abs(H$lo-H$hi) + abs(H$lo - H$mid)<.1)
H$hi[sames] <- H$mid[sames]*1.25
H$lo[sames] <- H$mid[sames]*0.75

## missing countries
load('paedHIV2015.Rdata') #2015 estimates up to 2014
load('../whokey.Rdata')

## join in iso3
H15 <- merge(paedHIV,WHOkey,by.x='Area',by.y = 'country',all.x = FALSE)
H15 <- reshape2::dcast(H15,iso3+year~quantity,value=value)
H16 <- merge(H,WHOkey[,c(1,3)],by.x='country',by.y='country')

head(H16)
head(H15)

H15 <- as.data.table(H15)
H16 <- as.data.table(H16)

## look at those not in 2016 but in 2015...
cn15 <- as.character(unique(H15[year==2014,iso3]))
cn16 <- as.character(unique(H16[year==2015,iso3]))
Cin15not16 <- cn15[!cn15 %in% cn16]

tot15 <- H15[year==2014]
tot15[is.na(mid),mid:=(hi+lo)/2]
tot15 <- tot15[,sum(mid)]

props <- H15[year==2014 & iso3 %in% Cin15not16,list(year,iso3,
                                                    mid = mid/tot15,
                                                    lo = lo/tot15,
                                                    hi = hi/tot15)]
tot16 <- H16[year==2015,sum(mid)]
H1615 <- H16[year==2015]
props <- merge(props,WHOkey[,c(1,3)],by='iso3')
props[,mid:=mid*tot16]
props[,lo:=lo*tot16]
props[,hi:=hi*tot16]

H1615c <- rbind(H1615,props)
H1615c[,sum(mid)]/1.8e6

save(H1615c,file='paedHIV2016.Rdata')



## ===================== ART =======================

A <- read.csv('KAcleaner.csv',header=FALSE)
for(i in 2:ncol(A))
    A[,i] <- as.numeric(gsub('[<>]','',as.character(A[,i])))
names(A)[1] <- 'country'
A <- merge(A,WHOkey[,c(1,3)],by='country')

## missing?
A$iso3[!as.character(A$iso3) %in% as.character(H1615c$iso3)] #no ART not in HIV
missing <- H1615c$iso3[!  as.character(H1615c$iso3) %in% as.character(A$iso3)] #from ART
A16 <- A

load('paedART2015.Rdata') #2015 data on 2014

paedART <- merge(paedART,WHOkey[,c(1,3)],by.x='Area.ID',by.y='iso3')
names(paedART)[1] <- 'iso3'


miscov <- data.frame(iso3=missing,cov=NA)
for(i in 1:nrow(miscov))
    if(as.character(miscov$iso3[i]) %in% paedART$iso3 & as.character(miscov$iso3[i]) %in% H15$iso3) miscov$cov[i] <- 1e2*paedART[iso3==as.character(miscov$iso3[i]) & year == 2014,value] / H15[iso3==as.character(miscov$iso3[i]) & year==2014,mid]

miscov <- na.omit(miscov)

A <- A[,c(20,19)]
names(A)[2] <- 'cov'
A <- rbind(A,miscov)

save(A,file='paedART2016.Rdata')


## ============== India ========================
## from AIDSinfo for 2015
## 50976 kids on ART
## 868165 adults on ART
## 919141 total on ART
## ART coverage adults 44 [36-55]
## ART coverage total = 43 [36-54]

tothiv <- 919141 / (1e-2*c(43,36,54)) 
adhiv <- 868165 /  (1e-2*c(44,36,55))
kidhiv <- tothiv*(1-mean(adhiv/tothiv))
kidart <- 50976

load(file='paedHIV2016.Rdata')
load(file='paedART2016.Rdata')

H1615c <- rbind(H1615c,
                data.frame(country='India',year=2015,
                           hi=kidhiv[2],lo=kidhiv[3],mid=kidhiv[1],iso3='IND'))

A <- rbind(A,data.frame(iso3='IND',cov=1e2*kidart/kidhiv[1]))


save(H1615c,file='paedHIV2016.Rdata')
save(A,file='paedART2016.Rdata')
