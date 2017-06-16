## US HIV/ART data
hadf <- data.frame(epoch=rep(c(0,1),each=4),   #1 is late
                   HIV=rep(c(0,1),4),          #1 is HIV+
                   ART=c(rep(0,4),c(0,1,0,1)), #1 is ART+
                   outcome=rep(0:1,each=2),    #1 is death
                   count=c(671,66,8,11,3814,85,17,3)
                   )
hadfL <- hadf[rep(seq_len(nrow(hadf)), hadf$count),1:4] ## long form

## logistic regression
mod1 <- glm(outcome ~ epoch + HIV + ART ,family=binomial,data=hadfL)
V <- vcov(mod1)                         #extract variance-covariance matrix
ha <- c('HIV','ART')                    #factors of interest

## computed quantites for onward use
V[ha,ha]                                #variance-covariance for HIV/ART
coef(mod1)[ha]                          #means for HIV/ART (log ORs)
confint(mod1)[ha,]                      #95% CIs via profile likelihood
## output values:
## mean =  2.6375681 -0.5683867
## V = [[0.2325509 -0.2325509],[-0.2325509  0.6367345]]

## CIs for reporting in table 1
## NB ART effect doesn't include HIV - compute by simulation
Z <- MASS::mvrnorm(n=5e3,mu=coef(mod1)[ha],Sigma = V[ha,ha])
quantile(Z[,1],probs=c(.025,.975))#very close to above profile-likelihood version for HIV 
quantile(exp(Z[,1]),probs=c(.025,.5,.975)) #median, 95% for HIV OR
quantile(exp(rowSums(Z)),probs=c(.025,.5,.975)) #median, 95% for HIV + ART OR
## output values (subject to MC error):
## HIV =  13.92270 (5.37483  36.19513 )
## HIV+ART = 7.764393 (2.304652   26.997024)
