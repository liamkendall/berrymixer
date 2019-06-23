##BUMBLE

#load packages
library(stringr)
library(plyr)
library(dplyr)
library(glmmTMB)
library(lme4)
library(ggplot2)

#load data
bumble <- read.csv("data/BLUEBERRY_MIXED_V5.csv", header=T)
bumble <- bumble[bumble$Species%in%"BR"&bumble$TREATMENT%in%"V",] %>%droplevels()
bumble$duration=bumble_time_time$duration
bumble$V1duration=bumble_time_time$V1duration

#CLEAN DATA
table(bumble$SPEC.COM)
bumble[bumble$SPEC.COM%in%"MX ",c("SPEC.COM")]="MX"
bumble[bumble$SPEC.COM%in%"NX",c("SPEC.COM")]="MX"
bumble[bumble$SPEC.COM%in%"HB   ",c("SPEC.COM")]="HB"
bumble[bumble$SPEC.COM%in%"HB  ",c("SPEC.COM")]="HB"
bumble[bumble$SPEC.COM%in%"HB ",c("SPEC.COM")]="HB"



#subset dataframe to just visitation treatment
#convert all visitor columns to character
bumble_updated <- bumble[bumble$TVN > 1,] %>% mutate_at(vars(19:33), as.character) %>% droplevels()
#retain only first letter of the string for visitors IDs
bumble_updated2 <- bumble_updated %>% mutate_at(vars(19:33), word)

#sum the number of visits from each taxa
bumble_updated2$honey_bee <- rowSums(bumble_updated2[, c(19:33)] == "H")
bumble_updated2$bumble_bee <- rowSums(bumble_updated2[, c(19:33)] == "B")

#sum the total number of visits
bumble_updated2$sumvisits <- rowSums(bumble_updated2[, c(40:41)])

#calculate the percent visits from each taxa
bumble_updated2$p_honey_bee <- bumble_updated2$honey_bee/bumble_updated2$sumvisits
bumble_updated2$p_bumble_bee <- bumble_updated2$bumble_bee/bumble_updated2$sumvisits

bumble_updated2$Visitor1=as.factor(bumble_updated2$Visitor1)



bumble_updated3=bumble_updated2[bumble_updated2$SPEC.COM%in%"MX",]

bumble_updated2$RP=paste0(bumble_updated2$Row,bumble_updated2$Plant.number)

bumble_updated3$RP=paste0(bumble_updated3$Row,bumble_updated3$Plant.number)


#priority effects models - fruit set
m1 <- glm(FS~Visitor1*sumvisits,
              family="binomial",
              data = bumble_updated3)
m2 <- glmmTMB(FS~Visitor1+sumvisits+offset(log(duration)),
          family="binomial",
          data = bumble_updated3)
summary(m2)

m2 <- glm(FS~SPEC.COM*sumvisits,
          family="binomial",
          data = bumble_updated2)

summary(m2)
rsq::rsq.n(m2)

#priority effects models - fruit weight
m3 <- glmmTMB(Fresh.wgt~Visitor1*sumvisits+(1|RP),
          family="gaussian",
          data = bumble_updated3)

m4 <- glmmTMB(Fresh.wgt~Visitor1*sumvisits+(1|RP)+(1|Date),
              family="gaussian",
              data = bumble_updated3)
m5 <- glmmTMB(Fresh.wgt~Visitor1*sumvisits+(1|Date/RP),
              family="gaussian",
              data = bumble_updated3)

m6 <- glmmTMB(Fresh.wgt~Visitor1+sumvisits+(1|RP)+(1|Date),
              family="gaussian",
              data = bumble_updated3)
AIC(m3,m4,m5,m6)

#duration
m7 <- glmmTMB(Fresh.wgt~Visitor1*duration*sumvisits+(1|RP)+(1|Date),
              family="gaussian",
              data = bumble_updated3)
MuMIn::dredge(m7)
m7 <- glmmTMB(Fresh.wgt~Visitor1+sumvisits+(1|RP)+(1|Date),
              family="gaussian",
              data = bumble_updated3)
AIC(m3,m4,m5,m6)
summary(m6)

#Species composition

bumble_updated2$lfw=log(bumble_updated2$Fresh.wgt)

bumble_updated4=bumble_updated2[!bumble_updated2$SPEC.COM%in%"MX",]
m7 <- glmmTMB(lfw~SPEC.COM+sumvisits+(1|RP)+(1|Date),
          family="gaussian",
          data = bumble_updated4)
MuMIn::dredge(m7)

m8 <- glmmTMB(Fresh.wgt~SPEC.COM*sumvisits+(1|RP)+(1|Date),
              family="gaussian",
              data = bumble_updated2)
m9 <- glmmTMB(Fresh.wgt~SPEC.COM+sumvisits+(1|RP)+(1|Date),
              family="gaussian",
              data = bumble_updated2)
m10 <- glmmTMB(log(Fresh.wgt)~SPEC.COM*sumvisits+(1|RP)+(1|Date),
              family="gaussian",
              data = bumble_updated2)
m11 <- glmmTMB(log(Fresh.wgt)~SPEC.COM+sumvisits+(1|RP)+(1|Date),
               family="gaussian",
               data = bumble_updated2)
m12 <- glmmTMB(log(Fresh.wgt)~SPEC.COM+(1|RP)+(1|Date),
               family="gaussian",
               data = bumble_updated2)
AIC(m7,m8,m9,m10,m11,m12)            

summary(m10)
