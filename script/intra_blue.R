#load packages
library(stringr)
library(plyr)
library(dplyr)
library(glmmTMB)
library(lme4)
library(ggplot2)

#load data
options(stringsAsFactors = FALSE)
berry_single <- read.csv("data/BLUEBERRY_MIXED_V3.csv", header=T)
berry_single$Year=as.factor(berry_single$Year)

#SUBSET TO JUST EVERGREEN
berry_single=berry_single[berry_single$Species%in%"EG",]%>%droplevels()

#subset dataframe to just visitation treatment
#convert all visitor columns to character
berry_single_updated <- berry_single[berry_single$TREATMENT%in%"V" & berry_single$TVN > 1,] %>%
                        mutate_at(vars(19:33), as.character) %>% 
                        droplevels()

berry_single_updated <- berry_single_updated[!berry_single_updated$SPEC.COM%in%"MX",] %>% 
                        droplevels()


#retain only first letter of the string for visitors IDs
berry_single_updated2 <- berry_single_updated %>% mutate_at(vars(19:33), function (x) word(x,2))


#sum the number of visits from each taxa
berry_single_updated2$N <- rowSums(berry_single_updated2[, c(19:33)] == "N",na.rm=T)
berry_single_updated2$P <- rowSums(berry_single_updated2[, c(19:33)] == "P",na.rm=T)

#sum the total number of visits
berry_single_updated2$sumvisits <- rowSums(berry_single_updated2[, c(37:38)])

#calculate the percent visits from each taxa
berry_single_updated2$p_N <- berry_single_updated2$N/berry_single_updated2$sumvisits
berry_single_updated2$p_P <- berry_single_updated2$P/berry_single_updated2$sumvisits

##Weight dataframe
berry_single_updated3=berry_single_updated2[!is.na(berry_single_updated2$Fresh.wgt)==TRUE,]
berry_single_updated3=berry_single_updated3[berry_single_updated3$sumvisits <6,]
table(berry_single_updated3$sumvisits)

#priority effects models
m1 <- glmmTMB(FS~SPEC.COM*Visitor1+sumvisits+(1|Year),
              family="binomial",
              data = berry_single_updated2)
summary(m1)

m2 <- glm(FS~SPEC.COM+Visitor1+sumvisits,
              family="binomial",
              data = berry_single_updated2)
summary(m2)
m2res=simulateResiduals(m2)
plot(m2res)


m2 <- glmmTMB(Fresh.wgt~SPEC.COM*Visitor1+sumvisits+(1|Year),
              family="gaussian",
              data = berry_single_updated3)
summary(m2)




m1 <- glmmTMB(Fresh.wgt~Visitor1*TVN+(1|Year),
              family="gaussian",
              data = berry_single_updated3)
summary(m1)










