##BLUE SPECIES COMPOSITION MODELS

#load packages
library(stringr)
library(plyr)
library(dplyr)
library(glmmTMB)
library(lme4)
library(ggplot2)

#load data
berry <- read.csv("data/BLUEBERRY_MIXED_V3.csv", header=T)

#subset dataframe to just visitation treatment
#convert all visitor columns to character
berry_updated <- berry[berry$TREATMENT%in%"V" & berry$TVN > 1,] %>% mutate_at(vars(19:33), as.character) %>% droplevels()

#retain only first letter of the string for visitors IDs
berry_updated2 <- berry_updated %>% mutate_at(vars(19:33), word)

#sum the number of visits from each taxa
berry_updated2$stingless_bee <- rowSums(berry_updated2[, c(19:33)] == "S")
berry_updated2$honey_bee <- rowSums(berry_updated2[, c(19:33)] == "H")
berry_updated2$bumble_bee <- rowSums(berry_updated2[, c(19:33)] == "B")

#sum the total number of visits
berry_updated2$sumvisits <- rowSums(berry_updated2[, c(37:39)])

#calculate the percent visits from each taxa
berry_updated2$p_stingless_bee <- berry_updated2$stingless_bee/berry_updated2$sumvisits
berry_updated2$p_honey_bee <- berry_updated2$honey_bee/berry_updated2$sumvisits
berry_updated2$p_bumble_bee <- berry_updated2$bumble_bee/berry_updated2$sumvisits

#remove Tasmania BB data and RE data 
#due to very low fresh weight sample size (10 Mixed visits)
table(berry_updated2[berry_updated2$Species%in%"RE" &
                       is.na(berry_updated2$Fresh.wgt)==FALSE,]$SPEC.COM)

berry_updated3 <- berry_updated2[!berry_updated2$Species%in%"BR",]%>%droplevels()
berry_updated3 <- berry_updated3[!berry_updated3$Species%in%"RE",]%>%droplevels()

##Subset to 2 or more visits and only mixed visits
berry_updated_comp=berry_updated3[berry_updated3$sumvisits >1,]%>%droplevels()


##SPEC COM FRUIT SET MODELS

blue.sc.fs.m1=glmmTMB(FS~SPEC.COM*sumvisits+(1|Year/Block),
                       family=binomial,
                       data=berry_updated_comp)

blue.sc.m1res=simulateResiduals(blue.sc.fs.m1)
plot(blue.sc.m1res)
testResiduals(blue.sc.m1res)

summary(blue.sc.fs.m1)

#simple visit model
blue.sc.fs.m2=glmmTMB(FS~sumvisits+(1|Year/Block),
                       family=binomial,
                       data=berry_updated_comp)
blue.sc.m2res=simulateResiduals(blue.sc.fs.m2)
plot(blue.sc.m2res)
testResiduals(blue.sc.m2res)

#compare with AIC
AIC(blue.sc.fs.m1,blue.sc.fs.m2)

#           df   AIC
#blue.sc.fs.m1  8 169.922
#blue.sc.fs.m2  4 163.801

#Compare pollinators
emtrends(blue.sc.fs.m1,pairwise~SPEC.COM,var="sumvisits")


##SPEC COM WEIGHT MODELS
blue.sc.wgt.m1=glmmTMB(Fresh.wgt~SPEC.COM*sumvisits+(1|Year/Block),
                   family=gaussian,
                   data=berry_updated_comp)

blue.sc.m1res=simulateResiduals(blue.sc.wgt.m1)
plot(blue.sc.m1res)
testResiduals(blue.sc.m1res)


#simple visit model
blue.sc.wgt.m2=glmmTMB(Fresh.wgt~sumvisits+(1|Year/Block),
                   family=gaussian,
                   data=berry_updated_comp)

blue.sc.m2res=simulateResiduals(blue.sc.wgt.m2)
plot(blue.sc.m2res)
testResiduals(blue.sc.m2res)


#compare with AIC
AIC(blue.sc.wgt.m1,blue.sc.wgt.m2)

#           df   AIC
#blue.sc.wgt.m1  9 211.5126
#blue.sc.wgt.m2  5 209.0676

#Compare pollinators
emtrends(blue.sc.wgt.m1,pairwise~SPEC.COM,var="sumvisits")

#no difference between MX, SB or HB

ggplot(berry_updated_comp,aes(x=sumvisits,y=Fresh.wgt,col=SPEC.COM))+
  geom_point()+
  stat_smooth(method=lm,se=FALSE)+
  facet_wrap(~SPEC.COM)+
  scale_colour_brewer(palette = "Set1")+
  theme_bw()
