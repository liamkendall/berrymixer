##RASP SPECIES COMPOSITION MODELS

#libraries
library(plyr)
library(dplyr)
library(DHARMa)
library(emmeans)
library(ggplot2)
library(glmmTMB)

#left function
left = function (string,char){
  substr(string,1,char)
}
#load data
rasp <- read.csv("data/raspberry.csv", header=T)


#remove zeros
rasp=rasp[!rasp$Weight==0,]
rasp=rasp[,-which(names(rasp) %in% c("Comments","COMMENTS","POLLEN_DONOR"))]

#subset dataframe to just visitation treatment
#convert all visitor columns to character

rasp_updated <- rasp[!rasp$POLLINATORS%in%"BC",]
rasp_updated <- rasp_updated[!rasp_updated$POLLINATORS%in%"OC",]%>%droplevels()

rasp_updated <- rasp_updated %>% mutate_at(vars(16:44), as.character) %>% droplevels()

#retain only first letter of the string for visitors IDs
rasp_updated2 <- rasp_updated %>% mutate_at(vars(16:44), function (x) left(x,char=1))

#sum the number of visits from each taxa
rasp_updated2$stingless_bee <- rowSums(rasp_updated2[, c(16:44)] == "S")
rasp_updated2$honey_bee <- rowSums(rasp_updated2[, c(16:44)] == "H")

#sum the total number of visits
rasp_updated2$sumvisits <- rowSums(rasp_updated2[, c(47:48)])

#calculate the percent visits from each taxa
rasp_updated2$p_stingless_bee <- rasp_updated2$stingless_bee/rasp_updated2$sumvisits
rasp_updated2$p_honey_bee <- rasp_updated2$honey_bee/rasp_updated2$sumvisits


##CORRELATION BETWEEN QUALITY VARIABLES

#pairs plot
plot(rasp_updated2[,c("Weight","Druplets","Crumbliness")])

#pearson correlation
cor(rasp_updated2[,c("Weight","Druplets","Crumbliness")],method="pearson")

#            Weight    Druplets    Crumbliness
#Weight      1.0000000 0.7721321   0.5788748
#Druplets    0.7721321 1.0000000   0.7297822
#Crumbliness 0.5788748 0.7297822   1.0000000

##SPEC COM MODELS

rasp.sc.m1=glmmTMB(Weight~POLLINATORS*sumvisits+(1|BLOCK),
                   family=gaussian,
                   data=rasp_updated2)

rasp.sc.m1res=simulateResiduals(rasp.sc.m1)
plot(rasp.sc.m1res)
testResiduals(rasp.sc.m1res)


#simple visit model
rasp.sc.m2=glmmTMB(Weight~sumvisits+(1|BLOCK),rasp_updated2)
rasp.sc.m2res=simulateResiduals(rasp.sc.m2)
plot(rasp.sc.m2res)
testResiduals(rasp.sc.m2res)

summary(rasp.sc.m2)

#compare with AIC
AIC(rasp.sc.m1,rasp.sc.m2)

#           df   AIC
#rasp.sc.m1  8 853.1449
#rasp.sc.m2  4 845.3692

#Compare pollinators
emtrends(rasp.sc.m1,pairwise~POLLINATORS,var="sumvisits")

#no difference between MX, SB or HB

ggplot(rasp_updated2,aes(x=sumvisits,y=Weight,col=POLLINATORS))+
  geom_point()+
  stat_smooth(method=lm,se=FALSE)+
  facet_wrap(~POLLINATORS)+
  scale_colour_brewer(palette = "Set1")+
  theme_bw()
