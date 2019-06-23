#####RASPY
####Single species duration models
##Single visits or multiple

#load packages
library(stringr)
library(stringi)
library(plyr)
library(dplyr)
library(glmmTMB)
library(lme4)
library(ggplot2)
library(MuMIn)
library(DHARMa)
library(emmeans)
library(sjPlot)
library(scales)

#load data
rasp <- read.csv("data/raspberry.csv", header=T)

#remove zeros
rasp <- rasp[!rasp$Weight==0,]
rasp <- rasp[,-which(names(rasp) %in% c("Comments","COMMENTS","POLLEN_DONOR"))]

#subset dataframe to just visitation treatment
#convert all visitor columns to character
rasp_updated <- rasp[!rasp$POLLINATORS%in%"BC",]
rasp_updated <- rasp_updated[!rasp_updated$POLLINATORS%in%"OC",]%>%droplevels()
rasp_updated <- rasp_updated %>% mutate_at(vars(16:44), as.character) %>% droplevels()

#retain only first letter of the string for visitors IDs
#left function
left = function (string,char){
  substr(string,1,char)
}
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

rasp_updated2 <- rasp_updated %>% 
                 mutate_at(vars(16:44), function (x) left(x,char=1))

rasp_updated_time <- rasp_updated %>%
                     mutate_at(vars(16:44), function (x) stri_sub(x,3))

rasp_updated_time$duration <- rowSums(apply(rasp_updated_time[,16:44],
                                          MARGIN=2,FUN = unlist(as.numeric)),na.rm=TRUE)

rasp_updated_time$V1duration <- as.numeric(rasp_updated_time$VISIT1)
rasp_updated_time$sub_duration <- rasp_updated_time$duration -rasp_updated_time$V1duration

#sum the number of visits from each taxa
rasp_updated2$stingless_bee <- rowSums(rasp_updated2[, c(16:44)] == "S")
rasp_updated2$honey_bee <- rowSums(rasp_updated2[, c(16:44)] == "H")

#sum the total number of visits
rasp_updated2$sumvisits <- rowSums(rasp_updated2[, c(47:48)])

#calculate the percent visits from each taxa
rasp_updated2$p_stingless_bee <- rasp_updated2$stingless_bee/rasp_updated2$sumvisits
rasp_updated2$p_honey_bee <- rasp_updated2$honey_bee/rasp_updated2$sumvisits

#add duration variables (duration code)
rasp_updated2$duration <- rasp_updated_time$duration
rasp_updated2$V1D <- rasp_updated_time$V1duration
rasp_updated2$sub.d <- rasp_updated_time$sub_duration

#remove sum visits greater than 20
rasp_updated2 <- rasp_updated2[!rasp_updated2$sumvisits >20,]%>%droplevels()

#remove sb or hb visits > 10
rasp_updated2 <- rasp_updated2[!(rasp_updated2$POLLINATORS%in%"SB" & rasp_updated2$sumvisits > 10),] %>% droplevels()
rasp_updated2 <- rasp_updated2[!(rasp_updated2$POLLINATORS%in%"HB" & rasp_updated2$sumvisits > 10),] %>% droplevels()

##Subset to 2 or more visits and only mixed visits
rasp.intra <- rasp_updated2[!rasp_updated2$POLLINATORS%in%"Mixed",] %>% droplevels()

#Scale duration variables
rasp.intra$scale_v1=scale(rasp.intra$V1D)
rasp.intra$scale_d=scale(rasp.intra$duration)
rasp.intra$scale_sd=scale(rasp.intra$sub.d)

#log fruit weight
rasp.intra$log.wgt <- log(rasp.intra$Weight)

#single visit dataset
rasp.intra.1 <- rasp.intra[rasp.intra$sumvisits ==1,]%>%droplevels()

table(rasp.intra.1$POLLINATORS) #11 HB 11 SB - small..

#multi-visit dataset
rasp.intra.2 <- rasp.intra[rasp.intra$sumvisits >1,]%>%droplevels()

table(rasp.intra.2$POLLINATORS) #86 HB 12 SB - small SBs

#################MODELS
###RUN THE MODELS

#priority effects for fruit weight
#################################

#V1D = Inital visitor duration
#Visitor1 = Inital visitor identity
#sumvisits = number of total visits
#scale_v1 = scaled first visit duration


#single visit models
#log wgt
rasp.int.s.mod <- glmmTMB(log.wgt~POLLINATORS*scale_v1+
                            (1|BLOCK),
                          family=gaussian,
                          data = rasp.intra.1)

rasp.int.s.mod.2 <- glmmTMB(log.wgt~POLLINATORS+scale_v1+
                              (1|BLOCK),
                            family=gaussian,
                            data = rasp.intra.1)

#fresh wgt
rasp.int.s.mod.3 <- glmmTMB(Weight~POLLINATORS*scale_v1+
                              (1|BLOCK),
                            family=gaussian,
                            data = rasp.intra.1)

rasp.int.s.mod.4 <- glmmTMB(Weight~POLLINATORS+scale_v1+
                              (1|BLOCK),
                            family=gaussian,
                            data = rasp.intra.1)



#AIC.BIC
AICc(rasp.int.s.mod,rasp.int.s.mod.2,rasp.int.s.mod.3,rasp.int.s.mod.4,rasp.int.s.mod.5)

#log wgt model without interaction best fit

#check residuals
rasp.int.s.mod.2.res=simulateResiduals(rasp.int.s.mod.2)
plot(rasp.int.s.mod.2.res)
testResiduals(rasp.int.s.mod.2.res)#looks real bad

#Run without block term

rasp.int.s.mod.5 <- glmmTMB(log.wgt~POLLINATORS+scale_v1,
                            family=gaussian,
                            data = rasp.intra.1)
summary(rasp.int.s.mod.5)

#check residuals
rasp.int.s.mod.5.res=simulateResiduals(rasp.int.s.mod.5)
plot(rasp.int.s.mod.5.res)
testResiduals(rasp.int.s.mod.5.res)#looks real bad

#multiple visit models

#run with first visit duration and total duration
#can only be run on 2 or more visit

####visitor identity * visit number * duration first visit
#log wgt
rasp.intm.mod <- glmmTMB(log.wgt~POLLINATORS*scale_v1*sumvisits+
                           (1|BLOCK),
                         family=gaussian,
                         data = rasp.intra.2)

#fresh wgt
rasp.intm.mod.2 <- glmmTMB(Weight~POLLINATORS*scale_v1*sumvisits+
                             (1|BLOCK),
                           family=gaussian,
                           data = rasp.intra.2)

####visitor identity * total duration * duration first visit
#log wgt
rasp.intm.mod.3 <- glmmTMB(log.wgt~POLLINATORS*scale_v1*scale_sd+
                             (1|BLOCK),
                           family=gaussian,
                           data = rasp.intra.2)

#fresh wgt
rasp.intm.mod.4 <- glmmTMB(Weight~POLLINATORS*scale_v1*scale_sd+
                             (1|BLOCK),
                           family=gaussian,
                           data = rasp.intra.2)

AICc(rasp.intm.mod,rasp.intm.mod.2,rasp.intm.mod.3,rasp.intm.mod.4)
#log wgt better fit than fresh weight

#                df      AIC
#rasp.intm.mod   10 255.0741
#rasp.intm.mod.2 10 339.7675
#rasp.intm.mod.3 10 250.4430
#rasp.intm.mod.4 10 332.9159

#duration model much better
summary(rasp.intm.mod.3)


#check residuals
#visit number model
rasp.intm.mod.res=simulateResiduals(rasp.intm.mod)
plot(rasp.intm.mod.res)
testResiduals(rasp.intm.mod.res)#could be better

#full duration model
rasp.intm.mod.3.res=simulateResiduals(rasp.intm.mod.3)
plot(rasp.intm.mod.3.res)
testResiduals(rasp.intm.mod.3.res)#could be better

#dredge it
rasp.intm.d <- dredge(rasp.intm.mod,rank="AICc") #visit number
rasp.intm.3.d <- dredge(rasp.intm.mod.3,rank="AICc") #full duration

r.intm.mods <-get.models(rasp.intm.d,subset=TRUE)
r.intm3.mods <-get.models(rasp.intm.3.d,subset=TRUE)

#summarise top mods
summary(r.intm.mods[[1]])
summary(r.intm3.mods[[1]])

r.squaredGLMM(r.intm.mods[[1]])
r.squaredGLMM(r.intm3.mods[[1]])

summary(model.avg(rasp.intm.d,subset = delta <4))
summary(model.avg(rasp.intm.3.d,subset = delta <4))

#nothing going on here except duration overall


