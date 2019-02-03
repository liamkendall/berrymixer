##RASP

library(plyr)
library(dplyr)
library(DHARMa)

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

#priority effect model
rasp_m1 <- glmmTMB(Weight~VISIT1*(sumvisits/p_honey_bee)+(1|BLOCK),
              family="gaussian",
              data = rasp_updated2)
summary(rasp_m1)
rasp_m1res=simulateResiduals(rasp_m1)
plot(rasp_m1res)

rasp_m2 <- glmmTMB(Weight~sumvisits*POLLINATORS+(1|BLOCK),
            family="gaussian",
            data = rasp_updated2)

summary(rasp_m2)
rasp_m2res=simulateResiduals(rasp_m2)
plot(rasp_m2res)
