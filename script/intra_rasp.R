##Intraspecific raspberry
library(plyr)
library(dplyr)
library(DHARMa)

#load data
rasp_single <- read.csv("data/raspberry.csv", header=T)

#remove zeros
rasp_single=rasp_single[!rasp_single$Weight==0,]
rasp_single=rasp_single[,-which(names(rasp_single) %in% c("Comments","COMMENTS","POLLEN_DONOR"))]

#subset dataframe to just visitation treatment
#convert all visitor columns to character
rasp_single_updated <- rasp_single[!rasp_single$POLLINATORS%in%"BC",]
rasp_single_updated <- rasp_single_updated[!rasp_single_updated$POLLINATORS%in%"OC",]%>%droplevels()
rasp_single_updated <- rasp_single_updated[!rasp_single_updated$POLLINATORS%in%"Mixed",]%>%droplevels()

table(rasp_single_updated2$Total.visits)

rasp_single_updated <- rasp_single_updated %>% mutate_at(vars(16:44), as.character) %>% droplevels()

#retain only first letter of the string for visitors IDs
rasp_single_updated2 <- rasp_single_updated %>% mutate_at(vars(16:44), function (x) substr(x,2,2))

#sum the number of visits from each taxa
rasp_single_updated2$P <- rowSums(rasp_single_updated2[, c(16:44)] == "N",na.rm=TRUE)
rasp_single_updated2$N <- rowSums(rasp_single_updated2[, c(16:44)] == "P",na.rm=TRUE)

#sum the total number of visits
rasp_single_updated2$sumvisits <- rowSums(rasp_single_updated2[, c(47:48)])

#calculate the percent visits from each taxa
rasp_single_updated2$p_P <- rasp_single_updated2$P/rasp_single_updated2$sumvisits
rasp_single_updated2$p_N <- rasp_single_updated2$N/rasp_single_updated2$sumvisits

#drop either honeybees or stingless bees
rasp_single_HB <- rasp_single_updated2[!rasp_single_updated2$POLLINATORS%in%"SB",]%>%droplevels()
rasp_single_SB <- rasp_single_updated2[!rasp_single_updated2$POLLINATORS%in%"HB",]%>%droplevels()


##Weight model
intra_rasp_m1 <- glmmTMB(Weight~VISIT1*sumvisits+(1|BLOCK),
              family="gaussian",
              data = rasp_single_HB)

summary(intra_rasp_m1)
