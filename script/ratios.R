#load packages
library(stringr)
library(plyr)
library(dplyr)
library(glmmTMB)

#load data
berry <- read.csv("data/BLUEBERRY_MIXED_V2.csv", header=T)

#subset dataframe to just visitation treatment
#convert all visitor columns to character
berry_updated <- berry[berry$TREATMENT%in%"V" & berry$TVN > 0,] %>% mutate_at(vars(19:33), as.character) %>% droplevels()

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

#priority effects models
m1 <- glmmTMB(FS~Visitor1*(TVN/p_honey_bee)+(1|Block/Plant.number)+(1|Year),
                        family="binomial",
                        data = berry_updated2)

m2 <- glmmTMB(Fresh.wgt~Visitor1*(TVN/p_honey_bee)+(1|Block/Plant.number)+(1|Year),
              family="gaussian",
              data = berry_updated2)

#model cases where HD:SB ratios are 0:1 or 1:1
#no cases where >1 visits include only singless bess
berry_updated3 <- berry_updated2[!berry_updated2$Species%in%"BR",]
keep <- c(0, 0.5, 1)
equal_ratio <- berry_updated3[berry_updated3$p_stingless_bee %in% keep, ]

equal_ratio$p_honey_bee <- as.factor(as.character(equal_ratio$p_honey_bee))
m3 <- glmmTMB(FS~p_honey_bee+(1|Block/Plant.number)+(1|Year),
              family="binomial",
              data = equal_ratio)

m4 <- glmmTMB(Fresh.wgt~p_honey_bee+(1|Block/Plant.number)+(1|Year),
              family="gaussian",
              data = equal_ratio)

#######################################
#END
#######################################