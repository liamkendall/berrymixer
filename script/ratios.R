#load packages
library(stringr)
library(plyr)
library(dplyr)

#load data
berry <- read.csv("data/BLUEBERRY_MIXED.csv", header=T)

#convert all visitor columns to character
berry_updated <- berry %>% mutate_at(vars(19:33), as.character)

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

#######################################
#END
#######################################