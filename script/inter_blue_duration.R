#load packages
library(stringr)
library(plyr)
library(dplyr)
library(glmmTMB)
library(lme4)
library(ggplot2)
library(MuMIn)
library(DHARMa)
library(emmeans)
library(sjPlot)

#load data
berry <- read.csv("data/BLUEBERRY_MIXED_V5.csv", header=T)

#subset dataframe to just visitation treatment
#convert all visitor columns to character
berry_updated <- berry[berry$TREATMENT%in%"V",] %>% mutate_at(vars(19:33), as.character) %>% droplevels()

#retain only first letter of the string for visitors IDs
berry_updated_time <- berry_updated %>% mutate_at(vars(19:33), function (x) word(x,3))

#sum the number of visits from each taxa
berry_updated_time$duration<-rowSums(apply(berry_updated_time[,19:33],
                         MARGIN=2,FUN = unlist(as.numeric)),na.rm=TRUE)
berry_updated_time$V1duration<-as.numeric(berry_updated_time$Visitor1)

#berry_updated2