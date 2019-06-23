####Single species duration models
##Single visits or multiple
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
library(ggplot2)

#load data
berry <- read.csv("data/BLUEBERRY_MIXED_V5.csv", header=T)

#subset dataframe to just visitation treatment
#convert all visitor columns to character
berry_updated <- berry[berry$TREATMENT%in%"V",] %>% mutate_at(vars(19:33), as.character) %>% droplevels()

#retain only first letter of the string for visitors IDs
berry_updated2 <- berry_updated %>% mutate_at(vars(19:33), word)

#sum the number of visits from each taxa
berry_updated2$stingless_bee <- rowSums(berry_updated2[, c(19:33)] == "S")
berry_updated2$honey_bee <- rowSums(berry_updated2[, c(19:33)] == "H")
berry_updated2$bumble_bee <- rowSums(berry_updated2[, c(19:33)] == "B")

length(berry_updated2)
#sum the total number of visits
berry_updated2$sumvisits <- rowSums(berry_updated2[, c(38:40)])

#calculate the percent visits from each taxa
berry_updated2$p_stingless_bee <- berry_updated2$stingless_bee/berry_updated2$sumvisits
berry_updated2$p_honey_bee <- berry_updated2$honey_bee/berry_updated2$sumvisits
berry_updated2$p_bumble_bee <- berry_updated2$bumble_bee/berry_updated2$sumvisits
berry_updated2$RP=paste0(berry_updated2$Row,berry_updated2$Plant.number)

#add duration variables (duration code)
str(berry_updated2)
str(berry_updated_time)

berry_updated2$duration=berry_updated_time$duration
berry_updated2$V1D=berry_updated_time$V1duration

berry_updated2[is.na(berry_updated2$V1D)==TRUE,]

#remove duplicated records
berry_updated2[duplicated(berry_updated2$Tag),]
berry_updated2 <- berry_updated2[!duplicated(berry_updated2$Tag),]%>%droplevels()

berry_updated2 <- berry_updated2[!(berry_updated2$SPEC.COM%in%"HB" &
                                   berry_updated2$sumvisits ==15),]%>%droplevels()

berry_updated2 <- berry_updated2[!(berry_updated2$SPEC.COM%in%"HB" &
                                     berry_updated2$sumvisits ==7),]%>%droplevels()

#remove Tasmania BB data and RE data 
#due to very low fresh weight sample size (10 Mixed visits)
table(berry_updated2[berry_updated2$Species%in%"RE" &
                       is.na(berry_updated2$Fresh.wgt)==FALSE,]$SPEC.COM)
berry_updated3 <- berry_updated2[!berry_updated2$Species%in%"BR",]%>%droplevels()
berry_updated3 <- berry_updated3[!berry_updated3$Species%in%"RE",]%>%droplevels()

##Subset to 2 or more visits and only mixed visits
berry.intra <- berry_updated3[!berry_updated3$SPEC.COM%in%c("MX"),]%>%droplevels()

#calculate log of berry weight
berry.intra$log.wgt <- log(berry.intra$Fresh.wgt)

#scale duration - first visit and full
berry.intra$V1D <- as.integer(berry.intra$V1D)
berry.intra$scale_v1 <- scale(berry.intra$V1D)
berry.intra$scale_d <- scale(berry.intra$duration)


###RUN THE MODELS

#priority effects for fruit weight
#################################

#V1D = Inital visitor duration
#Visitor1 = Inital visitor identity
#sumvisits = number of total visits
#scale_v1 = scaled first visit duration

table(berry.intra$Year,berry.intra$Block)

#remove single 2018 record (same as line that removes HB 15 visit'a)
berry.intra <- berry.intra[!berry.intra$Year%in%2018,]%>%droplevels()

#subset to single visits and 2 or more visits
berry.intra.1 <- berry.intra[berry.intra$sumvisits==1,]%>%droplevels()
table(berry.intra.1$SPEC.COM) #HB 39, SB 21

berry.intra.2 <- berry.intra[berry.intra$sumvisits>1,]%>%droplevels()
table(berry.intra.2$SPEC.COM,berry.intra.2$sumvisits) #HB 39, SB 21

#    2  3  4  5 #visit number
#HB 21 26 11 13
#SB 10  6  2  6

#single visit models
#log wgt
blue.int.s.mod <- glmmTMB(log.wgt~SPEC.COM*scale_v1+
                          (1|Block),
                        family=gaussian,
                        data = berry.intra.1)

blue.int.s.mod.2 <- glmmTMB(log.wgt~SPEC.COM+scale_v1+
                            (1|Block),
                          family=gaussian,
                          data = berry.intra.1)

#fresh wgt
blue.int.s.mod.3 <- glmmTMB(Fresh.wgt~SPEC.COM*scale_v1+
                            (1|Block),
                          family=gaussian,
                          data = berry.intra.1)

blue.int.s.mod.4 <- glmmTMB(Fresh.wgt~SPEC.COM+scale_v1+
                              (1|Block),
                            family=gaussian,
                            data = berry.intra.1)



#AIC.BIC
AIC(blue.int.s.mod,blue.int.s.mod.2,blue.int.s.mod.3,blue.int.s.mod.4)

#fresh wgt model with interaction best fit

#check residuals
blue.int.s.mod.3.res=simulateResiduals(blue.int.s.mod.3)
plot(blue.int.s.mod.3.res)
testResiduals(blue.int.s.mod.3.res)#looks good

#multiple visit models

#run with first visit duration and total duration
#can only be run on 2 or more visit

####visitor identity * visit number * duration first visit
#log wgt
blue.intm.mod <- glmmTMB(log.wgt~SPEC.COM*scale_v1*sumvisits+
                         (1|Block),
                       family=gaussian,
                       data = berry.intra.2)

#fresh wgt
blue.intm.mod.2 <- glmmTMB(Fresh.wgt~SPEC.COM*scale_v1*sumvisits+
                          (1|Block),
                        family=gaussian,
                        data = berry.intra.2)

####visitor identity * total duration * duration first visit
#log wgt
blue.intm.mod.3 <- glmmTMB(log.wgt~SPEC.COM*scale_v1*scale_d+
                            (1|Block),
                          family=gaussian,
                          data = berry.intra.2)

#fresh wgt
blue.intm.mod.4 <- glmmTMB(Fresh.wgt~SPEC.COM*scale_v1*scale_d+
                             (1|Block),
                           family=gaussian,
                           data = berry.intra.2)

AIC(blue.intm.mod,blue.intm.mod.2,blue.intm.mod.3,blue.intm.mod.4)
#Fresh wgt better fit than log weight

#df      AIC
#blue.intm.mod   10 83.69887
#blue.intm.mod.2 10 76.01434
#blue.intm.mod.3 10 69.19154
#blue.intm.mod.4 10 57.36443

#duration model much better
summary(blue.intm.mod.4)


#check residuals
#fresh wgt model
blue.intm.mod.2.res=simulateResiduals(blue.intm.mod.2)
plot(blue.intm.mod.2.res)
testResiduals(blue.intm.mod.2.res)#looks good

#full duration model
blue.intm.mod.4.res=simulateResiduals(blue.intm.mod.4)
plot(blue.intm.mod.4.res)
testResiduals(blue.intm.mod.4.res)#looks good

#dredge it and print csv
blue.intm.2.d <- dredge(blue.intm.mod.2,rank="AICc") #visit number
blue.intm.4.d <- dredge(blue.intm.mod.4,rank="AICc") #full duration

b.intm2.mods <-get.models(blue.intm.2.d,subset=TRUE)
b.intm4.mods <-get.models(blue.intm.4.d,subset=TRUE)

#summarise best priority effects model
summary(b.intm2.mods[[1]])
summary(b.intm4.mods[[1]])

r.squaredGLMM(b.intm2.mods[[1]])
r.squaredGLMM(b.intm4.mods[[1]])

summary(model.avg(blue.intm.2.d,subset = delta <5))
summary(model.avg(blue.intm.4.d,subset = delta <5))

#model with total duration outperforms any priority effect model
#negative interaction between total duration and stingless bees
#

###Model predictions and plots == multi.visit | visit duration dataset 
#extract estimates etc
b.d.pred <- emmip(b.intm4.mods[[1]], SPEC.COM ~ scale_d , mult.name = "SPEC.COM", cov.reduce = FALSE, CIs = T, plotit = TRUE)

b.d.pred <- b.d.pred$data



b.d.pred$xvar <- as.numeric(as.character(b.d.pred$xvar))

colnames(b.d.pred)=c("Taxa","Duration","Fruit weight","SE","df","LCL","UCL","Taxa2","xvar")

b.d.pred$Taxa <- revalue(b.d.pred$Taxa,c("HB" = "Honeybee",
                                                   "SB" = "Stingless bee"))

##Remove honeybee predictions beyond extent of dataset
b.d.pred <- b.d.pred[!(b.d.pred$Taxa%in%c("Honeybee") &
           b.d.pred$Duration >1),]

berry.intra.3 <- berry.intra.2

berry.intra.3$Taxa <- revalue(berry.intra.3$SPEC.COM,c("HB" = "Honeybee",
                                         "SB" = "Stingless bee"))

ggplot(b.d.pred,aes(x=Duration,y=`Fruit weight`)) + 
  xlab("Duration") + ylab("Fruit weight (g)") + 
  geom_ribbon(aes(ymin=LCL, ymax=UCL, x=Duration, fill=Taxa), alpha = 0.5,
              show.legend = FALSE)+
  geom_line(aes(Duration,`Fruit weight`, colour=Taxa), size=1,show.legend = FALSE)  +
  geom_jitter(data = berry.intra.3, aes(x = scale_d, y = Fresh.wgt, colour=Taxa),
              size=2.5, shape = 21, width=0, height =0.04,show.legend = FALSE)+
  theme_bw() + 
  theme(plot.title = element_text(vjust = 0,hjust = 0,face = "bold",size=16,family="Helvetica"),
        legend.box.background = element_rect(colour = "black"),
        aspect.ratio = 1,
        strip.background = element_blank(),
        text=element_text(),
        axis.ticks.length = unit(0,"mm"),
        legend.title.align=0.5,
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.text = element_text(size=15,face = "bold"),
        panel.spacing = unit(0.25,"lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.4),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm"))+
  scale_fill_brewer(palette="Set1") + 
  scale_colour_brewer(palette="Set1") + 
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
  facet_wrap(~Taxa,scales="free_x")
       
       
#Single visit predictions and graph

blue.int.s.mod.3

b.d.s.pred <- emmip(blue.int.s.mod.3, SPEC.COM ~ scale_v1 , mult.name = "SPEC.COM", cov.reduce = FALSE, CIs = T, plotit = TRUE)

b.d.s.pred <- b.d.s.pred$data

b.d.s.pred$xvar <- as.numeric(as.character(b.d.s.pred$xvar))

colnames(b.d.s.pred)=c("Taxa","Duration","Fruit weight","SE","df","LCL","UCL","Taxa2","xvar")

b.d.s.pred$Taxa <- revalue(b.d.s.pred$Taxa,c("HB" = "Honeybee",
                                         "SB" = "Stingless bee"))

##Remove honeybee predictions beyond extent of dataset
b.d.s.pred <- b.d.s.pred[!(b.d.s.pred$Taxa%in%c("Honeybee") &
                         b.d.s.pred$Duration >1),]

b.d.s.pred <- b.d.s.pred[!(b.d.s.pred$Taxa%in%c("Stingless bee") &
                             b.d.s.pred$Duration >3),]

berry.intra.4 <- berry.intra.1

berry.intra.4$Taxa <- revalue(berry.intra.4$SPEC.COM,c("HB" = "Honeybee",
                                                       "SB" = "Stingless bee"))


ggplot(b.d.s.pred,aes(x=Duration,y=`Fruit weight`)) + 
  xlab("Duration") + ylab("Fruit weight (g)") + 
  geom_ribbon(aes(ymin=LCL, ymax=UCL, x=Duration, fill=Taxa), alpha = 0.5,
              show.legend = FALSE)+
  geom_line(aes(Duration,`Fruit weight`, colour=Taxa), size=1,show.legend = FALSE)  +
  geom_jitter(data = berry.intra.4, aes(x = scale_d, y = Fresh.wgt, colour=Taxa),
              size=2.5, shape = 21, width=0, height =0.04,show.legend = FALSE)+
  theme_bw() + 
  theme(plot.title = element_text(vjust = 0,hjust = 0,face = "bold",size=16,family="Helvetica"),
        legend.box.background = element_rect(colour = "black"),
        aspect.ratio = 1,
        strip.background = element_blank(),
        text=element_text(),
        axis.ticks.length = unit(0,"mm"),
        legend.title.align=0.5,
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.text = element_text(size=15,face = "bold"),
        panel.spacing = unit(0.25,"lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.4),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm"))+
  scale_fill_brewer(palette="Set1") + 
  scale_colour_brewer(palette="Set1") + 
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
  facet_wrap(~Taxa,scales="free_x")

#what a ugly one

r.squaredGLMM(blue.int.s.mod.3)

summary(lm(Fresh.wgt~SPEC.COM*scale_v1, data = berry.intra.1))
summary(glm(FS~SPEC.COM*scale_v1,
            family="binomial",
            data = berry.intra.1))
