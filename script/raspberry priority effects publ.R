#Raspberry fruit weight models of priority effects and species differences

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
rasp_priority <- read.csv("data/rasp_priority.csv")

##Subset to 2 or more visits and only mixed visits
rasp_priority2 =rasp_priority[rasp_priority$POLLINATORS%in%"Mixed"
                            & rasp_priority$sumvisits >1,]%>%droplevels()

#remove mixed visits greater than 15 for priority effects model
rasp_priority2 <- rasp_priority2[!rasp_priority2$sumvisits >15,]%>%droplevels()

#log fruit weight
rasp_priority2$log.wgt <- log(rasp_priority2$Weight)

#####################################
#priority effects----
#####################################

#null
rasp_null.d <- glmmTMB(Weight~cent_d+(1|BLOCK),REML = FALSE,
                       family="gaussian",
                       data = rasp_priority2)

#run full priority effect model
rasp_m1.d <- glmmTMB(Weight~VISIT1*cent_sd+(1|BLOCK),REML = FALSE,
                     family="gaussian",
                     data = rasp_priority2)

#test model residuals
plot(simulateResiduals(rasp_m1.d))
testResiduals(simulateResiduals(rasp_m1.d))

rasp_m1.db <- glmmTMB(log(Weight)~VISIT1*cent_sd+(1|BLOCK),REML = FALSE,
                      family="gaussian",
                      data = rasp_priority2)

#test model residuals
plot(simulateResiduals(rasp_m1.db)) # log bad residuals

#go with raw values
rasp_m2.d <- glmmTMB(Weight~cent_v1*cent_sd+(1|BLOCK),REML = FALSE,
                     family="gaussian",
                     data = rasp_priority2)

plot(simulateResiduals(rasp_m2.d))
testResiduals(simulateResiduals(rasp_m2.d))

#check VIF
performance::check_collinearity(rasp_m1.d)
performance::check_collinearity(rasp_m2.d)

#likelihood ratio tests
#identity
anova(rasp_null.d,rasp_m1.d)

#duration
anova(rasp_null.d,rasp_m2.d)

#REML FOR PRESENTATION
rasp_m1.d.reml <- glmmTMB(Weight~VISIT1*cent_sd+(1|BLOCK),REML = TRUE,
                          family="gaussian",
                          data = rasp_priority2)

summary(rasp_m1.d.reml)

rasp_m2.d.reml <- glmmTMB(Weight~cent_v1*cent_sd+(1|BLOCK),REML = TRUE,
                          family="gaussian",
                          data = rasp_priority2)

summary(rasp_m2.d.reml)

#######################################
#species composition models
#######################################

#run the model
rasp_m2 <- glmmTMB(Weight~cent_d*POLLINATORS+(1|BLOCK),
                   family="gaussian",
                   data = rasp_priority)

rasp_m3 <- glmmTMB(log.wgt~cent_d2*POLLINATORS+(1|BLOCK),
                   family="gaussian",
                   data = rasp_priority)

rasp_m4 <- glmmTMB(Weight~cent_d2*POLLINATORS+(1|BLOCK),
                   family=Gamma(link="log"),
                   data = rasp_priority)

rasp_m5 <- glmmTMB(Weight~scale(sumvisits,center=TRUE)*POLLINATORS+(1|BLOCK),
                   family="gaussian",
                   data = rasp_priority)

AIC(rasp_m2,rasp_m3,rasp_m4,rasp_m5)

#df      AIC
#rasp_m2  8 808.6742
#rasp_m3  8 574.9701
#rasp_m4  8 900.9198

#additive
performance::check_collinearity(rasp_m2)
performance::check_collinearity(rasp_m3)
performance::check_collinearity(rasp_m4)
performance::check_collinearity(rasp_m5)

summary(rasp_m2)

#log  fail residual test
#check residuals
rasp_m2.res=simulateResiduals(rasp_m2)
plot(rasp_m2.res)
testResiduals(rasp_m2.res)

rasp_m3.res=simulateResiduals(rasp_m3)
plot(rasp_m3.res)
testResiduals(rasp_m3.res)

rasp_m4.res=simulateResiduals(rasp_m4)
plot(rasp_m4.res)
testResiduals(rasp_m4.res)

#so use un-logged as better residuals

#test that slopes for each taxa are different from 0
pairs(emmeans(rasp_m2,~POLLINATORS))
rasp.emm <- emtrends(rasp_m2, pairwise~POLLINATORS,var="cent_d")
test(rasp.emm, adjust="none")

#extract estimates etc
ras.comp.wght <- emmip(rasp_m2, POLLINATORS ~ cent_d, mult.name = "POLLINATORS", at = list(cent_d=unique(rasp_updated2$cent_d)),
                       cov.reduce = FALSE, CIs = T, plotit = TRUE)
ras.comp.wght <- ras.comp.wght$data

ras.comp.wght$xvar <- as.numeric(as.character(ras.comp.wght$xvar))

#rescale
ras.comp.wght$cent_d <- ras.comp.wght$cent_d*attr(rasp_updated2$cent_d, 'scaled:scale') +
  attr(rasp_updated2$cent_d, 'scaled:center')

ras.comp.wght <- ras.comp.wght[!(ras.comp.wght$POLLINATORS%in%"SB" & ras.comp.wght$cent_d > 175),] %>% droplevels()
ras.comp.wght <- ras.comp.wght[!(ras.comp.wght$POLLINATORS%in%"HB" & ras.comp.wght$cent_d > 250),] %>% droplevels()

#plot the data
rasp_priority4 <- rasp_priority

colnames(ras.comp.wght)[1] = "Taxa"
colnames(rasp_priority4)[16] = "Taxa"
ras.comp.wght$Taxa <- revalue(ras.comp.wght$Taxa,c("HB" = "Honeybee",
                                                   "SB" = "Stingless bee",
                                                   "Mixed" = "Mixed"))

ras.comp.wght$Taxa <- factor(ras.comp.wght$Taxa,levels=c("Honeybee","Stingless bee","Mixed"))

rasp_priority4$Taxa <- revalue(rasp_priority4$Taxa,c("HB" = "Honeybee",
                                                   "SB" = "Stingless bee",
                                                   "Mixed" = "Mixed"))

rasp_priority4$Taxa <- factor(rasp_priority4$Taxa,levels=c("Honeybee","Stingless bee","Mixed"))

rasp_comp <- ggplot() + 
  xlab("Total visit duration (s)") + 
  ylab("Fruit weight")+ 
  geom_point(data = rasp_priority4,aes(x = cent_d*attr(rasp_updated2$cent_d, 'scaled:scale') +
                                        attr(rasp_updated2$cent_d, 'scaled:center'), y = Weight, colour=Taxa), 
             size=3, alpha=0.5,
             shape = 16,show.legend = FALSE) +
  geom_ribbon(data=ras.comp.wght,
              aes(ymin=LCL, ymax=UCL, x=cent_d, fill=Taxa), 
              alpha = 0.5,show.legend = FALSE) +
  geom_line(data=ras.comp.wght, aes(x = cent_d,
                                    y = yvar,
                                    colour=Taxa), size=1,show.legend = FALSE) +
  theme_bw() + 
  ggtitle("B) Raspberry")+
  theme(plot.title = element_text(hjust = 0,face="bold",size=14),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=11),
        aspect.ratio = 1,
        strip.background = element_blank(),
        text=element_text(),
        axis.ticks.length = unit(0,"mm"),
        legend.title.align=0.5,
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        strip.text = element_blank(),
        panel.spacing = unit(0.25,"lines"),
        legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.4))+
       # plot.margin = unit(c(0.01,0.01,0.01,0.01), "lines"))+
  scale_fill_manual(values=c("#e41a1c","#377eb8","#984ea3"))+
  scale_colour_manual(values=c("#e41a1c","#377eb8","#984ea3"))+
  #scale_fill_brewer(palette="Set1") + 
  #scale_colour_brewer(palette="Set1") + 
  labs(fill="Pollinator taxa",col="Pollinator taxa")

rasp_comp
