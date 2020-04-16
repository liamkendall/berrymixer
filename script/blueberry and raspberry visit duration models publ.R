###Models of floral visit duration in relation to number of preceding visits in blueberry and raspberry

#libraries
library(reshape2)
library(tidyr)
library(cowplot)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(MuMIn)
library(DHARMa)
library(emmeans)
library(plyr)
library(dplyr)
library(stringi)
library(stringr)


#raspberry dataframe
rasp.length <- read.csv("data/rasp_visit_length.csv", header=T)

#blueberry visit duration dataframe
berry.length <- read.csv("data/blue_visit_length.csv", header=T)

#####
#Blueberry models
#####

#gaussian
blue.ord.vg.m=glmmTMB(v.d~poll.sp*ord.v+(1|Year/Block/Tag),
                      family="gaussian",berry.length)

#Gamma
blue.ord.gam.m=glmmTMB(v.d~poll.sp*ord.v+(1|Year/Block/Tag),
                       family=Gamma(link="log"),berry.length)

#lognormal
blue.ord.ln.m=glmmTMB(log.vd~poll.sp*ord.v+(1|Year/Block/Tag),
                      family="gaussian",berry.length)

AIC(blue.ord.vg.m, #gaussian
    blue.ord.gam.m, #gamma
    blue.ord.ln.m) #log-gaus

#blue.ord.ln.m  best

testResiduals(simulateResiduals(blue.ord.vg.m))
testResiduals(simulateResiduals(blue.ord.gam.m))
testResiduals(simulateResiduals(blue.ord.ln.m))

#check collinearity
performance::check_model(blue.ord.ln.m)

#marginal effects
pairs(emmeans(blue.ord.ln.m,~poll.sp))
emtrends(blue.ord.ln.m, pairwise~poll.sp,var="ord.v") #no differences between taxa

#####
#Raspberry models
#####

rasp.ord.v.m=glmmTMB(v.d~poll.sp*ord.v+(1|BLOCK/Tag),
                     family="gaussian",rasp.length)

rasp.ord.vg.m=glmmTMB(v.d~poll.sp*ord.v+(1|BLOCK/Tag),
                      family=Gamma(link="log"),rasp.length)

rasp.ord.ln.m=glmmTMB(log.vd~poll.sp*ord.v+(1|BLOCK/Tag),
                      family="gaussian",rasp.length)


AIC(rasp.ord.v.m,#gaussian
    rasp.ord.vg.m, #gamma
    rasp.ord.ln.m)#lognormal

#test residuals
testResiduals(simulateResiduals(rasp.ord.ln.m)) #D = 0.058245, p-value = 3.156e-06
testResiduals(simulateResiduals(rasp.ord.vg.m))
testResiduals(simulateResiduals(rasp.ord.v.m))

#collinearity
performance::check_collinearity(rasp.ord.ln.m)

##GRAPH DURATION VS VISIT NUMBER

#blueberry
b.ord.v.pred <- emmip(blue.ord.ln.m, poll.sp ~ ord.v,
                      mult.name = "poll.sp", 
                      cov.reduce = FALSE,
                      CIs = T, plotit = TRUE)

b.ord.v.pred <- b.ord.v.pred$data

#plot dataframe
berry.length3 <- berry.length

b.ord.v.pred$poll.sp<- revalue(b.ord.v.pred$poll.sp,c("H" = "Honeybee",
                                                      "S" = "Stingless bee"))
berry.length3$poll.sp <- revalue(berry.length3$poll.sp,c("H" = "Honeybee",
                                                     "S" = "Stingless bee"))

#plot the data
blue.dur.number.plot <- ggplot() + xlab("Floral visit number") + ylab("ln visit duration (s)") + 
  
  geom_jitter(data = berry.length3,aes(x = ord.v, y = log(v.d), colour=poll.sp), 
              size=2.5, shape = 16, width=0.2, height = 0.2,show.legend = TRUE,alpha=0.25)+ 
  geom_ribbon(data=b.ord.v.pred,
              aes(ymin=LCL, ymax=UCL, x=ord.v, fill=poll.sp), alpha = 0.5,show.legend = TRUE) + 
  geom_line(data=b.ord.v.pred, aes(ord.v,yvar, colour=poll.sp), size=1,show.legend = TRUE)+ 
  theme_bw() + 
  ggtitle("A) Blueberry")+
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
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(0.5,"lines"),
        legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.4))+
  scale_fill_brewer(palette="Set1") + 
  scale_colour_brewer(palette="Set1") + 
  scale_x_continuous(breaks = c(0,5,10,15))+
  #ylim(0,50)+
  #xlim(1,16)+
  labs(fill="Pollinator taxa",col="Pollinator taxa")
  
blue.dur.number.plot

##raspberry

r.ord.v.pred <- emmip(rasp.ord.ln.m, poll.sp ~ ord.v , mult.name = "poll.sp",
                      cov.reduce = FALSE,
                      CIs = T, plotit = TRUE)

r.ord.v.pred <- r.ord.v.pred$data

#plot dataframe
rasp.length3 <- rasp.length

r.ord.v.pred$poll.sp<- revalue(r.ord.v.pred$poll.sp,c("H" = "Honeybee",
                                                      "S" = "Stingless bee"))

rasp.length3$poll.sp <- revalue(rasp.length3$poll.sp,c("H" = "Honeybee",
                                                   "S" = "Stingless bee"))

#plot the data
rasp.dur.number.plot <- ggplot() + xlab("Floral visit number") + ylab("ln visit duration (s)") + 
  geom_jitter(data = rasp.length3%>% 
                drop_na(),aes(x = ord.v, y = log(v.d), colour=poll.sp), 
              size=2.5, shape = 16, width=0.2, height = 0.2,show.legend = FALSE,alpha=0.25)+ 
  geom_ribbon(data=r.ord.v.pred,
              aes(ymin=LCL, ymax=UCL, x=ord.v, fill=poll.sp), alpha = 0.5,show.legend = FALSE) + 
  geom_line(data=r.ord.v.pred, aes(ord.v,yvar, colour=poll.sp), size=1,show.legend = FALSE)+ 
  theme_bw() + 
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
        panel.spacing = unit(0.5,"lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.4))+
  scale_fill_brewer(palette="Set1") + 
  scale_colour_brewer(palette="Set1") + 
  labs(fill="Pollinator taxa",col="Pollinator taxa")+
  ggtitle("B) Raspberry")
  
rasp.dur.number.plot

##SLOPE CONTRASTS

#blueberry
emmeans(blue.ord.ln.m, ~poll.sp,transf = "exp")
b.ord.v.slopes <- emtrends(blue.ord.ln.m, pairwise~poll.sp,var="ord.v")
test(b.ord.v.slopes,adjust ="fdr")

#raspberry
r.ord.v.slopes <- emtrends(rasp.ord.ln.m2, pairwise~poll.sp,var="ord.v")
test(r.ord.v.slopes,adjust ="fdr")

#species diffs
blue.sp.dif <- emmeans(blue.ord.ln.m, pairwise ~ poll.sp)
rasp.sp.dif <- emmeans(rasp.ord.ln.m2, ~ poll.sp)

blue.sp.dif <- as.data.frame(blue.sp.dif)
rasp.sp.dif <- as.data.frame(rasp.sp.dif)

#back-transform
blue.sp.dif[,c("emmean","SE","lower.CL","upper.CL")] <- exp(blue.sp.dif[,c("emmean","SE","lower.CL","upper.CL")])
rasp.sp.dif[,c("emmean","SE","lower.CL","upper.CL")] <- exp(rasp.sp.dif[,c("emmean","SE","lower.CL","upper.CL")])

#GRID GRAPHS
Fig1=align_plots(blue.dur.number.plot, rasp.dur.number.plot,align="hv", axis="tblr")
Fig1A <- ggdraw(Fig1[[1]])
Fig1B <- ggdraw(Fig1[[2]])

Fig1_duration_number <- ggarrange(Fig1A,Fig1B,ncol=1,nrow=2,common.legend = TRUE)

ggsave(Fig1_duration_number,file="graphs/Figure 1A-B ln visit duration.pdf", height = 7, width = 5,dpi = 300)


