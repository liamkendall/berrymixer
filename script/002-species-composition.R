################################ 
##blueberry species composition models
################################
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
library(performance)
library(patchwork)


#####
#palette
comp.cols <- c("Honeybee" = "#2C3B75",
               "Stingless bee" = "#B8321A",
               "Mixed" = "#565052")

################################ 
##visit duration
################################

##dataframe
berry <- read.csv("blue_priority.csv")

#Fresh weight ~ species composition
b.sc.simp.mod <- glmmTMB(Fresh.wgt~SPEC.COM+(1|Year/Block/RP),
                     family="gaussian",
                     data = berry)

summary(b.sc.simp.mod)

plot(simulateResiduals(b.sc.simp.mod)) #good

#emmeans
test(emmeans(b.sc.simp.mod,pairwise ~SPEC.COM,adjust="none"))

##plot
berry.box <- berry.max
berry.box$SPEC.COM <- revalue(berry.box$SPEC.COM,c("HB" = "Honeybee",
                                                   "SB" = "Stingless bee",
                                                   "MX" = "Mixed"))

berry.box$SPEC.COM <- factor(berry.box$SPEC.COM,levels=c("Honeybee","Stingless bee","Mixed"))

berry.boxplot <- ggplot(berry.box,aes(x=SPEC.COM,
                                      y=Fresh.wgt,
                                      col=SPEC.COM,
                                      fill=SPEC.COM))+
  geom_boxplot(show.legend = F,fill=NA)+
  geom_jitter(width=0.2, height = 0,col="white",
              show.legend = F,alpha=0.5,shape=21,size=4)+
  theme_bw()+
  ylim(0,2.9)+
  #labs(tag = 'A')+
  geom_signif(comparisons = list(c("Mixed", "Honeybee"),
                                 c("Mixed","Stingless bee")),
              annotation = c("***", "*"),
              col="black",textsize=4,
              y_position = c(2.75, 2.55))+
  theme(plot.title = element_text(hjust = 0,
                                  face="bold",
                                  size=14),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=11),
        aspect.ratio = 1,
        strip.background = element_blank(),
        text=element_text(),
        axis.ticks.length = unit(0,"mm"),
        legend.title.align=0.5,
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        strip.text = element_blank(),
        panel.spacing = unit(0.25,
                             "lines"),
        legend.background = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.4))+
  scale_fill_manual(values=comp.cols)+
  scale_colour_manual(values=comp.cols)+
  labs(fill="Pollinator taxon",
       col="Pollinator taxon")+
  xlab("Pollinator taxon")+
  ylab("Fruit weight (g)")
berry.boxplot


#separate models per species
berry.sp.split <- split(berry,berry$SPEC.COM)

#####DURATION MODELS FOR EACH SPECIES OR MIXED

#DURATION
hb.mod <- glmmTMB::glmmTMB(Fresh.wgt~duration+(1|Block/RP),data=berry.sp.split[[1]])
sb.mod <- glmmTMB::glmmTMB(Fresh.wgt~duration+(1|Block/RP),data=berry.sp.split[[3]])

berry.sp.split[[2]]$YRP <- paste(berry.sp.split[[2]]$Year,berry.sp.split[[2]]$RP)
mx.mod <- glmmTMB::glmmTMB(Fresh.wgt~duration+(1|YRP),data=berry.sp.split[[2]])

hb.res1=simulateResiduals(hb.mod)
plot(hb.res1)
testResiduals(hb.res1)#good

sb.res1=simulateResiduals(sb.mod)
plot(sb.res1)
testResiduals(sb.res1)#good

mx.res1=simulateResiduals(mx.mod)
plot(mx.res1)
testResiduals(mx.res1)#good

summary(hb.mod)
summary(sb.mod)
summary(mx.mod)

###DESCRIPTION OF RESULTS/PLOTTING
#extract estimates for plotting etc
hb.pred <- emmip(hb.mod,  ~ duration,  cov.reduce = FALSE, CIs = T, plotit = TRUE)
hb.pred <- hb.pred$data
hb.pred$SPEC.COM <- "Honeybee"

sb.pred <- emmip(sb.mod,  ~ duration,  cov.reduce = FALSE, CIs = T, plotit = TRUE)
sb.pred <- sb.pred$data
sb.pred$SPEC.COM <- "Stingless bee"

mx.pred <- emmip(mx.mod,  ~ duration,  cov.reduce = FALSE, CIs = T, plotit = TRUE)
mx.pred <- mx.pred$data
mx.pred$SPEC.COM <- "Mixed"

comp.preds <- rbind(hb.pred,sb.pred,mx.pred)

berry.comp.plot <- berry
berry.comp.plot$SPEC.COM <- revalue(berry.comp.plot$SPEC.COM,c("HB" = "Honeybee",
                                                               "SB" = "Stingless bee",
                                                               "MX" = "Mixed"))

berry.comp.plot$SPEC.COM <- factor(berry.comp.plot$SPEC.COM,levels=c("Stingless bee","Mixed","Honeybee"))


bomp.plot <- ggplot(berry.comp.plot,aes(x=duration,
                                        y=Fresh.wgt,
                                        col=SPEC.COM,
                                        fill=SPEC.COM))+
  geom_point(show.legend = F,alpha=0.5,
             shape = 21,col="white",
             size=4)+
  theme_bw()+
  geom_ribbon(data=comp.preds,
                aes(ymin=LCL,
                    ymax=UCL,
                    x=duration,
                    y=yvar),show.legend=F,
              alpha=0.3,
              colour = NA)+
    geom_line(data=comp.preds,
              aes(x = duration,
                  linetype=SPEC.COM,
                  y = yvar),
              show.legend=F,
              size=1)+
  xlab("Total visit duration (s)")+
  ylab(NULL)+
  ylim(0,2.9)+
  theme(plot.title = element_text(hjust = 0,
                                  face="bold",
                                  size=14),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=11),
        aspect.ratio = 1,
        strip.background = element_blank(),
        text=element_text(),
        axis.ticks.length = unit(0,"mm"),
        legend.title.align=0.5,
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        strip.text = element_blank(),
        panel.spacing = unit(0.25,
                             "lines"),
        legend.background = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.4))+
  scale_fill_manual(values=comp.cols)+
  scale_colour_manual(values=comp.cols)+
  scale_linetype_manual(values=c("solid","dashed","dashed"))+
  labs(fill="Pollinator taxon",
       col="Pollinator taxon",
       linetype="Pollinator taxon")

berry.comp.duration <- berry.boxplot+bomp.plot

ggsave(berry.comp.duration,
       file="graphs/Fig 2.pdf",
       height = 4,
       width = 6,
       dpi = 600)

ggsave(berry.comp.duration,
       file="graphs/Fig 2.jpg",
       device = "jpg",
       height = 4,
       width = 6,
       dpi = 600)

################################ 
##number of visits
################################

#####DURATION MODELS FOR EACH SPECIES OR MIXED

#DURATION
hb.sv.mod <- glmmTMB::glmmTMB(Fresh.wgt~sumvisits+(1|Block/RP),data=berry.sp.split[[1]])
sb.sv.mod <- glmmTMB::glmmTMB(Fresh.wgt~sumvisits+(1|Block/RP),data=berry.sp.split[[3]])
mx.sv.mod <- glmmTMB::glmmTMB(Fresh.wgt~sumvisits+(1|YRP),data=berry.sp.split[[2]])

hb.sv.res1=simulateResiduals(hb.sv.mod)
plot(hb.sv.res1)
testResiduals(hb.sv.res1)#good

sb.sv.res1=simulateResiduals(sb.sv.mod)
plot(sb.sv.res1)
testResiduals(sb.sv.res1)#good

mx.sv.res1=simulateResiduals(mx.sv.mod)
plot(mx.sv.res1)
testResiduals(mx.sv.res1)#good

summary(hb.sv.mod)
summary(sb.sv.mod)
summary(mx.sv.mod)

###########
#PLOTTING
##########

###DESCRIPTION OF RESULTS/PLOTTING
#extract estimates for plotting etc
hb.sv.pred <- emmip(hb.sv.mod,  ~ sumvisits,  cov.reduce = FALSE, CIs = T, plotit = TRUE)
hb.sv.pred <- hb.sv.pred$data
hb.sv.pred$SPEC.COM <- "Honeybee"

sb.sv.pred <- emmip(sb.sv.mod,  ~ sumvisits,
                    cov.reduce = FALSE, CIs = T, plotit = TRUE)
sb.sv.pred <- sb.sv.pred$data
sb.sv.pred$SPEC.COM <- "Stingless bee"

mx.sv.pred <- emmip(mx.sv.mod,  ~ sumvisits,  cov.reduce = FALSE, CIs = T, plotit = TRUE)
mx.sv.pred <- mx.sv.pred$data
mx.sv.pred$SPEC.COM <- "Mixed"

comp.sv.preds <- rbind(hb.sv.pred,sb.sv.pred,mx.sv.pred)

berry.comp.plot <- berry
berry.comp.plot$SPEC.COM <- revalue(berry.comp.plot$SPEC.COM,c("HB" = "Honeybee",
                                                               "SB" = "Stingless bee",
                                                               "MX" = "Mixed"))

berry.comp.plot$SPEC.COM <- factor(berry.comp.plot$SPEC.COM,levels=c("Honeybee","Stingless bee","Mixed"))

bomp.sv.plot <- ggplot(berry.comp.plot,aes(x=sumvisits,
                                        y=Fresh.wgt,
                                        col=SPEC.COM,
                                        fill=SPEC.COM))+
  geom_point(show.legend = F,
             alpha=0.5,
             shape = 21,col="white",
             size=4)+
  theme_bw()+
  geom_ribbon(data=comp.sv.preds,
              aes(ymin=LCL,
                  ymax=UCL,
                  x=sumvisits,
                  y=yvar),show.legend=F,
              alpha=0.3,
              colour = NA)+
  geom_line(data=comp.sv.preds,
            aes(x = sumvisits,
                linetype=SPEC.COM,
                y = yvar),
            show.legend=F,
            size=1)+
  xlab("Number of visits")+
  ylab(NULL)+
  ylim(0,2.9)+
  #xlim(1,15)+
  scale_x_continuous(breaks=seq(1,15,2))+
  theme(plot.title = element_text(hjust = 0,
                                  face="bold",
                                  size=14),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=11),
        aspect.ratio = 1,
        strip.background = element_blank(),
        text=element_text(),
        axis.ticks.length = unit(0,"mm"),
        legend.title.align=0.5,
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        strip.text = element_blank(),
        panel.spacing = unit(0.25,
                             "lines"),
        legend.background = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 0.4))+
  scale_fill_manual(values=comp.cols)+
  scale_colour_manual(values=comp.cols)+
  scale_linetype_manual(values=c("solid","dashed","dashed"))+
  labs(fill="Pollinator taxon",
       col="Pollinator taxon",
       linetype="Pollinator taxon")
bomp.sv.plot

berry.comp.sv <- berry.boxplot+bomp.sv.plot
berry.comp.sv

ggsave(berry.comp.sv,
       file="graphs/Fig S1.pdf",
       height = 4,
       width = 6,
       dpi = 600)

ggsave(berry.comp.sv,
       file="graphs/Fig S1.jpg",
       device = "jpg",
       height =4,
       width = 6,
       dpi = 600)

