###Models of floral visit duration in relation to number of preceding visits in blueberry

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

##Dataframe
write.csv(berry.length,file="data/blueberry duration dataset.csv")

berry.length <- read.csv("data/blueberry duration dataset.csv",
                         header=T,stringsAsFactors = FALSE)

#species composition split - solo and mixed species visit
berry.length$solo <- ifelse(berry.length$SPEC.COM=="MX","MX","solo")

#####
#Model choice
#####
berry.length$solo <- factor(berry.length$solo,levels=c("solo","MX"))

berry.length$log.ord.v <- log(log.ord.v$ord.v)

blue.visit.duration.all.m <- glmmTMB(v.d~poll.sp*ord.v+solo+
                                     (1|Year/Block/RP/Tag),
                                  family=nbinom2(link="log"),
                                  data=berry.length)

plot(simulateResiduals(blue.visit.duration.all.m))

blue.visit.duration.all.m2=glmmTMB(v.d~poll.sp*log.ord.v+solo+
                                     (1|Year/Block/RP/Tag),
                                   dispformula = ~0+poll.sp,
                                   family=truncated_nbinom1(link="log"),
                                   data=berry.length)

plot(simulateResiduals(blue.visit.duration.all.m2))

AIC(blue.visit.duration.all.m,
    blue.visit.duration.all.m2)

summary(blue.visit.duration.all.m2)

##compare trends from solo and mixed species models
duration.pred <- emmip(blue.visit.duration.all.m2, poll.sp ~log.ord.v,
                       transform="response",
                       mult.name = "poll.sp", 
                       cov.reduce = FALSE,
                       CIs = T, plotit = TRUE)

#stingless bees visit for a lot longer in mixed species visit than solo species visits

#plot the data
berry.length.plot <- berry.length

berry.length.plot$poll.sp <- revalue(berry.length.plot$poll.sp,c("H" = "Honeybee",
                                                         "S" = "Stingless bee"))
duration.pred$data$poll.sp<- revalue(duration.pred$data$poll.sp,c("H" = "Honeybee",
                                                         "S" = "Stingless bee"))

blue.dur.number.plot <- 
  ggplot() + 
  xlab("Floral visit number") + 
  ylab("Visit duration (s)") + 
  scale_y_log10(breaks=c(1,10,100))+
  geom_jitter(data = berry.length.plot,
              aes(x = ord.v, y = v.d,
                  colour=poll.sp,
                  fill=poll.sp), 
              size=4, shape = 21,col="white",
              width=0.2, height = 0.05,
              show.legend = TRUE,alpha=0.15)+
  geom_ribbon(data=duration.pred$data,
              aes(ymin=(LCL),
                  ymax=(UCL),
                  x=exp(log.ord.v),
                 # linetype=solo,
                  fill=poll.sp),
              alpha = 0.5,show.legend = T)+
  geom_line(data=duration.pred$data,
            aes(x=exp(log.ord.v),
                yvar,
                #linetype=solo,
                colour=poll.sp),
            size=1,show.legend = F)+ 
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
        legend.background = element_blank(),
        panel.border = element_rect(color = "black",
                                    fill = NA, size = 0.4))+
  scale_fill_manual(values = c("#2C3B75","#B8321A"))+
  scale_colour_manual(values = c("#2C3B75","#B8321A"))+
  scale_linetype_manual(values=c("solid","dashed"),guide = "none") +
  scale_x_continuous(breaks = c(0,5,10,15))+
  labs(name="Pollinator taxa",
       fill="Pollinator taxa",
       col="Pollinator taxa")
blue.dur.number.plot


ggsave(blue.dur.number.plot,file="graphs/Fig 1.pdf",
       height = 4, width = 5, dpi = 300)
ggsave(blue.dur.number.plot,file="graphs/Fig 1.jpg",
       height = 4, width = 5, dpi = 300)


