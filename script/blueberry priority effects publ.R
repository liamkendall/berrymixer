#Blueberry fruit weight models of priority effects and species differences

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
library(cowplot)
library(ggpubr)

#load data
berry <- read.csv("blue_priority.csv")

#scale visitation variables
berry$cent_v1 <- scale(berry$V1D,center=TRUE,scale=TRUE)
berry$cent_d  <- scale(berry$duration,center=TRUE,scale=TRUE)
berry$cent_sd  <- scale((berry$duration-berry$V1D),center=TRUE,scale=TRUE)

#log berry weight if needed
berry$log.wgt <- log(berry$Fresh.wgt)

##Subset to 2 or more visits and only mixed visits
berry_pe=berry[berry$SPEC.COM%in%"MX"&berry$sumvisits >1,]%>%droplevels()

#########
##priority effects models
#########

#V1D = Inital visitor duration
#Visitor1 = Inital visitor identity
#sumvisits = number of total visits
#scale_v1 = scaled first visit duration
#scale_d = scaled total duration

#####
#null models
#####

#check distributions

blue.null1.mod <- glmmTMB(Fresh.wgt~cent_d+
                            (1|Year),REML=FALSE,
                          family="gaussian",
                          data = berry_pe)

blue.null2.mod <- glmmTMB(log.wgt~cent_d+
                           (1|Year),REML=FALSE,
                         family=gaussian,
                         data = berry_pe)

blue.null3.mod <- glmmTMB(Fresh.wgt~cent_d+
                            (1|Year),REML=FALSE,
                          family="Gamma",
                          data = berry_pe)
#test residuals

testResiduals(simulateResiduals(blue.null1.mod)) #good
testResiduals(simulateResiduals(blue.null2.mod)) #bad
testResiduals(simulateResiduals(blue.null3.mod)) #bad

#visitor identity * duration
#gaussian untransformmed best residuals

blue.sp.mod <- glmmTMB(Fresh.wgt~Visitor1*cent_sd+
                         (1|Year),REML=FALSE,
                       family=gaussian,
                       data = berry_pe)

testResiduals(simulateResiduals(blue.sp.mod)) #good

#inital duration * total duration
blue.dr.mod <- glmmTMB(Fresh.wgt~cent_v1*cent_sd+
                         (1|Year),REML=FALSE,
                       family=gaussian,
                       data = berry_pe)

testResiduals(simulateResiduals(blue.dr.mod)) #good

#VIF
performance::check_collinearity(blue.sp.mod)
performance::check_collinearity(blue.dr.mod)

#AIC
AICc(blue.sp.mod,blue.dr.mod,blue.null2.mod)

#species model test
anova(blue.null1.mod,blue.sp.mod)
AIC(blue.null1.mod)-AIC(blue.sp.mod)

#duration model test
anova(blue.null1.mod,blue.dr.mod)
AIC(blue.null1.mod)-AIC(blue.dr.mod)

#difference between species identity and visit duration
anova(blue.dr.mod,blue.sp.mod)
AIC(blue.sp.mod)-AIC(blue.dr.mod)

#Re-fit with REML

#visitor identity * duration
blue.sp.mod.reml <- glmmTMB(Fresh.wgt~Visitor1*cent_sd+
                         (1|Year),REML=TRUE,
                       family=gaussian,
                       data = berry_pe)

#inital duration * total duration
blue.dr.mod.reml <- glmmTMB(Fresh.wgt~cent_v1*cent_sd+
                         (1|Year),REML=TRUE,
                       family=gaussian,
                       data = berry_pe)

summary(blue.sp.mod.reml)
summary(blue.dr.mod.reml)

#test slopes dependent on first visitor
test(emtrends(blue.sp.mod.reml, pairwise~Visitor1,var="cent_sd"))

#for description
#calculate percent difference at two visits
#average from two visits 

#average duration of 1 minute of visitation
scale_low <- scale((60-mean(berry_pe[,c("V1D")])),
      attr(berry$cent_sd, "scaled:center"), attr(berry$cent_sd, "scaled:scale"))

scale_high <- scale((480-mean(berry_pe[,c("V1D")])),
                    attr(berry$cent_sd, "scaled:center"), attr(berry$cent_sd, "scaled:scale"))

#With stingless bee - low visitation
sb_low=blue.sp.mod.reml$fit$parfull[1]+blue.sp.mod.reml$fit$parfull[2]+
         blue.sp.mod.reml$fit$parfull[3]*scale_low[1]+blue.sp.mod.reml$fit$parfull[4]*scale_low[1]


hb_low=blue.sp.mod.reml$fit$parfull[1]+blue.sp.mod.reml$fit$parfull[3]*scale_low[1]

((sb_low-hb_low)/hb_low)*100

#57.28904 %

sb_high=blue.sp.mod.reml$fit$parfull[1]+blue.sp.mod.reml$fit$parfull[2]+
            blue.sp.mod.reml$fit$parfull[3]*scale_high[1]+blue.sp.mod.reml$fit$parfull[4]*scale_high[1]

hb_high=blue.sp.mod.reml$fit$parfull[1]+blue.sp.mod.reml$fit$parfull[3]*scale_high[1]

((hb_high-sb_high)/sb_high)*100

#-24.23534 

#extract and format model predictions
wght.pred <- emmip(blue.sp.mod.reml, Visitor1 ~ cent_sd, mult.name = "Visitor1",
                   cov.reduce = FALSE, CIs = T, plotit = TRUE)

wght.pred <- wght.pred$data
wght.pred$xvar <- as.numeric(as.character(wght.pred$xvar))
wght.pred$xvar <- wght.pred$xvar * attr(berry$cent_sd, 'scaled:scale') + attr(berry$cent_sd, 'scaled:center')
wght.pred$cent_sd <- wght.pred$cent_sd* attr(berry$cent_sd, 'scaled:scale') + attr(berry$cent_sd, 'scaled:center')
wght.pred$Visitor1 <- revalue(wght.pred$Visitor1,c("H" = "Honeybee",
                                                   "S" = "Stingless bee"))

#plot dataframe
berry_plot <- berry_pe

berry_plot$Visitor1 <- revalue(berry_plot$Visitor1,c("H" = "Honeybee",
                                                             "S" = "Stingless bee"))

blue.identity.plot <- ggplot() + 
  xlab("Post-initial visit duration (s)") + 
  ylab("Blueberry fruit weight (g)") + 
  geom_jitter(data = berry_plot,aes(x = cent_sd * attr(berry$cent_sd, 'scaled:scale') + attr(berry$cent_sd, 'scaled:center'), y = Fresh.wgt, colour=Visitor1), 
              size=3, shape = 16, width=0, height =0) + 
  geom_ribbon(data=wght.pred,
              aes(ymin=LCL, ymax=UCL, x=xvar, fill=Visitor1), alpha = 0.35) + 
  geom_line(data=wght.pred, aes(x = cent_sd,
                                y = yvar, colour=Visitor1), size=1) + 

  theme_bw()+ 
  #ylim(0,3.2)+
  #xlim(600,657)+
  ggtitle("B) Initial species identity")+
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
               strip.text = element_text(face = "bold",size=14,family="Helvetica"),
               panel.spacing = unit(0.5,"lines"),
               panel.border = element_rect(color = "black", fill = NA, size = 0.4))+
   scale_fill_brewer(palette="Set1")+
   scale_colour_brewer(palette="Set1")+
   labs(fill="Initial visitor",col="Initial visitor")

blue.identity.plot

#Save plot
#ggsave(blue.identity.plot,file="graphs/blueberry_priority_effect.pdf", height =6, width=8,dpi=300)

###PLOT - initial visit duration
blue.predict.grid <- data.frame(cent_v1 =rep(c(10,30,60,90,120),26),
                                 cent_sd =rep(c(25,25,25,25,25,
                                                50,50,50,50,50,
                                                100,100,100,100,100,
                                                150,150,150,150,150,
                                                200,200,200,200,200,
                                                250,250,250,250,250,
                                                300,300,300,300,300,
                                                350,350,350,350,350,
                                                400,400,400,400,400,
                                                450,450,450,450,450,
                                                500,500,500,500,500,
                                                650,650,650,650,650,
                                                600,600,600,600,600),2),
                                Year=rep(c(2017,2018),130))
                                
berry$rem_d <- berry$duration-berry$V1D

blue.predict.grid$cent_sd <- scale(blue.predict.grid$cent_sd, attr(berry$cent_sd, "scaled:center"), attr(berry$cent_sd, "scaled:scale"))

blue.predict.grid$cent_v1 <- scale(blue.predict.grid$cent_v1 , attr(berry$cent_v1, "scaled:center"), attr(berry$cent_v1, "scaled:scale"))

blue.predicted <- cbind(predict(blue.dr.mod.reml,newdata = blue.predict.grid,se.fit = TRUE),blue.predict.grid)

colnames(blue.predicted)[1] = "pred"

#Calculate 95% C.I's
blue.predicted$LCL <- blue.predicted$pred - (1.98*blue.predicted$se.fit)
blue.predicted$UCL <- blue.predicted$pred + (1.98*blue.predicted$se.fit)

blue.predicted$cent_v1 <-  blue.predicted$cent_v1* attr(berry$cent_v1, 'scaled:scale') +
  attr(berry$cent_v1, 'scaled:center')

blue.predicted$cent_sd <- blue.predicted$cent_sd* attr(berry$cent_sd, 'scaled:scale') +
                          attr(berry$cent_sd, 'scaled:center')

blue.predicted$cent_v1 <- as.factor(blue.predicted$cent_v1)

#extract estimates etc
dur.w.pred <- emmip(blue.dr.mod.reml, cent_v1 ~ cent_sd, mult.name = "cent_v1", 
                    cov.reduce = FALSE, CIs = T, plotit = TRUE)
dur.w.pred <- dur.w.pred$data
dur.w.pred$V1D <- dur.w.pred$cent_v1*attr(berry$cent_v1, 'scaled:scale') +
  attr(berry$cent_v1, 'scaled:center')

#plot dataframe
berry_updated7 <- berry_pe
berry_updated7$rem_d <- berry_updated7$duration-berry_updated7$V1D

#plot the data
dur.w.pred$cent_v1 <- as.factor(dur.w.pred$cent_v1*attr(berry$cent_v1, 'scaled:scale') +
                                   attr(berry$cent_v1, 'scaled:center'))
dur.w.pred$cent_v1

blue.predicted2 <- blue.predicted[blue.predicted$cent_v1%in%c(10,120),]

head(blue.predicted2)
str(blue.predicted2)
blue.duration.priority.plot <- ggplot(blue.predicted2,aes(x=cent_sd,y=pred,col=cent_v1))+
  xlab("Post-initial visit duration (s)")+
  ylab("Blueberry fruit weight (g)")+
  geom_point(data = berry_updated7,aes(x = rem_d, y = Fresh.wgt),col="black",alpha=0.5,
             size=3, shape = 16) +
  geom_ribbon(data=blue.predicted2,
              aes(ymin=LCL, ymax=UCL, x=cent_sd, fill=cent_v1), alpha = 0.35,
              colour=NA) + 
  geom_line(data=blue.predicted2, aes(x = cent_sd,
                                y = pred, fill=cent_v1,col=cent_v1), size=1) +
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
        strip.text = element_text(face = "bold",size=14,family="Helvetica"),
        panel.spacing = unit(0.5,"lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.4))+
  #ylim(0,3.2)+
  ggtitle("A) Initial visit duration")+
  scale_fill_brewer(palette="Set1")+
  scale_colour_brewer(palette="Set1")+
  labs(fill="Initial visit duration",col="Initial visit duration")

blue.duration.priority.plot

ggsave(blue.duration.priority.plot,file="graphs/blueberry_priority_effect.pdf", height = 6, width = 8,dpi = 300)


################################
##species composition models
################################

#check distribution families
b.sc.mod1 <- glmmTMB(Fresh.wgt~SPEC.COM*cent_d+(1|Year/Block),
                     family="gaussian",
                     data = berry)

b.sc.mod2 <- glmmTMB(log.wgt~SPEC.COM*cent_d+(1|Year/Block),
                      family="gaussian",
                      data = berry)

b.sc.mod3 <- glmmTMB(Fresh.wgt~SPEC.COM*cent_d+(1|Year/Block),
                     family=Gamma(link="log"),
                     data = berry)


AIC(b.sc.mod1,b.sc.mod2,b.sc.mod3)

#gaussian best AIC

#check residuals
b.sc.res1=simulateResiduals(b.sc.mod1)
plot(b.sc.res1)
testResiduals(b.sc.res1)#good

b.sc.res2=simulateResiduals(b.sc.mod2)
plot(b.sc.res2)
testResiduals(b.sc.res2)#good

b.sc.res3=simulateResiduals(b.sc.mod3)
plot(b.sc.res3)
testResiduals(b.sc.res3)#good

#all fine so go with Gaussian

#Check collinearity
performance::check_collinearity(b.sc.mod1)

#Percentage differences between 10 and 100 seconds of visitation

(predict(b.sc.mod1,newdata=list(SPEC.COM=c("HB","SB","MX"),
                               cent_d=c(scale(100,attr(berry$cent_d, 'scaled:center'),
                                              attr(berry$cent_d, 'scaled:scale')),
                                        scale(100,attr(berry$cent_d, 'scaled:center'),
                                              attr(berry$cent_d, 'scaled:scale')),
                                        scale(100,attr(berry$cent_d, 'scaled:center'),
                                              attr(berry$cent_d, 'scaled:scale'))),
                               Year=c(NA,NA,NA),
                               Block=c(NA,NA,NA)),re.form = NA)-

predict(b.sc.mod1,newdata=list(SPEC.COM=c("HB","SB","MX"),
                               cent_d=c(scale(10,attr(berry$cent_d, 'scaled:center'),
                                              attr(berry$cent_d, 'scaled:scale')),
                                        scale(10,attr(berry$cent_d, 'scaled:center'),
                                              attr(berry$cent_d, 'scaled:scale')),
                                        scale(10,attr(berry$cent_d, 'scaled:center'),
                                              attr(berry$cent_d, 'scaled:scale'))),
                               Year=c(NA,NA,NA),
                               Block=c(NA,NA,NA)),re.form = NA))/
  predict(b.sc.mod1,newdata=list(SPEC.COM=c("HB","SB","MX"),
                                 cent_d=c(scale(100,attr(berry$cent_d, 'scaled:center'),
                                                attr(berry$cent_d, 'scaled:scale')),
                                          scale(100,attr(berry$cent_d, 'scaled:center'),
                                                attr(berry$cent_d, 'scaled:scale')),
                                          scale(100,attr(berry$cent_d, 'scaled:center'),
                                                attr(berry$cent_d, 'scaled:scale'))),
                                 Year=c(NA,NA,NA),
                                 Block=c(NA,NA,NA)),re.form = NA)

#check trends in intercept/slope among species
emtrends(b.sc.mod1, pairwise~SPEC.COM,var="cent_d") #no differences between taxa
bb.emm <- emtrends(b.sc.mod1, pairwise~SPEC.COM,var="cent_d")#no differences between taxa
test(bb.emm, null = 0,adjust="none")


#extract estimates for plotting etc
comp.wght.pred <- emmip(b.sc.mod1, SPEC.COM ~ cent_d, mult.name = "SPEC.COM", cov.reduce = FALSE, CIs = T, plotit = TRUE)
comp.wght.pred <- comp.wght.pred$data

#remove model estimates for more than 5 visits for SB and HB 
comp.wght.pred <- comp.wght.pred[!(comp.wght.pred$yvar > 2.6),] %>% droplevels()

#rescale
comp.wght.pred$cent_d <- comp.wght.pred$cent_d*attr(berry$cent_d, 'scaled:scale') +
  attr(berry$cent_d, 'scaled:center')

#
comp.wght.pred <- comp.wght.pred[!(comp.wght.pred$SPEC.COM%in%"SB" & comp.wght.pred$cent_d > 400),] %>% droplevels()
comp.wght.pred <- comp.wght.pred[!(comp.wght.pred$SPEC.COM%in%"HB" & comp.wght.pred$cent_d > 100),] %>% droplevels()


berry_plot2  <- berry
colnames(comp.wght.pred)[1] = "Taxa"
colnames(berry_plot2)[13] = "Taxa"
comp.wght.pred$Taxa <- revalue(comp.wght.pred$Taxa,c("HB" = "Honeybee",
                                                     "SB" = "Stingless bee",
                                                     "MX" = "Mixed"))

comp.wght.pred$Taxa <- factor(comp.wght.pred$Taxa,levels=c("Honeybee","Stingless bee","Mixed"))

berry_plot2$Taxa <- revalue(berry_plot2$Taxa,c("HB" = "Honeybee",
                                                     "SB" = "Stingless bee",
                                                     "MX" = "Mixed"))

berry_plot2$Taxa <- factor(berry_plot2$Taxa,levels=c("Honeybee","Stingless bee","Mixed"))

#back-transform duration (2*sd(berry$duration))

#plot the data
blue_comp <- ggplot() + 
     xlab("Total visit duration (s)") + ylab("Fruit weight (g)") + 
     geom_point(data = berry_plot2, aes(x = duration,
                                           y = Fresh.wgt, colour=Taxa),
                alpha=0.5,size=3,
                shape = 16,show.legend = TRUE) +
  geom_ribbon(data=comp.wght.pred,
              aes(ymin=LCL, ymax=UCL, x=cent_d, fill=Taxa), alpha = 0.5,
              show.legend = TRUE)  +
  geom_line(data=comp.wght.pred, aes(x=cent_d,y=yvar, colour=Taxa), 
            size=1,show.legend = TRUE) +
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
        axis.title.x = element_text(size=14),
        strip.text = element_blank(),
        panel.spacing = unit(0.25,"lines"),
        legend.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.4))+
       # plot.margin = unit(c(0.01,0.01,0.01,0.01), "lines"))
     scale_fill_manual(values=c("#e41a1c","#377eb8","#984ea3"))+
     scale_colour_manual(values=c("#e41a1c","#377eb8","#984ea3"))+
     #scale_fill_brewer(palette="Set1") + 
     #scale_colour_brewer(palette="Set1") + 
     labs(fill="Pollinator taxa",col="Pollinator taxa")
     #scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) + 
     facet_wrap(~Taxa,ncol=3,scales = "free_x")

blue_comp
ggsave(p,file="graphs/Fig2A Blueberry.pdf", height = 5, width = 10, dpi = 300)

#arrange plots for publication

#priority effects plots
prior.plots=align_plots(blue.duration.priority.plot,blue.identity.plot,align="hv", axis="tblr")
prior.plots.a <- ggdraw(prior.plots[[1]])
prior.plots.b <- ggdraw(prior.plots[[2]])
priority.plots <- ggarrange(prior.plots.a,prior.plots.b,nrow=2)

#species composiion plots (ESM)
comp.fig=align_plots(blue_comp,rasp_comp,align="hv", axis="tblr") #rasp_comp from other code file
comp.fig.a <- ggdraw(comp.fig[[1]])
comp.fig.b <- ggdraw(comp.fig[[2]])
comp.arrange <- ggarrange(comp.fig.a,comp.fig.b,nrow=2,common.legend = TRUE)

ggsave(priority.plots,file="graphs/Figure 2 -  priority effects.pdf", height = 7, width = 5, dpi = 300)
ggsave(comp.arrange,file="graphs/Figure S1 - speciescomposition.pdf", height = 7, width = 5,dpi=300)
