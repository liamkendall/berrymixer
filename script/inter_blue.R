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
berry_updated2$duration=berry_updated_time$duration
berry_updated2$V1D=berry_updated_time$V1duration

berry_updated2[is.na(berry_updated2$V1D)==TRUE,]


berry_updated2 <- berry_updated2[!(berry_updated2$SPEC.COM%in%"HB" & berry_updated2$sumvisits ==15),]%>%droplevels()

#remove duplicated records
berry_updated2[duplicated(berry_updated2$Tag),]
berry_updated2 <- berry_updated2[!duplicated(berry_updated2$Tag),]%>%droplevels()

#remove Tasmania BB data and RE data 
#due to very low fresh weight sample size (10 Mixed visits)
table(berry_updated2[berry_updated2$Species%in%"RE" &
        is.na(berry_updated2$Fresh.wgt)==FALSE,]$SPEC.COM)
berry_updated3 <- berry_updated2[!berry_updated2$Species%in%"BR",]%>%droplevels()
berry_updated3 <- berry_updated3[!berry_updated3$Species%in%"RE",]%>%droplevels()

##Subset to 2 or more visits and only mixed visits
berry_updated4=berry_updated3[berry_updated3$SPEC.COM%in%"MX"
                              & berry_updated3$sumvisits >1,]%>%droplevels()

#calculate log of berry weight
berry_updated4$log.wgt <- log(berry_updated4$Fresh.wgt)

#add missing V1 Duration
berry_updated4[berry_updated4$Tag%in%"Juan-11",c("V1D")]=58
berry_updated4$V1D <- as.numeric(berry_updated4$V1D)
berry_updated4$scale_v1 <- scale(berry_updated4$V1D)
berry_updated4$scale_d <- scale(berry_updated4$duration)
berry_updated4$scale_sd <- scale(berry_updated4$duration-berry_updated4$V1D)

#################################
#priority effects----
#################################

#too many 1's in the data to run a fruit set model
m1 <- glmmTMB(FS~Visitor1*sumvisits+(1|Year),
                        family="binomial",
                        data = berry_updated3)

summary(m1)

plot_model(m1, type = "pred",terms = "Visitor1")
plot_model(m1, type = "pred",terms = "sumvisits")
plot_model(m1,type="int")

#priority effects for fruit weight
#################################

#V1D = Inital visitor duration
#Visitor1 = Inital visitor identity
#sumvisits = number of total visits
#scale_v1 = scaled first visit duration

#visitor identity
blue.sp.mod <- glmmTMB(log.wgt~Visitor1*sumvisits+
                (1|Year),
              family=gaussian,
              data = berry_updated4)

#duration
blue.dr.mod <- glmmTMB(log.wgt~scale_v1*sumvisits+
                (1|Year),
              family=gaussian,
              data = berry_updated4)

#both
blue.both.mod <- glmmTMB(log.wgt~scale_v1*sumvisits*Visitor1+
                         (1|Year),
                       family=gaussian,
                       data = berry_updated4)

#full duration model
blue.cm.dr.mod <- glmmTMB(log.wgt~scale_v1*scale_d*Visitor1+
                           (1|Year),
                         family=gaussian,
                         data = berry_updated4)

#full duration model
blue.cm.dr.mod <- glmmTMB(log.wgt~scale_v1*scale_d*Visitor1+
                             (1|Year),
                           family=gaussian,
                           data = berry_updated4)

#full duration but total duration = total - v1
blue.sb.dr.mod <- glmmTMB(log.wgt~scale_v1*scale_sd*Visitor1+
                            (1|Year),
                          family=gaussian,
                          data = berry_updated4)

#AIC.BIC
AICc(blue.sp.mod,blue.dr.mod,blue.cm.dr.mod,blue.sb.dr.mod)



#check residuals
#species identity model
blue.sp.mod.res=simulateResiduals(blue.sp.mod)
plot(blue.sp.mod.res)
testResiduals(blue.sp.mod.res)#looks good

#visit duration model
blue.dr.mod.res=simulateResiduals(blue.dr.mod)
plot(blue.dr.mod.res)
testResiduals(blue.dr.mod.res)#looks good


#visit full duration model
blue.cm.dr.mod.res=simulateResiduals(blue.cm.dr.mod)
plot(blue.cm.dr.mod.res)
testResiduals(blue.cm.dr.mod.res)#looks good

#dredge it and print csv
blue.sp.d <- dredge(blue.sp.mod,rank="AICc") 
blue.dr.d <- dredge(blue.dr.mod,rank="AICc") 
blue.both.d <- dredge(blue.both.mod,rank="AICc") 
blue.cm.dr.d <- dredge(blue.cm.dr.mod,rank="AICc") 
blue.sb.dr.d <- dredge(blue.sb.dr.mod,rank="AICc") 

write.csv(blue.sp.d, "dredge_blue_species.csv")
write.csv(blue.dr.d, "dredge_blue_duration.csv")

b.sp.mods <-get.models(blue.sp.d,subset=TRUE)
b.dr.mods <-get.models(blue.dr.d,subset=TRUE)
b.both.mods <-get.models(blue.both.d,subset=TRUE)
b.cm.dr.mods <-get.models(blue.cm.dr.d,subset=TRUE)
b.sb.dr.mods <-get.models(blue.sb.dr.d,subset=TRUE)

#summarise best priority effects model
summary(b.sp.mods[[1]])
summary(b.dr.mods[[1]])
summary(b.both.mods[[1]])
summary(b.cm.dr.mods[[1]])
summary(b.sb.dr.mods[[1]])

importance(b.sp.mods)
importance(b.dr.mods)
importance(b.both.mods)
importance(b.cm.dr.mods)

AICc(b.cm.dr.mods[[1]],b.sb.dr.mods[[1]]) #almost identical

#check residuals of best duration or identity models

sp.fit.res=simulateResiduals(b.sp.mods[[1]])
plot(sp.fit.res)
testResiduals(sp.fit.res)#looks good

dr.fit.res=simulateResiduals(b.dr.mods[[1]])
plot(dr.fit.res)
testResiduals(dr.fit.res)#looks good


AIC(blue.sp.d)


confint(b.dr.mods[[5]])
confint(b.sp.mods[[1]])
confint(b.sp.mods[[2]])
confint(b.sp.mods[[3]])
confint(b.sp.mods[[4]])
confint(b.sp.mods[[5]])


#
#
#
#
#
#
#calculate percent difference at two visits

#stingless bee initial with two visits
sb_2=exp(d.mods[[1]]$fit$par[1]+d.mods[[1]]$fit$par[2]*2+
  d.mods[[1]]$fit$par[3]+d.mods[[1]]$fit$par[4]*2)

hb_2=exp(d.mods[[1]]$fit$par[1]+d.mods[[1]]$fit$par[2]*2)

((sb_2-hb_2)/hb_2)*100
#84.72958%

sb_15=exp(d.mods[[1]]$fit$par[1]+d.mods[[1]]$fit$par[2]*15+
          d.mods[[1]]$fit$par[3]+d.mods[[1]]$fit$par[4]*15)

hb_15=exp(d.mods[[1]]$fit$par[1]+d.mods[[1]]$fit$par[2]*15)

((sb_15-hb_15)/hb_15)*100
#-30.18957%

#compare slopes
slopes.m3 <- emtrends(d.mods[[1]], pairwise~Visitor1,var="sumvisits")
test(slopes.m3)
sjPlot::plot_model(blue.m3, type = "int")

#$emtrends
#Visitor1 sumvisits.trend     SE df t.ratio p.value
#H                 0.0507 0.0174 45  2.923  0.0054 
#S                -0.0241 0.0238 45 -1.011  0.3172 


#$contrasts
#contrast estimate     SE df t.ratio p.value
#H - S      0.0749 0.0295 45 2.538   0.0147

#extract estimates etc
wght.pred <- emmip(d.mods[[1]], Visitor1 ~ sumvisits, mult.name = "Visitor1", cov.reduce = FALSE, CIs = T, plotit = TRUE)
wght.pred <- wght.pred$data
wght.pred$yvar <- exp(wght.pred$yvar)
wght.pred$UCL <- exp(wght.pred$UCL)
wght.pred$LCL <- exp(wght.pred$LCL)
wght.pred$xvar <- as.numeric(as.character(wght.pred$xvar))

#plot dataframe
berry_updated6 <- berry_updated4
#colnames(berry_updated6)[19] = "Initial visitor"
#colnames(wght.pred)[1] = "Initial visitor"

berry_updated6$Visitor1 <- revalue(berry_updated6$Visitor1,c("H" = "Honeybee",
                                                             "S" = "Stingless bee"))

wght.pred$Visitor1 <- revalue(wght.pred$Visitor1,c("H" = "Honeybee",
                                                             "S" = "Stingless bee"))
#plot the data
p <- ggplot()
p <- p + xlab("Number of visits") + ylab("Blueberry fruit weight (g)")
p <- p + geom_ribbon(data=wght.pred,
                     aes(ymin=LCL, ymax=UCL, x=xvar, fill=Visitor1), alpha = 0.35)
p <- p + geom_line(data=wght.pred, aes(xvar,yvar, colour=Visitor1), size=1)
p <- p + scale_x_continuous(breaks=seq(2,15,2))
p <- p + scale_y_continuous(breaks=seq(0.4,2.6,0.4))
p <- p + geom_jitter(data = berry_updated6,aes(x = sumvisits, y = Fresh.wgt, colour=Visitor1), 
                     size=2.5, shape = 21, width=0, height =0.04)
p <- p + theme_bw()
p <- p + theme(plot.title = element_text(hjust = 0.5),
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
               panel.spacing = unit(0.5,"lines"))
p <- p + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.4))
p <- p + scale_fill_brewer(palette="Set1")
p <- p + scale_colour_brewer(palette="Set1")
p <- p + labs(fill="Initial visitor",col="Initial visitor")
p
ggsave(p,file="graphs/blueberry_priority_effect.pdf", height =6, width=8,dpi=300)

################################
#species composition----
################################

#log berry weight is better model
berry_updated3$log.wgt <- log(berry_updated3$Fresh.wgt)
m5 <- glmmTMB(log.wgt~SPEC.COM*sumvisits+(1|Block)+(1|Year),
            family="gaussian",
            data = berry_updated3)

summary(m5)
emtrends(m5, pairwise~SPEC.COM,var="sumvisits") #no differences between taxa
bb.emm <- emtrends(m5, pairwise~SPEC.COM,var="sumvisits")#no differences between taxa
#test(bb.emm, null = 0, side = ">")
test(bb.emm, null = 0,adjust="fdr")

#check residuals
m5res=simulateResiduals(m5)
plot(m5res)
testResiduals(m5res)#good



#extract estimates etc
comp.wght.pred <- emmip(m5, SPEC.COM ~ sumvisits, mult.name = "SPEC.COM", cov.reduce = FALSE, CIs = T, plotit = TRUE)
comp.wght.pred <- comp.wght.pred$data
comp.wght.pred$yvar <- exp(comp.wght.pred$yvar)
comp.wght.pred$UCL <- exp(comp.wght.pred$UCL)
comp.wght.pred$LCL <- exp(comp.wght.pred$LCL)
comp.wght.pred$SE <- exp(comp.wght.pred$SE)
comp.wght.pred$xvar <- as.numeric(as.character(comp.wght.pred$xvar))

#remove model estimates for more than 5 visits for SB and HB 
comp.wght.pred <- comp.wght.pred[!(comp.wght.pred$SPEC.COM%in%"SB" & comp.wght.pred$sumvisits > 5),] %>% droplevels()
comp.wght.pred <- comp.wght.pred[!(comp.wght.pred$SPEC.COM%in%"HB" & comp.wght.pred$sumvisits > 5),] %>% droplevels()
keep <- c("HB", "MX", "SB")
berry_updated3 <- berry_updated3[berry_updated3$SPEC.COM %in% keep,]
berry_updated5  <- berry_updated3[!(berry_updated3$SPEC.COM%in%"HB" & berry_updated3$sumvisits>5),]

colnames(comp.wght.pred)[1] = "Taxa"
colnames(berry_updated5)[12] = "Taxa"
comp.wght.pred$Taxa <- revalue(comp.wght.pred$Taxa,c("HB" = "Honeybee",
                                                     "SB" = "Stingless bee",
                                                     "MX" = "Mixed"))

comp.wght.pred$Taxa <- factor(comp.wght.pred$Taxa,levels=c("Honeybee","Stingless bee","Mixed"))

berry_updated5$Taxa <- revalue(berry_updated5$Taxa,c("HB" = "Honeybee",
                                                     "SB" = "Stingless bee",
                                                     "MX" = "Mixed"))

berry_updated5$Taxa <- factor(berry_updated5$Taxa,levels=c("Honeybee","Stingless bee","Mixed"))

#plot the data
#might be worth plotting just 1-5 visits
p <- ggplot() + 
     xlab("Number of visits") + ylab("Fruit weight (g)") + 
     geom_ribbon(data=comp.wght.pred,
              aes(ymin=LCL, ymax=UCL, x=xvar, fill=Taxa), alpha = 0.5,
                  show.legend = TRUE) + 
     geom_line(data=comp.wght.pred, aes(xvar,yvar, colour=Taxa), size=1,show.legend = TRUE) + 
     geom_jitter(data = berry_updated5, aes(x = sumvisits, y = Fresh.wgt, colour=Taxa),
                 size=2.5, shape = 21, width=0, height =0.04,show.legend = TRUE) + 
     theme_bw() + 
     ggtitle("Blueberry")+
     theme(plot.title = element_text(vjust = 0,hjust = 0,face = "bold",size=16,family="Helvetica"),
          legend.box.background = element_rect(colour = "black"),
          aspect.ratio = 1,
          strip.background = element_blank(),
          text=element_text(),
          axis.ticks.length = unit(0,"mm"),
          legend.title.align=0.5,
          axis.text.y = element_text(size=11),
          axis.text.x = element_text(size=11),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=14),
          strip.text = element_blank(),
          panel.spacing = unit(0.25,"lines"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.4),
          plot.margin = unit(c(0, 0.5, 0, 0), "cm"))+
     scale_fill_brewer(palette="Set1") + 
     scale_colour_brewer(palette="Set1") + 
     scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) + 
     facet_wrap(~Taxa,scales = "free_x")
p
ggsave(p,file="graphs/blueberry_composition.pdf", height =5, width=10,dpi=300)


library(ggpubr)


comp.fig=align_plots(p,rasp_comp,align="hv", axis="tblr")
comp.fig.a <- ggdraw(comp.fig[[1]])
comp.fig.b <- ggdraw(comp.fig[[2]])

comp.arrange <- ggarrange(comp.fig.a,comp.fig.b,nrow=2,common.legend = TRUE)
comp.arrange
ggsave(comp.arrange,file="graphs/arranged_comp.pdf", height = 7, width = 10,dpi=300)

################################
#liam finished here (24/6/19)
################################



