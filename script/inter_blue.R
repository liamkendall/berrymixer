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
berry <- read.csv("data/BLUEBERRY_MIXED_V3.csv", header=T)

#subset dataframe to just visitation treatment
#convert all visitor columns to character
berry_updated <- berry[berry$TREATMENT%in%"V",] %>% mutate_at(vars(19:33), as.character) %>% droplevels()

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
berry_updated2$RP=paste0(berry_updated2$Row,berry_updated2$Plant.number)

#add duration variables (duration code)
berry_updated2$duration=berry_updated_time$duration
berry_updated2$V1D=berry_updated_time$V1duration

berry_updated2 <- berry_updated2[!(berry_updated2$SPEC.COM%in%"HB" & berry_updated2$sumvisits ==15),]%>%droplevels()

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

#################################
#priority effects----
#################################

#too many 1's in the data to run the fruitset model
m1 <- glmmTMB(FS~Visitor1*sumvisits+(1|Year),
                        family="binomial",
                        data = berry_updated4)
summary(m1)
plot_model(m1, type = "pred",terms = "Visitor1")
plot_model(m1, type = "pred",terms = "sumvisits")

#priority effects for fruit weight
m2 <- glmmTMB(log.wgt~Visitor1*sumvisits*p_honey_bee+(1|Year),
              family=gaussian,
              data = berry_updated4)
summary(m2)
emtrends(m2, pairwise~Visitor1,var="p_honey_bee")

#check residuals
m2res=simulateResiduals(m2)
plot(m2res)
testResiduals(m2res)#looks good

#dredge it and print csv
d.m2 <- dredge(m2) # model excluding ratio is the best
write.csv(d.m2, "/Users/macuser/Library/Mobile Documents/com~apple~CloudDocs/H_drive_DT/berrymixer/dredge.fruit.weight.csv")

#run best priority effects model as determined by dregde above 
m3 <- glmmTMB(log.wgt~Visitor1*sumvisits+(1|Year),#log fruit weight gives lower AIC
              family=gaussian,
              data = berry_updated4)
summary(m3)

#check residuals
m3res=simulateResiduals(m3)
plot(m3res)
testResiduals(m3res)#looks good

#compare slopes
slopes.m3 <- emtrends(m3, pairwise~Visitor1,var="sumvisits")
test(slopes.m3)
#$contrasts
#contrast  estimate       SE df t.ratio p.value
#H - S    0.0783168 0.027268 78   2.872  0.0052

#extract estimates etc
wght.pred <- emmip(m3, Visitor1 ~ sumvisits, mult.name = "Visitor1", cov.reduce = FALSE, CIs = T, plotit = TRUE)
wght.pred <- wght.pred$data
wght.pred$yvar <- exp(wght.pred$yvar)
wght.pred$UCL <- exp(wght.pred$UCL)
wght.pred$LCL <- exp(wght.pred$LCL)
wght.pred$xvar <- as.numeric(as.character(wght.pred$xvar))

#plot the data
p <- ggplot()
p <- p + xlab("Number of visits") + ylab("Blueberry fruit weight (g)")
p <- p + geom_ribbon(data=wght.pred,
                     aes(ymin=LCL, ymax=UCL, x=xvar, fill=Visitor1), alpha = 0.5)
p <- p + geom_line(data=wght.pred, aes(xvar,yvar, colour=Visitor1), size=1)
p <- p + scale_x_continuous(breaks=seq(2,15,2))
p <- p + scale_y_continuous(breaks=seq(0.4,2.6,0.4))
p <- p + geom_jitter(data = berry_updated4,aes(x = sumvisits, y = Fresh.wgt, colour=Visitor1), 
                     size=2.5, shape = 21, width=0, height =0.04)
p <- p + theme(axis.line.x = element_blank(),
                   axis.line.y = element_blank(),
                   panel.grid.major = element_line(size=.4, colour = "#d3d3d3"),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank()) +
  theme(axis.text.x=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =14),
        axis.title.x=element_text(size=20, vjust = 1),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =14),
        axis.title.y=element_text(size=20, vjust = 1),
        axis.text=element_text(colour = "black"))+
  theme(axis.ticks.length = unit(2, "mm"),
        axis.ticks = element_line(colour = 'black', size = 0.4))+
  theme(strip.background = element_rect(colour="NA", fill=NA),
        strip.text = element_text(size=16))
p <- p + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.4))
p <- p + scale_fill_brewer(palette="Set2")
p <- p + scale_colour_brewer(palette="Set2")
p

################################
#species composition----
################################

#log berry weight is better model
berry_updated3$log.wgt <- log(berry_updated3$Fresh.wgt)
m5 <- glmmTMB(log.wgt~SPEC.COM*sumvisits+(1|Block/RP)+(1|Year),
            family="gaussian",
            data = berry_updated3)
summary(m5)
emtrends(m5, pairwise~SPEC.COM,var="sumvisits")#no differences between taxa
bb.emm <- emtrends(m5, pairwise~SPEC.COM,var="sumvisits")#no differences between taxa
test(bb.emm, null = 0, side = ">")

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

#plot the data
p <- ggplot()
p <- p + xlab("Number of visits") + ylab("Blueberry fruit weight (g)")
p <- p + geom_ribbon(data=comp.wght.pred,
                     aes(ymin=LCL, ymax=UCL, x=xvar, fill=SPEC.COM), alpha = 0.5)
p <- p + geom_line(data=comp.wght.pred, aes(xvar,yvar, colour=SPEC.COM), size=1)
p <- p + scale_x_continuous(breaks=seq(1,15,2))
p <- p + scale_y_continuous(breaks=seq(0.22,2.6,0.4))
p <- p + geom_jitter(data = berry_updated3, aes(x = sumvisits, y = Fresh.wgt, colour=SPEC.COM), size=2.5, shape = 21, width=0, height =0.04)
p <- p + theme(axis.line.x = element_blank(),
               axis.line.y = element_blank(),
               panel.grid.major = element_line(size=.4, colour = "#d3d3d3"),
               panel.grid.minor = element_blank(),
               panel.background = element_blank()) +
  theme(axis.text.x=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =14),
        axis.title.x=element_text(size=20, vjust = 1),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =14),
        axis.title.y=element_text(size=20, vjust = 1),
        axis.text=element_text(colour = "black"))+
  theme(axis.ticks.length = unit(2, "mm"),
        axis.ticks = element_line(colour = 'black', size = 0.4))+
  theme(strip.background = element_rect(colour="NA", fill=NA),
        strip.text = element_text(size=12))
p <- p + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.4))
p <- p + scale_fill_brewer(palette="Set2")
p <- p + scale_colour_brewer(palette="Set2")
p

################################
#jamie finished here (11/2/19)
################################