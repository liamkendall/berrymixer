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

#remove Tasmania BB data and RE data 
#due to very low fresh weight sample size (10 Mixed visits)
table(berry_updated2[berry_updated2$Species%in%"RE" &
        is.na(berry_updated2$Fresh.wgt)==FALSE,]$SPEC.COM)

berry_updated3 <- berry_updated2[!berry_updated2$Species%in%"BR",]%>%droplevels()
berry_updated3 <- berry_updated3[!berry_updated3$Species%in%"RE",]%>%droplevels()

##Subset to 2 or more visits and only mixed visits
berry_updated4=berry_updated3[berry_updated3$SPEC.COM%in%"MX"
                              & berry_updated3$sumvisits >1,]%>%droplevels()
berry_updated4$log.wgt <- log(berry_updated4$Fresh.wgt)

#################################
#model priority effects
#################################

#too many 1's in the data to run the fruitset model
m1 <- glmmTMB(FS~Visitor1*sumvisits+(1|Year),
                        family="binomial",
                        data = berry_updated4)
summary(m1)
plot_model(m1, type = "pred",terms = "Visitor1")
plot_model(m1, type = "pred",terms = "sumvisits")

#run full model for fruit weight
m2 <- glmmTMB(log.wgt~Visitor1*sumvisits*p_honey_bee+(1|Year),
              family=gaussian,
              data = berry_updated4)
summary(m2)
#check residuals
m2res=simulateResiduals(m2)
plot(m2res)
testResiduals(m2res)#looks good

#dredge it and print csv
d.m2 <- dredge(m2) # model excluding ratio is the best
write.csv(d.m2, "/Users/macuser/Library/Mobile Documents/com~apple~CloudDocs/H_drive_DT/berrymixer/dredge.fruit.weight.csv")

#run best model as determined by dregde above 
m3 <- glmmTMB(log.wgt~Visitor1*sumvisits+(1|Year),#log fruit weight gives lower AIC
              family=gaussian,
              data = berry_updated4)
summary(m3)

#check residuals
m3res=simulateResiduals(m3)
plot(m3res)
testResiduals(m3res)
#looks good

#quickely look at estimates
sjPlot::plot_model(m3, type = "pred",terms = "Visitor1")
sjPlot::plot_model(m3, type = "pred",terms = "sumvisits")

#compare slopes
emtrends(m3, pairwise~Visitor1,var="sumvisits")
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
p <- p + xlab("Number of visits") + ylab("Fruit weight (g)")
p <- p + theme(text = element_text(size=18))
p <- p + geom_ribbon(data=wght.pred,
                     aes(ymin=LCL, ymax=UCL, x=xvar, fill=Visitor1), alpha = 0.5)
p <- p + geom_line(data=wght.pred, aes(xvar,yvar, colour=Visitor1), size=1)
p <- p + scale_x_continuous(breaks=seq(2,15,2))
p <- p + scale_y_continuous(breaks=seq(0.4,2.6,0.4))
p <- p + geom_jitter(data = berry_updated4,aes(x = sumvisits, y = Fresh.wgt, colour=Visitor1), size=1.5, shape = 21, width=0, height =0.08)
p <- p + theme(axis.line.x = element_blank(),
                   axis.line.y = element_blank(),
                   panel.grid.major = element_line(size=.4, colour = "#d3d3d3"),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank()) +
  theme(axis.text.x=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =10),
        axis.title.x=element_text(size=16, vjust = 1),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =10),
        axis.title.y=element_text(size=16, vjust = 1),
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
#model species composition
################################

#log berry weight is better model
berry_updated3$log.wgt <- log(berry_updated3$Fresh.wgt)
m5 <- glmmTMB(log.wgt~SPEC.COM*sumvisits+(1|Block/Row)+(1|Year),
            family="gaussian",
            data = berry_updated3)

emtrends(m5, pairwise~SPEC.COM,var="sumvisits")#no differences between taxa

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
p <- p + xlab("Number of visits") + ylab("Fruit weight (g)")
p <- p + theme(text = element_text(size=18))
p <- p + geom_ribbon(data=comp.wght.pred,
                     aes(ymin=LCL, ymax=UCL, x=xvar, fill=SPEC.COM), alpha = 0.5)
p <- p + geom_line(data=comp.wght.pred, aes(xvar,yvar, colour=SPEC.COM), size=1)
p <- p + scale_x_continuous(breaks=seq(1,15,2))
p <- p + scale_y_continuous(breaks=seq(0.22,2.6,0.4))
p <- p + geom_jitter(data = berry_updated3, aes(x = sumvisits, y = Fresh.wgt, colour=SPEC.COM), size=1.5, shape = 21, width=0, height =0.08)
p <- p + theme(axis.line.x = element_blank(),
               axis.line.y = element_blank(),
               panel.grid.major = element_line(size=.4, colour = "#d3d3d3"),
               panel.grid.minor = element_blank(),
               panel.background = element_blank()) +
  theme(axis.text.x=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =10),
        axis.title.x=element_text(size=16, vjust = 1),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =10),
        axis.title.y=element_text(size=16, vjust = 1),
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

















###LIAM STOPPED HERE - SPECIES COMP MODELS IN NEW SCRIPT


#create new dataframe
pred <- expand.grid(SPEC.COM=unique(berry_updated3$SPEC.COM),
                    TVN = unique(berry_updated3$TVN))
pred <- na.omit(pred)

pred$Fresh.wgt = predict(m4, newdata=pred, type="response", re.form = NA)
mm <- model.matrix(terms(m4), pred)

## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(m4),mm))
tvar1 <- pvar1+VarCorr(m1m)$plant_species[1]  ## must be adapted for more complex models
tvar2 <- tvar1+VarCorr(m1m)$transect_number[1]
cmult <- 1.96 ## could use 1.96
links.pred <- data.frame(
  links.pred
  , plo = links.pred$number_fruits_consumed-cmult*sqrt(pvar1)
  , phi = links.pred$number_fruits_consumed+cmult*sqrt(pvar1)
  , tlo = links.pred$number_fruits_consumed-cmult*sqrt(tvar2)
  , thi = links.pred$number_fruits_consumed+cmult*sqrt(tvar2)
)


#model cases where HD:SB ratios are 0:1 or 1:1
#no cases where >1 visits include only singless bess
keep <- c(0, 0.5, 1)
equal_ratio <- berry_updated3[berry_updated3$p_stingless_bee %in% keep, ]
equal_ratio$p_honey_bee <- as.factor(as.character(equal_ratio$p_honey_bee))
keep <- c(2,4,6,8)
equal_ratio <- equal_ratio[equal_ratio$TVN %in% keep, ]
equal_ratio <- na.omit(equal_ratio)

m3 <- glmmTMB(FS~p_honey_bee+(1|Block/RP)+(1|Year),
              family="binomial",
              data = equal_ratio)

m4 <- glmmTMB(Fresh.wgt~p_honey_bee+(1|Block/RP)+(1|Year),
              family="gaussian",
              data = equal_ratio)

mx <- glmer(Fresh.wgt~SPEC.COM*TVN+(1|Block/RP)+(1|Year),
              family="gaussian",
              data = equal_ratio)

#create new dataframe
pred.eq <- expand.grid(SPEC.COM=unique(equal_ratio$SPEC.COM),
                    TVN = unique(equal_ratio$TVN))
pred.eq <- na.omit(pred.eq)

pred.eq$Fresh.wgt = predict(mx, newdata=pred.eq, type="response", re.form = NA)

#plot the data
p <- ggplot()
p <- p + xlab("Metric value") + ylab("Number of fruits consumed")
p <- p + theme(text = element_text(size=18))
#p <- p + geom_ribbon(data=pred.links.comp,
#                     aes(ymin=tlo, ymax=thi, x=value), alpha = 0.5, fill = "#CC6666")
#p <- p + geom_ribbon(data=pred.links.comp,
#                     aes(ymin=plo, ymax=phi, x=value), alpha = 0.7, fill = "#9999CC")
p <- p + geom_line(data=pred.eq, aes(TVN,Fresh.wgt, colour=SPEC.COM), size=0.8, linetype="dashed")
#p <- p + geom_rug(data=rug.data,aes(x=value, y=NULL))
#p <- p + scale_y_continuous(breaks=seq(2,10,2))
#p <- p + geom_jitter(data = data.points,aes(x = value, y = number_fruits_consumed), width = 0.005)
#p <- p + facet_wrap(~metric, scales="free_x")
p <- p + theme(panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black", size=.5)) +
  theme(panel.border=element_rect(colour = "black", fill = "NA", size = .5)) +
  theme(axis.text.x=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =18),
        axis.title.x=element_text(size=30, vjust = 1),
        axis.text.y=element_text(angle= 360, hjust = 0.5, vjust = 0.5, size =18),
        axis.title.y=element_text(size=30, vjust = 1),
        axis.text=element_text(colour = "black"))+
  theme(strip.background = element_rect(colour="NA", fill=NA),
        strip.text = element_text(size=20))+
  theme(axis.ticks.length = unit(2, "mm"))
p <- p + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
p

#######################################
#END
#######################################