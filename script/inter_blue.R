#load packages
library(stringr)
library(plyr)
library(dplyr)
library(glmmTMB)
library(lme4)
library(ggplot2)

#load data
berry <- read.csv("data/BLUEBERRY_MIXED_V3.csv", header=T)

#subset dataframe to just visitation treatment
#convert all visitor columns to character
berry_updated <- berry[berry$TREATMENT%in%"V" & berry$TVN > 1,] %>% mutate_at(vars(19:33), as.character) %>% droplevels()

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

#remove Tasmania BB data
berry_updated3 <- berry_updated2[!berry_updated2$Species%in%"BR",]%>%droplevels()

#priority effects models
m1 <- glmmTMB(FS~Visitor1*sumvisits+(1|Block)+(1|Year),
                        family="binomial",
                        data = berry_updated3)
summary(m1)



m2 <- glmmTMB(Fresh.wgt~Visitor1*sumvisits+(1|Block),
              family="gaussian",
              data = berry_updated3[berry_updated3$SPEC.COM%in%"MX",])
summary(m2)

plot_model(m2, type = "pred",terms = "Visitor1")

emmeans::emtrends(m2, pairwise~Visitor1,var="sumvisits")

emmeans::emmip(m2, Visitor1 ~ sumvisits,mult.name = "Visitor1", cov.reduce = FALSE)


m3 <- glmmTMB(FS~SPEC.COM*sumvisits+(1|Block)+(1|Year),
              family="binomial",
              data = berry_updated3)
summary(m3)

m4 <- glmmTMB(Fresh.wgt~SPEC.COM*sumvisits+(1|Block)+(1|Year),
              family="gaussian",
              data = berry_updated3)
summary(m4)
emmeans::emtrends(m4, pairwise~SPEC.COM,var="sumvisits")
table(berry_updated3$SPEC.COM)
emmeans::emmip(m4, SPEC.COM ~ sumvisits,mult.name = "SPEC.COM", cov.reduce = FALSE)

#compare fruit set models
AIC(m1,m3)

#compare fruit weight models
AIC(m2,m4)

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

#plot the data
p <- ggplot()
p <- p + xlab("Metric value") + ylab("Number of fruits consumed")
p <- p + theme(text = element_text(size=18))
#p <- p + geom_ribbon(data=pred.links.comp,
#                     aes(ymin=tlo, ymax=thi, x=value), alpha = 0.5, fill = "#CC6666")
#p <- p + geom_ribbon(data=pred.links.comp,
#                     aes(ymin=plo, ymax=phi, x=value), alpha = 0.7, fill = "#9999CC")
p <- p + geom_line(data=pred, aes(TVN,Fresh.wgt, colour=SPEC.COM), size=0.8, linetype="dashed")
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