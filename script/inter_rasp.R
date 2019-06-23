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
rasp <- read.csv("data/raspberry.csv", header=T)

#remove zeros
rasp <- rasp[!rasp$Weight==0,]
rasp <- rasp[,-which(names(rasp) %in% c("Comments","COMMENTS","POLLEN_DONOR"))]

#subset dataframe to just visitation treatment
#convert all visitor columns to character
rasp_updated <- rasp[!rasp$POLLINATORS%in%"BC",]
rasp_updated <- rasp_updated[!rasp_updated$POLLINATORS%in%"OC",]%>%droplevels()
rasp_updated <- rasp_updated %>% mutate_at(vars(16:44), as.character) %>% droplevels()

#retain only first letter of the string for visitors IDs
#left function
left = function (string,char){
  substr(string,1,char)
}
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

rasp_updated2 <- rasp_updated %>% mutate_at(vars(16:44), function (x) left(x,char=1))
rasp_updated_time <- rasp_updated %>% mutate_at(vars(16:44), function (x) stri_sub(x,3))
rasp_updated_time$duration<-rowSums(apply(rasp_updated_time[,16:44],
                                          MARGIN=2,FUN = unlist(as.numeric)),na.rm=TRUE)
rasp_updated_time$V1duration<-as.numeric(rasp_updated_time$VISIT1)

#sum the number of visits from each taxa
rasp_updated2$stingless_bee <- rowSums(rasp_updated2[, c(16:44)] == "S")
rasp_updated2$honey_bee <- rowSums(rasp_updated2[, c(16:44)] == "H")

#sum the total number of visits
rasp_updated2$sumvisits <- rowSums(rasp_updated2[, c(47:48)])

#calculate the percent visits from each taxa
rasp_updated2$p_stingless_bee <- rasp_updated2$stingless_bee/rasp_updated2$sumvisits
rasp_updated2$p_honey_bee <- rasp_updated2$honey_bee/rasp_updated2$sumvisits

#add duration variables (duration code)
rasp_updated2$duration=rasp_updated_time$duration
rasp_updated2$V1D=rasp_updated_time$V1duration

#remove sum visits greater than 20
rasp_updated2 <- rasp_updated2[!rasp_updated2$sumvisits >20,]%>%droplevels()

#remove sb or hb visits only greater than 10
rasp_updated2 <- rasp_updated2[!(rasp_updated2$POLLINATORS%in%"SB" & rasp_updated2$sumvisits > 10),] %>% droplevels()
rasp_updated2 <- rasp_updated2[!(rasp_updated2$POLLINATORS%in%"HB" & rasp_updated2$sumvisits > 10),] %>% droplevels()

##Subset to 2 or more visits and only mixed visits
rasp_updated3=rasp_updated2[rasp_updated2$POLLINATORS%in%"Mixed"
                              & rasp_updated2$sumvisits >1,]%>%droplevels()

#remove mixed visits greater than 15 for priority effects model
rasp_updated3 <- rasp_updated3[!rasp_updated3$sumvisits >15,]%>%droplevels()

#count number of reps per number of visits for each taxa group
table(rasp_updated3$VISIT1, rasp_updated3$sumvisits)

#Scaled Visitor 1 duration
rasp_updated3$scale_v1=scale(rasp_updated3$V1D)

rasp_updated3$log.wgt <- log(rasp_updated3$Weight)
hist(rasp_updated3$Weight)

#####################################
#priority effects----
#####################################

#run full priority effect model
rasp_m1.d <- glmmTMB(Weight~scale_v1*VISIT1*sumvisits+(1|BLOCK),REML = FALSE,
                     family="gaussian",
                     data = rasp_updated3)

rasp_m2.d <- glmmTMB(log.wgt~scale_v1*VISIT1+sumvisits+(1|BLOCK),REML = FALSE,
                     family="gaussian",
                     data = rasp_updated3)

#compare AIC of raw/logged response
AIC(rasp_m1.d,rasp_m2.d)
#df      AIC
#rasp_m1.d 10 361.8604 #raw
#rasp_m2.d 10 197.5314 #log

#check residuals
rasp_m1res=simulateResiduals(rasp_m1.d)
plot(rasp_m1res)

rasp_m2res=simulateResiduals(rasp_m2.d)
plot(rasp_m2res)

#dredge it and print csv
dredge.rasp_m1 <- dredge(rasp_m1.d,rank=AICc,beta = "partial.sd") #
dredge.rasp_m2 <- dredge(rasp_m2.d,rank=AICc) #

head(dredge.rasp_m1)
head(dredge.rasp_m2)

write.csv(dredge.rasp_m1, "dredge.raspberryfruit.weight.csv")

#get models from dredge
rasp.mods.1 <- get.models(dredge.rasp_m1,subset=TRUE)
rasp.mods.2 <- get.models(dredge.rasp_m2,subset=TRUE)

#run priority effect model for reporting # top model
summary(rasp.mods.1[[1]])
summary(rasp.mods.2[[1]])

rasp.delta.2 <- model.avg(dredge.rasp_m1, subset = delta < 2)
summary(rasp.delta.2)

importance(dredge.rasp_m2) 

ma_df <- model.avg(dredge.rasp_m1, subset = delta < 2)
ma_coefs <- coefTable(ma_df, full=TRUE, adjust.se = TRUE)
coefnames <- row.names(ma_coefs)
ma_coefs <- as.data.frame(ma_coefs)
ma_coefs <- mutate(ma_coefs,
                   coefficient = coefnames,
                   t = Estimate / `Std. Error`,
                   p = pt(abs(t), df = 195, lower.tail = FALSE),
                   lower95 = Estimate - 1.96 * `Std. Error`,
                   upper95 = Estimate + 1.96 * `Std. Error`) %>% select(coefficient, Estimate, `Std. Error`, t, p, lower95, upper95)
knitr::kable(ma_coefs, digits = 2)


rasp.avg.95p <- model.avg(dredge.rasp_m2, cumsum(weight) <= .95)

summary(model.avg(dredge.rasp_m1, subset = delta < 5))

summary(rasp.delta.5)


#First visit duration
#plot the data
#extract estimates etc
ras.wght.pred <- emmip(dredge.rasp_mods[[1]], ~scale_v1, cov.reduce = FALSE, CIs = T, plotit = TRUE)

ras.wght.pred <- ras.wght.pred$data

ras.wght.pred$xvar <- as.numeric(as.character(ras.wght.pred$xvar))


p <- ggplot()
p <- p + xlab("First visit duration") + ylab("Raspberry fruit weight (g)")
p <- p + geom_ribbon(data=ras.wght.pred,
                     aes(ymin=LCL, ymax=UCL, x=xvar), alpha = 0.5)
p <- p + geom_line(data=ras.wght.pred, aes(xvar,yvar), size=1)
#p <- p + scale_x_continuous(breaks=seq(2,15,2))
#p <- p + scale_y_continuous(breaks=seq(0,6,1))
p <- p + geom_jitter(data = rasp_updated3,aes(x = scale_v1, y = Weight),
                     size=2, shape = 21, width=0, height =0.04)
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
p <- p + scale_fill_brewer(palette="Set1")
p <- p + scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
p <- p + scale_colour_brewer(palette="Set1")
p


#No. of visits
#plot the data
#extract estimates etc
ras.wght.pred2 <- emmip(dredge.rasp_mods[[1]], ~sumvisits, cov.reduce = FALSE, CIs = T, plotit = TRUE)

ras.wght.pred2 <- ras.wght.pred2$data

ras.wght.pred2$xvar <- as.numeric(as.character(ras.wght.pred2$xvar))

p <- ggplot()
p <- p + xlab("No. of visits") + ylab("Raspberry fruit weight (g)")
p <- p + geom_ribbon(data=ras.wght.pred2,
                     aes(ymin=LCL, ymax=UCL, x=xvar), alpha = 0.5)
p <- p + geom_line(data=ras.wght.pred2, aes(xvar,yvar), size=1)
#p <- p + scale_x_continuous(breaks=seq(2,15,2))
p <- p + scale_y_continuous(breaks=seq(0,6,1))
p <- p + geom_jitter(data = rasp_updated3,aes(x = sumvisits, y = Weight),
                     size=2, shape = 21, width=0, height =0.04)
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
p <- p + scale_fill_brewer(palette="Set1")
p <- p + scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
p <- p + scale_colour_brewer(palette="Set1")
p






#Species identity model for comparison
rasp.modID <- glmmTMB(Weight~(sumvisits*VISIT1)+(1|BLOCK),
                      family="gaussian",
                      data = rasp_updated3)
summary(rasp.modID)

#######################################
#species composition----
#######################################

#look at number of reps per number of visits
table(rasp_updated2$POLLINATORS, rasp_updated2$sumvisits)

#run the model
rasp_m2 <- glmmTMB(Weight~sumvisits*POLLINATORS+(1|BLOCK),
            family="gaussian",
            data = rasp_updated2)

summary(rasp_m2)

#test that slopes for each taxa are different from 0
rasp.emm <- emtrends(rasp_m2, pairwise~POLLINATORS,var="sumvisits")
test(rasp.emm, null = 0)

#check residuals
rasp_m2res=simulateResiduals(rasp_m2)
plot(rasp_m2res)

#extract estimates etc
ras.comp.wght <- emmip(rasp_m2, POLLINATORS ~ sumvisits, mult.name = "POLLINATORS", cov.reduce = FALSE, CIs = T, plotit = TRUE)
ras.comp.wght <- ras.comp.wght$data
ras.comp.wght$xvar <- as.numeric(as.character(ras.comp.wght$xvar))
ras.comp.wght <- ras.comp.wght[!(ras.comp.wght$POLLINATORS%in%"SB" & ras.comp.wght$sumvisits > 10),] %>% droplevels()
ras.comp.wght <- ras.comp.wght[!(ras.comp.wght$POLLINATORS%in%"HB" & ras.comp.wght$sumvisits > 10),] %>% droplevels()

#plot the data
colnames(ras.comp.wght)[1] = "Taxa"
colnames(rasp_updated2)[15] = "Taxa"
ras.comp.wght$Taxa <- revalue(ras.comp.wght$Taxa,c("HB" = "Honeybee",
                                                     "SB" = "Stingless bee",
                                                     "Mixed" = "Mixed"))

ras.comp.wght$Taxa <- factor(ras.comp.wght$Taxa,levels=c("Honeybee","Stingless bee","Mixed"))

rasp_updated2$Taxa <- revalue(rasp_updated2$Taxa,c("HB" = "Honeybee",
                                                     "SB" = "Stingless bee",
                                                     "Mixed" = "Mixed"))

rasp_updated2$Taxa <- factor(rasp_updated2$Taxa,levels=c("Honeybee","Stingless bee","Mixed"))

rasp_comp <- ggplot() + 
  xlab("Number of visits") + 
  ylab("Fruit weight (g)")+ 
  geom_ribbon(data=ras.comp.wght,
              aes(ymin=LCL, ymax=UCL, x=xvar, fill=Taxa), 
              alpha = 0.5,show.legend = FALSE) +
  geom_line(data=ras.comp.wght, aes(xvar,yvar, colour=Taxa), size=1,show.legend = FALSE) +
  scale_y_continuous(breaks=seq(0,6,1)) + 
  geom_jitter(data = rasp_updated2, 
              aes(x = sumvisits, y = Weight, colour=Taxa), 
              size=2.5, shape = 21, width=0, height =0.04,show.legend = FALSE) +
  theme_bw() + 
  ggtitle("Raspberry")+
  theme(plot.title = element_text(hjust = 0,face = "bold",size=16,family="Helvetica"),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text(),
        aspect.ratio = 1,
        strip.background = element_blank(),
        text=element_text(),
        axis.ticks.length = unit(0,"mm"),
        legend.title.align=0.5,
        axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=11),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.text = element_blank(),
        panel.spacing = unit(0.25,"lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.4)) +
  scale_fill_brewer(palette="Set1") +
  scale_colour_brewer(palette="Set1") +
  scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))+
  facet_wrap(~Taxa, scales = "free_x")
rasp_comp
#
ggsave(p1,file="graphs/raspberry_composition.pdf", height = 5, width = 10,dpi=300)


################################
#jamie finished here (12/2/19)
################################

##RASP EXCLUSION TEST FOR SANITY

rasp.exc <- rbind.fill(rasp[rasp$POLLINATORS%in%"BC",],rasp[rasp$TREATMENT%in%"V",])%>%droplevels()
table(rasp.exc$TREATMENT)

#run full priority effect model
rasp.exc.mod <- glmmTMB(Weight~TREATMENT+(1|BLOCK),
                     family="gaussian",
                     data = rasp.exc)

summary(rasp.exc.mod)

r.squaredGLMM(rasp.exc.mod)
