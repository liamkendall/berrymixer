##RASP

library(plyr)
library(dplyr)
library(DHARMa)
library(glmmTMB)
library(ggplot2)

#left function
left = function (string,char){
  substr(string,1,char)
}

#load data
rasp <- read.csv("data/raspberry.csv", header=T)

rasp$MONTH
#remove zeros
rasp=rasp[!rasp$Weight==0,]
rasp=rasp[,-which(names(rasp) %in% c("Comments","COMMENTS","POLLEN_DONOR"))]

#subset dataframe to just visitation treatment
#convert all visitor columns to character

rasp_updated <- rasp[!rasp$POLLINATORS%in%"BC",]
rasp_updated <- rasp_updated[!rasp_updated$POLLINATORS%in%"OC",]%>%droplevels()

rasp_updated <- rasp_updated %>% mutate_at(vars(16:44), as.character) %>% droplevels()

#retain only first letter of the string for visitors IDs
rasp_updated2 <- rasp_updated %>% mutate_at(vars(16:44), function (x) left(x,char=1))

#sum the number of visits from each taxa
rasp_updated2$stingless_bee <- rowSums(rasp_updated2[, c(16:44)] == "S")
rasp_updated2$honey_bee <- rowSums(rasp_updated2[, c(16:44)] == "H")

#sum the total number of visits
rasp_updated2$sumvisits <- rowSums(rasp_updated2[, c(47:48)])

#calculate the percent visits from each taxa
rasp_updated2$p_stingless_bee <- rasp_updated2$stingless_bee/rasp_updated2$sumvisits
rasp_updated2$p_honey_bee <- rasp_updated2$honey_bee/rasp_updated2$sumvisits

rasp_updated2 <- rasp_updated2[!rasp_updated2$sumvisits >20,]%>%droplevels()
rasp_updated2 <- rasp_updated2[!(rasp_updated2$POLLINATORS%in%"SB" & rasp_updated2$sumvisits > 10),] %>% droplevels()
rasp_updated2 <- rasp_updated2[!(rasp_updated2$POLLINATORS%in%"HB" & rasp_updated2$sumvisits > 10),] %>% droplevels()
rasp_updated2$RP <- paste0(rasp_updated2$ROW, rasp_updated2$PLANT)

##Subset to 2 or more visits and only mixed visits
rasp_updated3=rasp_updated2[rasp_updated2$POLLINATORS%in%"Mixed"
                              & rasp_updated2$sumvisits >1,]%>%droplevels()
rasp_updated3 <- rasp_updated3[!rasp_updated3$sumvisits >15,]%>%droplevels()

#count number of reps per number of visits for each taxa group
table(rasp_updated3$VISIT1, rasp_updated3$sumvisits)

#run priority effect model
rasp_m1.d <- glmmTMB(Weight~VISIT1*sumvisits*p_honey_bee+(1|BLOCK),
                   family="gaussian",
                   data = rasp_updated3)
summary(rasp_m1.d)
        
#dredge it and print csv
dredge.rasp_m1 <- dredge(rasp_m1.d) # model excluding ratio is the best
write.csv(dredge.rasp_m1, "/Users/macuser/Library/Mobile Documents/com~apple~CloudDocs/H_drive_DT/berrymixer/dredge.raspberryfruit.weight.csv")

#run priority effect model
rasp_m1 <- glmmTMB(Weight~VISIT1*sumvisits+(1|BLOCK),
              family="gaussian",
              data = rasp_updated3)

summary(rasp_m1)
emtrends(rasp_m1, pairwise~VISIT1,var="sumvisits")
#contrast   estimate         SE  df t.ratio p.value
#H - S    0.01052346 0.05573458 114   0.189  0.8506

rasp_m1res=simulateResiduals(rasp_m1)
plot(rasp_m1res)

#extract estimates etc
ras.wght.pred <- emmip(rasp_m1, VISIT1 ~ sumvisits, mult.name = "VISIT1", cov.reduce = FALSE, CIs = T, plotit = TRUE)
ras.wght.pred <- ras.wght.pred$data
ras.wght.pred$xvar <- as.numeric(as.character(ras.wght.pred$xvar))

#plot the data
p <- ggplot()
p <- p + xlab("Number of visits") + ylab("Raspberry fruit weight (g)")
p <- p + geom_ribbon(data=ras.wght.pred,
                     aes(ymin=LCL, ymax=UCL, x=xvar, fill=VISIT1), alpha = 0.5)
p <- p + geom_line(data=ras.wght.pred, aes(xvar,yvar, colour=VISIT1), size=1)
p <- p + scale_x_continuous(breaks=seq(2,15,2))
p <- p + scale_y_continuous(breaks=seq(0,6,1))
p <- p + geom_jitter(data = rasp_updated3,aes(x = sumvisits, y = Weight, colour=VISIT1), 
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
        strip.text = element_text(size=12))
p <- p + theme(axis.title.y=element_text(margin=margin(0,20,0,0)))
p <- p + theme(panel.border = element_rect(color = "black", fill = NA, size = 0.4))
p <- p + scale_fill_brewer(palette="Set2")
p <- p + scale_colour_brewer(palette="Set2")
p

#look at number of reps per number of visits
table(rasp_updated2$POLLINATORS, rasp_updated2$sumvisits)

#run the model
rasp_m2 <- glmmTMB(Weight~sumvisits*POLLINATORS+(1|BLOCK/RP),
            family="gaussian",
            data = rasp_updated2)

summary(rasp_m2)
rasp.emm <- emtrends(rasp_m2, pairwise~POLLINATORS,var="sumvisits")#no differences between taxa
test(rasp.emm, null = 0, side = ">")

rasp_m2res=simulateResiduals(rasp_m2)
plot(rasp_m2res)

#extract estimates etc
ras.comp.wght <- emmip(rasp_m2, POLLINATORS ~ sumvisits, mult.name = "POLLINATORS", cov.reduce = FALSE, CIs = T, plotit = TRUE)
ras.comp.wght <- ras.comp.wght$data
ras.comp.wght$xvar <- as.numeric(as.character(ras.comp.wght$xvar))
ras.comp.wght <- ras.comp.wght[!(ras.comp.wght$POLLINATORS%in%"SB" & ras.comp.wght$sumvisits > 10),] %>% droplevels()
ras.comp.wght <- ras.comp.wght[!(ras.comp.wght$POLLINATORS%in%"HB" & ras.comp.wght$sumvisits > 10),] %>% droplevels()

#plot the data
p <- ggplot()
p <- p + xlab("Number of visits") + ylab("Raspberry fruit weight (g)")
p <- p + geom_ribbon(data=ras.comp.wght,
                     aes(ymin=LCL, ymax=UCL, x=xvar, fill=POLLINATORS), alpha = 0.5)
p <- p + geom_line(data=ras.comp.wght, aes(xvar,yvar, colour=POLLINATORS), size=1)
p <- p + scale_x_continuous(breaks=seq(2,20,2))
p <- p + scale_y_continuous(breaks=seq(0,6,1))
p <- p + geom_jitter(data = rasp_updated2, aes(x = sumvisits, y = Weight, colour=POLLINATORS), size=2.5, shape = 21, width=0, height =0.04)
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
#jamie finished here (12/2/19)
################################