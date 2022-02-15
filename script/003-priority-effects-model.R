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
berry <- read.csv("data/blue_priority.csv")

#####
#SET-UP
#####

#scale visitation variables
berry$cent_v1 <- scale(berry$V1D,center=TRUE,scale=TRUE)
berry$cent_sd  <- scale((berry$duration-berry$V1D),center=TRUE,scale=TRUE)

#log berry weight if needed
berry$log.wgt <- log(berry$Fresh.wgt)

##Subset to 2 or more visits and only mixed visits
berry_pe <- berry[berry$SPEC.COM%in%"MX"&berry$sumvisits >1,]%>%droplevels()
berry_pe <- berry_pe[!is.na(berry_pe$Fresh.wgt)==TRUE,] #remove NAs for fruit weight

#PLANT ID
berry_pe$RP <- as.factor(berry_pe$RP)

#########
##priority effects models
#########

#V1D = Inital visitor duration
#Visitor1 = Inital visitor identity
#sumvisits = number of total visits
#scale_v1 = scaled first visit duration
#scale_sd = scaled remaining duration

#####
#null model - with only total visit duration
#####
#Plant as random as across two years, mixed visits are all from one block, but variance very low in year so better

berry_pe$YRP <- paste(berry_pe$Year,berry_pe$RP)

#visit duration model
blue.null1.mod <- glmmTMB(Fresh.wgt~cent_d+
                            (1|YRP),REML=FALSE,
                          family="gaussian",
                          data = berry_pe)
#test residuals
plot(simulateResiduals(blue.null1.mod)) #good


######
#priority effect model 1
#visitor identity * duration
######
blue.sp.mod <- glmmTMB(Fresh.wgt~Visitor1*cent_sd+
                         (1|YRP),REML=FALSE,
                       family=gaussian,
                       data = berry_pe)

blue.sp.res <- simulateResiduals(blue.sp.mod)
plot(blue.sp.res)#good-ish  - top quantile

blue.sp.mod2 <- glmmTMB(Fresh.wgt~Visitor1*cent_sd+
                         (1|YRP),REML=FALSE,
                       family=gaussian,
                       data = berry_pe)

blue.sp.res2 <- simulateResiduals(blue.sp.mod2)
plot(blue.sp.res2) #good

######
#priority effect model 2
#inital duration * remaining duration
######

blue.dr.mod <- glmmTMB(Fresh.wgt~cent_v1*cent_sd+
                         (1|YRP),REML=FALSE,
                       family=gaussian,
                       data = berry_pe)

plot(simulateResiduals(blue.dr.mod)) #good-ish - top quantile again

#VIF
performance::check_collinearity(blue.sp.mod) #good
performance::check_collinearity(blue.dr.mod) #good

#AIC
AICc(blue.sp.mod,blue.dr.mod,blue.null1.mod)

#species model test
anova(blue.null1.mod,blue.sp.mod)
AIC(blue.null1.mod)-AIC(blue.sp.mod)

#duration model test
anova(blue.null1.mod,blue.dr.mod)
AIC(blue.null1.mod)-AIC(blue.dr.mod)

#difference between species identity and visit duration
AIC(blue.sp.mod)-AIC(blue.dr.mod)

#Re-fit with REML

#visitor identity * duration
blue.sp.mod.reml <- glmmTMB(Fresh.wgt~Visitor1*cent_sd+
                         (1|YRP),REML=TRUE,
                       family=gaussian,
                       data = berry_pe)

plot(simulateResiduals(blue.sp.mod.reml))

#inital duration * total duration
blue.dr.mod.reml <- glmmTMB(Fresh.wgt~cent_v1*cent_sd+
                         (1|YRP),REML=TRUE,
                       family=gaussian,
                       data = berry_pe)

plot(simulateResiduals(blue.dr.mod.reml))

summary(blue.sp.mod.reml)
summary(blue.dr.mod.reml)

#test slopes dependent on first visitor
test(emtrends(blue.sp.mod.reml, pairwise~Visitor1,var="cent_sd"))

######BELOW####
####for description of results
####calculate percent difference at two visits
####average from two visits 

#average duration of 1 minute of visitation

#average first visit
mean(berry_pe[,c("V1D")])
     
#60 seconds minus average of first visit

scale_low <- scale((60-mean(berry_pe[,c("V1D")])),
      attr(berry$cent_sd, "scaled:center"), attr(berry$cent_sd, "scaled:scale"))

#480 seconds (8 minutes) minus average of first visit

scale_high <- scale((480-mean(berry_pe[,c("V1D")])),
                    attr(berry$cent_sd, "scaled:center"), attr(berry$cent_sd, "scaled:scale"))

#With stingless bee - low visitation
sb_low=blue.sp.mod.reml$fit$parfull[1]+blue.sp.mod.reml$fit$parfull[2]+
         blue.sp.mod.reml$fit$parfull[3]*scale_low[1]+blue.sp.mod.reml$fit$parfull[4]*scale_low[1]

hb_low=blue.sp.mod.reml$fit$parfull[1]+blue.sp.mod.reml$fit$parfull[3]*scale_low[1]

((sb_low-hb_low)/hb_low)*100

#59.81272%

sb_high=blue.sp.mod.reml$fit$parfull[1]+blue.sp.mod.reml$fit$parfull[2]+
            blue.sp.mod.reml$fit$parfull[3]*scale_high[1]+blue.sp.mod.reml$fit$parfull[4]*scale_high[1]

hb_high=blue.sp.mod.reml$fit$parfull[1]+blue.sp.mod.reml$fit$parfull[3]*scale_high[1]

((hb_high-sb_high)/sb_high)*100

#-24.27155%

#####plot trends

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
  xlab("Total visit duration (s)") + 
  ylab("Blueberry fruit weight (g)") + 
  geom_point(data = berry_plot,aes(x = cent_sd * attr(berry$cent_sd, 'scaled:scale') + attr(berry$cent_sd, 'scaled:center'),
                                    y = Fresh.wgt,
                                    colour=Visitor1,
                                    fill=Visitor1), alpha=0.5,col="white",
              size=4, shape = 21) + 
  geom_ribbon(data=wght.pred,
              aes(ymin=LCL, ymax=UCL, x=xvar, fill=Visitor1), alpha = 0.35) + 
  geom_line(data=wght.pred, aes(x = cent_sd,
                                y = yvar, colour=Visitor1),show.legend = F, size=1) + 

  theme_bw()+ 
  #ylim(0,3.2)+
  #xlim(600,657)+
  ggtitle("(B) Initial species identity")+
  theme(plot.title = element_text(hjust = 0,face="bold",size=14),
               legend.box.background = element_rect(colour = "black"),
               legend.title = element_text(size=12),
               legend.text = element_text(size=10),
               aspect.ratio = 1,
               strip.background = element_blank(),
               text=element_text(),
               axis.ticks.length = unit(0,"mm"),
               legend.title.align=0.5,
               axis.text.y = element_text(size=10),
               axis.text.x = element_text(size=10),
               axis.title.y = element_text(size=12),
               axis.title.x = element_text(size=12),
               strip.text = element_text(face = "bold",size=14,family="Helvetica"),
               panel.spacing = unit(0.5,"lines"),
               panel.border = element_rect(color = "black", fill = NA, size = 0.4))+
   scale_fill_manual(values = c("#2C3B75","#B8321A"))+
   scale_colour_manual(values = c("#2C3B75","#B8321A"))+
   labs(fill="Initial visitor",col="Initial visitor")

###Plot - initial visit duration - predictions
prepplot <- as.data.frame(matrix(ncol = 3, nrow = 6000))
colnames(prepplot) <- c("cent_v1", "cent_sd", "pred")

prepplot$cent_v1 <- rep(seq(10,600,10), 100)
prepplot <- prepplot[order(prepplot$cent_v1),]
prepplot$cent_sd <- rep(seq(10,600,10), 100)
#prepplot$RP <- NA
prepplot$YRP <- NA

prepplot$cent_sd <- scale(prepplot$cent_sd, 
                          attr(berry$cent_sd, "scaled:center"), 
                          attr(berry$cent_sd, "scaled:scale"))

prepplot$cent_v1 <- scale(prepplot$cent_v1 , 
                          attr(berry$cent_v1, "scaled:center"), 
                          attr(berry$cent_v1, "scaled:scale"))

prepplot$pred <- predict(blue.dr.mod.reml,
                         newdata = prepplot,
                         re.form = NA)

prepplot$cent_v1 <-  prepplot$cent_v1* attr(berry$cent_v1, 'scaled:scale') +
  attr(berry$cent_v1, 'scaled:center')
prepplot$cent_sd <- prepplot$cent_sd* attr(berry$cent_sd, 'scaled:scale') +
  attr(berry$cent_sd, 'scaled:center')

prepplot$pred <- ifelse(prepplot$pred<0,0,prepplot$pred)

prepplot2 <- prepplot%>%
  # filter(pred>0)%>%
  filter(cent_v1<range(berry$V1D,na.rm = T)[2])%>%
  filter(cent_sd<range((berry$duration-berry$V1D),na.rm = T)[2])

mid<-mean(berry_pe$Fresh.wgt,na.rm = T)

raster.surface <- ggplot(prepplot2, aes(y = cent_v1, 
                                        x = cent_sd,
                                        z=pred,
                                        fill = pred)) + 
  geom_raster() +
  stat_contour(col="black",
               size=0.4,
               bins=8) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0,face="bold",size=14),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        aspect.ratio = 1,
        strip.background = element_blank(),
        text=element_text(),
        axis.ticks.length = unit(0,"mm"),
        legend.title.align=0.5,
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        strip.text = element_text(face = "bold",size=14,family="Helvetica"),
        panel.spacing = unit(0.5,"lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.4))+
  ylab("Initial visit duration (s)") + 
  xlab("Total visit duration (s)") +
  scale_fill_gradient2(midpoint=1.5,
                       low="#2C3B75",
                       mid="white",
                       high="#B8321A",
                       space ="Lab",
                       guide="colourbar")+
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  ggtitle("(A) Initial visit duration")+
  labs(fill="Fruit weight (g)")

raster.surface
#arrange plots for publication

#priority effects plots
prior.plots.a <- ggplotGrob(raster.surface)
prior.plots.b <- ggplotGrob(blue.identity.plot)
priority.plots <- rbind(prior.plots.a,prior.plots.b)

ggsave(priority.plots,file="graphs/Fig 3.jpg",
       device="jpg",
       height = 7, 
       width = 5,
       dpi = 600)


ggsave(priority.plots,
       file="graphs/Fig 3.pdf",
       height = 7,
       width = 5, 
       dpi = 600)
