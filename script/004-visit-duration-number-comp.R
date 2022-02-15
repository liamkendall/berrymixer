####relationship between speceis, visit number and visit duration

#COR.TEST
cor.test(berry.sp.split[[1]]$duration,berry.sp.split[[1]]$sumvisits,
         method="spearman")
cor.test(berry.sp.split[[2]]$duration,berry.sp.split[[2]]$sumvisits,
         method="spearman")
cor.test(berry.sp.split[[3]]$duration,berry.sp.split[[3]]$sumvisits,
         method="spearman")

###Model comparisons
AIC(hb.mod)-AIC(hb.sv.mod)
AIC(sb.mod)-AIC(sb.sv.mod)
AIC(mx.mod)-AIC(mx.sv.mod)

summary(hb.sv.mod)
summary(sb.sv.mod)
summary(mx.sv.mod)
#########
##priority effects sumvisits models
#########

#####
#null model - with only total visit duration
#####
berry_pe$subvisits <- berry_pe$sumvisits-1

#visit duration model
bsv.null1.mod <- glmmTMB(Fresh.wgt~sumvisits+
                            (1|YRP),REML=FALSE,
                          family="gaussian",
                          data = berry_pe)
#test residuals
plot(simulateResiduals(bsv.null1.mod)) #good


######
#priority effect model 1
#visitor identity * duration
######
bsv.sp.mod <- glmmTMB(Fresh.wgt~Visitor1*subvisits+
                         (1|YRP),REML=FALSE,
                       family=gaussian,
                       data = berry_pe)

bsv.sp.res <- simulateResiduals(bsv.sp.mod)
plot(bsv.sp.res)#good-ish  - top quantile

bsv.sp.mod2 <- glmmTMB(Fresh.wgt~Visitor1*subvisits+
                          (1|YRP),REML=FALSE,
                        family=gaussian,
                        data = berry_pe)

bsv.sp.res2 <- simulateResiduals(bsv.sp.mod2)
plot(bsv.sp.res2) #good

######
#priority effect model 2
#inital duration * remaining duration
######

bsv.dr.mod <- glmmTMB(Fresh.wgt~cent_v1*subvisits+
                         (1|YRP),REML=FALSE,
                       family=gaussian,
                       data = berry_pe)

plot(simulateResiduals(bsv.dr.mod)) #good-ish - top quantile again

#VIF
performance::check_collinearity(bsv.sp.mod) #good
performance::check_collinearity(bsv.dr.mod) #good

#AIC
AICc(bsv.sp.mod,bsv.dr.mod,bsv.null1.mod)

#species model test
anova(bsv.null1.mod,bsv.sp.mod)
AIC(bsv.null1.mod)-AIC(bsv.sp.mod)

#duration model test
anova(bsv.null1.mod,bsv.dr.mod)
AIC(bsv.null1.mod)-AIC(bsv.dr.mod)

#difference between species identity and visit duration
AIC(bsv.sp.mod)-AIC(bsv.dr.mod)

AIC(bsv.dr.mod)
AIC(bsv.sp.mod)
AIC(bsv.null1.mod)

AIC(bsv.sp.mod)-AIC(blue.sp.mod)
AIC(bsv.dr.mod)-AIC(blue.dr.mod)

#Re-fit with REML

#visitor identity * duration
bsv.sp.mod.reml <- glmmTMB(Fresh.wgt~Visitor1*subvisits+
                              (1|YRP),REML=TRUE,
                            family=gaussian,
                            data = berry_pe)

plot(simulateResiduals(bsv.sp.mod.reml))

#inital duration * total duration
bsv.dr.mod.reml <- glmmTMB(Fresh.wgt~cent_v1*subvisits+
                              (1|YRP),REML=TRUE,
                            family=gaussian,
                            data = berry_pe)

plot(simulateResiduals(bsv.dr.mod.reml))

summary(bsv.sp.mod.reml)
summary(bsv.dr.mod.reml)

AIC(bsv.dr.mod.reml)
AIC(bsv.sp.mod.reml)

#test slopes dependent on first visitor
test(emtrends(bsv.sp.mod.reml, pairwise~Visitor1,var="subvisits"))

######BELOW####
####for description of results
####calculate percent difference at two visits
####average from two visits 

#average duration of 1 minute of visitation

#average first visit
mean(berry_pe[,c("V1D")])

#60 seconds minus average of first visit

sv_low <-60-mean(berry_pe[,c("V1D")])

#480 seconds (8 minutes) minus average of first visit

sv_high <- 480-mean(berry_pe[,c("V1D")])

#2 visits and ten visits

#With stingless bee - low visitation
sv_sb_low=bsv.sp.mod.reml$fit$parfull[1]+bsv.sp.mod.reml$fit$parfull[2]+
  bsv.sp.mod.reml$fit$parfull[3]*1+bsv.sp.mod.reml$fit$parfull[4]*1

sv_hb_low=bsv.sp.mod.reml$fit$parfull[1]+bsv.sp.mod.reml$fit$parfull[3]*1

((sb_low-hb_low)/hb_low)*100

#59.81272%

sv_sb_high=bsv.sp.mod.reml$fit$parfull[1]+bsv.sp.mod.reml$fit$parfull[2]+
  bsv.sp.mod.reml$fit$parfull[3]*10+bsv.sp.mod.reml$fit$parfull[4]*10

sv_hb_high=bsv.sp.mod.reml$fit$parfull[1]+bsv.sp.mod.reml$fit$parfull[3]*10

((sv_hb_high-sv_sb_high)/sv_sb_high)*100

#-24.27155%

sv_sb_high2=bsv.sp.mod.reml$fit$parfull[1]+bsv.sp.mod.reml$fit$parfull[2]+
  bsv.sp.mod.reml$fit$parfull[3]*15+bsv.sp.mod.reml$fit$parfull[4]*15

sv_hb_high2=bsv.sp.mod.reml$fit$parfull[1]+bsv.sp.mod.reml$fit$parfull[3]*15

((sv_hb_high2-sv_sb_high2)/sv_sb_high2)*100
#####plot trends

#extract and format model predictions
wght.sv.pred <- emmip(bsv.sp.mod.reml, Visitor1 ~ subvisits, mult.name = "Visitor1",
                   cov.reduce = FALSE, CIs = T, plotit = TRUE)

wght.sv.pred <- wght.sv.pred$data
wght.sv.pred$xvar <- as.numeric(as.character(wght.sv.pred$xvar))
#wght.pred$xvar <- wght.pred$xvar * attr(berry$subvisits, 'scaled:scale') + attr(berry$subvisits, 'scaled:center')
#wght.pred$subvisits <- wght.pred$subvisits* attr(berry$subvisits, 'scaled:scale') + attr(berry$subvisits, 'scaled:center')
wght.sv.pred$Visitor1 <- revalue(wght.sv.pred$Visitor1,c("H" = "Honeybee",
                                                   "S" = "Stingless bee"))

#plot dataframe
berry_plot <- berry_pe

berry_plot$Visitor1 <- revalue(berry_plot$Visitor1,c("H" = "Honeybee",
                                                     "S" = "Stingless bee"))

bsv.identity.plot <- ggplot() + 
  xlab("Number of visits") + 
  ylab("Blueberry fruit weight (g)") + 
  geom_point(data = berry_plot,aes(x = subvisits+1,
                                   y = Fresh.wgt,
                                   colour=Visitor1,
                                   fill=Visitor1), alpha=0.5,col="white",
             size=4, shape = 21) + 
  geom_ribbon(data=wght.sv.pred,
              aes(ymin=LCL, ymax=UCL, x=xvar+1, fill=Visitor1), alpha = 0.35) + 
  geom_line(data=wght.sv.pred, aes(x = subvisits+1,
                                y = yvar, colour=Visitor1), size=1) + 
  
  theme_bw()+ 
  scale_x_continuous(breaks=seq(1,15,2))+
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
sv.prepplot <- as.data.frame(matrix(ncol = 3, nrow = 6000))
colnames(sv.prepplot) <- c("cent_v1", "subvisits", "pred")

sv.prepplot$cent_v1 <- rep(seq(10,300,20), 100)
sv.prepplot <- sv.prepplot[order(sv.prepplot$cent_v1),]
sv.prepplot$subvisits <- rep(seq(1,15,1), 100)
sv.prepplot$YRP <- NA
sv.prepplot$subvisits <- sv.prepplot$subvisits-1
sv.prepplot$cent_v1 <- scale(sv.prepplot$cent_v1 , 
                          attr(berry$cent_v1, "scaled:center"), 
                          attr(berry$cent_v1, "scaled:scale"))
sv.prepplot <- sv.prepplot%>%
  filter(!subvisits==0)

sv.prepplot$pred <- predict(bsv.dr.mod.reml,
                         newdata = sv.prepplot,
                         re.form = NA)

sv.prepplot$cent_v1 <-  sv.prepplot$cent_v1* attr(berry$cent_v1, 'scaled:scale') +
  attr(berry$cent_v1, 'scaled:center')

sv.prepplot$pred <- ifelse(sv.prepplot$pred<0,0,sv.prepplot$pred)

sv.prepplot2 <- sv.prepplot%>%
  filter(cent_v1<range(berry$V1D,na.rm = T)[2])%>%
  filter(!pred<0)

mid<-mean(berry_pe$Fresh.wgt,na.rm = T)

sv.raster.surface <- ggplot(sv.prepplot2, aes(y = cent_v1, 
                                        x = subvisits+1,
                                        z=pred,
                                        fill = pred)) + 
  geom_raster() +
  stat_contour(col="black",
               size=0.4,
               bins=8) +
  theme_bw() +
  #scale_x_continuous()+
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
  xlab("Number of visits") +
  scale_fill_gradient2(midpoint=1.5,
                       low="#2C3B75",
                       mid="white",
                       high="#B8321A",
                       space ="Lab",
                       guide="colourbar")+
  scale_x_continuous(breaks=seq(1,15,2),expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  ggtitle("(A) Initial visit duration")+
  labs(fill="Fruit weight (g)")

sv.raster.surface
#arrange plots for publication

#priority effects plots
sv.prior.plots.a <- ggplotGrob(sv.raster.surface)
sv.prior.plots.b <- ggplotGrob(bsv.identity.plot)

sv.priority.plots <- rbind(sv.prior.plots.a,sv.prior.plots.b)

ggsave(sv.priority.plots,file="graphs/Fig S2.jpg",
       device="jpg",
       height = 7, 
       width = 5,
       dpi = 600)

ggsave(sv.priority.plots,
       file="graphs/Fig S2.pdf",
       height = 7,
       width = 5, 
       dpi = 600)
