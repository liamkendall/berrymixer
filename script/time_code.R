#duration models
library(glmmTMB)

#priority effects for fruit weight
dur1 <- glmmTMB(log.wgt~(scale(V1D)*Visitor1)+(sumvisits*Visitor1)+
                  (1|Year),
              family=gaussian,
              data = berry_updated4)
sum(is.na(berry_updated4$V1D))
summary(dur1)
dur1 <- glmmTMB(log.wgt~scale(V1D)*sumvisits+
                  (1|Year),
                family=gaussian,
                data = berry_updated4)
DUR_DREDGE=MuMIn::dredge(dur1)
summary(dur1)


#Full duration raspberry
dur_rasp1 <- glmmTMB(Weight~scale(V1D)+sumvisits+VISIT1+
                  VISIT1:sumvisits+VISIT1:scale(V1D)+
                  (1|BLOCK),
                family=gaussian,
                data = rasp_updated2)
summary(dur_rasp1)
dredge(dur_rasp1)


##t-tests for species duration differences 
t.test(V1D~Visitor1,berry_updated4)
((81.628-16.16)/16.16)*100

t.test(V1D~VISIT1,rasp_updated2)

model.avg(DUR_DREDGE,cumsum(weight) <= .95)

boxplot()

confset.95p <- get.models(DUR_DREDGE, cumsum(weight) <= .95)
model.avg(confset.95p)
  

