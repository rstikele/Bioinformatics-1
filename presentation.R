

library(car)
library(ggplot2)
library(pwr)
library(dunn.test)
library(plyr)

TSA=read.csv("TSA.csv", header = T)
FMA=read.csv("FMA.csv", header = T)

as.factor(TSA$Capture.Locale)
as.factor(FMA$Capture.Locale)
captureTSA<-factor(TSA$Capture.Locale,levels=c("Estrella", "South", "Phoenix", "Tempe"))
captureFMA<-factor(FMA$Capture.Locale,levels=c("Estrella", "South", "Phoenix", "Tempe"))


sample.size.TSA= aggregate(.~captureTSA, data=TSA, length)
sample.size.FMA= aggregate(.~captureFMA, data=FMA, length)

TSAcfu=TSA$Total.CFU.ml
FMAcfu=FMA$Adj..FDB.CFU.mL
TSAaov=aov(TSAcfu~captureTSA)
FMAaov=aov(FMAcfu~captureFMA)
par(mfrow = c(2, 2)) 
plot(TSAaov)
plot(FMAaov)
dev.off()

qqPlot(FMAcfu)
qqPlot(TSAcfu)

boxplot(TSAcfu~captureTSA,xlab="Capture Locale",ylab="CFU/ml")
boxplot(FMAcfu~captureFMA,xlab="Capture Locale",ylab="CFU/ml")

TSAaov=aov(TSAcfu~captureTSA)
summary(TSAaov)
FMAaov=aov(FMAcfu~captureFMA)
summary(FMAaov)

qqPlot(TSAaov$resid,ylab="Z value")
shapiro.test(TSAaov$resid)#Not normal
qqPlot(FMAaov$resid,ylab="Z value")
shapiro.test(FMAaov$resid)#Not normal

par(mfrow = c(2, 2))  
plot(TSAaov)
plot(FMAaov)

bartlett.test(TSAcfu~captureTSA)#Unequal vairance
bartlett.test(FMAcfu~captureFMA)#Unequal Variance

kruskal.test(TSAcfu,as.factor(captureTSA))#Significant Differences among Capture locales
kruskal.test(FMAcfu,as.factor(captureFMA))#Significant Differences among Capture locales

dunn.test(TSAcfu,g=captureTSA, method = c("none"), kw=T, label = TRUE, wrap = F, alpha = 0.05)
dunn.test(FMAcfu,g=captureFMA, method = c("none"), kw=T, label = TRUE, wrap = F, alpha = 0.05)
#####Raw data analysis complete, Transformed data set starts now
rtTSA=sqrt(TSAcfu)
rtFMA=sqrt(FMAcfu)
lgTSA=log(TSAcfu)
lgFMA=log(FMAcfu)


qqPlot(rtTSA)
qqPlot(rtFMA)
qqPlot(lgTSA)
qqPlot(lgFMA)

boxplot(lgTSA~captureTSA)
boxplot(lgFMA~captureFMA)#Oh, no what does that mean? Keep going

aovlgTSA=aov(lgTSA~captureTSA)
aovlgFMA=aov(lgFMA~captureFMA)#Oh no, keep going
#Let's try to do Shapiro on the raw transformed data
shapiro.test(lgFMA)#That doesn't look good, let's investigate the transformed data
lgFMA#Notice the negative infinity? Let's remove that
#Notice that the -Inf occurs at 54
#Open the CSV and remove sample 54, and adjust for capture locale length
lgFMA=read.csv("FMAlg.csv", header = T)
as.factor(lgFMA$Capture.Locale)
captureFMAlg<-factor(lgFMA$Capture.Locale,levels=c("Estrella", "South", "Phoenix", "Tempe"))
sample.size.lgFMA= aggregate(.~captureFMAlg, data=lgFMA, length)
lgFMAcfu=lgFMA$Adj..FDB.CFU.mL
lgFMAcfu#Notice the -inf is gone, caused by a 0 in the data

aovlgFMA=aov(lgFMAcfu~captureFMAlg)


shapiro.test(aovlgTSA$residuals)#Improved, but still not normal
shapiro.test(aovlgFMA$residuals)#Unchanged, and not normal

bartlett.test(lgTSA~captureTSA)#Final test is still Kruskal
#But we were really close to making it an ANOVA, so let's do an ANOVA too
bartlett.test(lgFMAcfu~captureFMAlg)#Kruskal is still the test

kruskal.test(lgTSA,as.factor(captureTSA))
kruskal.test(lgFMAcfu,as.factor(captureFMAlg))
#Significant differences among the capture locations
#Let's see where the differences are with a Dunn's test, a ranked version of Tukey
#And let's compare with our originals
dunn.test(TSAcfu,g=captureTSA, method = c("none"), kw=T, label = TRUE, wrap = F, alpha = 0.05)
dunn.test(lgTSA,g=captureTSA, method = c("none"), kw=T, label = TRUE, wrap = F, alpha = 0.05)#Transformation did not change our final conclusions
dunn.test(FMAcfu,g=captureFMA, method = c("none"), kw=T, label = TRUE, wrap = F, alpha = 0.05)
dunn.test(lgFMAcfu,g=captureFMAlg, method = c("none"), kw=T, label = TRUE, wrap = F, alpha = 0.05)#Transformation did not change our final conclusions

#Let's do ANOVA on our transformed TSA data
summary(aov(lm(lgTSA~captureTSA)))
#Original TSA
summary(aov(lm(TSAcfu~captureTSA)))
#It appears that the raw data gives us a better chance of our results being sig
summary(aov(lm(FMAcfu~captureFMA)))
summary(aov(lm(lgFMAcfu~captureFMAlg)))

#Let's use a for loop to look at central limit therom
n=1000
data1=rnorm(n,mean=1000,sd=250)#Creates a normal distrubution with n samples, mean 100 and standard dev 250

nn=10
reps=1000
means=c()
for (i in 1:reps ) {
  sample1=data1[round(runif(nn, min=0, max=n),0)]
  addM=mean(sample1)
  means=c(means,addM)
}

hist(data1,breaks=20,xlim=c(0,3000)) #Plots the distrubution generated
#Edit the sample size or reps to see different distrubutions
#Prove that a distrubution is normal (histogram is first step)
qqPlot(data1)
shapiro.test(data1)#Quatatative normaility check

#Use ggplot to view the data differently, and observe how graphing can affect interpretation
qplot(data1, binwidth=100)
qplot(data1, binwidth=10)
hist(data1)

###This is where we begin an automatic data analysis 
######This is where you set up your data to be auto-analyzed
#Sample-size comparison permutations must be edited in the code before analysis
#If analysis returns 666, an error has occured. Damn demonic intrusion.
#Start off with the TSA data
independent=TSA$Total.CFU.ml
dependent<-factor(TSA$Capture.Locale,levels=c("Estrella", "South", "Phoenix", "Tempe"))
shapiro_in=aov(independent~dependent)#(aov(independent~dependent)
sample.size1=count(TSA, "Capture.Locale")$freq[1]
sample.size2=count(TSA, "Capture.Locale")$freq[2]
sample.size3=count(TSA, "Capture.Locale")$freq[3]
sample.size4=count(TSA, "Capture.Locale")$freq[4]
#######

(ifelse(shapiro.test(shapiro_in$residuals)$p.value>0.05,(ifelse(shapiro.test(shapiro_in$residuals)$p.value>0.05, (ifelse(bartlett.test(independent~dependent)$p.value>0.2, summary(aov(lm(independent~dependent))), 666)), 
                                                                (ifelse(bartlett.test(independent~dependent)$p.value>0.05 | bartlett.test(independent~dependent)$p.value<=0.2 | bartlett.test(independent~dependent)$p.value<0.05, (ifelse(sample.size1==sample.size2 && sample.size3==sample.size4 && sample.size1==sample.size3 && sample.size2==sample.size4 && sample.size1==sample.size4, 
                                                                                                                                                                                                                                           (ifelse(bartlett.test(independent~dependent)$p.value>=0.05 | bartlett.test(independent~dependent)$p.value<=0.2, summary(aov(lm(independent~dependent))))), 
                                                                                                                                                                                                                                           (ifelse(bartlett.test(independent~dependent)$p.value<=0.05, oneway.test(independent~dependent,var.equal=F), oneway.test(independent~dependent,var.equal=T))))))))),  
        (ifelse(shapiro.test(shapiro_in$residuals)$p.value<0.05, 
                (ifelse(bartlett.test(independent~dependent)$p.value>0.2, 
                        (ifelse(sample.size1==sample.size2 && sample.size3==sample.size4 && sample.size1==sample.size3 && sample.size2==sample.size4 && sample.size1==sample.size4, summary(aov(lm(independent~dependent))), kruskal.test(independent,as.factor(dependent))$p.value)), 
                        (ifelse(bartlett.test(independent~dependent)$p.value>0.05 | bartlett.test(independent~dependent)$p.value<=0.2 | bartlett.test(independent~dependent)$p.value<0.05, kruskal.test(independent,as.factor(dependent)), 666))))))))
####Now yo know to do the Kruskal-Wallis test
#For some reason the automation only returns the chi-squared value
#This reason is why you manually repeat the appropriate test
kruskal.test(independent,as.factor(dependent))

#Now do the FMA data
######This is where you set up your data to be auto-analyzed
#Sample-size comparison permutations must be edited in the code before analysis
#If analysis returns 666, an error has occured. Damn demonic intrusion.
independent=FMAcfu
dependent<-factor(FMA$Capture.Locale,levels=c("Estrella", "South", "Phoenix", "Tempe"))
shapiro_in=aov(independent~dependent)#(aov(independent~dependent)
sample.size1=count(FMA, "Capture.Locale")$freq[1]
sample.size2=count(FMA, "Capture.Locale")$freq[2]
sample.size3=count(FMA, "Capture.Locale")$freq[3]
sample.size4=count(FMA, "Capture.Locale")$freq[4]
#######

(ifelse(shapiro.test(shapiro_in$residuals)$p.value>0.05,(ifelse(shapiro.test(shapiro_in$residuals)$p.value>0.05, (ifelse(bartlett.test(independent~dependent)$p.value>0.2, summary(aov(lm(independent~dependent))), 666)), 
                                                                (ifelse(bartlett.test(independent~dependent)$p.value>0.05 | bartlett.test(independent~dependent)$p.value<=0.2 | bartlett.test(independent~dependent)$p.value<0.05, (ifelse(sample.size1==sample.size2 && sample.size3==sample.size4 && sample.size1==sample.size3 && sample.size2==sample.size4 && sample.size1==sample.size4, 
                                                                                                                                                                                                                                           (ifelse(bartlett.test(independent~dependent)$p.value>=0.05 | bartlett.test(independent~dependent)$p.value<=0.2, summary(aov(lm(independent~dependent))))), 
                                                                                                                                                                                                                                           (ifelse(bartlett.test(independent~dependent)$p.value<=0.05, oneway.test(independent~dependent,var.equal=F), oneway.test(independent~dependent,var.equal=T))))))))),  
        (ifelse(shapiro.test(shapiro_in$residuals)$p.value<0.05, 
                (ifelse(bartlett.test(independent~dependent)$p.value>0.2, 
                        (ifelse(sample.size1==sample.size2 && sample.size3==sample.size4 && sample.size1==sample.size3 && sample.size2==sample.size4 && sample.size1==sample.size4, summary(aov(lm(independent~dependent))), kruskal.test(independent,as.factor(dependent))$p.value)), 
                        (ifelse(bartlett.test(independent~dependent)$p.value>0.05 | bartlett.test(independent~dependent)$p.value<=0.2 | bartlett.test(independent~dependent)$p.value<0.05, kruskal.test(independent,as.factor(dependent)), 666))))))))
###Again, the appropriate test is Kruskal
#YOu can see that the chi-squared is different, thus the automation worked twice
kruskal.test(independent,as.factor(dependent))
###Lastly, let's do the transformed TSA data
######This is where you set up your data to be auto-analyzed
#Sample-size comparison permutations must be edited in the code before analysis
#If analysis returns 666, an error has occured. Damn demonic intrusion.
independent=lgTSA
dependent<-factor(TSA$Capture.Locale,levels=c("Estrella", "South", "Phoenix", "Tempe"))
shapiro_in=aov(independent~dependent)#(aov(independent~dependent)
sample.size1=count(TSA, "Capture.Locale")$freq[1]
sample.size2=count(TSA, "Capture.Locale")$freq[2]
sample.size3=count(TSA, "Capture.Locale")$freq[3]
sample.size4=count(TSA, "Capture.Locale")$freq[4]
#######

(ifelse(shapiro.test(shapiro_in$residuals)$p.value>0.05,(ifelse(shapiro.test(shapiro_in$residuals)$p.value>0.05, (ifelse(bartlett.test(independent~dependent)$p.value>0.2, summary(aov(lm(independent~dependent))), 666)), 
                                                                (ifelse(bartlett.test(independent~dependent)$p.value>0.05 | bartlett.test(independent~dependent)$p.value<=0.2 | bartlett.test(independent~dependent)$p.value<0.05, (ifelse(sample.size1==sample.size2 && sample.size3==sample.size4 && sample.size1==sample.size3 && sample.size2==sample.size4 && sample.size1==sample.size4, 
                                                                                                                                                                                                                                           (ifelse(bartlett.test(independent~dependent)$p.value>=0.05 | bartlett.test(independent~dependent)$p.value<=0.2, summary(aov(lm(independent~dependent))))), 
                                                                                                                                                                                                                                           (ifelse(bartlett.test(independent~dependent)$p.value<=0.05, oneway.test(independent~dependent,var.equal=F), oneway.test(independent~dependent,var.equal=T))))))))),  
        (ifelse(shapiro.test(shapiro_in$residuals)$p.value<0.05, 
                (ifelse(bartlett.test(independent~dependent)$p.value>0.2, 
                        (ifelse(sample.size1==sample.size2 && sample.size3==sample.size4 && sample.size1==sample.size3 && sample.size2==sample.size4 && sample.size1==sample.size4, summary(aov(lm(independent~dependent))), kruskal.test(independent,as.factor(dependent))$p.value)), 
                        (ifelse(bartlett.test(independent~dependent)$p.value>0.05 | bartlett.test(independent~dependent)$p.value<=0.2 | bartlett.test(independent~dependent)$p.value<0.05, kruskal.test(independent,as.factor(dependent)), 666))))))))

##We got the exact same Kruskal value on the transformed data as on the raw
##Let's prove that is correct
kruskal.test(TSAcfu, as.factor(captureTSA))
kruskal.test(lgTSA, as.factor(captureTSA))






