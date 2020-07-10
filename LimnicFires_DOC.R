library(reshape2)
library(ggplot2)
library(outliers)
library(gridExtra)
library(data.table)
library(ggsignif)
library(rstatix)

setwd("C:/Users/Thuile/PhD/LimnicFires/DOC")


#####read in and orgnaize data#####
data<-read.csv("DOC.csv",sep=",", header=F, check.names = F)

info<-data[,1]
DOC<-data[,2]
treat<-substr(info, start=1,stop=1)
flume<-substr(info, start=3,stop=3)
time<-substr(info, start=4,stop=4)
info2<-cbind(treat,flume,time)


data1<-cbind.data.frame(info2,DOC)
colnames(data1)[4]<-"DOC"
rownames(data1) <- NULL
mean_C <- mean(data1[c(1:5,11:15,21:25,31:35,41:45,51:55),4])
mean_T <- mean(data1[c(6:10,16:20,26:30,36:40,46:50,56:60),4])

# melt df
df_m<-melt(data1)
df_m<-as.data.table(df_m)
ggplot(df_m, aes(x=time,  y=value, fill=treat)) + geom_boxplot() 

#####plots#####
# boxplot with molten df; per time point
ggplot(df_m, aes(x=time,  y=value, fill=treat)) + geom_boxplot() + scale_fill_manual(values=c("white", "grey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  xlab("Time after charcoal addition (h)") + ylab("DOC (mg/l)")  


#geom_signif(comparisons = list(c("1", "2", "3", "4", "6", "8")), map_signif_level=TRUE, na.rm = TRUE) #tool to include significance in the plot


# boxplot with molten df; general C vs T
ggplot(df_m, aes(x=treat,  y=value, fill=treat)) + geom_boxplot() + scale_fill_manual(values=c("white", "grey")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  ylab("DOC (mg/l)")



# de-cast while computing mean
doc_mean<-dcast(df_m, treat+time ~ variable, fun.aggregate = mean)
doc_sd<-dcast(df_m, treat+time ~ variable, fun.aggregate = sd)
doc <- cbind(doc_mean, doc_sd[,3])
names(doc)[4] <- "sd"
str(doc)
df_m1<-melt(doc[,1:3])
df_m1<-as.data.table(df_m1)
doc<-as.data.table(doc)



# point-plot with data.table
ggplot(doc, aes(x=time, y=DOC, color=treat, fill=treat)) + scale_color_manual(values=c("black", "black")) + scale_fill_manual(values=c("white", "darkgrey")) + 
  geom_errorbar(aes(ymin=DOC-sd, ymax=DOC+sd, color=treat), width=0.2) + geom_point(shape=21, size=3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5))


# point plot with base R; mean
plot(as.numeric(doc[doc$treat=="C",doc$time]), doc[doc$treat=="C",doc$DOC], type="p", ylim=c(2,6), pch=21, bg="grey80", cex=1.4, xlab="Time (h)")
points(as.numeric(doc$time[7:12]), doc$DOC[7:12], type="p", ylim=c(2,6), pch=21, bg="grey20", cex=1.4)

mod_contr <- lm(doc$DOC[1:6] ~ as.numeric(doc$time[1:6]))
abline(mod_contr, col="grey80")
mod_treat <- lm(doc$DOC[7:12] ~ as.numeric(doc$time[7:12]))
abline(mod_treat, col="grey20")


# point plot with base R; singel data points
par(mar=c(5,5,2,2))
plot(as.character(df_m$time[df_m$treat=="C"]), df_m$value[df_m$treat=="C"], type="p", ylim=c(2.5,7), pch=21, bg="grey80", cex=2, xlab="Time (h)", ylab= "DOC (mg/l)", cex.lab=2, cex.axis=2)
points(as.character(df_m$time[df_m$treat=="T"]), df_m$value[df_m$treat=="T"], type="p", ylim=c(2,7), pch=21, bg="grey20", cex=2)

mod_contr <- lm(df_m$value[df_m$treat=="C"] ~ as.numeric(as.character(df_m$time[df_m$treat=="C"])))
abline(mod_contr, col="grey80", lwd=2)
mod_treat <- lm(df_m$value[df_m$treat=="T"] ~ as.numeric(as.character(df_m$time[df_m$treat=="T"])))
abline(mod_treat, col="grey20", lwd=2)

legend("topleft", legend=c("Treatment", "Control"), pch=c(21,21), pt.bg=c("grey20", "grey80"), cex=2, y.intersp = 0.75, bty="n")





#####simple statistic#####
t.test(df_m[,value]~df_m[,treat], paired=T, p.adjust.methods="bonferroni") #paired=T because of dependent samples; bonferroni adjustments also because of dependent samples

# per time point
t.test(df_m[1:10,value]~df_m[1:10,treat]) #not
t.test(df_m[11:20,value]~df_m[11:20,treat]) #not
t.test(df_m[21:30,value]~df_m[21:30,treat]) #not
t.test(df_m[31:40,value]~df_m[31:40,treat]) #not
t.test(df_m[41:50,value]~df_m[41:50,treat]) #not
t.test(df_m[51:60,value]~df_m[51:60,treat]) #not


#####including dependendness; NOT WORKING#####
# normal tow way anova; statistically problematic!!!
mod<-aov(DOC~treat*time,data=df_m) #higher significance but probably statistcally not correct because of repeated measures!
summary(mod)

TukeyHSD(mod) #significance comes from all T vs all C; singel C vs T comparisons per timepoint are not significant; 


# -> repeated measures is old school and also did not work
#write.table(df_m, "C:/Users/Thuile/PhD/LimnicFires/DOC/DOC_molten.txt", sep=",")

