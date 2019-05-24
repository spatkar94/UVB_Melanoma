#Read in data
TCGA_data <- read_csv("./TCGA_data.csv")
clinical_TCGA <- read_csv("./clinical_TCGA.csv")
TCGA_data = as.data.frame(TCGA_data)
rownames(TCGA_data) = TCGA_data$X1
clinical_TCGA = as.data.frame(clinical_TCGA)
rownames(clinical_TCGA) = clinical_TCGA$X1

#divide into low and high groups
TCGA_data$ITH_high_low = c("low","high")[as.numeric(TCGA_data$ITH > median(TCGA_data$ITH,na.rm = T))+1]
TCGA_data$ML_high_low = c("low","high")[as.numeric(TCGA_data$mutation_load > median(TCGA_data$mutation_load,na.rm = T))+1]
TCGA_data$CNV_high_low = c("low","high")[as.numeric(TCGA_data$cnv_load > median(TCGA_data$cnv_load,na.rm = T))+1]
TCGA_data$diversity_high_low = c("low","high")[as.numeric(TCGA_data$diversity > median(TCGA_data$diversity,na.rm = T))+1]
TCGA_data$sig7_high_low = c("low","high")[as.numeric(TCGA_data$sig7 > median(TCGA_data$sig7,na.rm = T))+1]


#TCGA survival analysis

library(GGally)
library(survival)
library(ggplot2)
library(ggpubr)
common = intersect(clinical$X1, TCGA_data$X1)

#CYT score transformed to log-scale for visualization
survdat = data.frame(time = clinical_TCGA[common,5],                                          
                     status = clinical_TCGA[common,6],
                     age = clinical_TCGA[common,2],
                     stage = as.numeric(as.factor(clinical_TCGA[common,7])),
                     purity = TCGA_data[common,6],
                     ITH = TCGA_data[common,9],
                     diversity = TCGA_data[common,12],
                     MUTATION_LOAD = TCGA_data[common,10],
                     CNV_LOAD = TCGA_data[common,11],
                     SIG7 = TCGA_data[common,13],
                     CYT = log(TCGA_data[common,8]))

surv = survfit(Surv(time,status) ~ ITH, data = survdat)
diff = survdiff(Surv(time,status) ~ ITH, data = survdat)

p1 = ggsurv(surv) + labs(subtitle = paste("Log rank test p-value: ",format(1 - pchisq(diff$chisq,length(diff$n) - 1),digits = 2),sep = "")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + guides(linetype = F) + scale_color_manual(name = "ITH", values = c("blue", "red"), labels = c("blue" = "low","red" = "high")) + scale_y_continuous(labels = scales::percent, limits = c(0,1)) +  scale_x_continuous(limits = c(0,10000)) + xlab("Time (Days)") + ylab("Survival(%)")
print(p1)




surv = survfit(Surv(time,status) ~ MUTATION_LOAD, data = survdat)
diff = survdiff(Surv(time,status) ~ MUTATION_LOAD, data = survdat)

p2 = ggsurv(surv) + labs(subtitle = paste("Log rank test p-value: ",format(1 - pchisq(diff$chisq,length(diff$n) - 1),digits = 2),sep = "")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + guides(linetype = F) + scale_color_manual(name = "Mutation load", values = c("blue", "red"), labels = c("blue" = "low","red" =  "high")) + scale_y_continuous(labels = scales::percent, limits = c(0,1)) +  scale_x_continuous(limits = c(0,10000)) + xlab("Time (Days)") + ylab("Survival(%)")
print(p2)

surv = survfit(Surv(time,status) ~ CNV_LOAD, data = survdat)
diff = survdiff(Surv(time,status) ~ CNV_LOAD, data = survdat)

p2 = ggsurv(surv) + labs(subtitle = paste("Log rank test p-value: ",format(1 - pchisq(diff$chisq,length(diff$n) - 1),digits = 2),sep = "")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + guides(linetype = F) + scale_color_manual(name = "CNV load", values = c("blue", "red"), labels = c("blue" = "low","red" =  "high")) + scale_y_continuous(labels = scales::percent, limits = c(0,1)) +  scale_x_continuous(limits = c(0,10000)) + xlab("Time (Days)") + ylab("Survival(%)")
print(p2)


survdat$ITHvsMutationLoad = strata(survdat$ITH,survdat$MUTATION_LOAD)
surv = survfit(Surv(time,status)~ITHvsMutationLoad,data = survdat)
diff = survdiff(Surv(time,status)~ITHvsMutationLoad,data = survdat)


p3 = ggsurv(surv) + labs(subtitle = paste("Log rank test p-value: ",format(1 - pchisq(diff$chisq,length(diff$n) - 1),digits = 3),sep = "")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + guides(linetype = F) + scale_color_discrete( name = "ITH, Mutation Load") +  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +  scale_x_continuous(limits = c(0,10000)) + xlab("Time (Days)") + ylab("Survival(%)")
print(p3)

survdat$ITHvsCNVLoad = strata(survdat$ITH,survdat$CNV_LOAD)
surv = survfit(Surv(time,status)~ITHvsCNVLoad,data = survdat)
diff = survdiff(Surv(time,status)~ITHvsCNVLoad,data = survdat)


p4 = ggsurv(surv) + labs(subtitle = paste("Log rank test p-value: ",format(1 - pchisq(diff$chisq,length(diff$n) - 1),digits = 3),sep = "")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + guides(linetype = F) + scale_color_manual( name = "ITH, CNV Load", values = c("#7CAE00","#F8766D", "#C77CFF","#00BFC4"), labels = c("#7CAE00" = "low, low", "#F8766D" = "low, high", "#C77CFF" = "high, low", "#00BFC4" = "high, high")) +  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +  scale_x_continuous(limits = c(0,10000)) + xlab("Time (Days)") + ylab("Survival(%)")
print(p4)

survdat$clones = as.character(TCGA_data[common,2])

surv = survfit(Surv(time,status)~clones,data = survdat)
diff = survdiff(Surv(time,status)~clones,data = survdat)

cols = c("#42a7f4", "#3863d8", "#31f3f9","#9e32fc","#ff30e6","#fc1919")
vals = c("2", "3","1","4","5","6")
p5 = ggsurv(surv) + labs(subtitle = paste("Log rank test p-value: ",format(1 - pchisq(diff$chisq,length(diff$n) - 1),digits = 3),sep = "")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + guides(linetype = F) + scale_color_manual(name = "# clones", values = cols, labels = vals) +  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +  scale_x_continuous(limits = c(0,10000)) + xlab("Time (Days)") + ylab("Survival(%)")
print(p5)

cmpr = list(c("low","high"))
print(ggplot(data = survdat[!is.na(survdat$ITH),], aes(x = ITH, y = CYT, fill = ITH)) + geom_dotplot(binaxis='y', stackdir='center',dotsize = 0.7) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + scale_fill_manual(name = "ITH", values = c("red","blue"), labels = c("high","low")) + xlab("ITH") + ylab("CYT") + stat_compare_means(label.y = 10))


print(ggplot(data = survdat[!is.na(survdat$clones),], aes(x = clones, y = CYT, fill = clones)) + geom_dotplot(binaxis='y', stackdir='center',dotsize = 0.7) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + scale_fill_manual(name = "ITH", values = cols[order(as.numeric(vals))], labels = vals[order(as.numeric(vals))]) + xlab("# clones") + ylab("CYT") + stat_compare_means(label.y = 10))

