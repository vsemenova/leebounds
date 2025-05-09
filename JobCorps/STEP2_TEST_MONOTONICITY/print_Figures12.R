rm(list=ls())


# install.packages(c("latex2exp"),repo="http://cran.rstudio.com/",lib="~/apps/R_3.5.1:")
library(ggplot2)
library(latex2exp)
# install.packages(c("ggforce","latex2exp"))

name<-"Vira"
if (name =="Vira") {
  my_path<-"/Users/virasemenova/Dropbox (MIT)/leebounds/"
}
if (name == "Mengsi") {
  my_path<-"C:/Users/mengs/Dropbox/leebounds"
  # my_path<-"/Users/mac/Dropbox/leebounds"
}
aspect_ratio<-1.5

#### REPLICATE FIGURES 1 and 2
name="Vira"
Figure1_dataset<-read.csv(paste0(my_path,"/JobCorps/Tables/STEP2_Monotonicity/Vira/Figure1.csv"))
Figure1_dataset$group<-as.factor(Figure1_dataset$group)
  
Figure2_dataset<-read.csv(paste0(my_path,"/JobCorps/Tables/STEP2_Monotonicity/Vira/Figure2.csv"))
test_result<-read.csv(paste0(my_path,"/JobCorps/Tables/STEP2_Monotonicity/Vira/test_result.csv"))

Figure1_dataset$fraction<-c(Figure2_dataset$fraction,1-Figure2_dataset$fraction)

ggplot(data=Figure1_dataset[setdiff(4*(3:104),c(209:219)),])+aes(x=weeks,y=delta,group=group,size=fraction/2)+
  geom_line(data=Figure1_dataset[setdiff(4*(3:104),c(209:219)),],
            aes(x=weeks,y=emp_rate_diff),
            lwd=0.5)+
  geom_point(aes(color=group))+xlab("Weeks since random assignment")+
  geom_vline(xintercept=90,lwd=0.25)+
  ylab("Treatment-control difference in employment rate")+
  theme_bw()+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank(),
        axis.line = element_line(colour = "black"))+
  geom_hline(yintercept=0)+
  scale_color_manual(values=c("gray","black"))+ theme(legend.position = "none")+ 
  annotate("text", x=120, y=0.1, label='Predicted Positive Effect',size=6)+ 
  annotate("text", x=120, y=-0.1, label='Predicted Negative Effect',size=6)+
  scale_x_continuous(breaks = c(10,40,80,90,120,160,200),expand=c(0,0),limits=c(0,210))+
  scale_y_continuous(breaks=c(-0.1,0.,0.1))
ggsave(paste0(my_path,"/JobCorps/Figures/Figure1.png"),height=7,width=aspect_ratio*7)

sum(test_result$pvalue<=0.01)
sum(test_result$pvalue>0.01 & test_result$pvalue<=0.05)
weeks<-1:208
weeks_001<-(1:208)[test_result$pvalue<=0.01]
weeks_05<-setdiff((1:208)[test_result$pvalue<=0.05], weeks_001)

Figure2_dataset$pvalue<-test_result$pvalue
Figure2_dataset$group<-"pvalue > 0.1"
Figure2_dataset$group[Figure2_dataset$pvalue <=0.01]<-"pvalue < 0.01"
Figure2_dataset$group[Figure2_dataset$pvalue >0.01 & Figure2_dataset$pvalue<=0.05 ]<-"pvalue in [0.01, 0.1]"

ggplot(data=Figure2_dataset)+
  geom_point( aes(x=weeks,y=pvalue,color=group))+
  geom_point( aes(x=weeks,y=pvalue),color="white")+
  scale_color_manual(values=c("gray100","gray90","gray75"),name="Monotonicity test p-value",
                     labels=c("> 0.05", "[0.01, 0.05]", "< 0.01" ))+
  geom_rect(data = data.frame(xmin = 80,
                              xmax = 91,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray75")+
  geom_rect(data = data.frame(xmin = 104,
                              xmax = 114,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray75")+
  geom_rect(data = data.frame(xmin = 124,
                              xmax = 139,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray75")+
  geom_rect(data = data.frame(xmin = 155,
                              xmax = 186,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray75")+
  geom_rect(data = data.frame(xmin = 190,
                              xmax = 192,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray75")+
  geom_rect(data = data.frame(xmin = 63,
                              xmax = 68,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_rect(data = data.frame(xmin = 71,
                              xmax = 80,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_rect(data = data.frame(xmin = 90,
                              xmax = 105,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_rect(data = data.frame(xmin = 114,
                              xmax = 124,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_rect(data = data.frame(xmin = 139,
                              xmax = 148,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_rect(data = data.frame(xmin = 153,
                              xmax = 154,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_rect(data = data.frame(xmin = 193,
                              xmax = 208,
                              ymin = -Inf,
                              ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray90")+
  geom_line(data=Figure2_dataset,
            aes(x=weeks,y=fraction),
            lwd=1)+xlab("Weeks since random assignment")+
  ylab(TeX("Fraction of applicants with positive employment effect"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank(),
        axis.line = element_line(colour = "black"))+theme(legend.position="bottom",legend.text=element_text(size=15),legend.title = element_text(size=15))+
  scale_x_continuous(breaks = c(0,40,80,120,160,200),expand=c(0,0),limits=c(0,210))+
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,0.9),expand=c(0,0),limits=c(0,0.9))

ggsave(paste0(my_path,"/JobCorps/Figures/Figure2.png"),height=7,width=aspect_ratio*7)

