rm(list=ls())


# install.packages(c("latex2exp"),repo="http://cran.rstudio.com/",lib="~/apps/R_3.5.1:")
library(ggplot2)
library(latex2exp)
library(geomnet)
library(ggforce)
# install.packages(c("ggforce","latex2exp"))

aspect_ratio<-1.5

name<-"Vira"
if (name =="Vira") {
  my_path<-"/Users/virasemenova/Dropbox (MIT)/leebounds/"
}
if (name == "Mengsi") {
  # my_path<-"C:/Users/mengs/Dropbox/leebounds"
  my_path<-"/Users/mac/Dropbox/leebounds"
}


Figure3_dataset<-read.csv(paste0(my_path,"/JobCorps/Tables/STEP5_Misc/", name, "/Figure3.csv"))
Figure3_dataset$group<-as.factor(Figure3_dataset$group)


ggplot(data=Figure3_dataset)+aes(x=weeks,y=bound,group=group)+
  xlab("Weeks since random assignment")+
  geom_point(aes(col=group),size=3,shape = 21,fill = "white")+
  geom_smooth(aes(col=group),se=F)+
  ylab("Expected  log wage in control status ")+
  theme_bw()+
  theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),
        axis.title=element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_color_manual(values=c("gray","black"),
                     labels = c("lower bound", "upper bound"))+ theme(legend.position = "none")+ 
  scale_x_continuous(breaks = c(0,40,80,120,160,200),expand=c(0,0),limits=c(0,210))+
  scale_y_continuous(breaks=c(1.2,1.6,1.8,2.0))
ggsave(paste0(my_path,"/JobCorps/Figures/Figure3.png"),height=7,width=aspect_ratio*7)

