#This script corresponds to the following article
#Anion exchange capacity explains deep soil nitrate accumulation in Brazilian Amazon croplands
#Alexandra Huddell, Christopher Neill, Cheryl A. Palm, Darlisson Nunes, and Duncan N. L. Menge


#AEC calculations
setwd('') #set to the appropriate working directory containing the dat file
dat<-read.csv('AEC_soil_nutrients_for_github.csv')

library(ggplot2)
library(dplyr)
library(patchwork)

#translating depth variable to factor
dat$depth_num<-as.factor(dat$depth_num)

# plots -------------------------------------------------------------------

#theme for plots
theme_default <- function(axis_text_size = 11) {
  theme(text=element_text(size=16, colour="black"),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text.x=element_text(size=axis_text_size, colour="black"),
        axis.text.y=element_text(size=axis_text_size, colour="black"),
        axis.title.y=element_text(angle=90, vjust=2),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.background=element_rect(fill="white", color="white"),
        legend.position = "right",
        legend.key=element_blank())
}

cbPalette <- c( "#D55E00","blue","#009E73")

dat$treatment<-factor(dat$treatment, levels=c("maize","soy","forest")) 

#calculating median NO3-N
NO3_median<-dat%>% 
  group_by(depth_num,treatment)%>%
  summarise(median=median(KCL_ug_NO3_N_g_soil, na.rm=T))



#plotting soil NO3-N concentrations by treatment and depth
p1<-ggplot()+
  geom_jitter(data=dat,
             aes(x=forcats::fct_rev(depth_num), y=KCL_ug_NO3_N_g_soil,
                 col=treatment)) +
  geom_line(data=NO3_median,aes(x =forcats::fct_rev(depth_num), y=median,
                                    col=treatment, group=treatment)) +
  geom_point(data=NO3_median,aes(x =forcats::fct_rev(depth_num), y=median, 
                                     col=treatment, group=treatment),shape=17,size=3 )  +
  coord_flip() +
  scale_x_discrete(labels=rev(c('0-1','1-2', '2-3', '3-4', '4-5',
                                '5-6', '6-7', '7-8')),
                   name='soil depth (m)') +
  scale_y_continuous(name=expression(paste('Soil ', NO[3]^'-',~
                                           '(µg N g ',soil^-1, ')')),
                     position= 'right') +
    geom_segment(aes(x=0,xend=8.5,y=0,yend=0), colour="black") +
  geom_segment(aes(x=8.5,xend=8.5,y=0,yend=50),colour="black") +
  theme_default()+
  theme(legend.position = "none")+
  scale_color_manual(name="Land use",
                     values = c(cbPalette[3],cbPalette[1],"blue"),
                     breaks=c("forest","soy","maize"),
                     labels = c("forest",  "soybean","soybean-maize")) 
p1


NH4_median<-dat%>% 
  group_by(depth_num,treatment)%>%
  summarise(median=median(KCL_ug_NH4_N_g_soil, na.rm=T))

#plot of soil NH4-N in KCl
p2<-   ggplot()+
  geom_jitter(data=dat,
             aes(x=forcats::fct_rev(depth_num), y=KCL_ug_NH4_N_g_soil,
                 col=treatment),alpha=0.7, height=0.1, width=0.1) +
  geom_line(data=NH4_median,aes(x =forcats::fct_rev(depth_num), y=median,
                                col=treatment, group=treatment)) +
  geom_point(data=NH4_median,aes(x =forcats::fct_rev(depth_num), y=median, 
                                 col=treatment, group=treatment),shape=17,size=3 )  +
  
  coord_flip() +
  scale_x_discrete(labels=rev(c('0-1','1-2', '2-3', '3-4', '4-5',
                                '5-6', '6-7', '7-8')),
                   name='soil depth (m)') +
  scale_y_continuous(name=expression(paste('Soil ', NH[4]^'+',~
                                             '(µg N g ',soil^-1, ')')),
                     position= 'right') +
  geom_segment(aes(x=0,xend=8.5,y=0,yend=0), colour="black") +
  geom_segment(aes(x=8.5,xend=8.5,y=0,yend=8),colour="black") +
  theme_default()+
  theme(legend.position = "bottom", legend.text=element_text(size=16))+
  scale_color_manual(name="Land use",
                     values = c(cbPalette[3],cbPalette[1],"blue"),
                     breaks=c("forest","soy","maize"),
                     labels = c("forest",  "soybean","soybean-maize")) 
p2


#arrange plots on grid
AEC_soil_nutrients<-p1+ p2 +plot_annotation(tag_levels = 'a')
ggsave('AEC_soil_nutrients.jpeg', plot = AEC_soil_nutrients, device = NULL,
                       width=8, height = 4, units = 'in', path = NULL,
                      scale=1.5, dpi = 300, limitsize = TRUE)


# soil N concentration to kg N/ha m soil ----------------------------------

#multiplying soil N concentrations by bulk density values and converting to
#kg N per m of soil
soil_N_to_kg_N_m_soil_ha<-function(soil_N,bd){soil_N*bd*10}
dat$kg_N_m_soil<-soil_N_to_kg_N_m_soil_ha(soil_N = dat$KCL_ug_NO3_N_g_soil,
                                          bd=dat$bulk_density)

#summing total NO3-N in top 8m by depth and treatment
N_summed<- dat %>%
  group_by(treatment, site) %>%
    summarize(kg_N_soil_total=sum(kg_N_m_soil))
N_summed

write.csv(N_summed,file = "Nitrate summed by site top 8m.csv")

N_summed$treatment<-factor(N_summed$treatment, levels=c("forest","soy","maize")) 

#plot of total soil NO3-N in top 8m
p3<- ggplot(N_summed,aes(x =treatment, y=log10(kg_N_soil_total) ,
                    col=treatment))+
  geom_point(alpha=.6,position = position_jitterdodge(jitter.width=.4),
             size=5) +
  scale_x_discrete(name= " ",labels=(c('forest','soybean', 'soybean-\nmaize'))) +
  geom_segment(aes(x=0,xend=0,y=1,yend=3.4),colour="black") +
  geom_segment(aes(x=0,xend=3,y=1,yend=1),colour="black") +
  theme_default()+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1,
                                                            size = 18))+
  scale_color_manual(name="Land use",
                     values = c(cbPalette[3],cbPalette[1],"blue"),
                     breaks=c("forest","soy","maize"),
                     labels = c("forest",  "soybean","soybean-maize")) +
  annotate("text", x = c(1,2,3), y = c(3.8,3.8,3.8), 
           label = c("a","b","b"),
           size=5 ,  fontface="bold")+
  scale_y_continuous(breaks=c(1,2,2.69897,3,3.30103),
                     labels=c(10,100,500,1000,2000),
                     name=expression(paste('Total 0-8m soil ',~NO[3]^'-',~'(kg N',~
                                             ha^-1,')')))

p3


#arrange plots on grid
fig1<-p1+ p2 + p3+plot_layout(guides = "collect")+plot_annotation(tag_levels = 'a')
fig1
ggsave('fig1.jpeg', plot = fig1, device = NULL,
       width=14, height = 7, units = 'in', path = NULL,
       scale=1, dpi = 300, limitsize = TRUE)

#ANOVA on total soil N in kg N per ha on untransformed data
summary(aov(lm(kg_N_soil_total ~ as.factor(treatment), data = N_summed)))
lm1<-(lm(kg_N_soil_total ~ as.factor(treatment), data = N_summed))
aov1<-aov(lm1)
summary(aov1)
TukeyHSD(aov1)

#ANOVA on total soil N in kg N per ha on log10-transformed data
summary(aov(lm(log10(kg_N_soil_total) ~ as.factor(treatment), data = N_summed)))
lm1<-(lm(log10(kg_N_soil_total) ~ as.factor(treatment), data = N_summed))
aov1<-aov(lm1)
TukeyHSD(aov1)



