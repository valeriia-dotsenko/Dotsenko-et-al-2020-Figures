# Figure 2. from "Genome-Wide Transcriptomic Analysis of Intestinal Mucosa in Celiac Disease Patients on a Gluten-Free Diet and Postgluten Challenge" (https://doi.org/10.1016/j.jcmgh.2020.07.010) 
#Gluten content (ng GIP/g stool) in patients’ stool samples taken before, during, and after gluten challenge, measured by iVYLISA GIP assay. 
#Gluten consumption per day was assessed by number of cookies eaten, taking into account the content of gluten in them, 
#eaten between the visits by patient divided by number of days between the visits (median 1.99 [interquartile range, 1.74 - 2.91] g/d). 
#LLoQ, lower limit of quantification; ULoQ, upper limit of quantification.

library(RColorBrewer)
library(ggpubr)

data<-read.table(file = "https://raw.githubusercontent.com/valeriia-dotsenko/Dotsenko-et-al-2020-Figures/main/data%20files/Patient's%20GIP%20measurements.txt",sep = '\t', header = TRUE)
#modify data for plotting
data$stool1m <- data$stool_1
data$stool1m[data$stool1m == '<LLoQ'] <- 100 # lower limit of GIP detection
data$stool1m[data$stool1m == '>ULoQ'] <- 5000 # upper limit of GIP detection
data$stool2m <- data$stool_2
data$stool2m[data$stool2m == '<LLoQ'] <- 100
data$stool2m[data$stool2m == '>ULoQ'] <- 5000
data %>% mutate(Visit=as.character(Visit),
                gluten_ingested=as.numeric(as.character(gluten_ingested)),
                stool1m=as.numeric(as.character(stool1m)),
                stool2m=as.numeric(as.character(stool2m)),Visit=as.numeric(as.character(Visit))) -> data
data$mean <- rowMeans(data[,8:9])
data$sd <- sqrt((data$stool1m-data$mean)^2 +(data$stool2m-data$mean)^2)
data$label <- paste(round(data$mean,1),"±",round(data$sd,1),sep = '')
data$label[data$label == '100±0'] <- '<LLoQ'
data$label[data$label == '5000±0'] <- '>ULoQ'
data$label[data$label == 'NA±NA'] <- NA

data$gluten_daily <- data$gluten_ingested/14 # 14 days in between visits

#build the heatmap
pal2=brewer.pal(8,"OrRd") #create a palette
heatmap <- ggplot(data, aes(as.character(Visit), Patient, fill= mean)) + 
  geom_tile()+
  geom_text(aes(label = label), size=6/.pt)+
  scale_fill_gradientn(colors = pal2, na.value="white")+
  scale_y_continuous(breaks = c(1:15))+
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_blank(), axis.line = element_line(colour = "black", size=1/.pt),
        axis.text.x = element_text(size=8, color = "black"),
        axis.text.y = element_text(size=8, color = "black"),
        axis.title.y = element_text(size=10, color = "black"),
        axis.title.x = element_text(size=10, color = "black"),
        legend.text = element_text( size = 6, color = "black"),plot.margin = unit(c(0,0,0,0), "cm"))+
  xlab("Visit")+border(color = "black", size = 1/.pt, linetype = NULL)


shading <- data.frame(min = c(0.5,2.4,7.75), #specify the areas on the plot to be highlighted
                       max = c(2.4,7.75,8.5),
                       col = c(0,1,0))
xplot<-ggplot()+
  geom_rect(data = shading,
            aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf,
                fill = factor(col), alpha = 0.1)) +
  scale_fill_manual(values = c("white", "coral1"))+
  geom_point(data=data, aes(x=Visit, y=gluten_daily))+xlim(c(-0.3,8.5))+
  theme(panel.background = element_rect(fill = "white"),
        #axis.line=element_line(colour = "black", size=1/.pt),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.5,0,0,0), "cm"), 
        axis.text.y = element_text(size=8, color = "black"),
        axis.title.y = element_text(size=6, color = "black")) +
  guides(fill = FALSE, alpha = FALSE)+
  ylab('Gluten consumption \ng/day')+
  border(color = "black", size = 1/.pt,linetype = NULL)+ 
  annotate(geom="text", x=5, y=6.1, label="Gluten challenge", size = 6/.pt)


# Arranging the plot
p <- ggarrange(xplot,heatmap, 
                ncol = 1, nrow = 2,  align = "v", 
                heights = c(0.5, 2), legend='none')
p

ggexport(p,filename = "Fig 2 ng GIP/g stool.tiff",
         width = 1800, # 6 inch
         height = 1500, # 5 inch
         res = 300) 

