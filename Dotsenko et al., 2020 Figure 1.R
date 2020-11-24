# Figure 1. from "Genome-Wide Transcriptomic Analysis of Intestinal Mucosa in Celiac Disease Patients on a Gluten-Free Diet and Postgluten Challenge"
#(https://doi.org/10.1016/j.jcmgh.2020.07.010) 
#Histomorphometric and molecular histomorphometric data of the DC group, celiac disease patients on a GFD, and celiac disease patients PGC. 
#(A) VH:CrD measurements. (B) Counts of IELs per 100 enterocytes (EC). (C) PCA plot for all samples exploring biological differences between the DC (blue),
#GFD (yellow), and PGC (orange) groups based on the gene expression profile with bigger spheres depicting the center of a distribution. 
#(D) Violin plots showing the log-transformed expression of differentially expressed genes #detected in PGC vs DC groups comparisons 
#(768 downregulated in the left panel and 580 upregulated genes in the right panel) in each group. 
#Three lines (from the bottom to the top) in each violin plot show the location of the lower quartile, the median, and the upper quartile, respectively. 
#The blue line is fitted regression model, demonstrating the linear character of mean of the log-transformed expression changes. 
#Student’s t test was used for comparing the mean values.

#load libraries needed 
library(stringr)
library(ggpubr)
#library(devtools)
#install_github('trinker/numform')
library(numform)
library(factoextra)
library(GEOquery)
library(reshape)
library(plyr)
library(gridExtra)
library(cowplot)





#===========================================Fig 1.; panels A&B; histomorphological parameters of Samples plots==========================================

var <- read.delim('https://raw.githubusercontent.com/valeriia-dotsenko/Dotsenko-et-al-2020-Figures/main/data%20files/Histomorphological%20parameters.txt',sep = '\t', header=T) #open file with histomorphological parameters of Samples
var$Number <- str_replace_all(as.character(var$SAMPLE),'[EOSNORMBL"]', "") #extract the Number of the Sample, to match paired Samples on the plot 
text<-"VH:CrD"
text2<-"IELs per 100 EC"

#1. Specify the stat comparison groups
P <- "P"
anno_df1 <- compare_means(VH.CrD~Type, data = var, method="t.test",paired=F) #calculate t-test for means as non-paired samples
anno_df1[1,4] <- t.test(var[var$Type=="GFD",]$VH.CrD,var[var$Type=="PGC",]$VH.CrD,paired = T)$p.value #in case PGC vs. GFD comparisons, samples are paired 
anno_df2 <- compare_means(IELsper100EC~Type, data = var, method="t.test",paired=F) 
anno_df2[1,4] <- t.test(var[var$Type=="GFD",]$VH.CrD,var[var$Type=="PGC",]$VH.CrD,paired = T)$p.value
anno_df <- rbind(anno_df1,anno_df2)
anno_df %>%  #format the p-value  
  mutate(p_new = ifelse(p > 0.01, c(paste("italic('P')~`=", f_num(p,2), "`")), p))%>% 
  mutate(p_new = ifelse(p < 0.01, c(paste("italic('P')~`=", f_num(p,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(p < 0.001, c(paste("italic('P')~`", "<.001", "`")), p_new))->anno_df

#2. Plot the boxplot for VH:CrD
p <- ggplot(var, aes(x=Type, y=VH.CrD, group=Type)) + 
  geom_boxplot(outlier.colour=NA)+
  geom_line(aes(x = Type, group = Number),linetype="longdash",
            alpha=0.5,size=1/.pt)+
  geom_signif(annotations = anno_df[anno_df$.y.=="VH.CrD",]$p_new,
              y_position = c(3.5, 4,4.5), xmin=anno_df[anno_df$.y.=="VH.CrD",]$group1, xmax=anno_df[anno_df$.y.=="VH.CrD",]$group2,textsize=6/.pt, manual= F, parse=T, size=0.3)+
  theme(axis.title.x = element_blank(),legend.position="none",
        panel.background = element_rect(fill = 'White'), axis.line = element_line(colour = "black", size=1/.pt),
        plot.margin = margin(0, 0, 0.1, 0.4, "cm")
        )+
  labs(y = text)+
  scale_x_discrete(limits = c("DC","GFD","PGC"))+
  scale_y_continuous(expand = c(.2,0))
p <- p + geom_point(aes(x = Type), size=0.4)+
  theme(axis.text.x = element_text(size=8, color = "black"),
        axis.text.y= element_text(size=8, color = "black"),
        axis.title.y = element_text(size=10, color = "black"))
p

#3. Plot the boxplot for IELs per 100 EC
p2 <- ggplot(var, aes(x=Type, y=IELsper100EC, group=Type)) + 
  geom_boxplot(outlier.colour=NA)+
  geom_line(aes(x = Type, group = Number),linetype="dashed",
            alpha=0.5,size=1/.pt)+
  geom_signif(annotations = anno_df[anno_df$.y.=="IELsper100EC",]$p_new,
              y_position = c(75, 85,95), xmin=anno_df[anno_df$.y.=="IELsper100EC",]$group1, xmax=anno_df[anno_df$.y.=="IELsper100EC",]$group2,textsize=6/.pt, manual= F, parse=T, size=0.3)+
  #geom_text_repel(parse = TRUE)
  theme(axis.title.x = element_blank(),legend.position="none",
        panel.background = element_rect(fill = 'White'), axis.line = element_line(colour = "black", size=1/.pt),
        plot.margin = margin(0, 0, 0.1, 0.2, "cm")
        )+
  labs(y = text2)+
  scale_x_discrete(limits = c("DC","GFD","PGC"))+
  scale_y_continuous(expand = c(.2,0))
p2 <- p2 + geom_point(aes(x = Type),size=0.4)

p2<-p2+
  theme(axis.text.x = element_text(size=8, color = "black"),
        axis.text.y= element_text(size=8, color = "black"),
        axis.title.y = element_text(size=10, color = "black"))
p2

#=======================================================Fig 1.; panel C; PCA analysis ==========================================
#open the transformed raw counts
trcounts<-read.table(file = 'https://raw.githubusercontent.com/valeriia-dotsenko/Dotsenko-et-al-2020-Figures/main/data%20files/DESeq2_vsd_transformed_counts.tsv',sep = '\t', header = TRUE) 

alldata_sub_t<-t(trcounts)
rownames(alldata_sub_t) <- colnames(trcounts)
samples<-colnames(trcounts)
var <- var[match(samples, var$SAMPLE),]
#Compute PCA
res.pca <- prcomp(alldata_sub_t, center =T,scale = F)
pca <- fviz_pca_ind(res.pca, geom.ind = "point", pointshape = 21,
                    pointsize = 2, 
                    fill.ind = var$Type, 
                    palette = c( "#6495ed","#ffd700", "#ffa500"))+
 theme_classic()+
 theme(plot.title = element_blank(), panel.background = element_rect(fill = 'White'), axis.line = element_line(colour = "black", size=1/.pt),
       axis.text.x = element_text(size=8, color = "black"),
       axis.text.y = element_text(size=8, color = "black"),
       axis.title.y = element_text(size=10, color = "black"),
       axis.title.x = element_text(size=10, color = "black"),
       
       legend.text = element_text( size = 6, color = "black"),
       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       legend.title = element_blank(), 
       legend.position = c(0.85, 0.9),
       legend.margin=margin(0,0,0,0),
       #legend.key.size = unit(3,"mm"),
       legend.key.height = unit(0.5, "mm"),legend.key.width = unit(1, "mm"),
       legend.background = element_blank(),
       legend.box.margin=margin(-5,0,0,0),
       legend.box.background = element_rect(colour = "black", size=1/.pt),
       plot.margin = margin(0, 0, 0, 0.1, "cm"))
pca

##========================================================Fig 1.; panel D;  violin plots====================================
#Download the UMI counts data
#download supp.files from GSE145358 paper https://doi.org/10.1016/j.jcmgh.2020.07.010
eList <- getGEOSuppFiles("GSE145358")
tarArchive <- rownames(eList)[1]
#open supp.files as table
data <- read.table(tarArchive, header=T)
PGC_DC <- read.table(file = 'https://raw.githubusercontent.com/valeriia-dotsenko/Dotsenko-et-al-2020-Figures/main/data%20files/challenge%20VS%20control_ALL_Differentially_regulated_genes.txt',sep = '\t', header = TRUE) 
  
#plot 1 for downregulated genes
PGC_DC_sdown <- na.omit(PGC_DC[PGC_DC$padj<.05 & PGC_DC$log2FoldChange<=(-0.5),]) 
datasub<-subset(data, SYMBOL %in% PGC_DC_sdown$gene)
#match the order of the columns
samples<-c("NORM27A","NORM28A","NORM29B","NORM5B","NORM7A","NORM8B", 
           "X101.029BL","X103.016BL","X101.002BL","X103.008BL",
           "X101.019BL","X101.023BL","X101.027BL",
           "X101.035BL","X102.001BL","X102.007BL","X102.021BL","X102.026BL",
           "X103.001BL", "X103.006BL","X103.017BL", 
           "X101.029EOS","X103.016EOS","X101.002EOS", "X103.008EOS",
           "X101.019EOS","X101.023EOS","X101.027EOS",
           "X101.035EOS", "X102.001EOS", "X102.007EOS", "X102.021EOS","X102.026EOS",
           "X103.001EOS","X103.006EOS", "X103.017EOS",
           "gene","ENSEMBL","SYMBOL")
datasub<-datasub[,match(samples, colnames(datasub))]
datasub_melt <- melt(datasub[,-37:-38], id.vars=c("SYMBOL"))
datasub_melt <- merge(datasub_melt, var, by.x="variable", by.y="SAMPLE")
datasub_mean <- ddply(datasub_melt, .(SYMBOL,Type), summarize,  `mean(value)`=mean(value))
datasub_mean$`log10_value` <- log10(datasub_mean$`mean(value)`)
datasub_mean$`log10_value`[datasub_mean$log10_value=="-Inf"] <- (-0.01)
anno_df3 <- compare_means(log10_value~Type, data = datasub_mean, method="t.test",paired=F) 
anno_df3 %>% 
  mutate(p_new = ifelse(p > 0.01, c(paste("italic('P')~`=", f_num(p,2), "`")), p))%>% 
  mutate(p_new = ifelse(p < 0.01, c(paste("italic('P')~`=", f_num(p,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(p < 0.001, c(paste("italic('P')~`", "<.001", "`")), p_new))->anno_df3



vp1 <- ggplot(datasub_mean, aes(x = factor(Type), y = log10_value, fill=Type)) +
  geom_violin() +
  geom_boxplot(width=0.1,outlier.shape=NA)+
  stat_smooth(aes(y = log10_value, group=1),
              method= lm,formula = y ~ poly(x, 2), se=F, size=2/.pt) + #tweak the signifance level suitable for your study
  scale_fill_manual(values=c("#6495ed", "#ffd700", "#ffa500"))+
  geom_jitter(size=0.5, fill="black", alpha=0.1)+
  theme_classic()+
  theme(axis.title.x = element_blank(), legend.title = element_blank(), axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size=1/.pt))+
  scale_y_continuous(breaks=c(0,1,2,3),limits=c(-0.1, 4))+
  annotate(geom="text", x=2, y=-0.1, label="downregulated genes",size=6/.pt)+
  geom_signif(annotations = anno_df3$p_new,
              y_position = c(3.4, 3.7,4.0), xmin=anno_df3$group1, xmax=anno_df3$group2,textsize=6/.pt, manual= F, parse = T, size=0.3)

vp1 <- vp1+
  theme(plot.margin = unit(c(0,0,0.5,0.3), "cm"),legend.position = "none",
        axis.text.x = element_text(size=8, color = "black"),
        axis.text.y= element_text(size=8, color = "black"))
vp1
#=======================================plot2===============================================
#plot 1 for upregulated genes
PGC_DC_sup <- na.omit(PGC_DC[PGC_DC$padj<.05 & PGC_DC$log2FoldChange>=(0.5),]) 
datasub2<-subset(data, SYMBOL %in% PGC_DC_sup$gene)
datasub2<-datasub2[,match(samples, colnames(datasub2))]
datasub2_melt <- melt(datasub2[,-37:-38], id.vars=c("SYMBOL"))
datasub2_melt <- merge(datasub2_melt, var, by.x="variable", by.y="SAMPLE")
datasub2_mean <- ddply(datasub2_melt, .(SYMBOL,Type), summarize,  `mean(value)`=mean(value))
datasub2_mean$log10_value <- log10(datasub2_mean$`mean(value)`)
anno_df4 <- compare_means(log10_value~Type, data = datasub2_mean, method="t.test",paired=F) 
anno_df4 %>% 
  mutate(p_new = ifelse(p > 0.01, c(paste("italic('P')~`=", f_num(p,2), "`")), p))%>% 
  mutate(p_new = ifelse(p < 0.01, c(paste("italic('P')~`=", f_num(p,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(p < 0.001, c(paste("italic('P')~`", "<.001", "`")), p_new))->anno_df4

vp2 <- ggplot(datasub2_mean, aes(x = factor(Type), y = log10_value, fill=Type)) +
   geom_violin() +
  geom_boxplot(width=0.1,outlier.shape=NA)+
    stat_smooth(aes(y = log10_value, group=1),
              method= lm,formula = y ~ poly(x, 2), se=F, size=2/.pt) +
  scale_fill_manual(values=c("#6495ed", "#ffd700", "#ffa500"))+
  geom_jitter(size=0.5, fill="black", alpha=0.1)+
  scale_y_continuous(breaks=c(0,1,2,3),limits=c(-0.1, 4) #,position = "right", 
                     )+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y=element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_line(colour = "black", size=1/.pt))+
  annotate(geom="text", x=2, y=-0.1, label="upregulated genes",size=6/.pt)+
  geom_signif(annotations = anno_df4$p_new,
              y_position = c(3.4, 3.7,4.0), xmin=anno_df4$group1, xmax=anno_df4$group2,textsize=6/.pt, manual= F, parse = T, size=0.3)



vp2 <- vp2+
  theme(plot.margin = unit(c(0,0,0.5,0), "cm"), legend.position = "none", 
        axis.text.x = element_text(size=8, color = "black"))
vp2

#make one plot
plots <- ggarrange(vp1,vp2,ncol = 2, nrow = 1, #font.label = list(size = 16),
                   common.legend = TRUE, legend="none")
plots <- annotate_figure(plots,
                         left = text_grob("log10(Expression, UMI counts)", size = 10, rot=90))
plots

#================================================================Arrange all 4 plots in one========================================== 
# Arrange plots using arrangeGrob; returns a gtable (gt)
gt <- arrangeGrob(p,
                  p2,
                  pca,
                  plots,                            
                  ncol = 2, nrow = 2, 
                  layout_matrix = rbind(c(1,2), c(3,4)))
# Add labels to the arranged plots
p3 <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C","D"), size = 16, #add plot labels
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5))
p3
ggexport(p3,filename = "Fig.1.tiff",
         width = 1500, # 5 inch
         height = 1500, # 5 inch
         res = 300) 
