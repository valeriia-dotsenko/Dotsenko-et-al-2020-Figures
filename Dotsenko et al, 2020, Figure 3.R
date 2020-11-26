#load all libraries
library(GEOquery)
library(readxl)
library(dplyr)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap) 
library(circlize)
library(GOplot)
library(ggpubr)
library(numform) ##
library(plyr)
library(reshape)
library(gridExtra)
library(cowplot)

#open all the data needed
#Download the UMI counts data
eList <- getGEOSuppFiles("GSE145358") #download supp.files from GSE145358 
tarArchive <- rownames(eList)[1]
#open supp.files as table
UMIdata <- read.table(tarArchive, header=T)

# open DE genes list 
#download it from paper supplemental file
URL <- 'https://www.cmghjournal.org/cms/10.1016/j.jcmgh.2020.07.010/attachment/66f15ab8-38a9-4871-98c0-23d649b4819d/mmc1.xlsx'
pf <- tempfile()
download.file(URL, pf, mode="wb")
DE_genes <- read_excel(path = pf) #open file

#===============================================Fig3. PGC/GFD comparisons=========================
#===================================================panelA heatmap===============================

#prepare the data for plotting
PGC_GFD_neg_s <- DE_genes[DE_genes$comparison=="PGC vs. GFD" & DE_genes$log2FoldChange <= -0.5,]
PGC_GFD_neg_s  <- PGC_GFD_neg_s [order(PGC_GFD_neg_s$log2FoldChange),][1:75,]
#arrange genes by log2FC
PGC_GFD_neg_s %>% arrange(log2FoldChange)->PGC_GFD_neg_s
dataA<-subset(UMIdata, ENSEMBL %in% PGC_GFD_neg_s$`ensembl id`) #data has some NAs in gene names, so lets add new genes names column
dataA <- merge(dataA[,-3],PGC_GFD_neg_s[,1:2], by.x = 'ENSEMBL', by.y ="ensembl id")


#define the order of rows if needed
dataA<-dataA[match(PGC_GFD_neg_s$`gene name`, dataA$`gene name`),]
#match the order of the colums
samples<-c("NORM27A","NORM28A","NORM29B","NORM5B","NORM7A","NORM8B",    
           "X101.002BL","X101.019BL","X101.023BL","X101.027BL","X101.029BL",
           "X101.035BL","X102.001BL","X102.007BL","X102.021BL","X102.026BL",
           "X103.001BL", "X103.006BL","X103.008BL","X103.016BL","X103.017BL", 
           "X101.002EOS", "X101.019EOS","X101.023EOS","X101.027EOS","X101.029EOS",
           "X101.035EOS", "X102.001EOS", "X102.007EOS", "X102.021EOS","X102.026EOS",
           "X103.001EOS","X103.006EOS", "X103.008EOS","X103.016EOS","X103.017EOS",
           "gene","ENSEMBL","gene name")
dataA<-dataA[,match(samples, colnames(dataA))]


#modify data as matrix
mdataA<-as.matrix(dataA[,c(-39,-38,-37)]) 
rownames(mdataA)<-as.character(dataA$`gene name`)
colnames(mdataA)<-paste(c(rep("DC",6),rep("GFD",15), rep("PGC",15)), c(1:36), sep="")
m.dataA_sc<-t(scale(t(mdataA)))

#create annotation 1(Sample type)
type1<-c(rep("DC",6),rep("GFD",15),rep("PGC",15))
type2 <- paste(c(rep("DC",6),rep("GFD",15), rep("PGC",15)), c(1:36), sep="")

top = HeatmapAnnotation(Samples = anno_block(gp = gpar(fill = c("cornflowerblue", "gold", "darkorange"), lty="blank"),
                                             labels = c("DC", "GFD", "PGC"), 
                                             labels_gp = gpar(col = "white", fontsize = 10)),
                        height = unit(0.3, "cm"))
#create annotation 2 (log2FC)
col_fun = colorRamp2(c(min(PGC_GFD_neg_s$log2FoldChange), 
                       mean(PGC_GFD_neg_s$log2FoldChange), max(PGC_GFD_neg_s$log2FoldChange)), 
                     c("#67001F", "#D6604D", "#FDDBC7"))
row <- rowAnnotation(log2FC = PGC_GFD_neg_s$log2FoldChange,
                     col = list(log2FC = col_fun),
                     show_annotation_name = F,
                     annotation_legend_param = list(log2FC = list(title_position = "topcenter",labels_gp = gpar(fontsize = 6),
                                                                  title_gp = gpar(fontsize = 6),direction = "horizontal",
                                                                  legend_width=unit(2.7,"cm"),grid_height = unit(0.2, "cm"))),
                     simple_anno_size = unit(2.5, "mm"))

#create color scheme for heatmap
newcolors <- c("#67001F",#"#B2182B","#D6604D", 
               "#F4A582",#"#FDDBC7", 
               "#F7F7F7", #white 
               "#D1E5F0", "#92C5DE", "#4393C3", 
               "#2166AC", "#053061")
ht_list1 = Heatmap(m.dataA_sc, name = "Row Z-Score of Expression", top_annotation = top, 
                   show_column_names = FALSE, col=newcolors,
                   column_title_gp = gpar(fontsize = 1, col="White"),
                   column_order= type2,
                   row_names_side = "left",
                   cluster_rows = FALSE, 
                   row_names_gp = gpar(fontsize = 6), 
                   column_split = factor(type1, levels = c("DC","GFD","PGC")),
                   column_gap = unit(0.5,"mm"),
                   heatmap_legend_param = list(
                     at = c(-1:4),
                     labels = c(-1:4),
                     title = "Row Z-Score",
                     labels_gp = gpar(fontsize = 6),
                     title_gp = gpar(fontsize = 6),direction = "horizontal",
                     legend_width=unit(2.1,"cm"),grid_height = unit(0.2, "cm"),
                     title_position = "topcenter"))+row

hm1 <- grid.grabExpr(draw(ht_list1, merge_legend = F, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
                          auto_adjust = FALSE))




#===================================================panel B heatmap===============================
#prepare the data for plotting
PGC_GFD_pos_s <- DE_genes[DE_genes$comparison=="PGC vs. GFD" & DE_genes$log2FoldChange >= 0.5,]
PGC_GFD_pos_s  <- PGC_GFD_pos_s [order(-PGC_GFD_pos_s$log2FoldChange),][1:75,]
#arrange genes by log2FC
PGC_GFD_pos_s %>% arrange(desc(log2FoldChange))->PGC_GFD_pos_s
dataB<-subset(UMIdata, ENSEMBL %in% PGC_GFD_pos_s$`ensembl id`) #data has some NAs in gene names, so lets add new genes names column
dataB <- merge(dataB[,-3],PGC_GFD_pos_s[,1:2], by.x = 'ENSEMBL', by.y ="ensembl id")

#define the order of rows if needed
dataB<-na.omit(dataB[match(PGC_GFD_pos_s$`gene name`, dataB$`gene name`),])
#match the order of the colums
samples<-c("NORM27A","NORM28A","NORM29B","NORM5B","NORM7A","NORM8B",    
           "X101.002BL","X101.019BL","X101.023BL","X101.027BL","X101.029BL",
           "X101.035BL","X102.001BL","X102.007BL","X102.021BL","X102.026BL",
           "X103.001BL", "X103.006BL","X103.008BL","X103.016BL","X103.017BL", 
           "X101.002EOS", "X101.019EOS","X101.023EOS","X101.027EOS","X101.029EOS",
           "X101.035EOS", "X102.001EOS", "X102.007EOS", "X102.021EOS","X102.026EOS",
           "X103.001EOS","X103.006EOS", "X103.008EOS","X103.016EOS","X103.017EOS",
           "gene","ENSEMBL","gene name")
dataB<-dataB[,match(samples, colnames(dataB))]


#modify data as matrix
mdataB<-as.matrix(dataB[,c(-39,-38,-37)]) 
rownames(mdataB)<-as.character(dataB$`gene name`)
colnames(mdataB)<-paste(c(rep("DC",6),rep("GFD",15), rep("PGC",15)), c(1:36), sep="")
#m.data <- as.matrix(scale(mdata)) 
m.data_sc2<-t(scale(t(mdataB)))


top2 = HeatmapAnnotation(Samples = anno_block(gp = gpar(fill = c("cornflowerblue", "gold", "darkorange"), lty="blank"),
                                              labels = c("DC", "GFD", "PGC"), 
                                              labels_gp = gpar(col = "white", fontsize = 10)),
                                              height = unit(0.3, "cm"),show_annotation_name = F, show_legend = FALSE)


#create annotation 2 (log2FC)
col_fun2 = colorRamp2(c(min(PGC_GFD_pos_s$log2FoldChange), 
                        mean(PGC_GFD_pos_s$log2FoldChange), max(PGC_GFD_pos_s$log2FoldChange)), 
                      c("#D1E5F0", "#4393C3", "#053061"))
row2 <- rowAnnotation(log2FC = PGC_GFD_pos_s$log2FoldChange,
                      col = list(log2FC = col_fun2),
                      show_annotation_name = F,
                      annotation_legend_param = list(log2FC = list(labels_gp = gpar(fontsize = 6),title_position = "topcenter",
                                                                   title_gp = gpar(fontsize = 6),direction = "horizontal",
                                                                   legend_width=unit(2.7,"cm"),grid_height = unit(0.2, "cm"))),
                      simple_anno_size = unit(2.5, "mm"))

#create color scheme for heatmap
newcolors2 <- c("#67001F",#"#B2182B","#D6604D", 
                "#F4A582",#"#FDDBC7", 
                "#F7F7F7", #white 
                "#D1E5F0", "#92C5DE", "#4393C3", 
                "#2166AC", "#053061")
ht_list2 = Heatmap(m.data_sc2, name = "Row Z-Score of Expression", top_annotation = top2, 
                   show_column_names = FALSE, col=newcolors2,
                   cluster_columns = F,
                   column_title_gp = gpar(fontsize = 1, col="White"),
                   column_order= type2,
                   row_names_side = "left",
                   cluster_rows = FALSE, 
                   row_names_gp = gpar(fontsize = 6), 
                   column_split = factor(type1, levels = c("DC","GFD","PGC")),
                   column_gap = unit(0.5,"mm"),
                   heatmap_legend_param = list(
                     at = c(-1:4),
                     labels = c(-1:4),
                     title = "Row Z-Score",
                     labels_gp = gpar(fontsize = 6),
                     title_gp = gpar(fontsize = 6),direction = "horizontal",
                     #legend_height= unit(, "cm"),
                     legend_width=unit(2.1,"cm"),grid_height = unit(0.2, "cm"),
                     title_position = "topcenter",outline = TRUE)
                   )+row2

hm2 <- grid.grabExpr(draw(ht_list2, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
                          auto_adjust = FALSE))


##=========================================================panel C barplot=========================================

# Load the dataset
# GO annotation results obtained from https://tools.dice-database.org/GOnet/
GOres <- read.table("https://raw.githubusercontent.com/valeriia-dotsenko/Dotsenko-et-al-2020-Figures/main/data%20files/GOnet%20results_PGC%20vs%20GFD%20DE%20genes.txt", sep="\t", header = T)  
genes <- DE_genes[DE_genes$comparison=="PGC vs. GFD",][,c("gene name","log2FoldChange")]
colnames(genes) <- c("ID","logFC") #circle_dat requires a data frame with columns for 'ID', 'logFC'
circ <- circle_dat(GOres, genes)
circ$direction = NA
circ$direction[circ$logFC>0] = 'up'
circ$direction[circ$logFC<0] = 'down'
circ$`-log10(adj.p-value)` <- (-log10(circ$adj_pval))

IDs <- c("GO:0002376", #specify GO terms to be shown on the plot
         "GO:0050896",
         "GO:0043900",
         "GO:0071345",
         "GO:0060337",
         "GO:0034341",
         "GO:0008152",
         "GO:0009056",
         "GO:0019752")

circ2<-subset(circ, ID %in% IDs)
circ2 %>% arrange(desc(`-log10(adj.p-value)`)) -> circ2
circ2$group <- rep(1, dim(circ2)[1])

# prepare data for plotting
plotting_df <-   circ2 %>% group_by(term, direction, `-log10(adj.p-value)`) %>% 
  dplyr::summarise(Freq = n()) %>% mutate(Freq = if_else(direction == "down", -Freq, Freq))
plotting_df <- na.omit(plotting_df)
plotting_df$term1 <- c("- Carboxylic acid \n metabolic process", "- Carboxylic acid \n metabolic process","- Catabolic process","- Catabolic process",                     
                       "- Cellular response to \n cytokine stimulus", "- Cellular response to \n cytokine stimulus","- Immune system \n process", 
                       "- Immune system \n process","- Metabolic process" ,                    
                       "- Metabolic process","- Regulation of multi- \n organism process",   "- Regulation of multi- \n organism process" , 
                       "- Response to \n interferon-gamma", "- Response to \n stimulus" , "- Response to \n stimulus",
                       "- Type I interferon \n signaling pathway")
circ2 <- merge(circ2, plotting_df[,c(1,5)], by.x = "term", by.y="term")
circ2 %>% group_by(term1, `-log10(adj.p-value)`, direction) %>% filter(direction == 'up') %>% dplyr::summarize(counts=mean(count),p=n()) -> circ2.counts


# side-by-side plot
p <- plotting_df %>% 
  ggplot(aes(x = reorder(term1,`-log10(adj.p-value)`), y = Freq, group = direction, fill = direction)) +
  geom_bar(stat = "identity", width = 0.75) +  coord_flip()+
  ylab("Gene counts")+
  scale_fill_manual(values= c("#B2182B","#2166AC"),
                    name="log2FC",
                    breaks=c("down", "up"),
                    labels=c("down", "up"))+
  geom_text(aes(label=abs(Freq),y=Freq/1.6, x=term1), vjust = 0.5, size=6/.pt)

f <- sum(abs(range(plotting_df$Freq)))/diff(range(plotting_df$`-log10(adj.p-value)`))*0.5 #factor for top y-axis (coordinates are flipped) scaling
#add scatter plot on top of barplot
bp <-  p + 
  geom_line(data=circ2, aes(x=reorder(term1,`-log10(adj.p-value)`),
                           y=f*(`-log10(adj.p-value)`-1.2*median(`-log10(adj.p-value)`)), group=group),
            color="#4daf4a",size = 0.5) +
  geom_point(data=circ2, aes(x=reorder(term1,`-log10(adj.p-value)`),
                            y=f*(`-log10(adj.p-value)`-1.2*median(`-log10(adj.p-value)`)), group=group),
             shape=21, color="black", fill="#69b3a2", size=2) +
  theme_bw()+ 
  geom_text(data=circ2.counts, aes(x = reorder(term1,`-log10(adj.p-value)`), y=135,label = counts),size=6/.pt) +
  scale_y_continuous(sec.axis = sec_axis(~.*1/f+1.2*median(circ$`-log10(adj.p-value)`), name = "-log10(adj.pvalue)"))+
  theme(plot.title = element_blank(), panel.background = element_rect(fill = 'White'), axis.line = element_line(colour = "black", size=1/.pt),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title.x.top = element_text(color='#4daf4a',size=10), 
      axis.text.x.top = element_text(color='#4daf4a',size=8),
      axis.text.x.bottom = element_text(size=8, color = "black"),
      axis.title.x.bottom = element_text(size=10, color = "black"), 
      axis.text.y.left = element_text(size=8,color = "black"),
      #legend.position = "bottom",
      legend.title = element_text( size = 6),
      legend.text = element_text( size = 6),
      axis.title.y = element_blank(),#text(size=16),
      strip.text.x = element_text(size = 6))+
  theme(plot.margin = unit(c(0,0.5,0,0), "cm"),legend.margin=margin(0,0,0,0),#legend.position="bottom", 
        legend.position = c(0.15, 0.9),
        legend.key.height = unit(0.1, "mm"),legend.key.width = unit(1, "mm"),legend.background = element_blank(),
        legend.box.margin=margin(1,2,1,1),
        legend.box.background = element_rect(colour = "black", size=1/.pt))

##========================================================D panel violin plots====================================
var <- read.delim('https://raw.githubusercontent.com/valeriia-dotsenko/Dotsenko-et-al-2020-Figures/main/data%20files/Histomorphological%20parameters.txt',sep = '\t', header=T) #open file with histomorphological parameters of Samples

dataA_melt <- melt(dataA[,-37:-38], id.vars=c("gene name"))
dataA_melt <- merge(dataA_melt, var, by.x="variable", by.y="SAMPLE")
dataA_mean <- ddply(dataA_melt, .(`gene name`,Type), summarize,  `mean(value)`=mean(value))
dataA_mean$`log10_value` <- log10(dataA_mean$`mean(value)`)
dataA_mean$`log10_value`[dataA_mean$log10_value=="-Inf"] <- (-0.01)
anno_df3 <- compare_means(log10_value~Type, data = dataA_mean, method="t.test",paired=F) 
anno_df3 %>% 
  mutate(p_new = ifelse(p > 0.01, c(paste("italic('P')~`=", f_num(p,2), "`")), p))%>% 
  mutate(p_new = ifelse(p < 0.01, c(paste("italic('P')~`=", f_num(p,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(p < 0.001, c(paste("italic('P')~`", "<.001", "`")), p_new))->anno_df3



vp1 <- ggplot(dataA_mean, aes(x = factor(Type), y = log10_value, fill=Type)) +
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
              y_position = c(3.4, 3.7,4.0), xmin=anno_df3$group1, xmax=anno_df3$group2,textsize=6/.pt, manual= F, parse=T, size=0.3)

vp1 <- vp1+
  theme(plot.margin = unit(c(0.2,0,0.6,0.2), "cm"),legend.position = "none",
        axis.text.x = element_text(size=8, color = "black"),
        axis.text.y= element_text(size=8, color = "black"))
vp1

#=======================================plot2===============================================
dataB_melt <- melt(dataB[,-37:-38], id.vars=c("gene name"))
dataB_melt <- merge(dataB_melt, var, by.x="variable", by.y="SAMPLE")
dataB_mean <- ddply(dataB_melt, .(`gene name`,Type), summarize,  `mean(value)`=mean(value))
dataB_mean$log10_value <- log10(dataB_mean$`mean(value)`)
dataB_mean$log10_value[dataB_mean$log10_value=="-Inf"] <- (-0.01)

anno_df4 <- compare_means(log10_value~Type, data = dataB_mean, method="t.test",paired=F) 
anno_df4 %>% 
  mutate(p_new = ifelse(p > 0.01, c(paste("italic('P')~`=", f_num(p,2), "`")), p))%>% 
  mutate(p_new = ifelse(p < 0.01, c(paste("italic('P')~`=", f_num(p,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(p < 0.001, c(paste("italic('P')~`", "<.001", "`")), p_new))->anno_df4

vp2 <- ggplot(dataB_mean, aes(x = factor(Type), y = log10_value, fill=Type)) +
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
              y_position = c(3.4, 3.7,4.0), xmin=anno_df4$group1, xmax=anno_df4$group2,textsize=6/.pt, manual= F, parse=T, size=0.3)



vp2 <- vp2+
  theme(plot.margin = unit(c(0.2,0.2,0.6,0), "cm"), legend.position = "none", 
        axis.text.x = element_text(size=8, color = "black"))
plots <- ggarrange(vp1,vp2,ncol = 2, nrow = 1, #font.label = list(size = 16),
                   common.legend = TRUE, legend="none")
plots <- annotate_figure(plots,
                         left = text_grob("log10(Expression, UMI counts)", size = 10, rot=90))

##=========================================================Arrange all panels=================================
# Arrange plots using arrangeGrob
# returns a gtable (gt)
gt <- arrangeGrob(hm1,                               # bar plot spaning two columns
                  hm2, bp, plots,                          # box plot and scatter plot
                  ncol = 2, nrow = 4, 
                  layout_matrix = rbind(c(1,2),c(1,2),c(3,4)),
                  padding=unit(0.5, 'line')
                  #heights=c(0.5,0.5,0.5,0.5)
                  )
# Add labels to the arranged plots
p<-as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C","D"), size = 16,
                  x = c(0, 0.5, 0,0.5), y = c(1, 1, 0.35,0.35))

ggexport(p,filename = "Fig.3 PGC_GFD with violin plots.tiff",
         width = 2100, # 7 inch
         height = 2700, # 9 inch
         res = 300) #max 2100 X 2700

