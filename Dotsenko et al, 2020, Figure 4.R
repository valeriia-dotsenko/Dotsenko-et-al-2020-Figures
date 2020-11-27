# Figure 4. from "Genome-Wide Transcriptomic Analysis of Intestinal Mucosa in Celiac Disease Patients on a Gluten-Free Diet and Postgluten Challenge" (https://doi.org/10.1016/j.jcmgh.2020.07.010) 
#Genes affected by gluten challenge shown as heatmaps organized to Gene Ontology (GO) terms.

# load packages
library(GEOquery)
library(readxl)
library(dplyr)
library(ggpubr)
library(GOplot)
library(circlize)
library(ComplexHeatmap) #BiocManager::install("ComplexHeatmap")
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



#====================================================panelA=========================================
##===================================================heatmap with GO===============================

#prepare the data for plotting
PGC_GFD_s <- DE_genes[DE_genes$comparison=="PGC vs. GFD",]
dataA<-subset(UMIdata, ENSEMBL %in% PGC_GFD_s$`ensembl id`) #data has some NAs in gene names, so lets add new genes names column
dataA <- merge(dataA[,-3],PGC_GFD_s[,1:2], by.x = 'ENSEMBL', by.y ="ensembl id")

#open GO data
# GO annotation results obtained from https://tools.dice-database.org/GOnet/
GOres<-read.table("https://raw.githubusercontent.com/valeriia-dotsenko/Dotsenko-et-al-2020-Figures/main/data%20files/GOnet%20results_PGC%20vs%20GFD%20DE%20genes.txt", sep="\t", header = T) 
genes <- DE_genes[DE_genes$comparison=="PGC vs. GFD",][,c("gene name","log2FoldChange")]
colnames(genes) <- c("ID","logFC") #circle_dat requires a data frame with columns for 'ID', 'logFC'
circ <- circle_dat(GOres, genes)
not <- as.data.frame(cbind(category= rep("no", 81),ID=rep("no enrichment", 81), term=rep("no enrichment", 81),
                           count=rep("0", 81), genes= as.character(PGC_GFD_s[!(PGC_GFD_s$`gene name` %in% circ$genes),]$`gene name`),
                           logFC=PGC_GFD_s[!(PGC_GFD_s$`gene name` %in% circ$genes),]$log2FoldChange,
                           adj_pval=PGC_GFD_s[!(PGC_GFD_s$`gene name` %in% circ$genes),]$padj,
                           zscore=rep("0", 81)))

circ <- rbind(circ, not)
IDs <- c("GO:0002376",
         "GO:0050896",
         "GO:0043900",
         "GO:0071345",
         "GO:0060337",
         "GO:0034341",
         "GO:0008152",
         "GO:0009056",
        "GO:0019752")

circ2 <- circ[circ$ID %in% IDs,]
circ3 <- circ2[order(match(circ2$ID,IDs)),]


data <- subset(dataA, `gene name` %in% circ3$genes)
circ3 <- circ3[circ3$genes %in% intersect(unique(circ3$genes),data$`gene name`),]
#define the order of rows if needed
data<-data[match(circ3$gene, data$`gene name`),]
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
           "gene","ENSEMBL","gene name")
data<-data[,match(samples, colnames(data))]


#modify data as matrix
mdata<-as.matrix(data[,c(-39,-38,-37)]) 
rownames(mdata)<-as.character(data$SYMBOL)
colnames(mdata)<-paste(c(rep("DC",6),rep("GFD",15),rep("PGC",15)), c(1:36), sep="")
m.data_sc<-t(scale(t(mdata)))

#create annotation 1(Sample type)
type1<-c(rep("DC",6),rep("GFD",15),rep("PGC",15))
type2 <- circ3$ID
type3 <- paste(c(rep("DC",6),rep("GFD",15),rep("PGC",15)), c(1:36), sep="")
top = HeatmapAnnotation(Samples = anno_block(gp = gpar(fill = c("cornflowerblue", "gold", "darkorange"), lty="blank"),
                                             labels = c("DC", "GFD", "PGC"), 
                                             labels_gp = gpar(col = "white", fontsize = 10)),
                        height = unit(0.3, "cm"))

#create annotation 2 (log2FC)
col_fun = colorRamp2(c(-1,-0.5,0,0.5,1), 
                     c("#cc6677", "#f4e0e3","#ffffff", "#b9dad4","#51a395"))
row <- rowAnnotation(Log2FC = as.numeric(as.character(circ3$logFC)),
                     col = list(Log2FC = col_fun),
                     show_annotation_name = F,
                     annotation_legend_param = list(Log2FC = list(title_position = "topcenter",labels_gp = gpar(fontsize = 9),
                                                                  title_gp = gpar(fontsize = 10),direction = "horizontal",
                                                                  legend_width=unit(2.7,"cm"),grid_height = unit(0.2, "cm"))),
                     simple_anno_size = unit(2.5, "mm"))

row2 = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#3b565a", "#729ca3","#36b8ea", "#52a68a","#cd8162","#ff6666"), lty="blank"),width = unit(0.3, "cm")))


#create color scheme for heatmap
newcolors <- c("#67001F","#F4A582","#F7F7F7", #white 
               "#D1E5F0", "#92C5DE", "#2166AC", "#053061")
ht_list = Heatmap(m.data_sc, name = "Row Z-Score of Expression", top_annotation = top,right_annotation = row,
                  left_annotation = row2,
                  show_column_names = F, col=newcolors,
                  column_order= type3,
                  column_title_gp = gpar(fontsize = 1, col="White"),
                  cluster_rows = T,
                  show_row_dend = FALSE,
                  column_split = factor(type1, levels = c("DC","GFD", "PGC")),
                  column_gap = unit(0.5,"mm"),
                  row_split = factor(type2, levels=unique(circ3$ID)),
                  cluster_row_slices = FALSE, 
                  show_row_names = F,
                  row_title_gp = gpar(fontsize = 6, direction = "horizontal"),
                  row_title_rot = 0,
                  heatmap_legend_param = list(
                    at = c(-1:4),
                    labels = c(-1:4),
                    title = "Row Z-Score",
                    labels_gp = gpar(fontsize = 9),
                    title_gp = gpar(fontsize = 10),direction = "horizontal",
                    legend_width=unit(2.1,"cm"),grid_height = unit(0.2, "cm"),
                    title_position = "topcenter"))

hm1 <- grid.grabExpr(draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", row_gap = unit(1, "mm"),auto_adjust = FALSE))
draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", row_gap = unit(0.5, "mm"),auto_adjust = FALSE)

#===========================================================panel B====================================================
#==========================================================terms table==================================================
termdf<-GOres[,c("ID","Term")]
termdf <- termdf[termdf$ID %in% IDs,]
termdf<- termdf[order(match(termdf$ID,IDs)),]

stable.p <- ggtexttable(termdf, rows = NULL,theme = ttheme(
  colnames.style = colnames_style(color = "black", size=10),
  tbody.style = tbody_style(color = "black", size=8)
))
stable.p

##=========================================================Arrange all panels=================================
gt <- arrangeGrob(hm1,                               
                  stable.p, 
                  ncol = 2, nrow = 2, 
                  layout_matrix = rbind(c(1,2), c(1,3)),
                  padding=unit(0.5, 'line'),
                  top=c(0,0,0,0,0,0),
                  bottom=c(0,0,0,0,0,0))
# Add labels to the arranged plots
p <- as_ggplot(gt) +                                # transform to a ggplot
  draw_plot_label(label = c("A", "B"), size = 16,
                  x = c(0, 0.5), y = c(1, 1))

ggexport(p,filename = "Fig.4 PGC_GFD_GO.tiff",
         width = 2100, # 7 inch
         height = 2700, # 9 inch
         res = 300) 
