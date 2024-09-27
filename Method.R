library(Seurat)
library(cetcolor)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
#####
#UMAP
#load data
all <- readRDS(file = "~/all.rds")
#add metadata
Idents(all) <- all$seurat_clusters
new.cluster.ids <- c("Monocytes", "Monocytes", "Monocytes", "Microglia", "T", "NK", 
                     "T", "MdC", "T", "cDC", "MDSC-like", "pDC", "B", "T", "Macrophage", "Oligodendrocyte")
names(new.cluster.ids) <- levels(all)
all <- RenameIdents(all, new.cluster.ids)
all <- AddMetaData(object = all,metadata = all@active.ident ,col.name = 'Cell')
#divide R mice (group1) and NR mice (group2) 
Idents(all) <- all$Group
R <- subset(all, idents = c("group1")) 
NR <- subset(all, idents = c("group2")) 
#make plot
Idents(R) <- R$Cell
Idents(NR) <- NR$Cell
pl1 <- DimPlot(R, reduction = "umap", label = F, pt.size = 0.05) + NoLegend()+ xlim (c(-15, 10))+ ylim (c(-15, 15))
pl2 <- DimPlot(NR, reduction = "umap", label = F, pt.size = 0.05) + NoLegend()+ xlim (c(-15, 10))+ ylim (c(-15, 15))
pl1+pl2
#make density contour plot
scale.col <- cet_pal(16, name = "blues")
plot1 <- pl1[[1]] & 
  stat_density_2d(aes_string(x = "UMAP_1", y = "UMAP_2", fill = "after_stat(level)"), 
                  linewidth = 0.15, geom = "density_2d_filled", 
                  colour = "black", alpha = 0.02, n = 150, h = c(1.2, 1.2)) & 
  scale_fill_gradientn(colours = scale.col)
plot2 <- pl2[[1]] & 
  stat_density_2d(aes_string(x = "UMAP_1", y = "UMAP_2", fill = "after_stat(level)"), 
                  linewidth = 0.15, geom = "density_2d_filled", 
                  colour = "black", alpha = 0.02, n = 150, h = c(1.2, 1.2)) & 
  scale_fill_gradientn(colours = scale.col)
plot1+plot2
#make dot plot and feature plot
DotPlot(all, features = c("Ly6c1","Cd14","Itgam","Cx3cr1","Cd3e","Cd4","Pdcd1","Cd8a","Klrb1c",
                          "Cd86","Trac","Trbc1","Trdc","Trgc1","Ccr7","Cxcr2","Cd274","Ly6g","Siglech","Cd19","Plp1"),
        cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) ),
        dot.min = 0,scale.min = 0,col.min=0) + RotatedAxis()
DotPlot(all, features = c("Cxcr2","S100a8","S100a9","Cd274","Clec2d","Acod1","Il1rn"),
        cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) ),
        dot.min = 0,scale.min = 0,col.min=0) + RotatedAxis()
FeaturePlot(all, features = c("Cxcr2"), order = T, cols = c(rgb(230,230,230,max=255), rgb(255,109,9,max=255) )) 
#find signiture gene of MDSC-like cells
Idents(all) <- all$seurat_clusters
MDSClike <- FindMarkers(all, ident.1 = 10,test.use = "MAST")

#####
#show cell number
count <- FetchData(all, 
                   vars = c("Cell", "Group")) %>%
  group_by(Cell) %>%
  dplyr::count(Group) 
write.csv(count, file = "~/scRNAseq_cell number.csv") #calculate in excel file
cell <- read.table(file = "~/scRNAseq_cell number.csv",sep=",",header = T)
ggplot(cell, aes(x = cell, y = number, fill = sample)) +
  geom_col(colour = "black", position = "dodge", orientation = "x")+
  xlab(expression("Cells")) + 
  ylab(expression("% of proportion")) +
  scale_x_discrete(limit = c("Monocytes","T","Microglia", "NK", "MdC", "cDC", "MDSC-like", "pDC", "B", "Macrophage", "Oligodendrocyte"))+
  scale_fill_manual(values = c(rgb(50,1,245,max=255), rgb(245,98,29,max=255)))+
  geom_text(aes(label = paste(format(number, nsmall = 1))), colour = "white", size = 2.5, 
            position = position_dodge(.9), vjust = 1.5, fontface = "bold")+ 
  theme(panel.background = element_rect(fill = "transparent", colour = "black", size = 1.2),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title =element_text(size=15, color="black"),
        axis.text.x =element_text(size=12, color="black",angle = 60, hjust = 1),
        axis.text.y =element_text(size=12, color="black"),
        axis.ticks = element_line(colour = "black", size = 0.7))

#####
#volcano plot of altered rates of body weight
wight <- read.table(file = "~/weight.csv",sep=",",header = T)
wight <- wight %>% 
  mutate(
    Expression = case_when(log2FC >= 0 & p >= 1.3 ~ "Up-regulated",
                           log2FC <= -0 & p >= 1.3 ~ "Down-regulated",
                           TRUE ~ "No-significant")
  )
rownames(wight) <- wight[,1]
names <- wight[,1]
ggplot(wight, aes(log2FC, p, label = names)) +
  geom_point(aes(color = Expression), size = 1.5) +
  xlab(expression("Died mice / Survived mice (log"[2]*"FC)")) + 
  ylab(expression("-log"[10]*italic("P")*"-value")) +
  scale_color_manual(values = c("dodgerblue3", "gray80", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1))) +
  theme_minimal()+
  geom_text_repel()+ 
  theme(panel.background = element_rect(fill = "transparent", colour = "black", size = 1.2),
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title =element_text(size=15,  family="sans", color="black"),
        axis.text.x =element_text(size=12,  family="sans", color="black"),
        axis.text.y =element_text(size=12,  family="sans", color="black"),
        axis.ticks = element_line(colour = "black", size = 0.7))














