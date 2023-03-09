setwd("/Users/shuntaro/Documents/R/CW/230302_scBehavSeq_12081/R_output")
getwd()

install.packages("factoextra")
install.packages("tidyverse")
install.packages("viridis")
install.packages("pheatmap")

library(factoextra)
library(tidyverse)
library(gplots)
library(viridis)
library(pheatmap)

df.normalized.1 = read_csv("12081_ATG_Ly6G_Normalized_M1.csv")
df.normalized.2 = read_csv("12081_ATG_Ly6G_Normalized_M2.csv")

# 相関マトリックス
CorMat1 = cor(df.normalized.1)
CorMat2 = cor(df.normalized.2)

## dendrogram
DenCor1 = as.dendrogram(hclust(as.dist(1 - CorMat1)))
DenCor2 = as.dendrogram(hclust(as.dist(1 - CorMat2)))

## heatmap
tiff("CorMat_1.tiff", width = 1200, height = 1200)
gplots::heatmap.2(
        x = as.matrix(CorMat1), 
        Rowv = DenCor1, 
        dendrogram = "none", 
        reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean), 
        symm = TRUE,
        revC = T,
        col = rev(heat.colors(256)), 
        scale = "none", 
        key = TRUE,
        keysize = 1,
        symkey = TRUE, 
        density.info = "none", 
        trace = "none",
        lwid = c(1, 5), 
        lhei = c(1, 5),
        key.title = NA, 
        key.xlab = "Correlation score",
        margin = c(10, 10),
        # cellnote = as.matrix(round(cor.mat, 1)),
        # notecol = "black", # cell label colour
        # notecex = 1, # cell label cex
        cexCol = 0.7, # Column label cex
        cexRow = 0.7, 
        offsetRow = 0, 
        offsetCol = 0
)

dev.off()

tiff("CorMat_2.tiff", width = 1200, height = 1200)
gplots::heatmap.2(
        x = as.matrix(CorMat2), 
        Rowv = DenCor2, 
        dendrogram = "none", 
        reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean), 
        symm = TRUE,
        revC = T,
        col = rev(heat.colors(256)), 
        scale = "none", 
        key = TRUE,
        keysize = 1,
        symkey = TRUE, 
        density.info = "none", 
        trace = "none",
        lwid = c(1, 5), 
        lhei = c(1, 5),
        key.title = NA, 
        key.xlab = "Correlation score",
        margin = c(10, 10),
        # cellnote = as.matrix(round(cor.mat, 1)),
        # notecol = "black", # cell label colour
        # notecex = 1, # cell label cex
        cexCol = 0.7, # Column label cex
        cexRow = 0.7, 
        offsetRow = 0, 
        offsetCol = 0
)

dev.off()

# 階層的クラスタリング
Cden1 = dist(df.normalized.1, method = "euclidean")
Cden2 = dist(df.normalized.2, method = "euclidean")

length(Cden1)
head(as.matrix(cden))

HcRes1 = hclust(Cden1, method = "ward.D2")
HcRes2 = hclust(Cden2, method = "ward.D2")

tiff(filename = "hc_1.tiff", width = 1200, height = 800, pointsize = 20)
plot(HcRes1, hang = -1, cex = 0.4, xlab = "", sub = "")
dev.off()

tiff(filename = "hc_2.tiff", width = 1200, height = 800, pointsize = 20)
plot(HcRes2, hang = -1, cex = 0.4, xlab = "", sub = "")
dev.off()

# まず、M2によるパラメータ抽出とクラスター数2を選択
HcRes2.2 = cutree(HcRes2, k = 2) 
head(HcRes2.2)
length(HcRes2.2)
nrow(HcRes2.2) # 102
table(HcRes2.2) # 74, 28

HcRes2.RowName = HcRes2$order

HcRes2.2.new = as_tibble(data.frame(HcRes2.2)) %>%
        mutate(Cluster = rep(c("2", "1"), c(74, 28)), 
               Cluster = factor(Cluster, levels = c("1", "2"))) %>%
        select(Cluster) %>%
        print()

# %>%
#         select(Cluster)

# heatmap, clustering
matrix.normalized.2 = as.matrix(df.normalized.2)
rownames(matrix.normalized.2) = paste(1:102)

matrix.normalized.2.t = t(matrix.normalized.2)

heatmap.1 = pheatmap(matrix.normalized.2.t, 
                     clustering_method = "ward.D2")

tiff(filename = "heatmap.1.tiff", width = 1200, height = 1000)
heatmap.1
dev.off()

# クラスタリングした結果変更された行名を取得する。
HeatOrderId = heatmap.1$tree_col$order
HeatOrderId

length(HeatOrderId)

# どの行がどのクラスターに入ったかを調べる。
ClusterAttr = cutree(heatmap.1$tree_col, k = 2)
table(ClusterAttr) # 74, 28

annotation_col = data.frame(
        Cluster = factor(rep(c("1", "2"), c(28, 74)))
)

rownames(annotation_col) = paste(HeatOrderId)
row.names(annotation_col)

head(annotation_col)

ann_colors = list(
        Cluster = c("1" = "salmon", "2" = "Cyan")
)

heatmap.2 = pheatmap(matrix.normalized.2.t, 
                     clustering_method = "ward.D2", 
                     legend_labels = "Scaled score",
                     annotation_col = annotation_col, 
                     annotation_colors = ann_colors,
                     show_colnames = FALSE)

tiff(filename = "heatmap.2.tiff", width = 800, height = 600)
heatmap.2
dev.off()

df.annotation = rownames_to_column(annotation_col, "Order") %>% 
        as_tibble(annotation_col) %>%
        mutate(Order = as.integer(Order)) %>%
        print()
head(df.annotation)

write_csv(df.annotation, "12081_ATG_Ly6G_Annotation.csv")

# それぞれのクラスターにおいて、他のクラスターに対して有意（<0.25: log scale）に変化しているパラメータを検索する。
dsp = df.hc4 %>%
        rowid_to_column("Cluster") %>%
        mutate(Cluster = str_c("Cluster", Cluster, sep = "_")) %>%
        pivot_longer(!Cluster, names_to = "variable", values_to = "value") %>%
        pivot_wider(names_from = Cluster)

View(dsp)

2^0.25
exp(0.25)
10^0.25

## Cluster1
dsp.c1 = dsp %>%
        rowwise() %>%
        mutate(Difference = Cluster_1 - mean(c(Cluster_2, Cluster_3, Cluster_4)), 
               Log_Difference = log10(abs(Difference))) %>%
        arrange(-Log_Difference)

View(dsp.c1)

list.dsp.c1 = dsp.c1 %>%
        filter(Log_Difference >= 0.25) %>%
        .$variable 
print(list.dsp.c1) #[1] "Perimeter_Mean"          "BoundingBoxArea_Std"     "EquivalentDiameter_Mean"

best.list.c1 = c("Perimeter_Mean", 
                 "BoundingBoxArea_Std", 
                 "EquivalentDiameter_Mean")

## Cluster2
dsp.c2 = dsp %>%
        rowwise() %>%
        mutate(Difference = Cluster_2 - mean(c(Cluster_1, Cluster_3, Cluster_4)), 
               Log_Difference = log10(abs(Difference))) %>%
        arrange(-Log_Difference)

View(dsp.c2) 

list.dsp.c2 = dsp.c2 %>%
        filter(Log_Difference >= 0.25) %>%
        .$variable %>%
        print() #character(0)

best.list.c2 = c("Solidity_Std", 
                 "MedianRadius_Mean", 
                 "Compactness_Std", 
                 "Linearity")

## Cluster3
dsp.c3 = dsp %>%
        rowwise() %>%
        mutate(Difference = Cluster_3 - mean(c(Cluster_1, Cluster_2, Cluster_4)), 
               Log_Difference = log10(abs(Difference))) %>%
        arrange(-Log_Difference)

View(dsp.c3) 

list.dsp.c3 = dsp.c3 %>%
        filter(Log_Difference >= 0.25) %>%
        .$variable %>%
        print() #character(0)

best.list.c3 = c("FormFactor_Std", 
                 "Eccentricity_Mean")

## Cluster4
dsp.c4 = dsp %>%
        rowwise() %>%
        mutate(Difference = Cluster_4 - mean(c(Cluster_1, Cluster_2, Cluster_3)), 
               Log_Difference = log10(abs(Difference))) %>%
        arrange(-Log_Difference)

View(dsp.c4)

list.dsp.c4 = dsp.c4 %>%
        filter(Log_Difference >= 0.25) %>%
        .$variable 
print(list.dsp.c4)
# [1] "CentralMoment_0_2_Std"           "CentralMoment_2_2_Std"          
# [3] "CentralMoment_2_1_Std"           "BoundingBoxArea_Std"            
# [5] "CentralMoment_1_2_Std"           "CentralMoment_1_1_Std"          
# [7] "Area_Std"                        "CentralMoment_0_0_Std"          
# [9] "CentralMoment_2_2_Mean"          "InertiaTensorEigenvalues_0_Std" 
# [11] "InertiaTensorEigenvalues_1_Std"  "InertiaTensor_0_0_Std"          
# [13] "InertiaTensor_1_1_Std"           "Perimeter_Std"                  
# [15] "CentralMoment_0_1_Std"           "MinFeretDiameter_Std"           
# [17] "InertiaTensor_0_1_Std"           "InertiaTensor_1_0_Std"          
# [19] "MinorAxisLength_Std"             "EquivalentDiameter_Std"         
# [21] "MaxFeretDiameter_Std"            "MajorAxisLength_Std"            
# [23] "MaximumRadius_Std"               "MeanRadius_Std"                 
# [25] "CentralMoment_0_2_Mean"          "InertiaTensor_1_1_Mean"         
# [27] "InertiaTensorEigenvalues_0_Mean" "BoundingBoxArea_Mean"           
# [29] "IntegratedDistance"              "Speed_Max"                      
# [31] "Area_Mean"                       "CentralMoment_0_0_Mean"  

best.list.c4 = c("Area_Mean", 
                 "Area_Std", 
                 "IntegratedDistance", 
                 "Speed_Max", 
                 "Speed_Std")

best.list = c(best.list.c1, best.list.c2, best.list.c3, best.list.c4)
print(best.list)
# [1] "Perimeter_Mean"          "BoundingBoxArea_Std"     "EquivalentDiameter_Mean"
# [4] "Solidity_Std"            "MedianRadius_Mean"       "Compactness_Std"        
# [7] "Linearity"               "FormFactor_Std"          "Eccentricity_Mean"      
# [10] "Area_Mean"               "Area_Std"                "IntegratedDistance"     
# [13] "Speed_Max"               "Speed_Std"  

df.hc4.best = df.hc4 %>%
        select(best.list)

# それぞれのクラスターにおけるDSPをまとめたパラメータで再度ヒートマップ
list.dsp = c(list.dsp.c1, list.dsp.c4)
df.hc4.dsp = df.hc4 %>%
        select(list.dsp)

write_csv(df.hc4.dsp, "Neutrophil_hc4_dsp.csv")

tiff(filename = "heatmap.hc4.dsp.tiff", width = 800, height = 600)
heatmap.2(
        as.matrix(df.hc4.dsp), 
        dendrogram = "col", 
        col = viridis,  
        key = TRUE,
        keysize = 1,
        key.title = NA,
        key.xlab = "scaled score",
        symkey = FALSE,
        density.info = "none",
        trace = "none",
        lwid = c(1, 5), 
        lhei = c(1, 5),
        margin = c(26,6),
        cexRow = 2,
        cexCol = 2
)
dev.off()

q()
