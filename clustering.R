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
        cexCol = 1, # Column label cex
        cexRow = 1, 
        offsetRow = 0, 
        offsetCol = 0
)

dev.off()

tiff("CorMat_2.tiff", width = 1400, height = 1400)
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
        margin = c(14, 14),
        # cellnote = as.matrix(round(cor.mat, 1)),
        # notecol = "black", # cell label colour
        # notecex = 1, # cell label cex
        cexCol = 1.2, # Column label cex
        cexRow = 1.2, 
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


