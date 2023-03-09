setwd("/Users/shuntaro/Documents/R/CW/230302_scBehavSeq_12081/R_output")
getwd()

install.packages("vctrs")
install.packages("Rtsne")

library(tidyverse)
library(cowplot)
library(scales)
library(Rtsne)
library(umap)

set.seed(2022)

mytheme1 <- theme(panel.background = element_blank(),
                  plot.margin = unit(c(1, 1, 0, 1), "lines"),
                  plot.title = element_text(hjust = 0.5, size=20, family = "Arial", face = "bold"),
                  axis.line = element_line(size = 0.7),
                  axis.ticks = element_line(size = 0.7), 
                  axis.title = element_text(size=20, family = "Arial", face = "bold"), 
                  axis.text = element_text(size=20, color = "black", family = "Arial", face = "bold"), 
                  legend.position = "top", 
                  legend.key = element_rect(fill = "white"), 
                  legend.title = element_text(size = 20, family = "Arial", face = "bold"),
                  legend.text = element_text(size = 20, family = "Arial", face = "bold"), 
                  strip.background = element_rect(colour = "black", fill = "white"), 
                  strip.text.x = element_text(size=18, family = "Arial", face = "bold"), 
                  strip.text.y = element_text(size=18, family = "Arial", face = "bold"))

mytheme2 = theme_bw() + 
        theme(panel.background = element_blank(),
              plot.margin = unit(c(1, 1, 0, 1), "lines"),
              plot.title = element_text(hjust = 0.5, size=20, family = "Arial", face = "bold"),
              axis.title = element_text(size=20, family = "Arial", face = "bold"), 
              axis.text = element_text(size=20, color = "black", family = "Arial", face = "bold"), 
              panel.grid = element_blank(),
              legend.position = "right", 
              legend.justification = "top",
              legend.key = element_rect(fill = "white"), 
              legend.title = element_text(size = 20, family = "Arial", face = "bold"),
              legend.text = element_text(size = 20, family = "Arial", face = "bold"), 
              strip.background = element_rect(colour = "black", fill = "white"), 
              strip.text.x = element_text(size=18, family = "Arial", face = "bold"), 
              strip.text.y = element_text(size=18, family = "Arial", face = "bold"))

df.meta.new = read_csv("12081_ATG_Ly6G_Metadata.csv")
df.joined.1.new = read_csv("12081_ATG_Ly6G_Parameters_M1.csv")
df.joined.2.new = read_csv("12081_ATG_Ly6G_Parameters_M2.csv")

# metadataから、Lifetime = 1のオブジェクトのみを抽出し、ナンバリング
df.meta.1 = df.meta.new %>%
        filter(Lifetime == 1) %>%
        mutate(Order = 1:nrow(.)) %>%
        relocate(Order, .before = ID)

nrow(df.meta.1)
View(df.meta.1)

# 標準化
Normalization = function(data = data){
        df.normalized = data %>%
                select(!c(Order, Run, Label)) %>%
                scale() %>%
                data.frame()
        
        return(df.normalized)
}

df.normalized.1 = Normalization(data = df.joined.1.new)
df.normalized.2 = Normalization(data = df.joined.2.new)

head(df.normalized)
View(df.normalized)

write_csv(df.normalized.1, "12081_ATG_Ly6G_Normalized_M1.csv")
write_csv(df.normalized.2, "12081_ATG_Ly6G_Normalized_M2.csv")

# PCA実行
RunPCA = function(df.normalized = df.normalized){
        PCA = prcomp(df.normalized)
        
        return(PCA)
}

PCA1 = RunPCA(df.normalized = df.normalized.1)
PCA2 = RunPCA(df.normalized = df.normalized.2)

PCA$x

PCA$sdev # 各種成分の固有値
PCA$sdev / sum(pr$sdev) # 各種成分の寄与率の計算
cumsum(PCA$sdev / sum(PCA$sdev)) # 累積寄与率の計算
head(PCA$x) # PCへのアクセス
summary(PCA.1)
summary(PCA.2)

eigen_value = PCA$sdev^2
plot(eigen_value)

# ggplotを用いたPCA結果の可視化
PlotPCA.pathology = function(PCA = PCA, 
                             metadata = df.meta.1){
        
        BindedPCA = as_tibble(as.data.frame(PCA$x)) %>%
                bind_cols(df.meta.1, .)
        
        Plot = ggplot(BindedPCA) +
                aes(x = PC1, y = PC2) +
                geom_point(aes(color = Pathology)) +
                mytheme1 + 
                theme(axis.ticks = element_blank(), 
                      axis.text = element_blank())
        
        return(Plot)
        
}

PlotPCA1 = PlotPCA.pathology(PCA = PCA1, 
                             metadata = df.meta.1)

PlotPCA2 = PlotPCA.pathology(PCA = PCA2, 
                             metadata = df.meta.1)

PlotPCA1
ggsave("pca_1.tiff", width = 7, height = 7)

PlotPCA2
ggsave("pca_2.tiff", width = 7, height = 7)

## Run毎に色分けして表示する
# ggplot(new.pr) +
#         aes(x = PC1, y = PC2) +
#         geom_point(aes(color = Run)) +
#         # facet_grid(cols = vars(Time)) + 
#         mytheme1 + 
#         theme(axis.ticks = element_blank(), 
#               axis.text = element_blank())
# 
# ggsave("pca.run.tiff", width = 8, height = 6)

# tSNE
RunTSNE = function(df.normalized = df.normalized, 
                   perplexity = perplexity){
        TSNE = Rtsne(as.matrix(df.normalized), 
                     perplexity = perplexity)
        
        return(TSNE)
}

TSNE1.40 = RunTSNE(df.normalized = df.normalized.1, 
                   perplexity = 40) 

TSNE1.20 = RunTSNE(df.normalized = df.normalized.1, 
                   perplexity = 20) 

TSNE1.10 = RunTSNE(df.normalized = df.normalized.1, 
                   perplexity = 10) 

TSNE2.40 = RunTSNE(df.normalized = df.normalized.2, 
                   perplexity = 40) 

TSNE2.20 = RunTSNE(df.normalized = df.normalized.2, 
                   perplexity = 20) 

TSNE2.10 = RunTSNE(df.normalized = df.normalized.2, 
                   perplexity = 10) 

TSNE1$Y
View(TSNE1)
summary(TSNE1)

## 可視化
PlotTSNE.pathology = function(TSNE = TSNE, 
                              metadata = metadata){
        
        BindedTSNE = as_tibble(as.data.frame(TSNE$Y)) %>%
                rename(tSNE1 = V1, tSNE2 = V2) %>%
                bind_cols(metadata, .)
        
        Plot = ggplot(BindedTSNE) +
                aes(x = tSNE1, y = tSNE2) +
                geom_point(aes(color = Pathology), size = 2) +
                mytheme1 + 
                theme(axis.ticks = element_blank(), 
                      axis.text = element_blank())
        
        return(Plot)
        
}

PlotTSNE1.40 = PlotTSNE.pathology(TSNE = TSNE1.40, 
                                  metadata = df.meta.1)
PlotTSNE1.40 

PlotTSNE1.20 = PlotTSNE.pathology(TSNE = TSNE1.20, 
                                  metadata = df.meta.1)
PlotTSNE1.20
ggsave("tSNE_1_pe20.tiff", width = 7, height = 7)

PlotTSNE1.10 = PlotTSNE.pathology(TSNE = TSNE1.10, 
                                  metadata = df.meta.1)
PlotTSNE1.10

PlotTSNE2.40 = PlotTSNE.pathology(TSNE = TSNE2.40, 
                                  metadata = df.meta.1)
PlotTSNE2.40

PlotTSNE2.20 = PlotTSNE.pathology(TSNE = TSNE2.20, 
                                  metadata = df.meta.1)
PlotTSNE2.20
ggsave("tSNE_2_pe20.tiff", width = 7, height = 7)

PlotTSNE2.10 = PlotTSNE.pathology(TSNE = TSNE2.10, 
                                  metadata = df.meta.1)
PlotTSNE2.10
ggsave("tSNE_2_pe10.tiff", width = 7, height = 7)


## Time毎に分割表示する
# ggplot(new.tsne) +
#         aes(x = tSNE1, y = tSNE2) +
#         geom_point(aes(color = Time)) +
#         facet_grid(cols = vars(Time)) + 
#         mytheme1 + 
#         theme(axis.ticks = element_blank(), 
#               axis.text = element_blank(), 
#               legend.position = "none")
# 
# ggsave("tSNE_facet_after.resize.tiff", width = 10, height = 6)

#UMAP
RunUMAP = function(df.normalized = df.normalized){
        UMAP = umap(df.normalized)
        
        return(UMAP) 
}

UMAP1 = umap(df.normalized.1)
UMAP2 = umap(df.normalized.2)


UMAP1$layout
str(umap)

# 可視化
PlotUMAP.pathology = function(UMAP = UMAP, 
                              metadata = metadata){
        
        BindedUMAP = as_tibble(as.data.frame(UMAP$layout)) %>%
                rename(UMAP1 = V1, UMAP2 = V2) %>%
                bind_cols(metadata, .)
        
        set.seed(2022)
        
        Plot = ggplot(BindedUMAP) +
                aes(x = UMAP1, y = UMAP2) +
                geom_point(aes(color = Pathology), size = 2) +
                mytheme1 + 
                theme(axis.ticks = element_blank(), 
                      axis.text = element_blank())
        
        return(Plot)
        
}

PlotUMAP1 = PlotUMAP.pathology(UMAP = UMAP1, 
                               metadata = df.meta.1)
PlotUMAP1
ggsave("UMAP_1.tiff", width = 7, height = 7)

PlotUMAP2 = PlotUMAP.pathology(UMAP = UMAP2, 
                     metadata = df.meta.1)
PlotUMAP2
ggsave("UMAP_2.tiff", width = 7, height = 7)

## Time毎に分割表示する
# ggplot(new.umap) +
#         aes(x = UMAP1, y = UMAP2) +
#         geom_point(aes(color = Time), 
#                    alpha = 0.5) +
#         facet_grid(cols = vars(Time)) + 
#         mytheme1 + 
#         theme(axis.ticks = element_blank(), 
#               axis.text = element_blank(), 
#               legend.position = "none")
# 
# ggsave("umap_facet_after.resize.tiff", width = 10, height = 6)
