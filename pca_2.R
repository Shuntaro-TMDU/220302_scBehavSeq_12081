setwd("/Users/shuntaro/Documents/R/CW/230302_scBehavSeq_12081/R_output")
getwd()

install.packages("vctrs")
install.packages("Rtsne")

library(tidyverse)
library(cowplot)
library(scales)
library(Rtsne)
library(tsne)
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
df.normalized.1 = read_csv("12081_ATG_Ly6G_Normalized_M1.csv")
df.normalized.2 = read_csv("12081_ATG_Ly6G_Normalized_M2.csv")

df.annotation = read_csv("12081_ATG_Ly6G_Annotation.csv") %>%
        mutate(Order = as.integer(Order), 
               Cluster = factor(Cluster, levels = c("1", "2"))) %>%
        print()

summary(df.annotation)

# metadataから、Lifetime = 1のオブジェクトのみを抽出
df.meta.1 = df.meta.new %>%
        filter(Lifetime == 1) %>%
        mutate(Order = 1:102) %>%
        relocate(Order, .after = Label) %>%
        print()

head(df.meta.1)

df.meta.hc2.2 = df.meta.1 %>%
        full_join(., df.annotation, by = "Order") %>%
        print()

View(df.meta.hc2.2)

# PCA実行
RunPCA = function(df.normalized = df.normalized){
        PCA = prcomp(df.normalized)
        
        return(PCA)
}

PCA1 = RunPCA(df.normalized = df.normalized.1)
PCA2 = RunPCA(df.normalized = df.normalized.2)

# ggplotを用いたPCA結果の可視化
PlotPCA.cluster = function(PCA = PCA, 
                           metadata = metadata){
        
        BindedPCA = as_tibble(as.data.frame(PCA$x)) %>%
                bind_cols(metadata, .)
        
        Plot = ggplot(BindedPCA) +
                aes(x = PC1, y = PC2) +
                geom_point(aes(color = Cluster)) +
                mytheme1 + 
                theme(axis.ticks = element_blank(), 
                      axis.text = element_blank())
        
        return(Plot)
        
}

PlotPCA.M2.Hc2 = PlotPCA.cluster(PCA = PCA2, 
                                 metadata = df.meta.hc2.2)

PlotPCA.M2.Hc2
ggsave("pca_M2_Hc2.tiff", width = 7, height = 7)

# tSNE
RunTSNE = function(df.normalized = df.normalized, 
                   perplexity = perplexity){
        TSNE = Rtsne(as.matrix(df.normalized), 
                     perplexity = perplexity)
        
        return(TSNE)
}

TSNE2.20 = RunTSNE(df.normalized = df.normalized.2, 
                   perplexity = 20) 

## 可視化
PlotTSNE.cluster = function(TSNE = TSNE, 
                            metadata = metadata){
        
        BindedTSNE = as_tibble(as.data.frame(TSNE$Y)) %>%
                rename(tSNE1 = V1, tSNE2 = V2) %>%
                bind_cols(metadata, .)
        
        Plot = ggplot(BindedTSNE) +
                aes(x = tSNE1, y = tSNE2) +
                geom_point(aes(color = Cluster), size = 2) +
                mytheme1 + 
                theme(axis.ticks = element_blank(), 
                      axis.text = element_blank())
        
        return(Plot)
        
}

PlotTSNE.M2.Hc2.Pe20 = PlotTSNE.cluster(TSNE = TSNE2.20, 
                                        metadata = df.meta.hc2.2)

PlotTSNE.M2.Hc2.Pe20
ggsave("tSNE_M2_Hc2_pe20.tiff", width = 7, height = 7)

#UMAP
RunUMAP = function(df.normalized = df.normalized){
        UMAP = umap(df.normalized)
        
        return(UMAP) 
}

UMAP1 = umap(df.normalized.1)
UMAP2 = umap(df.normalized.2)

# 可視化
PlotUMAP.cluster = function(UMAP = UMAP, 
                            metadata = metadata){
        
        BindedUMAP = as_tibble(as.data.frame(UMAP$layout)) %>%
                rename(UMAP1 = V1, UMAP2 = V2) %>%
                bind_cols(metadata, .)
        
        set.seed(2022)
        
        Plot = ggplot(BindedUMAP) +
                aes(x = UMAP1, y = UMAP2) +
                geom_point(aes(color = Cluster), size = 2) +
                mytheme1 + 
                theme(axis.ticks = element_blank(), 
                      axis.text = element_blank())
        
        return(Plot)
        
}

PlotUMAP.M2.Hc2 = PlotUMAP.cluster(UMAP = UMAP2, 
                                   metadata = df.meta.hc2.2)
PlotUMAP.M2.Hc2
ggsave("UMAP_M2_Hc2.tiff", width = 7, height = 7)
