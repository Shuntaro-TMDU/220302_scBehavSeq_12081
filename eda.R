setwd("/Users/shuntaro/Documents/R/CW/230302_scBehavSeq_12081/R_output")
getwd()

install.packages("rstatix")
install.packages("ggpubr")

library(tidyverse)
library(cowplot)
library(scales)
library(rstatix)
library(ggpubr)

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

df.joined.2.new = read_csv("12081_ATG_Ly6G_Parameters_M2.csv")

df.annotation = read_csv("12081_ATG_Ly6G_Annotation.csv") %>%
        mutate(Order = as.integer(Order), 
               Cluster = factor(Cluster, levels = c("1", "2"))) %>%
        print()

# クラスター間で有意（>0.25: log scale）に変化しているパラメータを検索する。
## df.joined.2.newとdf.annotationをOrderをアンカーにして結合
df.joined.2.hc = df.joined.2.new %>%
        full_join(df.annotation, by = "Order") %>%
        relocate(Cluster, .after = Label)

## クラスターごとのパラメータの平均値を計算
df.joined.2.hc.mean = df.joined.2.hc %>%
        group_by(Cluster) %>%
        summarise(across(Area_Mean:IntegratedDistance, mean)) %>%
        print()

df.joined.2.hc.mean.t = t(as.matrix(df.joined.2.hc.mean))
colnames(df.joined.2.hc.mean.t) = c("Cluster_1", "Cluster_2")


df.joined.2.hc.mean.t = df.joined.2.hc.mean %>%
        mutate(Cluster = str_c("Cluster", Cluster, sep = "_")) %>%
        pivot_longer(!Cluster, names_to = "Variable", values_to = "value") %>%
        pivot_wider(names_from = Cluster)

df.joined.2.hc.summary = df.joined.2.hc.mean.t %>%
        mutate(loge_fold = log(Cluster_2 / Cluster_1), 
               abs_loge_fold = abs(loge_fold)) %>%
        arrange(-abs_loge_fold) %>%
        print()

write_csv(df.joined.2.hc.summary, "12081_ATG_Ly6G_Hc2_Summary.csv")

Top10.param = c("Displacement", 
               "DistanceTraveled_Max", 
               "Speed_Max", 
               "DistanceTraveled_Std", 
               "Speed_Std", 
               "IntegratedDistance", 
               "DistanceTraveled_Mean", 
               "Speed_Mean", 
               "MajorAxisLength_Std", 
               "Area_Std")

# Top5.paramをクラスター間で比較する
## 指定したパラメータのtukey's t-test結果を表示する関数anova.pwcを作成する。
tukey = function(data = df.joined.2.hc, 
                     parameter){
        
        tukey.param = tukey_hsd(data, as.formula(paste0(parameter, "~ Cluster")))
        
        return(tukey.param)
}

tukey.Displacement = tukey(df.joined.2.hc, 
                                   parameter = "Displacement")
tukey.Displacement

## 指定したパラメータのavova結果を反映したviolin plotを表示する関数PlotViolinを作成する。
PlotViolin = function(data_full = df.joined.2.hc, 
                      data_mean = df.joined.2.hc.mean, 
                      parameter = parameter, 
                      pvalue.y.position = pvalue.y.position){
        
        pwc.param = tukey_hsd(data_full, as.formula(paste0(parameter, "~ Cluster")))
        
        p1 = ggplot(data_full) +
                mytheme1 +
                theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
        
        p2_1 = p1 +
                aes(x = Cluster, y = .data[[parameter]]) +
                geom_violin(aes(fill = Cluster), 
                            trim = FALSE, 
                            adjust = 1) +
                theme(axis.title.x = element_blank(), 
                      axis.ticks.x = element_blank(), 
                      axis.text.x = element_blank(), 
                      legend.position = "right", 
                      legend.justification = "top")
        
        p2_2 = p2_1 + 
                geom_point(data = data_mean,
                           aes(x = Cluster, y = .data[[parameter]]),
                           size = 2) +
                theme(axis.title.x = element_blank(), 
                      axis.ticks.x = element_blank(), 
                      axis.text.x = element_blank())
        
        ViolinPlot = p2_2 +
                stat_pvalue_manual(pwc.param[pwc.param$term == "Cluster", ], 
                                   label = "P = {p.adj}", 
                                   y.position = pvalue.y.position, 
                                   step.increase = 0.1, 
                                   tip.length = 0, 
                                   hide.ns = TRUE)
        
        return(ViolinPlot)
}

### Displacement
ViolinPlot.Displacement = PlotViolin(data_full = df.joined.2.hc, 
                                     data_mean = df.joined.2.hc.mean,
                                     parameter = "Displacement", 
                                     pvalue.y.position = 180)

ViolinPlot.Displacement
ggsave("Violin_Displacement.tiff", width = 5, height = 5)

### DistanceTraveled_Max
ViolinPlot.DistanceTraveled_Max = PlotViolin(parameter = "DistanceTraveled_Max", 
                                     pvalue.y.position = 60)

ViolinPlot.DistanceTraveled_Max
ggsave("Violin_DistanceTraveled_Max.tiff", width = 5, height = 5)

### Speed_Max
ViolinPlot.Speed_Max = PlotViolin(parameter = "Speed_Max", 
                                  pvalue.y.position = 20)

ViolinPlot.Speed_Max
ggsave("Violin_Speed_Max.tiff", width = 5, height = 5)

### DistanceTraveled_Std
ViolinPlot.DistanceTraveled_Std = PlotViolin(parameter = "DistanceTraveled_Std", 
                                             pvalue.y.position = 20)

ViolinPlot.DistanceTraveled_Std
ggsave("Violin_DistanceTraveled_Std.tiff", width = 5, height = 5)

### Speed_Std
ViolinPlot.Speed_Std = PlotViolin(parameter = "Speed_Std", 
                                             pvalue.y.position = 6)

ViolinPlot.Speed_Std
ggsave("Violin_Speed_Std.tiff", width = 5, height = 5)

### IntegratedDistance
ViolinPlot.IntegratedDistance = PlotViolin(parameter = "IntegratedDistance", 
                                             pvalue.y.position = 200)

ViolinPlot.IntegratedDistance
ggsave("Violin_IntegratedDistance.tiff", width = 5, height = 5)

### MajorAxisLength_Std
ViolinPlot.MajorAxisLength_Std = PlotViolin(parameter = "MajorAxisLength_Std", 
                                           pvalue.y.position = 10)

ViolinPlot.MajorAxisLength_Std
ggsave("Violin_MajorAxisLength_Std.tiff", width = 5, height = 5)

### EquivalentDiameter_Mean
summary(df.joined.2.hc$EquivalentDiameter_Mean)
summary(df.joined.2.hc$EquivalentDiameter_Std)
ViolinPlot.EquivalentDiameter_Mean = PlotViolin(parameter = "EquivalentDiameter_Mean", 
                                 pvalue.y.position = 120)

ViolinPlot.EquivalentDiameter_Mean
ggsave("Violin_EquivalentDiameter_Mean.tiff", width = 5, height = 5)

### Area_Std
ViolinPlot.Area_Std = PlotViolin(parameter = "Area_Std", 
                                      pvalue.y.position = 120)

ViolinPlot.Area_Std
ggsave("Violin_Area_Std.tiff", width = 5, height = 5)

### Perimeter_Std
ViolinPlot.Perimeter_Std = PlotViolin(parameter = "Perimeter_Std", 
                                            pvalue.y.position = 20)

ViolinPlot.Perimeter_Std
ggsave("Violin_Perimeter_Std.tiff", width = 5, height = 5)

### Eccentricity_Std
ViolinPlot.Eccentricity_Std = PlotViolin(parameter = "Eccentricity_Std", 
                                      pvalue.y.position = 0.4)

ViolinPlot.Eccentricity_Std
ggsave("Violin_Eccentricity_Std.tiff", width = 5, height = 5)

### Solidity_Std
ViolinPlot.Solidity_Std = PlotViolin(parameter = "Solidity_Std", 
                                         pvalue.y.position = 0.1)

ViolinPlot.Solidity_Std
ggsave("Violin_Solidity_Std.tiff", width = 5, height = 5)

# 全パラメータのviolin plotをまとめて表示
plot_grid(ViolinPlot.Speed_Max + theme(legend.position = "none"), 
          ViolinPlot.Speed_Std + theme(legend.position = "none"), 
          ViolinPlot.Area_Std + theme(legend.position = "none"), 
          ViolinPlot.Eccentricity_Std + theme(legend.position = "none"), 
          ViolinPlot.Solidity_Std + theme(legend.position = "right", 
                        legend.justification = "top"), 
          nrow = 2, 
          ncol = 3,
          rel_widths = c(1, 1.2, 1))
ggsave("Violin_Top.tiff", width = 15, height = 10)          
