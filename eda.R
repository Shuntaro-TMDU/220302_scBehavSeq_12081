setwd("/Users/shuntaro/Documents/R/CW/221229_sc-seq_rIgG,ATG_3241,9292/R_output")
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


df1 = read_csv("9292_rIgG_24hr_Neutrophil.csv")
df2 = read_csv("3241_ATG_3hr_Neutrophil.csv")
df3 = read_csv("3241_ATG_24hr_Neutrophil.csv")

df4 = read_csv("9292_rIgG_24hr_Metadata.csv")
df5 = read_csv("3241_ATG_3hr_Metadata.csv")
df6 = read_csv("3241_ATG_24hr_Metadata.csv")

hc.clusters.4 = read_csv("hc4.csv")

df.joined = rbind(df1, df2, df3)%>%
        mutate(Condition = factor(Condition, levels = (c("rIgG_24hr", "ATG_3hr", "ATG_24hr"))))

df.meta = rbind(df4, df5, df6) %>%
        mutate(Condition = factor(Condition, levels = (c("rIgG_24hr", "ATG_3hr", "ATG_24hr"))))

nrow(df.joined)
nrow(df.meta)

df.joined.hc4 = df.joined %>%
        cbind(Cluster = hc.clusters.4) %>%
        mutate(Cluster = factor(Cluster, levels = c(1, 2, 3, 4)))

str(df.joined.hc4)

write_csv(df.joined.hc4, "Neutrophil_hc4.csv")

# metadataから、Lifetime = 1のオブジェクトのみを抽出
df.meta.1 = df.meta.new %>%
        filter(Lifetime == 1)

df.meta.hc2.2 = df.meta.1 %>%
        bind_cols(., HcRes2.2)

# 各パラメータのクラスター毎の平均を求める
df.joined.hc4.mean = df.joined.hc4 %>%
        group_by(Cluster) %>%
        summarize_at(vars(Area_Mean:MeanderRatio), mean)

View(df.joined.hc4.mean)

# 可視化
print(best.list)
# [1] "Perimeter_Mean"          "BoundingBoxArea_Std"     "EquivalentDiameter_Mean"
# [4] "Solidity_Std"            "MedianRadius_Mean"       "Compactness_Std"        
# [7] "Linearity"               "FormFactor_Std"          "Eccentricity_Mean"      
# [10] "Area_Mean"               "Area_Std"                "IntegratedDistance"     
# [13] "Speed_Max"               "Speed_Std"  

## Area_Mean
p1 = ggplot(df.joined.hc4) +
        mytheme1 +
        theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
        
p2_1 = p1 +
        aes(x = Cluster, y = Area_Mean) +
        geom_violin(aes(fill = Cluster), 
                    trim = FALSE) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())

p2_1
p2_2 = p2_1 + 
        geom_point(data = df.joined.hc4.mean,
                   aes(x = Cluster, y = Area_Mean),
                   size = 2) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p2_2
aov.AreaMean =  anova_test(df.joined.hc4, Area_Mean ~ Cluster)
print(aov.AreaMean)

pwc.AreaMean = tukey_hsd(df.joined.hc4, Area_Mean ~ Cluster)
print(pwc.AreaMean)

p2_3 = p2_2 +
        stat_pvalue_manual(pwc.AreaMean[pwc.AreaMean$term == "Cluster", ], 
                           label = "P = {p.adj}", 
                           y.position = 200, 
                           step.increase = 0.1, 
                           tip.length = 0, 
                           hide.ns = TRUE)

p2_3

## Area_Std
p3_1 = p1 + aes(x = Cluster, y = Area_Std) +
        geom_violin(aes(fill = Cluster), 
                    trim = FALSE) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p3_1

p3_2 = p3_1 + 
        geom_point(data = df.joined.hc4.mean,
                       aes(x = Cluster, y = Area_Std),
                       size = 2) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p3_2

aov.AreaStd =  anova_test(df.joined.hc4, Area_Std ~ Cluster)
print(aov.AreaStd)

pwc.AreaStd = tukey_hsd(df.joined.hc4, Area_Std ~ Cluster)
print(pwc.AreaStd)

p3_3 = p3_2 +
        stat_pvalue_manual(pwc.AreaStd[pwc.AreaStd$term == "Cluster", ], 
                           label = "P = {p.adj}", 
                           y.position = 80, 
                           step.increase = 0.1, 
                           tip.length = 0, 
                           hide.ns = TRUE)

p3_3

# Solidity_Std
p4_1 = p1 + aes(x = Cluster, y = Solidity_Std) +
        geom_violin(aes(fill = Cluster), trim = FALSE) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p4_1

p4_2 = p4_1 + 
        geom_point(data = df.joined.hc4.mean,
                         aes(x = Cluster, y = Solidity_Std),
                         size = 2) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p4_2

aov.SolidityStd =  anova_test(df.joined.hc4, Solidity_Std ~ Cluster)
print(aov.SolidityStd)

pwc.SolidityStd = tukey_hsd(df.joined.hc4, Solidity_Std ~ Cluster)
print(pwc.SolidityStd)

p4_3 = p4_2 +
        stat_pvalue_manual(pwc.SolidityStd[pwc.SolidityStd$term == "Cluster", ], 
                           label = "P = {p.adj}", 
                           y.position = 0.1, 
                           step.increase = 0.1, 
                           tip.length = 0, 
                           hide.ns = TRUE)

p4_3

# Compactness_Std
p5_1 = p1 + aes(x = Cluster, y = Compactness_Std) +
        geom_violin(aes(fill = Cluster), trim = FALSE) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p5_1

p5_2 = p5_1 + 
        geom_point(data = df.joined.hc4.mean,
                   aes(x = Cluster, y = Compactness_Std),
                   size = 2) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p5_2

aov.CompactnessStd =  anova_test(df.joined.hc4, Compactness_Std ~ Cluster)
print(aov.CompactnessStd)

pwc.CompactnessStd = tukey_hsd(df.joined.hc4, Compactness_Std ~ Cluster)
print(pwc.CompactnessStd)

p5_3 = p5_2 +
        stat_pvalue_manual(pwc.CompactnessStd[pwc.CompactnessStd$term == "Cluster", ], 
                           label = "P = {p.adj}", 
                           y.position = 0.5, 
                           step.increase = 0.1, 
                           tip.length = 0, 
                           hide.ns = TRUE)

p5_3

# Eccentricity_Mean
p6_1 = p1 + aes(x = Cluster, y = Eccentricity_Mean) +
        geom_violin(aes(fill = Cluster), trim = FALSE) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p6_1

p6_2 = p6_1 + 
        geom_point(data = df.joined.hc4.mean,
                   aes(x = Cluster, y = Eccentricity_Mean),
                   size = 2) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p6_2

aov.EccentricityMean =  anova_test(df.joined.hc4, Eccentricity_Mean ~ Cluster)
print(aov.EccentricityMean)

pwc.EccentricityMean = tukey_hsd(df.joined.hc4, Eccentricity_Mean ~ Cluster)
print(pwc.EccentricityMean)

p6_3 = p6_2 +
        stat_pvalue_manual(pwc.EccentricityMean[pwc.EccentricityMean$term == "Cluster", ], 
                           label = "P = {p.adj}", 
                           y.position = 0.9, 
                           step.increase = 0.1, 
                           tip.length = 0, 
                           hide.ns = TRUE)

p6_3

# FormFactor_Std
p7_1 = p1 + aes(x = Cluster, y = FormFactor_Std) +
        geom_violin(aes(fill = Cluster), trim = FALSE) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p7_1

p7_2 = p7_1 + 
        geom_point(data = df.joined.hc4.mean,
                   aes(x = Cluster, y = FormFactor_Std),
                   size = 2) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p7_2

aov.FormFactorStd =  anova_test(df.joined.hc4, FormFactor_Std ~ Cluster)
print(aov.FormFactorStd)

pwc.FormFactorStd = tukey_hsd(df.joined.hc4, FormFactor_Std ~ Cluster)
print(pwc.FormFactorStd)

p7_3 = p7_2 +
        stat_pvalue_manual(pwc.FormFactorStd[pwc.FormFactorStd$term == "Cluster", ], 
                           label = "P = {p.adj}", 
                           y.position = 0.3, 
                           step.increase = 0.1, 
                           tip.length = 0, 
                           hide.ns = TRUE)

p7_3

# Speed_Max
p8_1 = p1 + aes(x = Cluster, y = Speed_Max) +
        geom_violin(aes(fill = Cluster), trim = FALSE) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p8_1

p8_2 = p8_1 + 
        geom_point(data = df.joined.hc4.mean,
                   aes(x = Cluster, y = Speed_Max),
                   size = 2) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p8_2

aov.SpeedMax =  anova_test(df.joined.hc4, Speed_Max ~ Cluster)
print(aov.SpeedMax)

pwc.SpeedMax = tukey_hsd(df.joined.hc4, Speed_Max ~ Cluster)
print(pwc.SpeedMax)

p8_3 = p8_2 +
        stat_pvalue_manual(pwc.SpeedMax[pwc.SpeedMax$term == "Cluster", ], 
                           label = "P = {p.adj}", 
                           y.position = 2, 
                           step.increase = 0.1, 
                           tip.length = 0, 
                           hide.ns = TRUE)

p8_3

# Speed_Std
p9_1 = p1 + aes(x = Cluster, y = Speed_Std) +
        geom_violin(aes(fill = Cluster), trim = FALSE) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p9_1

p9_2 = p9_1 + 
        geom_point(data = df.joined.hc4.mean,
                   aes(x = Cluster, y = Speed_Std),
                   size = 2) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p9_2

aov.SpeedStd =  anova_test(df.joined.hc4, Speed_Std ~ Cluster)
print(aov.SpeedStd)

pwc.SpeedStd = tukey_hsd(df.joined.hc4, Speed_Std ~ Cluster)
print(pwc.SpeedStd)

p9_3 = p9_2 +
        stat_pvalue_manual(pwc.SpeedStd[pwc.SpeedStd$term == "Cluster", ], 
                           label = "P = {p.adj}", 
                           y.position = 0.4, 
                           step.increase = 0.1, 
                           tip.length = 0, 
                           hide.ns = TRUE)

p9_3

# IntegratedDistance
p10_1 = p1 + aes(x = Cluster, y = IntegratedDistance) +
        geom_violin(aes(fill = Cluster), trim = FALSE) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p10_1

p10_2 = p10_1 + 
        geom_point(data = df.joined.hc4.mean,
                   aes(x = Cluster, y = IntegratedDistance),
                   size = 2) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p10_2

aov.IntegratedDistance =  anova_test(df.joined.hc4, IntegratedDistance ~ Cluster)
print(aov.IntegratedDistance)

pwc.IntegratedDistance = tukey_hsd(df.joined.hc4, IntegratedDistance ~ Cluster)
print(pwc.IntegratedDistance)

p10_3 = p10_2 +
        stat_pvalue_manual(pwc.IntegratedDistance[pwc.IntegratedDistance$term == "Cluster", ], 
                           label = "P = {p.adj}", 
                           y.position = 100, 
                           step.increase = 0.1, 
                           tip.length = 0, 
                           hide.ns = TRUE)

p10_3

# Linearity
p11_1 = p1 + aes(x = Cluster, y = Linearity) +
        geom_violin(aes(fill = Cluster), trim = FALSE) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p11_1

p11_2 = p11_1 + 
        geom_point(data = df.joined.hc4.mean,
                   aes(x = Cluster, y = Linearity),
                   size = 2) +
        theme(axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.x = element_blank())
p11_2

aov.Linearity =  anova_test(df.joined.hc4, Linearity ~ Cluster)
print(aov.Linearity)

pwc.Linearity = tukey_hsd(df.joined.hc4, Linearity ~ Cluster)
print(pwc.Linearity)

p11_3 = p11_2 +
        stat_pvalue_manual(pwc.Linearity[pwc.Linearity$term == "Cluster", ], 
                           label = "P = {p.adj}", 
                           y.position = 100, 
                           step.increase = 0.1, 
                           tip.length = 0, 
                           hide.ns = TRUE)

p11_3

# 全パラメータのviolin plotをまとめて表示
plot_grid(p2_3 + theme(legend.position = "none"), 
          p3_3 + theme(legend.position = "none"), 
          p4_3 + theme(legend.position = "none"), 
          p5_3 + theme(legend.position = "none"), 
          p6_3 + theme(legend.position = "none"), 
          p7_3 + theme(legend.position = "none"), 
          p8_3 + theme(legend.position = "none"), 
          p9_3 + theme(legend.position = "none"), 
          p10_3 + theme(legend.position = "none"), 
          p11_3 + theme(legend.position = "right", 
                        legend.justification = "top"), 
          nrow = 3, 
          ncol = 4,
          rel_widths = c(1, 1.2, 1, 1))
ggsave("violin.best.pwc.tiff", width = 20, height = 20)          
