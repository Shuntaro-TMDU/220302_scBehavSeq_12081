setwd("/Users/shuntaro/Documents/R/CW/230302_scBehavSeq_12081/R_output")
getwd()


library(tidyverse)
library(cowplot)
library(scales)

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

# df.meta.1にクラスター列df.annotationを結合し、Y座標を上下反転させたCenter_Yを作成
df.meta.hc2.2 = df.meta.1 %>%
        full_join(., df.annotation, by = "Order") %>%
        mutate(Center_Y = 512 - Center_Y) %>%
        print()

View(df.meta.hc2.2)

# apply cluster to XY field
## NonThrombosis, Run7, ImageNumber1
NonThro_7_0 = df.meta.hc2.2 %>%
        filter(Run == 7 & ImageNumber == 1) %>%
        print()

BackGate.NonThro_7_0 = ggplot() +
        geom_point(data = NonThro_7_0, 
                   aes(x = Center_X, y = Center_Y, color = Cluster), 
                   size = 2) +
        mytheme1 + 
        scale_x_continuous(limits = c(0, 512)) +
        scale_y_continuous(limits = c(0, 512)) +
        theme(axis.ticks = element_blank(), 
              axis.text = element_blank())
BackGate.NonThro_7_0

ggsave("backgate_nonthro_7_t1.tiff", width = 6, height = 6)

## NonThrombosis, Run7
NonThro_7 = df.meta.new %>%
        filter(Run == 7) %>%
        # mutate(Label = as.factor(Label)) %>%
        mutate(Center_Y = 512 - Center_Y) %>%
        print()

BackGate.NonThro_7 = BackGate.NonThro_7_0 + 
        geom_path(data = NonThro_7, 
                  aes(x = Center_X, y = Center_Y, group = Label),
                  color = "black") 

BackGate.NonThro_7

?geom_line
ggsave("backgate_nonthro_7.tiff", width = 6, height = 6)

## NonThrombosis, Run8, ImageNumber1
NonThro_8 = df.meta.hc2.2 %>%
        filter(Run == 8 & ImageNumber == 1) %>%
        print()

BackGate.NonThro_8 = ggplot(NonThro_8) +
        aes(x = Center_X, y = Center_Y) +
        geom_point(aes(color = Cluster), 
                   size = 2) +
        mytheme1 + 
        scale_x_continuous(limits = c(0, 512)) +
        scale_y_continuous(limits = c(0, 512)) +
        theme(axis.ticks = element_blank(), 
              axis.text = element_blank())
BackGate.NonThro_8

ggsave("backgate_nonthro_8_t1.tiff", width = 6, height = 6)

## Thrombosis, Run3, ImageNumber1
Thro_3_0 = df.meta.hc2.2 %>%
        filter(Run == 3 & ImageNumber == 1) %>%
        print()

BackGate.Thro_3_0 = ggplot() +
        geom_point(data = Thro_3_0, 
                   aes(x = Center_X, y = Center_Y, color = Cluster), 
                   size = 2) +
        mytheme1 + 
        scale_x_continuous(limits = c(0, 512)) +
        scale_y_continuous(limits = c(0, 512)) +
        theme(axis.ticks = element_blank(), 
              axis.text = element_blank())
BackGate.Thro_3_0

ggsave("backgate_thro_3_t1.tiff", width = 6, height = 6)

## Thrombosis, Run3
Thro_3 = df.meta.new %>%
        filter(Run == 3) %>%
        mutate(Center_Y = 512 - Center_Y) %>%
        print()

BackGate.Thro_3 = BackGate.Thro_3_0 + 
        geom_path(data = Thro_3, 
                  aes(x = Center_X, y = Center_Y, group = Label),
                  color = "black") 

BackGate.Thro_3

?geom_line
ggsave("backgate_thro_3.tiff", width = 6, height = 6)

# NonThro_7, Thro_3を横並びで表示
## タイトルあり
plot_grid(BackGate.NonThro_7 + 
                  labs(title = "Nonthrombosis") +
                  theme(legend.position = "none"), 
          BackGate.Thro_3 + 
                  labs(title = "Thrombosis") +
                  theme(legend.position = "right", 
                                          legend.justification = "top"), 
          nrow = 1, 
          ncol = 2,
          rel_widths = c(1, 1.2))
ggsave("backgate_nonthro_7_thro_3_withTitle.tiff", width = 12, height = 6)  

## タイトルなし
plot_grid(BackGate.NonThro_7 + 
                  theme(legend.position = "none"), 
          BackGate.Thro_3 + 
                  theme(legend.position = "right", 
                        legend.justification = "top"), 
          nrow = 1, 
          ncol = 2,
          rel_widths = c(1, 1.2))
ggsave("backgate_nonthro_7_thro_3.tiff", width = 12, height = 5.5)  

# Pathology毎のそれぞれのClusterの数を調べる
df.meta.hc2.2 %>%
        count(Pathology)

# A tibble: 2 × 2
# Pathology         n
# <chr>         <int>
#         1 NonThrombosis    51
# 2 Thrombosis       51

ggplot(df.meta.hc2.2, aes(Pathology)) +
        geom_bar(aes(fill = Cluster)) +
        mytheme1 + 
        theme(axis.title.x = element_blank())

ggplot(df.meta.hc2.2, aes(Pathology)) +
        geom_bar(aes(fill = Cluster),
                 width = 0.75,
                 position = "fill") +
        coord_flip() +
        mytheme1 + 
        scale_y_continuous(labels = percent) + # library(scales)の実行が必要
        theme(axis.text.x = element_text(hjust = 1), 
              axis.title.y = element_blank(), 
              legend.position = "right", 
              legend.justification = "top")

ggsave("percent_m2_hc2.tiff", width = 8, height = 4)

# p1 = ggplot(df.meta.hc4) +
#         mytheme1 
# # +
# #         theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
# 
# p2 = p1 + aes(x = Condition) +
#         geom_bar(aes(fill = Cluster)) +
#         theme(legend.position = "none",
#               axis.text.x = element_text(angle = 45, hjust = 1), 
#               axis.title.x = element_blank())
# p2
# 
# p3 = p1 + aes(x = Condition) +
#         geom_bar(aes(fill = Cluster), 
#                  position = "fill") +
#         scale_y_continuous(labels = percent) +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1), 
#               axis.title.x = element_blank(), 
#               axis.title.y = element_blank(), 
#               legend.position = "right", 
#               legend.justification = "top")
# p3
# 
# plot_grid(p2, p3, 
#           ncol = 2,
#           rel_widths = c(1, 1.5))
# 
# ggsave("Count.Percent.hc4.tiff", width = 8, height = 6)
