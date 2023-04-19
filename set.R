setwd("/Users/shuntaro/Documents/R/CW/230302_scBehavSeq_12081/R_input")
getwd()

sessionInfo()

library(tidyverse)
library(assertthat)

df1_1 = read_csv("12081_ATG_NonThrombosis_7_Ly6G.csv")
df1_2 = read_csv("12081_ATG_NonThrombosis_8_Ly6G.csv")
df1_3 = read_csv("12081_ATG_Thrombosis_1_Ly6G.csv")
df1_4 = read_csv("12081_ATG_Thrombosis_3_Ly6G.csv")
df1_5 = read_csv("12081_ATG_Thrombosis_4_Ly6G.csv")
df1_6 = read_csv("12081_ATG_Thrombosis_5_Ly6G.csv")

df1 = rbind(df1_1, df1_2, df1_3, df1_4, df1_5, df1_6) %>%
        print()

View(df1)

df2 = df1 %>%
        select(!c(Metadata_Run...10, Location_Center_X, Location_Center_Y, Location_Center_Z, TrackObjects_Area)) %>%
        rename_with( ~ str_extract(., '(?<=_).+'), .cols = contains('_')) %>%
        rename_with( ~ str_remove(., '...9')) %>%
        mutate(Species = factor(Pathology, levels = (c("NonThrombosis", "Thrombosis")))) %

# LifeTime < 2や細胞凝集体を取り除く
df3 = df2 %>% 
        filter(Label != "NaN", 
               !(Run == 3 & Label == 5)) %>%
        print()
View(df3)

# Speed, DistanceTraveled, Linearity, MeanderRatioを集計したdf.K.1を作成
ScanFreq = 3.25

df.K.1 = df3 %>%
        filter(Lifetime != 1) %>%
        group_by(Run, Label) %>%
        summarise(DistanceTraveled_Max = max(DistanceTraveled), 
                  DistanceTraveled_Mean = mean(DistanceTraveled), 
                  DistanceTraveled_Std = sd(DistanceTraveled), 
                  Speed_Max = max(DistanceTraveled / ScanFreq), 
                  Speed_Mean = mean(DistanceTraveled / ScanFreq), 
                  Speed_Std = sd(DistanceTraveled / ScanFreq),
                  Linearity_Mean = mean(Linearity), 
                  Linearity_Std = sd(Linearity), 
                  MeanderRatio_Mean = mean(1 / Linearity), 
                  MeanderRatio_Std = sd(1 / Linearity),
                  .groups = 'drop') 

View(df.K.1)

# Displacement, を集計したdf.Kを作成
## 練習
df.test = df3 %>%
        filter(Run == 7)

View(df.test)

df.test.K.2 = df.test %>%
        filter(Lifetime != 1, 
               FinalAge != "NaN") %>%
        select(Run, 
               Label, 
               FinalAge,
               Displacement, 
               IntegratedDistance) %>%
        mutate(Displacement = max(Displacement), 
               IntegratedDistance = max(IntegratedDistance)) # Displacement, IntegratedDistanceは複数あるため、最大値を採用

head(df.test.K.2)

## 本番
df.K.2 = df3 %>%
        filter(Lifetime != 1, 
               FinalAge != "NaN") %>%
        group_by(Run, Label) %>%
        select(Run, 
               Label, 
               FinalAge,
               Displacement, 
               IntegratedDistance) %>%
        mutate(Displacement = max(Displacement), 
               IntegratedDistance = max(IntegratedDistance)) %>%
        distinct(Run, Label, .keep_all = TRUE)

head(df.K.2)


## そのほかの統計量を集計したdf.Mを作成
## おおむね全ての形質学的パラメータを抽出したdf.M1を作成
M1 = c( 
        "Area",
        "BoundingBoxArea",
        "CentralMoment_0_0", 
        "CentralMoment_0_1", 
        "CentralMoment_0_2",
        "CentralMoment_1_1",
        "CentralMoment_1_2",
        "CentralMoment_2_1",
        "CentralMoment_2_2",
        "Compactness",
        "ConvexArea",
        "Eccentricity",
        "EquivalentDiameter",
        "Extent",
        "FormFactor",
        "HuMoment_0",
        "HuMoment_1",
        "HuMoment_2",
        "InertiaTensorEigenvalues_0",
        "InertiaTensorEigenvalues_1", 
        "InertiaTensor_0_0", 
        "InertiaTensor_0_1",
        "InertiaTensor_1_0",
        "InertiaTensor_1_1",
        "MajorAxisLength",
        "MaxFeretDiameter",
        "MaximumRadius",
        "MeanRadius",
        "MedianRadius",
        "MinFeretDiameter",
        "MinorAxisLength",
        "NormalizedMoment_0_2",
        "NormalizedMoment_1_1",
        "NormalizedMoment_1_2",
        "NormalizedMoment_2_0",
        "NormalizedMoment_2_1",
        "NormalizedMoment_2_2",
        "Perimeter",
        "Solidity",
        "Zernike_0_0",
        "Zernike_1_1",
        "Zernike_2_0",
        "Zernike_2_2")

M2 = c( 
        "Area",
        "BoundingBoxArea",
        "Compactness",
        "ConvexArea",
        "Eccentricity",
        "EquivalentDiameter",
        "Extent",
        "FormFactor",
        "MajorAxisLength",
        "MaxFeretDiameter",
        "MaximumRadius",
        "MeanRadius",
        "MedianRadius",
        "MinFeretDiameter",
        "MinorAxisLength",
        "Perimeter",
        "Solidity")

M1[1]

Make.df.M = function(data = df3, 
                     ParamM = M1){
        df.M = data %>%
                select(Run, Label, all_of(ParamM)) %>%
                group_by(Run, Label) %>%
                summarise(across(ParamM[1]:ParamM[length(ParamM)],
                                 mean,
                                 .names = "{.col}_Mean"),
                          across(ParamM[1]:ParamM[length(ParamM)],
                                 sd,
                                 .names = "{.col}_Std"), 
                          .groups = 'drop')
        
        return(df.M)
}

df.M1 = Make.df.M(data = df3, 
                  ParamM = M1)

View(df.M1)

df.M2 = Make.df.M(data = df3, 
                  ParamM = M2)

## データの結合
Join.df.M.K = function(df.M = df.M1, 
                       df.K.1, 
                       df.K.2){
        df.joined = inner_join(df.M, df.K.1, 
                               by = c("Run", "Label")) %>%
                inner_join(., df.K.2, 
                           by = c("Run", "Label"))
        
        return(df.joined)
}

df.joined.1 = Join.df.M.K(df.M = df.M1, 
                          df.K.1, 
                          df.K.2)

df.joined.2 = Join.df.M.K(df.M = df.M2, 
                          df.K.1, 
                          df.K.2)

View(df.joined)  

assert_that(noNA(df.joined.1))
subset(df.joined.1, complete.cases(df.joined.1)==F)

# 欠損値を含むオブジェクト抽出→Run: 3, Label: 17
missID.1 = df.joined.1 %>%
        filter(complete.cases(.)==F)

missID.2 = df.joined.2 %>%
        filter(complete.cases(.)==F)

View(missID.1)
View(missID.2)

# 欠損値を含むオブジェクトの行を削除
df.joined.1.new = df.joined.1 %>%
        drop_na(everything())

df.joined.2.new = df.joined.2 %>%
        drop_na(everything())

# Metadataを作成
df.meta = df3 %>%
        select(ID, Drug, Pathology, Run, Stain, ImageNumber, Label, Center_X, Center_Y, Lifetime) %>%
        arrange(Run, Label) %>% # Parameterデータフレーム作成時の並べ替え法則に則る
        print()

View(df.meta)

# 欠損値を含むオブジェクトの行を削除
df.meta.new = df.meta %>%
        filter(!(Run == 3 & Label == 17))

setwd("/Users/shuntaro/Documents/R/CW/230302_scBehavSeq_12081/R_output")
write_csv(df.meta.new, "12081_ATG_Ly6G_Metadata.csv")
write_csv(df.joined.1.new, "12081_ATG_Ly6G_Parameters_M1.csv")
write_csv(df.joined.2.new, "12081_ATG_Ly6G_Parameters_M2.csv")
