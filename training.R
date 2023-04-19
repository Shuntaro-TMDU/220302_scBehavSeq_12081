library(caret)
data(iris)
x <- iris[, 1:4]
y <- iris$Species

x
class(y)

# データをトレーニングデータとテストデータに分割する
set.seed(123) # 乱数のシードを設定
train_indices <- createDataPartition(iris$Species, p = 0.7, list = FALSE)
train_data <- iris[train_indices, ]
test_data <- iris[-train_indices, ]

# ロジスティック回帰モデルを作成する
model <- train(Species ~ ., data = train_data, method = "glmnet", trControl = trainControl(method = "cv"), 
               family = "binomial", tuneLength = 10, preProcess = c("center", "scale"), 
               metric = "Accuracy", alpha = 1)
