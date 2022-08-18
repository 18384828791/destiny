#读取 OTUs 丰度表
otu <- read.table('otu_table.txt', sep = '\t', row.names = 1, header = TRUE, fill = TRUE)

#过滤低丰度 OTUs 类群，它们对分类贡献度低，且影响计算效率
#120 个样本，就按 OTUs 丰度的行和不小于 120 为准吧
otu <- otu[which(rowSums(otu) >= 120), ]

#合并分组，得到能够被 randomForest 识别计算的格式
group <- read.table('group.txt', sep = '\t', row.names = 1, header = TRUE, fill = TRUE)
otu <- data.frame(t(otu))
otu_group <- cbind(otu, group)

#将总数据集分为训练集（占 70%）和测试集（占 30%）
set.seed(123)
select_train <- sample(120, 120*0.7)
otu_train <- otu_group[select_train, ]
otu_test <- otu_group[-select_train, ]

####################################################
#randomForest 包的随机森林
library(randomForest)

#随机森林计算（默认生成 500 棵决策树），详情 ?randomForest
set.seed(123)
otu_train.forest <- randomForest(groups ~ ., data = otu_train, importance = TRUE)
otu_train.forest

plot(margin(otu_train.forest, otu_train$groups), main = '观测值被判断正确的概率图')

#训练集自身测试
train_predict <- predict(otu_train.forest, otu_train)
compare_train <- table(train_predict, otu_train$groups)
compare_train
sum(diag(compare_train)/sum(compare_train))

#使用测试集评估
test_predict <- predict(otu_train.forest, otu_test)
compare_test <- table(otu_test$groups, test_predict, dnn = c('Actual', 'Predicted'))
compare_test

###关键 OTUs 识别
#查看表示每个变量（OTUs）重要性的得分
#summary(otu_train.forest)
importance_otu <- otu_train.forest$importance
head(importance_otu)

#或者使用函数 importance()
importance_otu <- data.frame(importance(otu_train.forest))
head(importance_otu)

#可以根据某种重要性的高低排个序，例如根据“Mean Decrease Accuracy”指标
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_otu)

#输出表格
#write.table(importance_otu, 'importance_otu.txt', sep = '\t', col.names = NA, quote = FALSE)

#作图展示 top30 重要的 OTUs
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)), main = 'Top 30 - variable importance')

###交叉验证帮助选择特定数量的 OTUs
#5 次重复十折交叉验证
set.seed(123)
otu_train.cv <- replicate(5, rfcv(otu_train[-ncol(otu_train)], otu_train$group, cv.fold = 10,step = 1.5), simplify = FALSE)
otu_train.cv

#提取验证结果绘图
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))

#拟合线图
library(ggplot2)
library(splines)  #用于在 geom_smooth() 中添加拟合线，或者使用 geom_line() 替代 geom_smooth() 绘制普通折线

p <- ggplot(otu_train.cv, aes(otus, value)) +
    geom_smooth(se = FALSE,	method = 'glm', formula = y~ns(x, 6)) +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +  
    labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')

p

#大约提取前 30 个重要的 OTUs
p + geom_vline(xintercept = 30)

#根据 OTUs 重要性排序后选择，例如根据“Mean Decrease Accuracy”指标
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_otu)

#输出表格
#write.table(importance_otu[1:30, ], 'importance_otu_top30.txt', sep = '\t', col.names = NA, quote = FALSE)

###简约分类器
#选择 top30 重要的 OTUs，例如上述已经根据“Mean Decrease Accuracy”排名获得
otu_select <- rownames(importance_otu)[1:30]

#数据子集的训练集和测试集
otu_train_top30 <- otu_train[ ,c(otu_select, 'groups')]
otu_test_top30 <- otu_test[ ,c(otu_select, 'groups')]

#随机森林计算（默认生成 500 棵决策树），详情 ?randomForest
set.seed(123)
otu_train.forest_30 <- randomForest(groups ~ ., data = otu_train_top30, importance = TRUE)
otu_train.forest_30

plot(margin(otu_train.forest_30, otu_test_top30$groups), main = '观测值被判断正确的概率图')

#训练集自身测试
train_predict <- predict(otu_train.forest_30, otu_train_top30)
compare_train <- table(train_predict, otu_train_top30$groups)
compare_train

#使用测试集评估
test_predict <- predict(otu_train.forest_30, otu_test_top30)
compare_test <- table(otu_test_top30$groups, test_predict, dnn = c('Actual', 'Predicted'))
compare_test

##NMDS 排序图中展示分类
#NMDS 降维
nmds <- vegan::metaMDS(otu, distance = 'bray') 
result <- nmds$points
result <- as.data.frame(cbind(result, rownames(result)))

#获得上述的分类预测结果
predict_group <- c(train_predict, test_predict)
predict_group <- as.character(predict_group[rownames(result)])

#作图
colnames(result)[1:3] <- c('NMDS1', 'NMDS2', 'samples')
result$NMDS1 <- as.numeric(as.character(result$NMDS1)) 
result$NMDS2 <- as.numeric(as.character(result$NMDS2))
result$samples <- as.character(result$samples)
result <- cbind(result, predict_group)
head(result)

ggplot(result, aes(NMDS1, NMDS2, color = predict_group)) +  
    geom_polygon(data = plyr::ddply(result, 'predict_group', function(df) df[chull(df[[1]], df[[2]]), ]), fill = NA) +
    geom_point()

