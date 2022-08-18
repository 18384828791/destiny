#��ȡ OTUs ��ȱ�
otu <- read.table('otu_table.txt', sep = '\t', row.names = 1, header = TRUE, fill = TRUE)

#���˵ͷ�� OTUs ��Ⱥ�����ǶԷ��๱�׶ȵͣ���Ӱ�����Ч��
#120 ���������Ͱ� OTUs ��ȵ��кͲ�С�� 120 Ϊ׼��
otu <- otu[which(rowSums(otu) >= 120), ]

#�ϲ����飬�õ��ܹ��� randomForest ʶ�����ĸ�ʽ
group <- read.table('group.txt', sep = '\t', row.names = 1, header = TRUE, fill = TRUE)
otu <- data.frame(t(otu))
otu_group <- cbind(otu, group)

#�������ݼ���Ϊѵ������ռ 70%���Ͳ��Լ���ռ 30%��
set.seed(123)
select_train <- sample(120, 120*0.7)
otu_train <- otu_group[select_train, ]
otu_test <- otu_group[-select_train, ]

####################################################
#randomForest �������ɭ��
library(randomForest)

#���ɭ�ּ��㣨Ĭ������ 500 �þ������������� ?randomForest
set.seed(123)
otu_train.forest <- randomForest(groups ~ ., data = otu_train, importance = TRUE)
otu_train.forest

plot(margin(otu_train.forest, otu_train$groups), main = '�۲�ֵ���ж���ȷ�ĸ���ͼ')

#ѵ�����������
train_predict <- predict(otu_train.forest, otu_train)
compare_train <- table(train_predict, otu_train$groups)
compare_train
sum(diag(compare_train)/sum(compare_train))

#ʹ�ò��Լ�����
test_predict <- predict(otu_train.forest, otu_test)
compare_test <- table(otu_test$groups, test_predict, dnn = c('Actual', 'Predicted'))
compare_test

###�ؼ� OTUs ʶ��
#�鿴��ʾÿ��������OTUs����Ҫ�Եĵ÷�
#summary(otu_train.forest)
importance_otu <- otu_train.forest$importance
head(importance_otu)

#����ʹ�ú��� importance()
importance_otu <- data.frame(importance(otu_train.forest))
head(importance_otu)

#���Ը���ĳ����Ҫ�Եĸߵ��Ÿ���������ݡ�Mean Decrease Accuracy��ָ��
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_otu)

#������
#write.table(importance_otu, 'importance_otu.txt', sep = '\t', col.names = NA, quote = FALSE)

#��ͼչʾ top30 ��Ҫ�� OTUs
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)), main = 'Top 30 - variable importance')

###������֤����ѡ���ض������� OTUs
#5 ���ظ�ʮ�۽�����֤
set.seed(123)
otu_train.cv <- replicate(5, rfcv(otu_train[-ncol(otu_train)], otu_train$group, cv.fold = 10,step = 1.5), simplify = FALSE)
otu_train.cv

#��ȡ��֤�����ͼ
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))

#�����ͼ
library(ggplot2)
library(splines)  #������ geom_smooth() ���������ߣ�����ʹ�� geom_line() ��� geom_smooth() ������ͨ����

p <- ggplot(otu_train.cv, aes(otus, value)) +
    geom_smooth(se = FALSE,	method = 'glm', formula = y~ns(x, 6)) +
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +  
    labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')

p

#��Լ��ȡǰ 30 ����Ҫ�� OTUs
p + geom_vline(xintercept = 30)

#���� OTUs ��Ҫ�������ѡ��������ݡ�Mean Decrease Accuracy��ָ��
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_otu)

#������
#write.table(importance_otu[1:30, ], 'importance_otu_top30.txt', sep = '\t', col.names = NA, quote = FALSE)

###��Լ������
#ѡ�� top30 ��Ҫ�� OTUs�����������Ѿ����ݡ�Mean Decrease Accuracy���������
otu_select <- rownames(importance_otu)[1:30]

#�����Ӽ���ѵ�����Ͳ��Լ�
otu_train_top30 <- otu_train[ ,c(otu_select, 'groups')]
otu_test_top30 <- otu_test[ ,c(otu_select, 'groups')]

#���ɭ�ּ��㣨Ĭ������ 500 �þ������������� ?randomForest
set.seed(123)
otu_train.forest_30 <- randomForest(groups ~ ., data = otu_train_top30, importance = TRUE)
otu_train.forest_30

plot(margin(otu_train.forest_30, otu_test_top30$groups), main = '�۲�ֵ���ж���ȷ�ĸ���ͼ')

#ѵ�����������
train_predict <- predict(otu_train.forest_30, otu_train_top30)
compare_train <- table(train_predict, otu_train_top30$groups)
compare_train

#ʹ�ò��Լ�����
test_predict <- predict(otu_train.forest_30, otu_test_top30)
compare_test <- table(otu_test_top30$groups, test_predict, dnn = c('Actual', 'Predicted'))
compare_test

##NMDS ����ͼ��չʾ����
#NMDS ��ά
nmds <- vegan::metaMDS(otu, distance = 'bray') 
result <- nmds$points
result <- as.data.frame(cbind(result, rownames(result)))

#��������ķ���Ԥ����
predict_group <- c(train_predict, test_predict)
predict_group <- as.character(predict_group[rownames(result)])

#��ͼ
colnames(result)[1:3] <- c('NMDS1', 'NMDS2', 'samples')
result$NMDS1 <- as.numeric(as.character(result$NMDS1)) 
result$NMDS2 <- as.numeric(as.character(result$NMDS2))
result$samples <- as.character(result$samples)
result <- cbind(result, predict_group)
head(result)

ggplot(result, aes(NMDS1, NMDS2, color = predict_group)) +  
    geom_polygon(data = plyr::ddply(result, 'predict_group', function(df) df[chull(df[[1]], df[[2]]), ]), fill = NA) +
    geom_point()

