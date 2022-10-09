 library(DESeq2)

# sugarcane RNA-seq data

#link to the paper : https://www.nature.com/articles/s41598-019-45184-1


count_matrix <- as.matrix(read.csv("https://reneshbedre.github.io/assets/posts/gexp/df_sc.csv", row.names = "gene"))

count_matrix

count_matrix <- count_matrix[, -7]

dim(count_matrix )
head(count_matrix)

colnames(count_matrix)

class( count_matrix)


ColData <- data.frame(condition= c( rep("ctr",3),rep("trt",3)))
ColData

ColData$condition <- factor( ColData$condition )


ddso <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = ColData,
                              design = ~ condition)

counts(ddso)

length(which(rowSums(counts(ddso)) < 10)) # 31 low count genes




