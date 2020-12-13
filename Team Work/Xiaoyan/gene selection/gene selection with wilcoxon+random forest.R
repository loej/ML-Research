# repeat p value = 0.05 5 times to get overlap
setwd("C:/Users/Lenovo/Desktop/Yu Lab/gene selection")
install.packages("readxl")
library(readxl)
data = read_excel("gene.xlsx") # the path to the file

## restructure data frame
df = data.frame(t(data))
df = df[2:13,2:length(df)]
df = lapply(df,as.character)
df = data.frame(lapply(df,as.numeric))
df$y = c(-1,-1,-1,-1,-1,-1,1,1,1,1,1,1)
df$y = as.factor(df$y)


# step one: wlcoxon test
# maybe the lack of certain gene will result in excessive drinkin????
pval = c()
index = c()
meand = c()
meann = c()
# conduct wilcoxon test for each gene (or variable)
for (i in 1:(length(df)-1)){
  a = wilcox.test(df[df$y == "1",i],df[df$y == "-1",i],alternative = "greater", )
  pval = c(pval,a$p.value)
  index = c(index,i)
  meand = c(meand,mean(df[df$y == "1",i]))
  meann = c(meann,mean(df[df$y == "-1",i]))
}
result = data.frame(cbind(index,pval,meand, meann))
re = result[result$pval<0.05,] # pick p-value as 0.05 to include more variables
dim(re)[1] # number of variables with significance difference
df1 = df[,re$index]

# continue with random forest
install.packages("randomForest")
library(randomForest)
df1$y = c(0,0,0,0,0,0,1,1,1,1,1,1) # organeze outcome variable
df1$y = as.factor(df1$y)
str(df1)
r = randomForest(y~.,data = df1, na.action = na.pass, ntree = 500, importance = TRUE) # run the model
importance = data.frame(r$importance) # get variable importance
View(importance)
im = importance[importance$MeanDecreaseGini>0 & importance$MeanDecreaseAccuracy>0,] # get variable with importance > 0
variable = rownames(im) # get variable index

