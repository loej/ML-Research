# repeat p value = 0.05 5 times to get overlap
library(readxl)
data = read_excel("gene.xlsx")

## restructure data frame
df = data.frame(t(data))
df = df[2:13,2:length(df)]
df = lapply(df,as.character)
df = data.frame(lapply(df,as.numeric))
df$y = c(-1,-1,-1,-1,-1,-1,1,1,1,1,1,1)
df$y = as.factor(df$y)
#df$Y = c(1,1,1,2,2,2,3,3,3,4,4,4)


# step one: wlcoxon test
# maybe the lack of certain gene will result in excessive drinkin????
pval = c()
index = c()
meand = c()
meann = c()
for (i in 1:(length(df)-1)){
  a = wilcox.test(df[df$y == "1",i],df[df$y == "-1",i],alternative = "greater")
  pval = c(pval,a$p.value)
  index = c(index,i)
  meand = c(meand,mean(df[df$y == "1",i]))
  meann = c(meann,mean(df[df$y == "-1",i]))
}
result = data.frame(cbind(index,pval,meand, meann))
re = result[result$pval<0.05,] # pick p-value as 0.05 to include more variables
d = dim(re)[1]
df1 = df[,re$index]

#scad svm
install.packages("penalizedSVM")
library(penalizedSVM)
model <- scadsvc(as.matrix(df1), y=df$y, lambda=0.01, seed = 123)
print(str(model))
ind = model$index
model1 <- scadsvc(as.matrix(df[1:(dim(df)[2]-1)]), y=df$y, lambda=0.01, seed = 123)
str(model1)
# continue with random forest
library(randomForest)
df2 = df1[,model$xind]
df2$y = c(0,0,0,0,0,0,1,1,1,1,1,1)
df2$y = as.factor(df2$y)
str(df2)
r = randomForest(y~.,data = df2, na.action = na.pass, ntree = 500, importance = TRUE)
importance = data.frame(r$importance)
importance
im = importance[importance$MeanDecreaseGini>0 & importance$MeanDecreaseAccuracy>0,]
variable = rownames(im)
##########################end########################
## step two: random forest
library(randomForest)
df1 = df[,re$index]
df1 = data.frame(scale(as.matrix(df1)))
df1$Y = c(0,0,0,0,0,0,1,1,1,1,1,1)
df1$Y = as.factor(df1$Y)
# repeat random forest 5 times to get overlap result since variables are selected randomly
# tke the first random forest as initialization
r = randomForest(Y~.,data = df1, na.action = na.pass, ntree = 4000, importance = TRUE)
importance = data.frame(r$importance)
im = importance[importance$MeanDecreaseGini>0 & importance$MeanDecreaseAccuracy>0,]
variable = rownames(im)
for (i in 2:5){
  r = randomForest(Y~.,data = df1, na.action = na.pass, ntree = 4000, importance = TRUE)
  importance = data.frame(r$importance)
  im = importance[importance$MeanDecreaseGini>0 & importance$MeanDecreaseAccuracy>0,]
  variable = intersect(variable, rownames(im)) # contain significant gene that are in every random forest model
  print(length(variable))
}

## step 3: group 3 mean/group 2 mean, group 4 mean/group 2 mean and find extreme ratio
df$X11573 = NULL # remove gene X11573 since all data is the same accross four groups
d1 = df[1:3,]
d2 = df[4:6,]
d3 = df[7:9,]
d4 = df[10:12,]
# mean ratio group 3 vs group 4
ratio = data.frame(name = colnames(d2),m1 = colMeans(d1),
                   m2=colMeans(d2),m3=colMeans(d3),m4=colMeans(d4))
ratio$r32 = ratio$m3/ratio$m2
ratio$r42 = ratio$m4/ratio$m2

# 20000+ genes have group 3 or 4 mean greater than group 2 mean
dim(ratio[ratio$r32>1,])
dim(ratio[ratio$r42>1,])

# ratio mean + 3SD and overlap with random forest result
# group 3
extreme32 = ratio[ratio$r32 > mean(ratio$r32)+3*sd(ratio$r32),]
inter32 = intersect(extreme32$name,variable)
df[,inter32] # individual gene expression
# group 4
extreme42 = ratio[ratio$r42 > mean(ratio$r42)+3*sd(ratio$r42),]
inter42 = intersect(extreme42$name,variable)
df[,inter42] # individual gene expression

# gene index in inter32 and inter42 is different from the index in the excel sheet
# how to find corrsponding gene (example): 
# "X2461" corresponds to gene on the 2461+1 = 2462th row in the excel sheet 
# or the 2461th row in the "data" dataframe that was created at the begenning of this script

# test gene functionality
g1 = c("X13064","X13434",	"X17788",	"X17861",	"X18982",	"X19293",	"X20783",	"X21118",	"X22933",	"X32100",	"X32818",	"X33184",	"X34811",	"X40897")
g2 = c("X2461",	"X3964",	"X13434",	"X17788",	"X17861",	"X20783",	"X21118",	"X22933",	"X33184",	"X34042",	"X34811")
g = union(g1,g2)
gene = df[,g]
gene$y = c(rep(0,6),rep(1,6))
gene$y = factor(gene$y)
cormat = cor(gene[,1:17])
fit = glm(y~., data = gene, family = "binomial") #X13064+X32818+X33184+X34811+X40897+X2461+X3964+X34042
pred = predict(fit, type = "response")

x = data.matrix(gene[,c(1:17)])
y = as.matrix(gene[18])

# check contigency table 
summary(fit) #- NA coefficient:check colinearity
step(fit) #X22933 


# glmnet with cross validation on lambda
library(glmnet)
glm = glmnet(x,y,alpha = 1, family = "binomial")
cvfit = cv.glmnet(x, y, alpha = 1, family = "binomial")
plot(cvfit)
coef(cvfit, s = "lambda.min")
pred = predict(cvfit, newx = x, s = "lambda.min", type = "response")

"""
X13064        0.2709328
X32818        0.2301437
X33184        0.7575431
X34811        2.7472402
X40897        1.7707731
X2461         1.1973496
X3964         1.8729009
X34042        2.3772605
"""
# simulate data and simulation method/algorithm
# put significant gene name to folder

#  wilcoxon tets for group 3 vs group 2, group 4 vs group 2, group 3 vs group 4
# step one: wlcoxon test
# maybe the lack of certain gene will result in excessive drinkin????
# 2 vs 3
pval = c()
index = c()
meand = c()
meann = c()
for (i in 1:(length(df)-1)){
  a = wilcox.test(df[df$Y %in% c(3),i],df[df$Y %in% c(2),i],alternative = "greater")
  pval = c(pval,a$p.value)
  index = c(index,i)
  meand = c(meand,mean(df[df$Y %in% c(3),i]))
  meann = c(meann,mean(df[df$Y %in% c(2),i]))
}
result = data.frame(cbind(index,pval,meand, meann))
re32 = result[result$pval<0.05,] 
# 2 vs 4
pval = c()
index = c()
meand = c()
meann = c()
for (i in 1:(length(df)-1)){
  a = wilcox.test(df[df$Y %in% c(4),i],df[df$Y %in% c(2),i],alternative = "greater")
  pval = c(pval,a$p.value)
  index = c(index,i)
  meand = c(meand,mean(df[df$Y %in% c(4),i]))
  meann = c(meann,mean(df[df$Y %in% c(2),i]))
}
result = data.frame(cbind(index,pval,meand, meann))
re42 = result[result$pval<0.05,] 
# 3 vs 4
pval = c()
index = c()
meand = c()
meann = c()
for (i in 1:(length(df)-1)){
  a = wilcox.test(df[df$Y %in% c(3),i],df[df$Y %in% c(4),i])
  pval = c(pval,a$p.value)
  index = c(index,i)
  meand = c(meand,mean(df[df$Y %in% c(3),i]))
  meann = c(meann,mean(df[df$Y %in% c(4),i]))
}
result = data.frame(cbind(index,pval,meand, meann))
re34 = result[result$pval<0.05,]  #result is null

names(df[re42$index])


