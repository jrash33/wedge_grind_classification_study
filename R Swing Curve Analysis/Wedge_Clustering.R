rm(list=ls())

#functions
LoadLibraries=function(){
  library(caret)
  library(AppliedPredictiveModeling)
  library(corrplot)
  library(rpart) 
  library(C50)
  library(sfsmisc)
  library(reshape2)
  library(stats)
  library(cluster)
  library(factoextra)
  library(timeDate)
  print("The Libraries have been loaded!")
}

LoadLibraries()


#--------------------- Load Data and configure for modeling -----------------------------------
#open csv file
data=read.csv("raw combined.csv")
#remove any missing data (NaN)
data = na.omit(data)

#--------------------------------------- Clean Repeats ----------------------------------------
duplicates=which(duplicated(data[,1]))
result=diff(duplicates)
splits=split(duplicates,cumsum(c(1,diff(duplicates)!=1)))

#delete players column
data=data[,-1]
col=ncol(data)

# #normalize data
# data=scale(data)
# data=as.data.frame(data)

#split up swings for each player into lists (splits data.frame)
j=length(splits)
for (i in 1:j) {
  x=splits[[i]][1]-1
  splits[[i]]=append(splits[[i]],x,0)
  
  #now fill each list with the actual data
  splits[[i]]=data[splits[[i]],]
  
  i=i+1
}

VARi=matrix(nrow=j,ncol=col)
AVGp=matrix(nrow=j,ncol=col)
STDV=matrix(nrow=j,ncol=col)
DataP=matrix(nrow=j,ncol=col)

#analyze any significant variability in the swing between tests
for (i in 1:j){
  u=as.data.frame(splits[[i]])
  u=data.matrix(u)
  for (k in 1:col){
    AVGp[i,k]=mean(u[,k])
    VARi[i,k]=var(u[,k])
    STDV[i,k]=sqrt(var(u[,k]))
    
    if (k==1){
    DataP[i,]=u[k,]
    }
    
    k=k+1
  }
  i=i+1
}

#Variance between test sessions for club var
boxplot.matrix(STDV[,7:col],use.cols=TRUE,ylab="STDEV", xlab="MoCap Measurement",main="STDEV for 9 Players in MoCap",
               names=c("Speed","Face","Loft","Path","Attack",
                       "EffLoft","EffFace","ImpactX","ImpactY","GripLean","GripLie"))

#Variance between test sessions for ball var
boxplot.matrix(STDV[,1:6],use.cols=TRUE,ylab="STDEV", xlab="TMAN Measurement",main="STDEV for 9 Players in MoCap",
               names=c("Ball Speed","L.Angle","L.Direction","Height",
                       "Carry","Side"))


#Take 2 data sets: 1 with the mean of the players, and 1 with just one of their swings
AVGdata=AVGp
Playerdata=DataP

#convert our avg data back to data frame
colnames(AVGdata)=names(data)
colnames(Playerdata)=names(data)
colnames(AVGp)=names(data)
AVGp=as.data.frame(AVGp)


#ANOVA
# test=lm(Ball.Speed ~ Launch.Angle + Launch.Direction+Back.Spin+Side.Spin+Height+Carry+Side,data=AVGp)
# summary(test)
# anova(test)

#analyze variability from player to player
p.sd=apply(AVGp,2,sd)
p.sd=t(p.sd)

#Variance between test sessions for ball var
plot(c(1:col),p.sd)
barplot(p.sd)




#------------------------- Pre-Processing: Normalize/Center --------------------------------
#skewness
skew=rep(NA,col)

for (i in 1:col){
  skew[i]=skewness(AVGdata[,i])
  i=i+1
}

#skewness isn't actually that bad
hist(skew)
plot(skew)

skewValues=apply(AVGdata,2,skewness)
head(skewValues)

#Center and Scale Data
trans=preProcess(AVGdata,method=c("BoxCox","center","scale","pca"))
trans=preProcess(AVGdata,method=c("BoxCox","center","scale"))
trans

#apply the transformations
transformed=predict(trans,AVGdata)
transformed
#different from before bc they were transformed prior to PCA
head(transformed[,1:5])


#Easier way to center and scale data
AVGdata=scale(AVGdata,center=TRUE,scale=TRUE)
Playerdata=scale(Playerdata,center=TRUE,scale=TRUE)

skewfixed=rep(NA,28)
skewfixed=apply(AVGdata,2,skewness)
hist(skewfixed)

#------------------------ Clustering Methods ------------------------------

#******HIERARCHIAL CLUSTERING******

#agglomerative hierarchical clustering (HC): computes all pairwise dissimilarities between the elements
# in cluster 1 and the elements in cluster 2, and considers the largest value of these dissimilarities
#as the distance between the 2 clusters

# #clean outliers
# AVGdata=AVGdata[-c(2,6),]
# Playerdata=Playerdata[-c(2,6),]
# Playerdata=Playerdata[-c(32,39,62,66),]


#delete unnecessary columns
#AVGdata=AVGdata[,-c(1:11,20:27)]
#Playerdata=Playerdata[,-c(1:11,20:27)]

#clean outliers
# AVGdata=AVGdata[-c(2,30),]
# Playerdata=Playerdata[-c(2,30),]



#dissimilarity matrix
d=dist(AVGdata,method="euclidean")
d2=dist(Playerdata,method="euclidean")

#hierarchical clustering using complete linkage with AVG data
HC1=hclust(d,method="complete")
plot(HC1)
HC2=hclust(d,method="average")
plot(HC2)
HC3=hclust(d,method="single")
plot(HC3)
HC4=hclust(d,method="ward.D")
plot(HC4)
HC5=hclust(d,method="ward.D2")
plot(HC5)

#hierarchical clustering using complete linkage with ind player data
HC6=hclust(d2,method="complete")
plot(HC6)
HC7=hclust(d2,method="average")
plot(HC7)
HC8=hclust(d2,method="single")
plot(HC8)
HC9=hclust(d2,method="ward.D")
plot(HC9)
HC10=hclust(d2,method="ward.D2")
plot(HC10)


#compute with agnes method to get the agglomerative function 
#which measures the amount of clustering structure found (values close to 1 are best)
HC11=agnes(AVGdata,method="complete")
HC11$ac
pltree(HC11,cex=.6,hang=-1,main="Dendrogram of agnes")

HC12=agnes(AVGdata,method="average")
HC12$ac
pltree(HC12,cex=.6,hang=-1,main="Dendrogram of agnes")

HC13=agnes(AVGdata,method="single")
HC13$ac
pltree(HC13,cex=.6,hang=-1,main="Dendrogram of agnes")

#best method based on agglomerative coefficient
HC14=agnes(AVGdata,method="ward")
HC14$ac
pltree(HC14,cex=.6,hang=-1,main="Dendrogram of agnes")


#Divisive Hierarchical Clustering
HC15=diana(AVGdata)
#divisive coefficient; amount of clustering structure found
HC15$dc
pltree(HC15,cex=.6,hang=-1,main="Dendrogram of Diana")


#identify clusters/sub-groups with cutree

#cut tree into 3 groups
sub_grp1=cutree(HC1,k=3)
sub_grp2=cutree(HC2,k=3)
sub_grp3=cutree(HC3,k=3)
sub_grp4=cutree(HC4,k=3)
sub_grp5=cutree(HC5,k=3)
sub_grp6=cutree(HC6,k=3)
sub_grp7=cutree(HC7,k=3)
sub_grp8=cutree(HC8,k=3)
sub_grp9=cutree(HC9,k=3)
sub_grp10=cutree(HC10,k=3)


#number of members in each cluster
table(sub_grp1)
table(sub_grp2)
table(sub_grp3)
#*
table(sub_grp4)
#*
table(sub_grp5)

table(sub_grp6)
table(sub_grp7)
table(sub_grp8)
#*
table(sub_grp9)
#*
table(sub_grp10)


#******Moving forward  we will analyze Sub Groups 1**,4,5,9 and 10 with the ward method******


#way to see all sub-groups!
counts=sapply(2:6,function(ncl)table(cutree(HC1,ncl)))
names(counts)=2:6
counts

counts=sapply(2:6,function(ncl)table(cutree(HC4,ncl)))
names(counts)=2:6
counts

counts2=sapply(2:6,function(ncl)table(cutree(HC5,ncl)))
names(counts2)=2:6
counts2

counts3=sapply(2:6,function(ncl)table(cutree(HC9,ncl)))
names(counts3)=2:6
counts3

counts4=sapply(2:6,function(ncl)table(cutree(HC10,ncl)))
names(counts4)=2:6
counts4

#we can add a cluster column to our original data
#AVGdata %>% mutate(cluster=sub_grp1) %>% head

#observe clusters
plot(HC4)
rect.hclust(HC4,k=3,border=2:5)
fviz_cluster(list(data=AVGdata,cluster=sub_grp4))
clusplot(AVGdata,sub_grp4,color=TRUE,shade=TRUE)

plot(HC5)
rect.hclust(HC5,k=3,border=2:5)
fviz_cluster(list(data=AVGdata,cluster=sub_grp5))
clusplot(AVGdata,sub_grp5,color=TRUE,shade=TRUE)

plot(HC9)
rect.hclust(HC9,k=3,border=2:5)
fviz_cluster(list(data=Playerdata,cluster=sub_grp9))
clusplot(AVGdata,sub_grp9,color=TRUE,shade=TRUE)

plot(HC10)
rect.hclust(HC10,k=3,border=2:5)
fviz_cluster(list(data=Playerdata,cluster=sub_grp10))
clusplot(AVGdata,sub_grp10,color=TRUE,shade=TRUE)


# #use cutree with agnes and diana
# HC_a=agnes(AVGdata,method="ward")
# cutree(as.hclust(HC_a),k=3)
# 
# #cut diana() tree into 3 groups
# HC_d=diana(AVGdata)
# cutree(as.hclust(HC_d),k=4)


#we can also compare 2 dendrograms: compare hierarchical clustering with complete linkage vs ward's method
#tanglegram plots 2 dendrograms side by side with their labels connected by lines

library(dendextend)
dend1=as.dendrogram(HC4)
dend2=as.dendrogram(HC5)
dend3=as.dendrogram(HC9)
dend4=as.dendrogram(HC10)

#tanglegram plot as visual
tanglegram(dend1,dend2)
tanglegram(dend3,dend4)
tanglegram(dend1,dend3)
tanglegram(dend2,dend4)


#we can measure the quality of the alignment of these 2 trees (1 is no match/0 is perfect match)
dend_list=dendlist(dend3,dend4)
tanglegram(dend3,dend4,
           highlight_distinct_edges=FALSE, #turn off dashed lines
           common_subtrees_color_lines = FALSE, #turn off line colors
           common_subtrees_color_branches = TRUE, #color common branches
           main=paste("entanglement=",round(entanglement(dend_list),2)))


#elbow method
fviz_nbclust(AVGdata,FUN=hcut,method="wss")
fviz_nbclust(Playerdata,FUN=hcut,method="wss")

#avg silhouette method
fviz_nbclust(AVGdata,FUN=hcut,method="silhouette")
fviz_nbclust(Playerdata,FUN=hcut,method="silhouette")

#gap statistic method
gap_stat = clusGap(AVGdata, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)



#look at each cluster group
AVGdata[sub_grp4==1,]


#sapply(unique(sub_grp5),function(g)AVGdata[sub_grp5==g,])

#cluster characterization

#look at the medians of each parameter for each cluster
Cluster1=aggregate(AVGdata,list(sub_grp4),mean)
Cluster2=aggregate(AVGdata,list(sub_grp5),mean)
Cluster3=aggregate(Playerdata,list(sub_grp9),mean)
Cluster4=aggregate(Playerdata,list(sub_grp10),mean)
Cluster1
Cluster2
Cluster3
Cluster4

#plot these clusters
Cluster1=as.matrix(Cluster1)
Cluster1=Cluster1[,-1]

barplot(Cluster1, col=colors()[c(125,258,365)] , border="white", font.axis=2, 
        beside=T, legend=rownames(Cluster1), xlab="group", font.lab=2)

  
  


#add in cluster and freq column to above matrix
#a3 = aggregate(AVGdata,list(sub_grp5),median)
#data.frame(Cluster=a3[,1],Freq=as.vector(table(sub_grp5)),a3[,-1])


#let's see how our decisive tree compares to K-means clustering/mediods partioning methods
library(cluster)
PAMa=pam(AVGdata,3)
PAMb=pam(Playerdata,3)
names(PAMa)
names(PAMb)

#compare the results to hclust 
table(sub_grp4,PAMa$clustering)
table(sub_grp5,PAMa$clustering)
table(sub_grp9,PAMb$clustering)
table(sub_grp10,PAMb$clustering)
#solutions seem to agree except for a couple 

#feature with PAM known as the silhouette plot

#.71-1 strong structure has been found
#.5-.7 a reasonable structure has been found
#.26-.5 The structure is weak and could be artificial
# No substantial structure has been found
plot(PAMa)
plot(PAMb)

#silhouette plot for heirarchical cluster analysis
plot(silhouette(cutree(HC4,4),d))
plot(silhouette(cutree(HC5,4),d))
plot(silhouette(cutree(HC9,4),d2))
plot(silhouette(cutree(HC10,4),d2))










#--------------------- Classification Trees and Rule-Based Models -----------------------------


#Classification and Regression Tree
AVGdata=as.data.frame(AVGdata)
CARTmodel=rpart(Side~.,data=AVGdata)
plot(CARTmodel)



trctrl <- trainControl(method = "repeatedcv", number = 10,repeats = 10)

set.seed(1)
rpartGrouped=train(Side~.,
                   data=AVGdata,
                   method="rpart",
                   tuneLength=30,
                   trControl=trctrl)
summary(rpartGrouped)
plot(rpartGrouped)


set.seed(1)
library(ipred)
bagged=bagging(AVGdata$Side~.,data=AVGdata)

rpartGrouped=train(Side~.,
                   data=AVGdata,
                   method="bag",
                   tuneLength=30,
                   trControl=trctrl)
summary(rpartGrouped)
plot(rpartGrouped)


set.seed(1)
library(randomForest)
RanForest=randomForest(AVGdata$Side~.,data=AVGdata)









##Tree based clustering






































