
# Set the working directory
setwd("/Users/data/EBV_genoms/")

# Load the library
library(CHNOSZ)

## Convert the protein sequences into dataframe of AA substitutions
aa = read.fasta("taget_protein_seq_aligned.fas", ret = "sequence", type="protein")
counts=read.fasta("taget_protein_seq_aligned.fas", ret = "count", type="protein")[,1]
names(aa)=counts

split=function(seq){
  str=strsplit(unlist(seq),split = "",fixed=T)
  return(str)
}
matrx=sapply(aa,split)
matrx = do.call(rbind,matrx)
matrx_conv=matrix(NA,nrow=nrow(matrx),ncol=ncol(matrx))

for(i in 1:ncol(matrx)){
  matrx_conv[,i] = matrx[,i]
}
rownames(matrx_conv)=names(aa)
colnames(matrx_conv)=1:ncol(matrx)
exclude=function(x){
  return(length(table(x)))
}
id=apply(matrx_conv,MARGIN=2,FUN=exclude)
sub_matrx=matrx_conv[,which(id>1)]

sub_matrx1 <- sub_matrx
sub_matrx1[sub_matrx1=="X"] <- NA
sub_matrx1[sub_matrx1=="-"] <- NA
write.csv(sub_matrx1,"target_protein_aa_substitution.csv",quote=F,row.names = T)

##############################################################
#################### PCA analysis  ###########################
##############################################################

# Load the library
library(adegenet)
library(ade4)

# Read AA substition data
aa_subs<-read.csv("target_protein_aa_substitution.csv",header=T,row.names = 1)

# Convert dataframe of AA substitution into a genind object
aa_subs_genind<-df2genind(aa_subs,ploidy=1,NA.char="X") 
aa_subs_genind 
X <- tab(aa_subs_genind, NA.method="mean")

# Make PCA
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE)

## Make basic plot
plot(pca1$li[,1],pca1$li[,2], pch=19)

# PCA visulization
# Load factorMineR and factoextra libraries
library("FactoMineR")
library("factoextra")

# Read the PCA results
eig.val <- get_eigenvalue(pca1)
# Extract the results for variables 
var <- get_pca_var(pca1)
head(var$coord)
head(var$cos2)
head(var$contrib)
eig.val

# Visualize variables
fviz_pca_var(pca1, col.var = "black")

# Read geographic information of each sequences
ann<-read.csv("9.seq list id_BALF2.csv",header=T, row.names=1)


# Plot for individuals
fviz_pca_ind(pca1, geom="point",repel=T)

# Assign colors to different geographic origin
col <- c("red","yellow","deepskyblue","lightskyblue","palegreen","orchid","indianred","slateblue")
names(col)=c("China","Africa","America","Australia","Europe","East Asia","Southeast Asia","Oceania")
ann$color[ann$Geographic_origin=="China"] <- "red"
ann$color[ann$Geographic_origin=="Africa"] <- "yellow"
ann$color[ann$Geographic_origin=="America"] <- "deepskyblue"
ann$color[ann$Geographic_origin=="Australia"] <- "lightskyblue"
ann$color[ann$Geographic_origin=="Europe"] <- "palegreen"
ann$color[ann$Geographic_origin=="East Asia"] <- "orchid"
ann$color[ann$Geographic_origin=="Southeast Asia"] <- "indianred"
ann$color[ann$Geographic_origin=="Oceania"] <- "slateblue"

# Color individuals by groups
sp <- fviz_pca_ind(pca1, col.ind=ann$Geographic_origin, palette=col,geom="point",mean.point=F, addEllipses = F,legend.title="Geographic_origin",pch=19)
sp

# Color individuals by groups, add concentration ellipses
sp1 <- fviz_pca_ind(pca1, col.ind=ann$Geographic_origin, palette=col,geom="point",mean.point=F, addEllipses = T,legend.title="Geographic_origin",pch=19)
sp1

# Color individuals by groups,jitter labels to avoid overplotting
sp2 <- fviz_pca_ind(pca1, col.ind=ann$Geographic_origin, palette=col,geom="point",mean.point=F, legend.title="Geographic_origin",pch=19) 
sp2 <-sp2 + geom_jitter(width = 0.15,height=0.15,color=ann$color)
sp2