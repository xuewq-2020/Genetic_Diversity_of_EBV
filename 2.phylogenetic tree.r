
# Set the working directory
setwd("/Users/data/EBV_genoms/treefile")

# Load the libraries
library("treeio")
library("ggtree")
library("ggplot2")
library("phangorn")

# read the treefile (from IQtree)
tree_ebv_protein <- read.newick("iqtree_target_protein.treefile",node.label="support")

# visualize the tree
p<-ggtree(tree_ebv_protein,color="black",size=0.2)+geom_tiplab(size=0.3)

# read annotation for the heatmap
# define labels for individual sequence:disease type/health condition and geographic origin
ann1<-read.csv("9.seq list id_BALF2.csv",header=TRUE,row.names = 1)

# add symbols colored by footstrap values
p1 <- p+geom_nodepoint(aes(colour=cut(support,c(0,40,70,100),include.lowest=TRUE)),shape=16,size=0.5)+theme(legend.position="right")
p1 <-p1 + scale_color_manual(values=c("grey","orange","red"),labels=c("<40","40<= & <70",">=70"),name='Bootstrap Percentage(BP)')

# visualize the tree with diseases/health conditions
p2 <- gheatmap(p1, ann1[,"Diseases",drop=FALSE],offset=0.0002,width=.05,colnames=FALSE)

# visualize the tree with geographic origion of the isolates
p3 <- gheatmap(p2,ann1[,"Geographic_origin",drop=FALSE],offset=0.0011,width=.05,colnames=FALSE)

# change colors for diseases/health conditions 
col <- c("seashell","red","yellow","deepskyblue","lightskyblue","cyan","palegreen","orangered","lightblue","ivory","purple","tomato","lightgreen","skyblue","darkred","orchid","lightyellow","beige","steelblue","indianred","black","white","ivory2","ivory3","slateblue","red","orange","grey")
# change colors for geographic origion of isolates
names(col)=c("CAEBV","China","Africa","America","Australia","BL","Europe","GC","HL","IM","LC","Lymphoepithelioma","NHL","NKTL","NPC","East Asia","PTLD","SOT","Lymphoma","Southeast Asia","TS","Healthy_LCLs_Saliva","RH","HIV","Oceania",c('(70,100]', '(40,70]', '[0,40]'))

# visualize the tree with both diseases/health condition and geographic origin info
pp <- p3 + scale_fill_manual(values=col)
p1x <- p1 + scale_color_manual(values=c("grey","orange","red"),labels=c("<40","40<= & <70",">=70"),name='Bootstrap Percentage(BP)')
p2x <- gheatmap(p, ann1[,"Diseases",drop=FALSE],offset=5, width=.1)+ scale_fill_manual(values=col,name='Disease')
p3x <- gheatmap(p, ann1[,"Geographic_origin",drop=FALSE],offset=5, width=.1)+ scale_fill_manual(values=col,name='Geographical origin')

# get the legends for bootstrap values, diseases/health condition and geographic origin
require(cowplot)
leg1 <- get_legend(p1x)
leg2 <- get_legend(p2x)
leg3 <- get_legend(p3x)

# combined the tree with the legends
pp <- pp + theme(legend.position="none")
right_legend <- plot_grid(leg1,leg2,leg3,ncol=1,rel_heights=c(1,3,2))
plot_grid(pp,right_legend,ncol=2,rel_widths=c(3,1))



