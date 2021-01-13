# set the working directory
setwd("/Users/data/EBV_genoms/")

library(adegenet)
library(ade4)

# Read AA substition data
aa_subs<-read.csv("target_protein_aa_substitution.csv",header=T,row.names = 1)

# Convert dataframe of AA substitution into a genind object
aa_subs_genind<-df2genind(aa_subs,ploidy=1,NA.char="X") 

# load the treeWAS library
library(treeWAS)

# exclude AA substitution minor frequency < 0.05
aa_subsnew<-aa_subs
aa_subname<-colnames(aa_subs)

dellist<-c()
for (i in (1:length(aa_subname))) 
{
if (any(table(aa_subsnew[,i])-(length(rownames(aa_subsnew))*0.05)<0)) 
{dellist<-c(dellist,i)}
}

for (i in (1:length(aa_subname))) 
{
if (any(table(aa_subsnew[,i])-(length(na.omit(aa_subsnew[,i])))==0)) 
{dellist<-c(dellist,i)}
}

aa_sub_maf5<-aa_subsnew[,-dellist]

write.csv(aa_sub_maf5,"aa_sub_maf5_exclude.csv")

# load data
aa_sub_maf5 <- read.csv("aa_sub_maf5_exclude.csv",row.names=1)
phen1<-read.csv("phenotypes.csv", row.names=1) #phenotypes can be healthy condition (NPC/health) or geographic areas (China/East Asia)
tree <-load("iqtree_target_protein.rda") #tree file from IQtree

is.null(tree$tip.label)
is.null(rownames(snps))
is.null(names(phen1))

## Cross‐check labels with each other:
all(tree$tip.label %in% rownames(aa_sub_maf5))
all(rownames(aa_sub_maf5) %in% tree$tip.label)
all(tree$tip.label %in% names(phen1))
all(rownames(phen1) %in% tree$tip.label)
all(rownames(phen1) %in% rownames(aa_sub_maf5))

# association analysis
out <- treeWAS(snps = aa_sub_maf5,
phen = phen1,
tree = tree,
n.subs = NULL,n.snps.sim = ncol(snps)*10,
test = c("terminal", "simultaneous", "subsequent"),
snps.reconstruction = "ML",
snps.sim.reconstruction = "parsimony",
phen.reconstruction = "parsimony",
na.rm = TRUE,
p.value = 0.01,
p.value.correct = "bonf",
p.value.by = "count",
dist.dna.model = "JC69",
plot.tree = FALSE,
plot.manhattan = TRUE,
plot.null.dist = TRUE,
plot.dist = FALSE,
snps.assoc = NULL,
filename.plot = NULL,
seed = 1)

print(out, sort.by.p=FALSE)

