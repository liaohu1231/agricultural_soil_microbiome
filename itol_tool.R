library(itol.toolkit)
library(data.table)
library(ape)
tree <- read.tree(file = "norQ_foldtree_struct_tree.nwk")
data=read.csv("bacteria_norQ_gtdb_itol.tsv",sep='\t',header = F)

#group=data[,c(5,4)]
hub <- create_hub(tree = tree)
unit_21 <- create_unit(data = data,
                       key = "E021_color_strip_1",
                       color = "wesanderson",
                       type = "TREE_COLORS", 
                       subtype = "range",
                       tree = tree)

hub <- hub + unit_21
write_hub(hub,getwd())


data2=read.table("../prokaryotic/Table_S_SOB_SRM_classification.csv",sep='\t',header = 1)

dsrB=read.table("DsrB_SOB_thisstudy2.txt",sep='\t')
dsrB=dsrB[!duplicated(dsrB$V3),]
table(dsrB$V2)

gtdbtk=read.table("../prokaryotic/gtdbtk.bac120.summary.tsv",row.names = 1,sep='\t',header = 1)
tmp=strsplit(gtdbtk$classification,split = ";")
tmp=do.call(rbind,tmp)
gtdbtk$class=tmp[,4]
dsrB$family=gtdbtk[as.character(dsrB$V3),20]

group=dsrB[,c(1,4)]

hub <- create_hub(tree = tree)

unit_21 <- create_unit(data = group,
                       key = "E021_color_strip_family_2",
                       type = "DATASET_COLORSTRIP",
                       tree = tree)
write_unit(unit_21)

data$study=rep(c("other study","this study"),c(406,28))
group=data[,c(1,22)]
unit_8 <- create_unit(data = group, 
                      key = "E008_tree_colors_2", 
                      type = "TREE_COLORS", 
                      subtype = "clade", 
                      size_factor = 1, 
                      tree = tree)
write_unit(unit_8)

## write template file
write_hub(hub,getwd())

