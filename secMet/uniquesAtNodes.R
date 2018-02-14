library(ape)
library(ggtree)



args<-commandArgs(TRUE)
print(args[1])
setwd(args[2])
filename = args[1]
treeFile = args[3]

# Can only be executed after cluster families have been generated

# getwd()
# setwd("/home/seth/asptoolbox/secMet/flavi_orgs")
aspDf <- read.table(filename, stringsAsFactors = FALSE, header = TRUE, sep = "\t")
## READ TREE IN HERE

# dummyTree <- rtree(length(tiplabs), tip.label = as.character(tiplabs))
dummyTree <- read.tree(treeFile) # Tree needs new organism ids
# plot(dummyTree)

p <- ggtree(dummyTree, branch.length = "none")
# p + geom_tiplab()

# viFams <- read.table(paste0(aspSetName,filename),stringsAsFactors = FALSE,header = TRUE)

getUniquesAtNodes <- function(p, tree, clusters){
  "function changed to return df, changed function to take name"

  # tree <- dummyTree
  # clusters <- viFams

  nodeContent <- lapply(p$data$node, function(x){
    get.offspring.tip(tree, p$data$node[x])
  })

  names(nodeContent) <- p$data$node

  orgsByFam <- split(clusters$name, clusters$clusterFam)

  # Comparing orgs at nodes and in families

  ndcUnique <- lapply(nodeContent, function(x){
    pattern <- sapply(orgsByFam, function(y){ all(x  %in% y) & all(y %in% x)})
    as.integer(names(orgsByFam)[pattern])
  })


  tmpUniques <- lapply(names(ndcUnique), function(x){
    tmpFams <- as.character(ndcUnique[[x]])
    if(length(tmpFams) == 0){ tmpFams <- "none"}
    data.frame(clusterFam = tmpFams, nodePosition = as.character(x) )
  })

  unNodeDf <- do.call(rbind, tmpUniques)

  return(unNodeDf)
}


head(aspDf)
uFamsAtNodes <- getUniquesAtNodes(p, dummyTree, aspDf)
head(uFamsAtNodes)

viFamsNodes <- merge(aspDf, uFamsAtNodes, by = "clusterFam", all.x = TRUE)

write.table(viFamsNodes, paste0(gsub(".tsv","", filename),"_n.tsv"), row.names = TRUE)

uNodes <- lapply(split(uFamsAtNodes$clusterFam, uFamsAtNodes$nodePosition),function(x){
    if(x == 'none'){res = 0
        }else{
            res <- length(x)
        }
        res
    })

uNodes <- data.frame(counts = do.call(rbind, uNodes), node = as.integer(names(uNodes)))

write.table(uNodes, sprintf("%s_nodePosition.tsv", gsub(".tsv","", filename)), row.names = FALSE)
# Remeber node can also be a tip.... anyway, shared families are now in uFamsAtNodes
# p$data$allgcCounts <- uNodes$counts
#
# p + geom_tiplab() + geom_text2(aes(label = allgcCounts, vjust = -.2, hjust = 2)) + xlim(0, 20)
