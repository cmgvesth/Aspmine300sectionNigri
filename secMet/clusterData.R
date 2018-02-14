#!/usr/bin/env Rscript
library(optparse)

library(igraph)
library(parallel)
library(ggtree)

getwd()
setwd("/home/seth/aspSM/secMet")
source('hConnectAspMine.R')



log_con <- file("run.log")

args<-commandArgs(TRUE)



print("Working on file")
print(args[1])
filename = args[1]

print("Changing directory to")
print(args[2])
setwd(args[2])


biblast = args[3]

print("Biblast table:")
print(biblast)

aspDf <- read.table(filename, header = TRUE,  stringsAsFactors = FALSE, sep = '\t', quote = "\n")
# head(aspDf)
query <- sprintf("
        SELECT *, ROUND(
                 ( COALESCE( ( pident_tailoring /(clust_size-q_max_bb)  )*0.35,0)+
                 COALESCE( (pident_bb/q_max_bb)*0.65, 0)
                 ),2) AS pident_score FROM (

                 SELECT bidir.q_org, bidir.q_clust_id, bidir.h_org, bidir.h_clust_id,
                 SUM(CASE WHEN sm_short != 'none' then pident else 0 end) AS pident_bb,
                 smurf.clust_size,
                 COALESCE(SUM(sm_short != 'none'),0) AS count_bb,
                 tqmax.q_max_bb,
                 SUM(CASE WHEN sm_short = 'none' then pident else 0 end) AS pident_tailoring,
                 COALESCE(SUM(sm_short = 'none'),0) AS count_tailoring

                 FROM (
                 SELECT bia.* FROM (
                 SELECT * FROM %s) bia
                 JOIN (
                 SELECT * FROM %s ) bib
                 ON bia.h_clust_id = bib.q_clust_id
                 WHERE bia.q_clust_id = bib.h_clust_id
                 GROUP BY bia.q_clust_id, bia.q_protein_id, bia.h_clust_id, bia.h_protein_id
                 ) AS bidir

                 JOIN smurf
                 ON bidir.q_org = smurf.org_id AND bidir.q_protein_id = smurf.sm_protein_id AND bidir.q_clust_id != bidir.h_clust_id

                 JOIN (SELECT CONCAT(org_id, '_' , clust_backbone,'_', clust_size) AS q_clust_id, SUM(sm_short != 'none') AS q_max_bb FROM smurf
                 GROUP BY q_clust_id) tqmax
                 ON bidir.q_clust_id = tqmax.q_clust_id

                 GROUP BY q_clust_id, h_clust_id ) ta;", biblast, biblast)

clusterBlastAll <- aspDbFetch(query)


clusterBlast <- clusterBlastAll

str(clusterBlast)
clusterBlast <- clusterBlast[clusterBlast$pident_score != 0, c('q_clust_id', 'h_clust_id', 'pident_score')]

names(clusterBlast)[3] <- 'weight'


head(aspDf)
head(subset(aspDf, subset = org_id==3))

head(clusterBlast)
print("Testing aspDf")
head(aspDf)
cat("Missing gene clusters in clusterblast query", file = log_con, append = TRUE)
cat(paste(setdiff(clusterBlast$q_clust_id, aspDf$cluster_id), collapse = ",") , file = log_con, append = TRUE)
setdiff(clusterBlast$q_clust_id, aspDf$cluster_id)
cat("Missing gene clusters in downloaded data frame", file = log_con, append = TRUE)
cat(paste(setdiff(aspDf$cluster_id, clusterBlast$q_clust_id), collapse = ",") , file = log_con, append = TRUE)
missingOnes <- setdiff(aspDf$cluster_id, clusterBlast$q_clust_id)
print(missingOnes)

# cat("Missing gene clusters in clusterblast query", file = log_con, append = TRUE)
cat(paste(setdiff(clusterBlast$h_clust_id, aspDf$cluster_id), collapse = ","), file = log_con, append = TRUE)

clusterBlast <- clusterBlast[(clusterBlast$q_clust_id %in% aspDf$cluster_id & clusterBlast$h_clust_id %in% aspDf$cluster_id),]
print("Checking clusterBlast")

str(clusterBlast)

clusterDups <- clusterBlast[duplicated(clusterBlast[,c('q_clust_id', 'h_clust_id')]),c('q_clust_id', 'h_clust_id')]
cat(paste(clusterDups, collapse = ","), file = log_con, append = TRUE)

sum(clusterBlast$q_clust_id == clusterBlast$h_clust_id)
clusterDups

write.csv( clusterBlast, file = "clusterBlast.csv", row.names = FALSE )
# }else{
#   clusterBlast <- read.csv("clusterBlast.csv")
# }
# setwd(file.path(mainDir, ))

head(clusterBlast)


# Creating a whole network, subsetting later with dataframe which is already created


set.seed(123)

cat("Info about clusterBLast", file = log_con, append = TRUE)

# clusterBlastSub <- clusterBlast[clusterBlast$weight >20,]
cat("Minimum clusterBlast weight", file = log_con, append = TRUE)
cat(min(clusterBlast$weight), file = log_con, append = TRUE)

cat("ClusterBlast weight over 100", file = log_con, append = TRUE)
# cat(clusterBlastAll[clusterBlastAll$pident_score > 100,], file = log_con, append = TRUE)


print("Creating Network")

# present_clusters <- unique(aspDf$cluster_id[aspDf$cluster_id %in% clusterBlast$q_clust_id | aspDf$cluster_id %in% clusterBlast$h_clust_id])

g.raw <- graph_from_data_frame(d=clusterBlast, vertices= present_clusters, directed=FALSE) # Decided on undirected graph because I it's a bidirectional relationship for all anyway
g <- igraph::simplify(g.raw, edge.attr.comb = list(weight = "mean")) #remove.multiple = TRUE)
walks <- cluster_walktrap(as.undirected(g), steps = 1)
walks
# sizes(walks)
# gw <- induced.subgraph(g, vids = walks[[16]])
# walkMore <- cluster_walktrap(gw)

# fg <- cluster_fast_greedy(g)

if(max(as.numeric(sizes(walks))) > length(unique(aspDf$org_id)) ){
    print("Network families are too large, needs more clustering")
    print("There are too many gene clusters per family, since they are exceeding the number of orgs (we assume that duplications will form another family, hence another round of clustering)")

    # print("Complete this section")
   superNetworks <- as.numeric(which(sizes(walks) > length(unique(aspDf$org_id))))

     firstNetworks <- lapply(names(groups(walks))[!(names(groups(walks)) %in% superNetworks)],function(x){
         sub <- induced.subgraph(g, vids = walks[[x]])
         })

     secondNetwork <- lapply(names(groups(walks))[names(groups(walks)) %in% superNetworks],function(x){
         sub <- induced.subgraph(g, vids = walks[[x]])
         walkMore <- cluster_walktrap(sub, steps = 1)
         lapply(names(groups(walkMore)), function(x){
             induced.subgraph(g, vids = walkMore[[x]])
             } )
         })

     testg <- induced_subgraph(g, vids = walks[[superNetworks[[1]]]])
     # plot(testg)
     # plot(secondNetwork[[1]][[4]])
     # secondNetwork

     secondNetwork <- unlist(secondNetwork, recursive = FALSE)
     cNetwork <- c(firstNetworks, secondNetwork) # complete network
     sizes <- sapply(cNetwork, function(x){length(V(x))})
          counter = 0
 clusterFamDf <- lapply(cNetwork, function(x){
     counter <<- counter + 1
     clusters <- V(x)$name
     data.frame(cluster_id = clusters, clusterFam = rep(counter, length(clusters)),
     stringsAsFactors = FALSE)
     })

 clusterFamDf

 clusterFamDf <- do.call(rbind, clusterFamDf)
 clusterFamDf$cluster_id <- as.character(clusterFamDf$cluster_id)

 # Filling up the ones with no hits
famLimit <-  max(clusterFamDf$clusterFam)
names(clusterFamDf)
missingClusters <- setdiff(unique(aspDf$cluster_id), clusterFamDf$cluster_id)
missingDf <- data.frame(cluster_id = missingClusters, clusterFam = seq(from = famLimit+1, to = famLimit+length(missingClusters), by = 1))

clusterFamDf <- rbind(clusterFamDf, missingDf)

viFams <- merge(aspDf, clusterFamDf, by = c("cluster_id"), all.x = TRUE)

viFams
# write.table(viFams, "smurfOrgIprMibigFams.tsv", row.names = TRUE)
# WRITE IT TO DISK
 # ggplot(as.data.frame(table(sizes))) + geom_bar(aes(sizes, Freq), stat = 'identity') + geom_text(aes(x = sizes, y = Freq, label = Freq) ,position=position_dodge(width=0.9), vjust=-0.25, size =  8 )+ labs( x = "Number of SMGC in family", y = "Family counts" , size = 2) + theme(axis.text = element_text(size = 14) , axis.title=element_text(size=12,face="bold"))

}else{

    print("No further clustering")
    coms <- communities(walks)
    tmpFam <- lapply(names(coms), function(x){
        data.frame(cluster_id = coms[[x]], clusterFam = (x), stringsAsFactors = FALSE)
        })

    famDf <- do.call(rbind, tmpFam)

    str(famDf)
    head(aspDf)
    viFams <- merge(aspDf, famDf, by = c("cluster_id"), all.x = TRUE)
    str(viFams)

    # write.table(viFams, "smurfOrgIprMibigFams.tsv", row.names = TRUE)
}

write.table(viFams, paste0(gsub(".tsv","", filename),"_c.tsv"), row.names = FALSE, sep = "\t")
