
input_dir = paste0('../OneDrive-GCCORP/Community distance/MoBC_data_code/input_file')
output_dir = paste0('../OneDrive-GCCORP/Community distance/MoBC_data_code/output_file')

# Module detection and Background network

# library
library(igraph)
library(dplyr)
library(readr)
library(ggraph)
library(ggplot2)
library(plyr)

## brca ppi network and module detection
# Human ppi network
human.alias <- read.delim(paste0(input_dir, '/9606.protein.aliases.v12.0.txt'), header=TRUE)
human.info <- read.delim(paste0(input_dir, '/9606.protein.info.v12.0.txt'), header=TRUE)
colnames(human.info) <- c('protein_id', 'Symbol', 'protein_size', 'annotation')
human.info.final <- human.info[,c('protein_id', 'Symbol')]
colnames(human.alias) <- c('protein_id', 'alias', 'source')

human.ppi.network <- read.table(paste0(input_dir, '/9606.protein.links.v12.0.txt'), header=TRUE)
human.info.final <- human.info[,c('protein_id', 'Symbol')]
human.alias.final <- human.alias[,c('protein_id', 'alias')]
human.ppi.network.700 <-  human.ppi.network[human.ppi.network$combined_score > 700,]
colnames(human.ppi.network.700) <- c('protein_id', 'protein2', 'combined_score')
human.ppi.network2 <- merge(human.ppi.network.700, human.info.final, by = 'protein_id')
colnames(human.ppi.network2) <- c('protein1', 'protein_id', 'combined_score', 'gene1')
human.ppi.network3 <- merge(human.ppi.network2, human.info.final, by = 'protein_id')
colnames(human.ppi.network3) <- c('protein1', 'protein2', 'combined_score', 'gene1', 'gene2')
human.ppi.network.final <- human.ppi.network3[,c('gene1', 'gene2')]

#save(human.ppi.network.final, file= paste0(output_dir, '/human.ppi.network.final.RData'))
load(paste0(output_dir, '/human.ppi.network.final.RData')) # human.ppi.network.final)
# brca deg results
load(paste0(output_dir, '/breast.deseq.results.RData')) # res2.sd (brast deg result)


# brca background network
human.background = unique(c(human.ppi.network.final$gene1, human.ppi.network.final$gene2))
human.itst = intersect(human.background, res2.sd$Symbol)
human.ppi.network.final.itst = human.ppi.network.final[(human.ppi.network.final$gene1 %in% human.itst) & (human.ppi.network.final$gene2 %in% human.itst),] # human background network
#save(human.ppi.network.final.itst, file= paste0(output_dir, '/human.ppi.background.network.RData'))



# brca deg selection
load(paste0(output_dir, '/breast.deg.result.final.2.RData')) 
fcut = 2
pcut = 0.01
brca.deg <- res.sd.2[abs(res.sd.2$log2FoldChange) > fcut & res.sd.2$padj < pcut, ]
brca.deg %>% dim
write.csv(brca.deg, file = paste0(output_dir, '/brca.deg.csv'), row.names=FALSE)


breast.network <- human.ppi.network.final[(human.ppi.network.final$gene1 %in% brca.deg$Symbol) & (human.ppi.network.final$gene2 %in% brca.deg$Symbol),] # brca deg network

# brca deg network
breast.g = graph_from_data_frame(breast.network, directed= F)
#save(breast.g, file= paste0(output_dir, '/breast.g.RData'))

load(file = paste0(output_dir, '/breast.g.RData'))
# brca module detection
mscore<- c()
for (i in 1:10) {
    walktrap_result <- cluster_walktrap(
  graph= breast.g,
#   weights = E(g.all.rb)$weight,
  steps = i,
  merges = TRUE,
  modularity = TRUE,
  membership = TRUE
)
mscore <- c(mscore, modularity(walktrap_result))
}

step.breast = 4

walktrap_result.breast <- cluster_walktrap(
  graph= breast.g,
#   weights = E(g)$weight,
  steps = step.breast,
  merges = TRUE,
  modularity = TRUE,
  membership = TRUE
)

comm.breast <- lapply(1:length(walktrap_result.breast), function(x) {
    re = walktrap_result.breast[[x]]
    if (length(re)> 20) {
        re
    }
    })

comm.breast <- compact(comm.breast)
comm.breast = comm.breast[order(lengths(comm.breast), decreasing=T)]
names(comm.breast) = 1:length(comm.breast)
#save(comm.breast, file = paste0(output_dir, '/comm.breast.RData'))


## cac ppi network and module detection
# Mouse ppi network
mouse.ppi = read.table(paste0(input_dir, '/10090.protein.links.v12.0.txt'), header=TRUE)
mouse.info = read.delim(paste0(input_dir, '/10090.protein.info.v12.0.txt'), header=TRUE)

colnames(mouse.info) <- c('protein_id', 'Symbol', 'protein_size', 'annontation')
mouse.info.final = mouse.info[,c('protein_id', 'Symbol')]
mouse.ppi.700 <-  mouse.ppi[mouse.ppi$combined_score > 700,]
colnames(mouse.ppi.700) <- c('protein_id', 'protein2', 'combined_score')
mouse.ppi.2 <- merge(mouse.ppi.700, mouse.info.final, by = 'protein_id')
colnames(mouse.ppi.2) <- c('protein1', 'protein_id', 'combined_score', 'gene1')
mouse.ppi.3 <- merge(mouse.ppi.2, mouse.info.final, by = 'protein_id')
colnames(mouse.ppi.3) <- c('protein1', 'protein2', 'combined_score', 'gene1', 'gene2')
mouse.ppi.final <- mouse.ppi.3[,c('gene1', 'gene2')]

#save(mouse.ppi.final, file = paste0(output_dir, '/mouse.ppi.final.RData'))

# cac deg results
load(paste0(input_dir, '/for_valid_cachexia_degs.Rdata')) # cac deg result - eml.cac (DEG results)
cac.res = eml.cac$GSE157251_cancer_LateFemale
#save(cac.res, file = paste0(output_dir, '/cac.deg.result.RData'))

# cac background network
mouse.background = c(mouse.ppi.final$gene1, mouse.ppi.final$gene2) %>% unique
mouse.itst = intersect(cac.res$gene_name, mouse.background)
mouse.ppi.final.itst = mouse.ppi.final[(mouse.ppi.final$gene1 %in% mouse.itst) & (mouse.ppi.final$gene2 %in% mouse.itst),]

#save(mouse.ppi.final.itst, file = paste0(output_dir, '/mouse.ppi.background.network.RData'))

# cac deg selection

pcut.cac = 0.05
fcut.cac = 1
cac.deg = subset(cac.res, abs(log2FoldChange) > fcut.cac & padj < pcut.cac)
cac.deg %>% dim
write.csv(cac.deg, file = paste0(output_dir, '/cac.deg.csv'), row.names=FALSE)

# cac deg network
cac.network<- mouse.ppi.final[(mouse.ppi.final$gene1 %in% cac.deg$gene_name) & (mouse.ppi.final$gene2 %in% cac.deg$gene_name),]

cac.g = graph_from_data_frame(cac.network, directed= F)
#save(cac.g, file= paste0(output_dir, '/cac.g.RData'))

# cac module detection
mscore.cac<- c()
for (i in 1:10) {
  walktrap_result <- cluster_walktrap(
  graph= cac.g,
#   weights = E(g.all.rb)$weight,
  steps = i,
  merges = TRUE,
  modularity = TRUE,
  membership = TRUE
)
mscore.cac <- c(mscore.cac, modularity(walktrap_result))
}

step.cac = 3

walktrap_result.cac <- cluster_walktrap(
  graph= cac.g,
#   weights = E(g)$weight,
  steps = step.cac,
  merges = TRUE,
  modularity = TRUE,
  membership = TRUE
)

comm.cac <- lapply(1:length(walktrap_result.cac), function(x) {
    re = walktrap_result.cac[[x]]
    if (length(re)> 10) {
        re
    }
    })


comm.cac <- compact(comm.cac)
comm.cac = comm.cac[order(lengths(comm.cac), decreasing=T)]
names(comm.cac) = 1:length(comm.cac)
#save(comm.cac, file = paste0(output_dir, '/comm.cac.RData'))

