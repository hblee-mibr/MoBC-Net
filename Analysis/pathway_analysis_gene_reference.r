input_dir = paste0('/Users/yoomibaek/Library/CloudStorage/OneDrive-GCCORP/Community distance/MoBC_data_code/input_file')
output_dir = paste0('/Users/yoomibaek/Library/CloudStorage/OneDrive-GCCORP/Community distance/MoBC_data_code/output_file')

# pathway analysis and gene reference

# library
library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)
library(purrr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(tibble)
library(tidyr)
library(rtracklayer)

# gene reference
# brca gene reference - gencode

z <- import('gencode.v46.annotation.gtf')
z1 = subset(z, gene_type=='protein_coding')
z1 = z1@elementMetadata %>% as.data.frame
z1 = z1[,c('gene_id','transcript_id','gene_name')] %>% unique
z1 = z1[which(!is.na(z1$transcript_id)),]
z1.new <- gsub('\\..*','',z1$gene_id)


z.entrez = read.table(paste0(input_dir, '/gencode.v46.metadata.EntrezGene.txt'), header=TRUE)
colnames(z.entrez) = c('transcript_id', 'Entrez')
z1.merge = merge(z1, z.entrez, by=c('transcript_id'),all.x=T)
z1.merge$transcript_id = gsub('\\..*','',z1.merge$transcript_id)
z1.merge$gene_id = gsub('\\..*','',z1.merge$gene_id)
z1.merge.unique = z1.merge[,2:4] %>% unique

save(z1.merge.unique, file = paste0(output_dir, '/human.ginfo.final.genecode.RData'))

# cac gene reference
load(paste0(output_dir, '/cac.deg.result.RData')) # cac.res
mouse.ginfo.final = cac.res[,c('GeneId', 'EntrezID', 'gene_name')] %>% unique
save(mouse.ginfo.final, file = paste0(output_dir, '/mouse.ginfo.final.RData'))

# background gene
# human background gene
load(paste0(output_dir, '/human.ppi.network.final.RData')) # human.ppi.network.final
load(paste0(output_dir, '/breast.deg.result.final.2.RData')) # res2.sd

human.background = unique(c(human.ppi.network.final$gene1, human.ppi.network.final$gene2))
human.itst = intersect(human.background, res.sd.2$Symbol)
human.itst.info =subset(z1.merge.unique, gene_name %in% human.itst)
human.info.entrez.genecode = human.itst.info$Entrez %>% unique %>% as.character %>% na.omit

# mouse background gene
load(paste0(output_dir, '/mouse.ppi.final.RData'))
load(paste0(output_dir, '/cac.deg.result.RData')) # cac.res

mouse.background = c(mouse.ppi.final$gene1, mouse.ppi.final$gene2) %>% unique
mouse.itst = intersect(cac.res$gene_name, mouse.background)
mouse.itst.info = subset(cac.res, gene_name %in% mouse.itst)
mouse.info.entrez.deg = mouse.itst.info$EntrezID %>% unique %>% na.omit

save(human.info.entrez.genecode, mouse.info.entrez.deg, file= paste0(output_dir, '/human.mouse.entrez.background.RData'))

# community gene mapping - Entrez
load(paste0(output_dir, '/comm.breast.RData'))
load(paste0(output_dir, '/comm.cac.RData'))


comm.cac.entrez<- lapply(comm.cac, function(x) {
  gene.list <- subset(cac.res, gene_name %in% x)$EntrezID
  gene.list = na.omit(gene.list) %>% unique %>% as.character
})

comm.breast.entrez<- lapply(comm.breast, function(x) {
  gene.list <- subset(z1.merge.unique, gene_name %in% x)$Entrez
  gene.list = na.omit(gene.list) %>% unique %>% as.character
})
save(comm.cac.entrez, comm.breast.entrez, file= paste0(output_dir, '/comm.breast.cac.RData'))


# pathway analysis - KEGG, GO
pathway.breast.go.narm.raw.final.all <- lapply(comm.breast.entrez, function(x) {
  path.res <- enrichGO(gene=x, 
                              OrgDb = org.Hs.eg.db,
                              ont = "BP",
                              universe = unique(human.info.entrez.genecode),
                              pvalueCutoff=Inf, pAdjustMethod = 'fdr', qvalueCutoff=Inf)

  path.res <- data.frame(path.res) %>% mutate(Description = sapply(str_split(data.frame(path.res)$Description, ' -'), function(x) x[1])) 
  # path.res <- path.res[path.res$p.adjust < 0.05,]
  return(path.res)
})

pathway.cac.go.narm.raw.final.all <- lapply(comm.cac.entrez, function(x) {
  path.res <- enrichGO(gene=x, 
                              OrgDb = org.Mm.eg.db,
                              ont = "BP",
                              universe = unique(mouse.info.entrez.deg),
                              pvalueCutoff=Inf, pAdjustMethod = 'fdr', qvalueCutoff=Inf)

  path.res <- data.frame(path.res) %>% mutate(Description = sapply(str_split(data.frame(path.res)$Description, ' -'), function(x) x[1])) 
  # path.res <- path.res[path.res$p.adjust < 0.05,]
  return(path.res)
})

pathway.breast.kegg.narm.raw.final.all <- lapply(comm.breast.entrez, function(x) {
  path.res <- enrichKEGG(gene=x ,
                            organism = 'hsa',
                            universe = unique(human.info.entrez.genecode),
                            pvalueCutoff=Inf, pAdjustMethod = 'fdr', qvalueCutoff=Inf)

  path.res <- data.frame(path.res) %>% mutate(Description = sapply(str_split(data.frame(path.res)$Description, ' -'), function(x) x[1])) 
  # path.res <- path.res[path.res$p.adjust < 0.05,]
  return(path.res)
})

pathway.cac.kegg.narm.raw.final.all <- lapply(comm.cac.entrez, function(x) {
  path.res <- enrichKEGG(gene=x ,
                            organism = 'mmu',
                            universe =unique(mouse.info.entrez.deg),
                            pvalueCutoff=Inf, pAdjustMethod = 'fdr', qvalueCutoff=Inf)

  path.res <- data.frame(path.res) %>% mutate(Description = sapply(str_split(data.frame(path.res)$Description, ' -'), function(x) x[1])) 
  # path.res <- path.res[path.res$p.adjust < 0.05,]
  return(path.res)
})

save(pathway.breast.go.narm.raw.final.all, pathway.cac.go.narm.raw.final.all, pathway.breast.kegg.narm.raw.final.all, pathway.cac.kegg.narm.raw.final.all, file =  paste0(output_dir, '/pathway.breast.cac.narm.raw.final.all.RData'))
