

input_dir = paste0('../OneDrive-GCCORP/Community distance/MoBC_data_code/input_file')
output_dir = paste0('../OneDrive-GCCORP/Community distance/MoBC_data_code/output_file')
dir.comm = paste0('../OneDrive-GCCORP/Community distance/community_distance_file')

# library
library(rtracklayer)
library(dplyr)
library(stringr)
library(DESeq2)


# brca input file
tcga.brca.tpm <- read.delim(paste0(input_dir, "/TcgaTargetGtex_rsem_gene_tpm.tsv"), header =TRUE, row.names = 1)
tcga.phe.info <- read.csv(paste0(input_dir, "/TcgaTargetGTEX_phenotype.csv"))
tcga.brca.count <- read.delim(paste0(input_dir, "/TcgaTargetGtex_gene_expected_count.tsv"), header=TRUE, row.names = 1)

# input data processing - info, tpm data processing
colnames(tcga.phe.info) <- c("ID", "detailed_category", "primary_disease", "primary_site", "sample_type", "gender", 'study')
tcga.phe.info.breast <- tcga.phe.info[tcga.phe.info$primary_site == 'Breast',]
tcga.phe.info.breast.normal <- tcga.phe.info.breast[tcga.phe.info.breast$sample_type %in% "Solid Tissue Normal",]
tcga.cancer.normal.id <- sapply(sapply(tcga.phe.info.breast.normal$ID, function(x) str_split(x, '-11')[[1]][1]), function(y) grep(y, tcga.phe.info.breast$ID)) %>% unlist
tcga.cancer.normal.id.results <- tcga.phe.info.breast$ID[tcga.cancer.normal.id] 

# data info table
tcga.cancer.normal.info <- tcga.phe.info.breast[tcga.phe.info.breast$ID %in% tcga.cancer.normal.id.results[!tcga.cancer.normal.id.results %in% c("TCGA-BH-A0AY-11", "TCGA-E2-A15K-06", "TCGA-BH-A18V-06", "TCGA-BH-A1FE-06")],]
tcga.cancer.normal.info$condition <- ifelse(tcga.cancer.normal.info$sample_type == 'Solid Tissue Normal', "Normal", "Tumor")
tcga.cancer.normal.info %>% dim

# tpm data processing
colnames(tcga.brca.tpm) <- lapply(colnames(tcga.brca.tpm), function(x) str_replace_all(x, '[.]', '-')) %>% unlist
# colnames(tcga.brca.tpm)
colnames(tcga.brca.tpm)[grep('-11', colnames(tcga.brca.tpm))]
tcga.cancer.normal.coln <- colnames(tcga.brca.tpm)[colnames(tcga.brca.tpm) %in% tcga.cancer.normal.info$ID]
tcga.cancer.normal.coln <- c(tcga.cancer.normal.coln[grep('-11', tcga.cancer.normal.coln)], tcga.cancer.normal.coln[grep('-01', tcga.cancer.normal.coln)])
tcga.cancer.normal.coln %>% length
# tcga.brca.tpm.paired %>% dim

tcga.brca.tpm.paired <- tcga.brca.tpm[,tcga.cancer.normal.coln]
tcga.brca.tpm.paired <- 2^(tcga.brca.tpm.paired) - 0.0001
id.paired <- gsub('\\..*','',rownames(tcga.brca.tpm.paired))
rownames(tcga.brca.tpm.paired) <- id.paired

# gene mapping - genecode
z <- import('gencode.v46.annotation.gtf')
# z <- import('gencode.v22.annotation.gtf')
z1 = subset(z, gene_type=='protein_coding' & type=='gene')
z1 = z1[,c('gene_id','gene_name')] %>% unique
z1 = z1@elementMetadata %>% as.data.frame
z1 %>% dim # 20062
z %>% length
z1 %>% head
z1.new <- gsub('\\..*','',z1$gene_id)
z1.new %>% length # 20062
save(z, z1, z1.new, file = paste0(output_dir, '/gencode.geneid.RData'))
# load(paste0(dir.comm, '/gencode.geneid.RData')) # z, z1, z1.new

tcga.brca.tpm.paired.2 <- tcga.brca.tpm.paired[rownames(tcga.brca.tpm.paired) %in% z1.new,]
tcga.brca.tpm.paired.2 <- tcga.brca.tpm.paired.2[!(tcga.brca.tpm.paired.2 %>% apply(1, sd) <= 0.1 ),]
tcga.brca.tpm.paired %>% dim
tcga.brca.tpm.paired.2 %>% dim

# save(tcga.brca.tpm.paired, tcga.brca.tpm.paired.2, file= paste0(output_dir, '/tcga.brca.tpm.RData'))


# count data preprocessing
colnames(tcga.brca.count) <- lapply(colnames(tcga.brca.count), function(x) str_replace_all(x, '[.]', '-')) %>% unlist
tcga.cancer.normal.coln.count <- colnames(tcga.brca.count)[colnames(tcga.brca.count) %in% tcga.cancer.normal.info$ID]
tcga.cancer.normal.coln.count <- c(tcga.cancer.normal.coln.count[grep('-11', tcga.cancer.normal.coln.count)], tcga.cancer.normal.coln.count[grep('-01', tcga.cancer.normal.coln.count)])
tcga.brca.count.paired <- tcga.brca.count[,tcga.cancer.normal.coln.count]

tcga.brca.count.paired <- 2^(tcga.brca.count.paired) - 1
tcga.brca.count.paired %>% dim
id.paired <- gsub('\\..*','',rownames(tcga.brca.count.paired))
rownames(tcga.brca.count.paired) <- id.paired

tcga.brca.count.paired <- tcga.brca.count.paired[rownames(tcga.brca.count.paired) %in% z1.new,]

# gene filtering - with tpm filtered data 
tcga.brca.count.paired  <- tcga.brca.count.paired[!(tcga.brca.count.paired %>% apply(1, sum) == 0 ),]
tcga.brca.count.paired.2 <- tcga.brca.count.paired[rownames(tcga.brca.count.paired) %in% rownames(tcga.brca.tpm.paired.2),]
tcga.brca.count.paired.2 %>% dim


save(tcga.cancer.normal.info, file = paste0(output_dir, '/tcga.brca.coldata.RData'))


# DEG analysis
dds.brca <- DESeqDataSetFromMatrix(countData = round(tcga.brca.count.paired.2)[,tcga.cancer.normal.info$ID],  
                      colData = tcga.cancer.normal.info, 
                      design = ~ condition)
dds.brca <- DESeq(dds.brca)

res.sd <- results(dds.brca, contrast=c('condition', 'Tumor', 'Normal'))
res.sd <- res.sd[!is.na(res.sd$stat),]

# ppi network gene mapping
human.alias <- read.delim(paste0(input_dir, '/9606.protein.aliases.v12.0.txt'), header=TRUE)
human.info <- read.delim(paste0(input_dir, '/9606.protein.info.v12.0.txt'), header=TRUE)
colnames(human.alias) <- c('protein_id', 'alias', 'source')
colnames(human.info) <- c('protein_id', 'Symbol', 'protein_size', 'annotation')


ids = gsub('\\..*','',rownames(res.sd))
res1 = res.sd %>% as.data.frame
res1$gene_id = ids
prid = subset(human.alias, alias %in% ids & source =='Ensembl_gene')
prid = prid[,c('protein_id','alias')] %>% unique %>% 'colnames<-'(c('protein_id','gene_id'))
gnid = human.info[,c(1:2)] %>% unique
allid = merge(gnid, prid, all=T)
res2.sd = merge(res1, allid)
res.sd.2 = merge(res1, allid, all.x=TRUE)
res.sd.2 %>% dim
# brca.deg %>% dim
# write.csv(brca.deg, file = paste0(output_dir , '/brca.deg.csv'), row.names=FALSE)

save(tcga.brca.tpm.paired, tcga.brca.tpm.paired.2, file= paste0(output_dir, '/tcga.brca.tpm.RData'))
save(tcga.brca.count.paired.2, file= paste0(output_dir, '/tcga.brca.count.RData'))

save(res.sd, file = paste0(output_dir, '/breast.deg.result.final.RData')) # 17520 (gene mapping 이전 deseq results)
save(res2.sd, file = paste0(output_dir, '/breast.deseq.results.RData')) # 17323 (gene mapping 이후 deseq results)
save(res.sd.2, file = paste0(output_dir, '/breast.deg.result.final.2.RData')) # 17520 (gene mapping 이후 deseq results)
