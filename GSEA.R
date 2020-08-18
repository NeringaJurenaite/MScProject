load("~/Documents/project/lumAHER2_C_N_CPE.Rdata")
load("~/Documents/project/lumABasal_C_N_CPE.Rdata")
load("~/Documents/project/lumANormal_C_N_CPE.Rdata")
load("~/Documents/project/lumAlumB_C_N_CPE.Rdata")
load("~/Documents/project/lumBBasal_C_N_CPE.Rdata")
load("~/Documents/project/lumBHER2_C_N_CPE.Rdata")
load("~/Documents/project/lumBNormal_C_N_CPE.Rdata")
load("~/Documents/project/HER2Basal_C_N_CPE.Rdata")
load("~/Documents/project/HER2Normal_C_N_CPE.Rdata")
load("~/Documents/project/BasalNormal_C_N_CPE.Rdata")


GSEA = function(gene_list, GO_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(gage)
  library(fgsea)
  library(ggplot2)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  fgRes <- fgsea(pathways = myGO, stats = gene_list, 
                 minSize=15, eps=0,
                 maxSize=Inf)
  fgRes <- fgRes[order(padj), ]


  gaRes = gage::gage(gene_list, gsets=myGO, rank.test=T, compare='unpaired', set.size =c(10,length(gene_list)))
  gse <- gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "SYMBOL", 
               minGSSize = 10, 
               maxGSSize = length(gene_list), 
               pvalueCutoff = 0.05, eps=0,
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "fdr")
  c5 <- read.gmt(GO_file)
  egmt <- enricher(names(gene_list), TERM2GENE=c5, pvalueCutoff =0.05, pAdjustMethod = "fdr", minGSSize = 10,  maxGSSize = length(gene_listlumAlumB),  
                   qvalueCutoff = 0.05)

  ups = as.data.frame(gaRes$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < 0.05 ) %>%
    dplyr::select("Pathway")
  
  downs = as.data.frame(gaRes$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < 0.05 ) %>%
    dplyr::select("Pathway")
  
  keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
  keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
  
  ## Collapse redundant pathways
  Up = collapsePathways(keepups[order(keepups$pval)][keepups$padj < 0.05], pathways = myGO, gene_list, pval.threshold = 0.05)
  Down = fgsea::collapsePathways(keepdowns[order(keepdowns$pval)][keepdowns$padj < 0.05], myGO, gene_list,  nperm = 500, pval.threshold = 0.05) 
  #fgRes = fgRes[ !is.na(match(fgRes$pathway, c( Up$mainPathways, Down$mainPathways))), ] %>% arrange(desc(NES))
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  
  fgRes$pathway <- gsub('_', ' ', fgRes$pathway, fixed = TRUE)
  filtRes = rbind(head(fgRes, n = 20),
                  tail(fgRes, n = 20 ))
  
  g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_segment( aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
    geom_point( size=5, aes( fill = Enrichment),
                shape=21, stroke=2) +
    scale_fill_manual(values = c("Down-regulated" = "darkgreen",
                                 "Up-regulated" = "purple4") ) +
    coord_flip() + labs(x="Pathway", y="Normalized Enrichment Score",title="Normal Cell Associated Significant GO Pathways") +
    theme_minimal()
  g
  
  p1 <- cnetplot(egmt, foldChange=gene_list)
  p2 <- cnetplot(egmt, foldChange=gene_list, categorySize="qvalue")
  p3 <- cnetplot(egmt, foldChange=gene_list, circular=T, colorEdge = TRUE, node_label="all", categorySize="qvalue")
  p5 <- heatplot(egmt, foldChange=gene_list)
  output = list("Results" = fgRes, "Plot" = g, "egmt" = egmt,"gse" = gse , "p1"=p1,"p2"=p2,"p3"=p3,"p5"=p5 )
  return(output)
}
?cnetplot

library(ensembldb)
library(EnsDb.Hsapiens.v86)

GO_file = "~/Downloads/c5.all.v7.1.symbols.gmt"

geneslumAlumB <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(normallumAlumB), keytype =   "GENEID", columns = c(  "GENEID", "SYMBOL", "ENTREZID"))
geneslumAlumB<- geneslumAlumB[!duplicated(geneslumAlumB$GENEID),]
geneslumAlumB$logFC <- normallumAlumB[which(rownames(normallumAlumB) %in% geneslumAlumB$GENEID),][geneslumAlumB$GENEID,]$logFC
gene_listlumAlumB = as.numeric(geneslumAlumB$logFC )
names(gene_listlumAlumB) = geneslumAlumB$SYMBOL
all.equal(geneslumAlumB$GENEID, rownames(normallumAlumB[which(rownames(normallumAlumB) %in% geneslumAlumB$GENEID),][geneslumAlumB$GENEID,]))
gene_listlumAlumB = sort(gene_listlumAlumB, decreasing = T)
gene_listlumAlumB = gene_listlumAlumB[!duplicated(names(gene_listlumAlumB))]
RESnormallumAlumB = GSEA(gene_listlumAlumB, GO_file, pval = 0.05)
dim(RESnormallumAlumB$Results)
lumAlumBgsea <- RESnormallumAlumB$Results ####################
RESnormallumAlumB$Plot#160
ggsave("RESnormallumAlumB$Plot.png")
RESnormallumAlumB$p1
RESnormallumAlumB$p2
RESnormallumAlumB$p3
RESnormallumAlumB$p5


geneslumANormal <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(normallumANormal), keytype =   "GENEID", columns = c(  "GENEID", "SYMBOL", "ENTREZID"))
geneslumANormal<- geneslumANormal[!duplicated(geneslumANormal$GENEID),]
geneslumANormal$logFC <- normallumANormal[which(rownames(normallumANormal) %in% geneslumANormal$GENEID),][geneslumANormal$GENEID,]$logFC
gene_listlumANormal = as.numeric(geneslumANormal$logFC )
names(gene_listlumANormal) = geneslumANormal$SYMBOL
all.equal(geneslumANormal$GENEID, rownames(normallumANormal[which(rownames(normallumANormal) %in% geneslumANormal$GENEID),][geneslumANormal$GENEID,]))
gene_listlumANormal = sort(gene_listlumANormal, decreasing = T)
gene_listlumANormal = gene_listlumANormal[!duplicated(names(gene_listlumANormal))]
RESnormallumANormal = GSEA(gene_listlumANormal, GO_file, pval = 0.05)
dim(RESnormallumANormal$Results)
lumANormalgsea <- RESnormallumANormal$Results
RESnormallumANormal$Plot #160
ggsave("RESnormallumANormal$Plot.png")
RESnormallumANormal$p1
RESnormallumANormal$p2
RESnormallumANormal$p3
RESnormallumANormal$p5


geneslumABasal <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(normallumABasal), keytype =   "GENEID", columns = c(  "GENEID", "SYMBOL", "ENTREZID"))
geneslumABasal<- geneslumABasal[!duplicated(geneslumABasal$GENEID),]
geneslumABasal$logFC <- normallumABasal[which(rownames(normallumABasal) %in% geneslumABasal$GENEID),][geneslumABasal$GENEID,]$logFC
gene_listlumABasal = as.numeric(geneslumABasal$logFC )
names(gene_listlumABasal) = geneslumABasal$SYMBOL
all.equal(geneslumABasal$GENEID, rownames(normallumABasal[which(rownames(normallumABasal) %in% geneslumABasal$GENEID),][geneslumABasal$GENEID,]))
gene_listlumABasal = sort(gene_listlumABasal, decreasing = T)
gene_listlumABasal = gene_listlumABasal[!duplicated(names(gene_listlumABasal))]
RESnormallumABasal = GSEA(gene_listlumABasal, GO_file, pval = 0.05)
dim(RESnormallumABasal$results)
lumABasalgsea <- RESnormallumABasal$Results ################################
RESnormallumABasal$Plot #160
ggsave("RESnormallumABasal$Plot.png")
RESnormallumABasal$p1
RESnormallumABasal$p2
RESnormallumABasal$p3
RESnormallumABasal$p5

geneslumAHER2 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(normallumAHER2), keytype =   "GENEID", columns = c(  "GENEID", "SYMBOL", "ENTREZID"))
geneslumAHER2<- geneslumAHER2[!duplicated(geneslumAHER2$GENEID),]
geneslumAHER2$logFC <- normallumAHER2[which(rownames(normallumAHER2) %in% geneslumAHER2$GENEID),][geneslumAHER2$GENEID,]$logFC
gene_listlumAHER2 = as.numeric(geneslumAHER2$logFC )
names(gene_listlumAHER2) = geneslumAHER2$SYMBOL
all.equal(geneslumAHER2$GENEID, rownames(normallumAHER2[which(rownames(normallumAHER2) %in% geneslumAHER2$GENEID),][geneslumAHER2$GENEID,]))
gene_listlumAHER2 = sort(gene_listlumAHER2, decreasing = T)
gene_listlumAHER2 = gene_listlumAHER2[!duplicated(names(gene_listlumAHER2))]
RESnormallumAHER2 = GSEA(gene_listlumAHER2, GO_file, pval = 0.05)
dim(RESnormallumAHER2$Results)
lumAHER2gsea <- RESnormallumAHER2$Results
RESnormallumAHER2$Plot #160
ggsave("RESnormallumAHER2$Plot.png")
RESnormallumAHER2$p1
RESnormallumAHER2$p2
RESnormallumAHER2$p3
RESnormallumAHER2$p5
gene_list <- gene_listlumAHER2

geneslumBHER2 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(normallumBHER2), keytype =   "GENEID", columns = c(  "GENEID", "SYMBOL", "ENTREZID"))
geneslumBHER2<- geneslumBHER2[!duplicated(geneslumBHER2$GENEID),]
geneslumBHER2$logFC <- normallumBHER2[which(rownames(normallumBHER2) %in% geneslumBHER2$GENEID),][geneslumBHER2$GENEID,]$logFC
gene_listlumBHER2 = as.numeric(geneslumBHER2$logFC )
names(gene_listlumBHER2) = geneslumBHER2$SYMBOL
all.equal(geneslumBHER2$GENEID, rownames(normallumBHER2[which(rownames(normallumBHER2) %in% geneslumBHER2$GENEID),][geneslumBHER2$GENEID,]))
gene_listlumBHER2 = sort(gene_listlumBHER2, decreasing = T)
gene_listlumBHER2 = gene_listlumBHER2[!duplicated(names(gene_listlumBHER2))]
RESnormallumBHER2 = GSEA(gene_listlumBHER2, GO_file, pval = 0.05)
dim(RESnormallumBHER2$Results)
lumBHER2gsea <- RESnormallumBHER2$Results
RESnormallumBHER2$Plot #160
ggsave("RESnormallumBHER2$Plot.png")
RESnormallumBHER2$p1
RESnormallumBHER2$p2
RESnormallumBHER2$p3
RESnormallumBHER2$p5

geneslumBBasal <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(normallumBBasal), keytype =   "GENEID", columns = c(  "GENEID", "SYMBOL", "ENTREZID"))
geneslumBBasal<- geneslumBBasal[!duplicated(geneslumBBasal$GENEID),]
geneslumBBasal$logFC <- normallumBBasal[which(rownames(normallumBBasal) %in% geneslumBBasal$GENEID),][geneslumBBasal$GENEID,]$logFC
gene_listlumBBasal = as.numeric(geneslumBBasal$logFC )
names(gene_listlumBBasal) = geneslumBBasal$SYMBOL
all.equal(geneslumBBasal$GENEID, rownames(normallumBBasal[which(rownames(normallumBBasal) %in% geneslumBBasal$GENEID),][geneslumBBasal$GENEID,]))
gene_listlumBBasal = sort(gene_listlumBBasal, decreasing = T)
gene_listlumBBasal = gene_listlumBBasal[!duplicated(names(gene_listlumBBasal))]
RESnormallumBBasal = GSEA(gene_listlumBBasal, GO_file, pval = 0.05)
dim(RESnormallumBBasal$Results)
lumBBasalgsea <- RESnormallumBBasal$Results ##############################
RESnormallumBBasal$Plot #160
ggsave("RESnormallumBBasal$Plot.png")
RESnormallumBBasal$p1
RESnormallumBBasal$p2
RESnormallumBBasal$p3
RESnormallumBBasal$p5


geneslumBNormal <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(normallumBNormal), keytype =   "GENEID", columns = c(  "GENEID", "SYMBOL", "ENTREZID"))
geneslumBNormal<- geneslumBNormal[!duplicated(geneslumBNormal$GENEID),]
geneslumBNormal$logFC <- normallumBNormal[which(rownames(normallumBNormal) %in% geneslumBNormal$GENEID),][geneslumBNormal$GENEID,]$logFC
gene_listlumBNormal = as.numeric(geneslumBNormal$logFC )
names(gene_listlumBNormal) = geneslumBNormal$SYMBOL
all.equal(geneslumBNormal$GENEID, rownames(normallumBNormal[which(rownames(normallumBNormal) %in% geneslumBNormal$GENEID),][geneslumBNormal$GENEID,]))
gene_listlumBNormal = sort(gene_listlumBNormal, decreasing = T)
gene_listlumBNormal = gene_listlumBNormal[!duplicated(names(gene_listlumBNormal))]
RESnormallumBNormal = GSEA(gene_listlumBNormal, GO_file, pval = 0.05)
dim(RESnormallumBNormal$Results)
lumBNormalgsea <- RESnormallumBNormal$Results
RESnormallumBNormal$Plot #160
ggsave("RESnormallumBNormal$Plot.png")
RESnormallumBNormal$p1
RESnormallumBNormal$p2
RESnormallumBNormal$p3
RESnormallumBNormal$p5



genesHER2Normal <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(normalHER2Normal), keytype =   "GENEID", columns = c(  "GENEID", "SYMBOL", "ENTREZID"))
genesHER2Normal<- genesHER2Normal[!duplicated(genesHER2Normal$GENEID),]
genesHER2Normal$logFC <- normalHER2Normal[which(rownames(normalHER2Normal) %in% genesHER2Normal$GENEID),][genesHER2Normal$GENEID,]$logFC
gene_listHER2Normal = as.numeric(genesHER2Normal$logFC )
names(gene_listHER2Normal) = genesHER2Normal$SYMBOL
all.equal(genesHER2Normal$GENEID, rownames(normalHER2Normal[which(rownames(normalHER2Normal) %in% genesHER2Normal$GENEID),][genesHER2Normal$GENEID,]))
gene_listHER2Normal = sort(gene_listHER2Normal, decreasing = T)
gene_listHER2Normal = gene_listHER2Normal[!duplicated(names(gene_listHER2Normal))]
RESnormalHER2Normal = GSEA(gene_listHER2Normal, GO_file, pval = 0.05)
dim(RESnormalHER2Normal$Results)
HER2Normalgsea <- RESnormalHER2Normal$Results
RESnormalHER2Normal$Plot #160
ggsave("RESnormalHER2Normal$Plot.png")
RESnormalHER2Normal$p1
RESnormalHER2Normal$p2
RESnormalHER2Normal$p3
RESnormalHER2Normal$p5


genesHER2Basal <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(normalHER2Basal), keytype =   "GENEID", columns = c(  "GENEID", "SYMBOL", "ENTREZID"))
genesHER2Basal<- genesHER2Basal[!duplicated(genesHER2Basal$GENEID),]
genesHER2Basal$logFC <- normalHER2Basal[which(rownames(normalHER2Basal) %in% genesHER2Basal$GENEID),][genesHER2Basal$GENEID,]$logFC
gene_listHER2Basal = as.numeric(genesHER2Basal$logFC )
names(gene_listHER2Basal) = genesHER2Basal$SYMBOL
all.equal(genesHER2Basal$GENEID, rownames(normalHER2Basal[which(rownames(normalHER2Basal) %in% genesHER2Basal$GENEID),][genesHER2Basal$GENEID,]))
gene_listHER2Basal = sort(gene_listHER2Basal, decreasing = T)
gene_listHER2Basal = gene_listHER2Basal[!duplicated(names(gene_listHER2Basal))]
RESnormalHER2Basal = GSEA(gene_listHER2Basal, GO_file, pval = 0.05)
dim(RESnormalHER2Basal$Results)
HER2Basalgsea <- RESnormalHER2Basal$Results
RESnormalHER2Basal$Plot #160
ggsave("RESnormalHER2Basal$Plot.png")
RESnormalHER2Basal$p1
RESnormalHER2Basal$p2
RESnormalHER2Basal$p3
RESnormalHER2Basal$p5


genesBasalNormal <- ensembldb::select(EnsDb.Hsapiens.v86, keys= rownames(normalBasalNormal), keytype =   "GENEID", columns = c(  "GENEID", "SYMBOL", "ENTREZID"))
genesBasalNormal<- genesBasalNormal[!duplicated(genesBasalNormal$GENEID),]
genesBasalNormal$logFC <- normalBasalNormal[which(rownames(normalBasalNormal) %in% genesBasalNormal$GENEID),][genesBasalNormal$GENEID,]$logFC
gene_listBasalNormal = as.numeric(genesBasalNormal$logFC )
names(gene_listBasalNormal) = genesBasalNormal$SYMBOL
all.equal(genesBasalNormal$GENEID, rownames(normalBasalNormal[which(rownames(normalBasalNormal) %in% genesBasalNormal$GENEID),][genesBasalNormal$GENEID,]))
gene_listBasalNormal = sort(gene_listBasalNormal, decreasing = T)
gene_listBasalNormal = gene_listBasalNormal[!duplicated(names(gene_listBasalNormal))]
RESnormalBasalNormal = GSEA(gene_listBasalNormal, GO_file, pval = 0.05)
dim(RESnormalBasalNormal$Results)
BasalNormalgsea <- RESnormalBasalNormal$Results
RESnormalBasalNormal$Plot #160
ggsave("RESnormalBasalNormal$Plot.png")
RESnormalBasalNormal$p1
RESnormalBasalNormal$p2
RESnormalBasalNormal$p3
RESnormalBasalNormal$p5

save(BasalNormalgsea,HER2Basalgsea, HER2Normalgsea, lumBNormalgsea , lumBHER2gsea , lumBBasalgsea, lumAlumBgsea , lumANormalgsea,lumAHER2gsea , lumABasalgsea, file="~/Documents/project/GSEARESnormal _C5_Results.Rdata")















