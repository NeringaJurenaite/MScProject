library(DESeq2)

load("~/Documents/project/unnormaliseddata.Rdata")
load("~/Documents/project/modelmatrices.Rdata")

DESeqres<-function(dataframe, model){
  model$subtype <- factor(model$subtype)
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(dataframe), colData = model, design = ~ subtype*CPE)
  dds <- DESeq(dds)
  return(dds)
}
dds1 <- DESeqres(matrix_temp_2_lumAlumB1, model_lumAlumB)
dds2 <- DESeqres(matrix_temp_2_BasalNormal1, model_BasalNormal)
dds3 <- DESeqres(matrix_temp_2_HER2Basal1, model_HER2Basal)
dds4 <- DESeqres(matrix_temp_2_HER2Normal1, model_HER2Normal)
dds5 <- DESeqres(matrix_temp_2_lumBHER21, model_lumBHER2)
dds6 <- DESeqres(matrix_temp_2_lumBNormal1, model_lumBNormal)
dds7 <- DESeqres(matrix_temp_2_lumBBasal1, model_lumBBasal)
dds8 <- DESeqres(matrix_temp_2_lumANormal1, model_lumANormal)
dds9 <- DESeqres(matrix_temp_2_lumAHER21, model_lumAHER2)
dds10 <- DESeqres(matrix_temp_2_lumABasal1, model_lumABasal)
save(dds1, dds2, dds3, dds4, dds5, dds6, dds7, dds8, dds9, dds10, file="~/Documents/project/dds.Rdata")
