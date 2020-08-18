#plot t-SNE plots of the normalised expression data to obtain a visualisation of subtype clustering:
load("~/Documents/project/normalisedALLandDATASUBTYPES_MY.Rdata")
load("~/Documents/project/datasubtypes_my.Rdata")
library(Rtsne)
library("ggsci")
library("ggplot2")
library("gridExtra")

tsne_results <- Rtsne(t(yCPE), perplexity=39, check_duplicates = FALSE, theta=0.0, normalize=F) # You can change the value of perplexity and see how the plot changes
## Generate the t_SNE plot
datasubtypes_my_new <- datasubtypes_my
datasubtypes_my_new$subtype <- as.character(datasubtypes_my_new$subtype)
datasubtypes_my_new[datasubtypes_my_new$subtype == "BRCA.LumA",]$subtype <- "Luminal A"
datasubtypes_my_new[datasubtypes_my_new$subtype == "BRCA.LumB",]$subtype <- "Luminal B"
datasubtypes_my_new[datasubtypes_my_new$subtype == "BRCA.Normal",]$subtype <- "Normal-like"
datasubtypes_my_new[datasubtypes_my_new$subtype == "BRCA.Basal",]$subtype <- "Basal-like"
datasubtypes_my_new[datasubtypes_my_new$subtype == "BRCA.Her2",]$subtype  <- "HER2-enriched"

df <- data.frame(x = tsne_results$Y[,1],
                 y = tsne_results$Y[,2],
                 colour = as.character(datasubtypes_my_new$subtype))

ggplot(df, aes(x, y, colour = colour)) +geom_point()+ scale_fill_d3() + theme_classic()
g <-ggplot(df, aes(x, y, colour = colour)) +
  geom_point(size=1) + 
   xlab(" ") + # new label for x axis
  ylab (" ") +
  scale_color_locuszoom()+
  ggtitle("Normalised Count t-SNE Plot") +
  theme_bw()  +
  theme (axis.text = element_blank(),
         axis.ticks = element_blank(),
        axis.title = element_text(size =12),
        legend.title = element_blank(),
        legend.text = element_text(face="italic",  size=10),
        plot.title = element_text(face="bold",  size=14))

g 
ggsave("t-SNEofALL.png") # do not use rstudio if this keeps crashing
