library("DEP")
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(colorspace)
library(gplots)
library(limma) 
library("SummarizedExperiment")
library("ComplexHeatmap")
library(remotes)
library("reshape2")
library(ggbiplot)
library(factoextra)
library(imguR)
library(photobiology)
library(pca3d)

PRO <- data.table::fread("C:/Users/k16206jp/OneDrive - The University of Manchester/OMD paper joint editing folder/6 OMD Paper/OMD Cells MS analysis/Proteome/proteinGroups filtered no missing values.txt")

names(PRO)[names(PRO) == "Gene names"] <- "Gene.names"
names(PRO)[names(PRO) == "Protein names"] <- "Protein.IDs"

#log2 transform data
intensityPRO <- PRO %>% 
  select(grep("0", colnames(.))) 
ColumnTitles <- PRO[,c(13:31)]

IntensityPROTitles <- cbind(ColumnTitles,intensityPRO)
IntensityPROTitles[IntensityPROTitles== -Inf] <- NA
IntensityPROTitles[IntensityPROTitles== "NaN"] <- NA

#import txt file with experimental design on label, condition, replicate
experimental_design <- data.table::fread("C:/Users/k16206jp/OneDrive - The University of Manchester/OMD paper joint editing folder/6 OMD Paper/OMD Cells MS analysis/Proteome/experimental design for DEP.txt")
data_unique <- make_unique(IntensityPROTitles, "Gene.names", "Protein.IDs", delim = ";")

#prepare data for plots
PRO_columns <- grep("0", colnames(data_unique)) 
# get LFQ column numbers
data_se <- make_se(data_unique, PRO_columns, experimental_design)

#plots
plot_frequency(data_se)
plot_numbers(data_se)
plot_coverage(data_se)

##############################################
#normaliseQuantiles
data_norm_quant <- data_se
assay(data_norm_quant)<-limma::normalizeQuantiles(assay(data_se))

plot_normalization(data_se, data_norm_quant)

#plot missing value presence of  data
plot_missval(data_se)
plot_missval(data_norm_quant)
plot_detect(data_se)

####################################################
#data imputation to the left of centre and plotting
data_imp <- impute(data_norm_quant, fun = "man", shift = 1.8, scale = 0.3)
plot_imputation(data_norm_quant, data_imp)

#save normalised dataframe pre imputation
#remove the word filtered if using unfiltered data
filtered.norm.data.limma.test <- as.data.frame(assays(data_norm_quant)[[1]])
write.table(filtered.norm.data.limma.test, file="C:/Users/k16206jp/OneDrive - The University of Manchester/OMD paper joint editing folder/6 OMD Paper/OMD Cells MS analysis/Proteome/DEP analysis no missing values/filtered normalised proteome pre imputation.txt", sep="\t", quote=F)
filtered.norm.data.limma.test.with.titles <- cbind(ColumnTitles,filtered.norm.data.limma.test)
write.table(filtered.norm.data.limma.test.with.titles, file="C:/Users/k16206jp/OneDrive - The University of Manchester/OMD paper joint editing folder/6 OMD Paper/OMD Cells MS analysis/Proteome/DEP analysis no missing values/filtered normalised proteome with titles.txt", sep="\t", quote=F)


########################################################
#correlation 
cors <- cor(filtered.norm.data.limma.test, use='complete.obs', method='pearson')
rownames(cors) <- colnames(cors) <- sub("LFQ intensity ", "", colnames(cors))

heatmap.2(cors, 
          Rowv = F, Colv=F, dendrogram = "none", # Dendogram
          trace = 'none', # remove trace in heatmap
          col=colorRampPalette(c('cadetblue1','white','red')), # find more colors: colors()
          margins = c(7,7)) 
##################################################################
#PCA
#go to perseus

#####################################################################################
#limma differential expression test

#here we change it depending what we want to compare. for example if we want to compare T(tumor) with All the rest we write 5 Ts and 8Fs
#take a look at the file
BOMBOO <- filtered.norm.data.limma.test[,c(7:12)]
head(BOMBOO)
m.design <- c("N","N","N","Y","Y","Y")
colnames(BOMBOO) <- m.design
head(BOMBOO)
design <- model.matrix(~0+factor(m.design, levels=unique(m.design)))
design
colnames(design) <- unique(m.design)
fit <- lmFit(BOMBOO, design)
cont.matrix <- makeContrasts("BOO vs BOM" =Y - N,
                             levels=unique(m.design))
cont.matrix
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)
vennDiagram(results, include=c("up", "down"), counts.col=c("red","green"))

dev.off();
topTableF(fit2, adjust="BH")

tT <- topTableF(fit2, adjust.method="fdr", number=nrow(fit2))
tT$gene <- rownames(tT)

#make a volcano plot with the results
v1 <-
  ggplot(tT, aes(y=-log10(P.Value), x="N vs Y")) +
  #  scale_colour_manual(values = c("red","springgreen3", "dodgerblue1","yellow"))+
  geom_point(alpha = 0.6, size=3) +
  geom_point(data=tT[tT$adj.P.Val<=0.1,], color="black", fill="red", size=3, shape=21) +
  scale_color_manual(values = c("FALSE" = "dodgerblue1", "TRUE" = "springgreen3")) +
  #geom_text_repel(data=tT[tT$adj.P.Val<=0.1,], aes(label=gene),
  #geom_text_repel(data=tT[which(tT$sig==TRUE|tT$gene ==bait.gene.name),], aes(label=gene),
  # segment.size = 0.2, segment.color = "grey50", size = 3,
  #  nudge_x= 0.2,nudge_y= 0.2 )+ 
  #ggtitle(paste0("Proteins detected: ",length(tT$id),"\n", "Significant interactors: ", length(selected$id)))+
  theme_bw()+
  theme(aspect.ratio=4/4)+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    #legend.text =element_text(size=16),
    #legend.title = element_text(size=16),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )+
  xlab("Log2FC") +
  ylab("-log10(P)") 
v1

#export results
write.table(tT, file="C:/Users/k16206jp/OneDrive - The University of Manchester/OMD paper joint editing folder/6 OMD Paper/OMD Cells MS analysis/Proteome/DEP analysis no missing values/Adjusted P value Data BOM vs BOO.txt", sep="\t",quote=F, row.names = F)
