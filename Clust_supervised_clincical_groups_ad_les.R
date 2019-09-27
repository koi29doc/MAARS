################################################################
###### Supervised clustering based on clinical groups #####
#### from the script Buit_sup_gp_for_hclust ######
################################################################


### Define working directory
setwd("/Users/alain_lu/Dropbox/These_de_science/Data/Processed/Clustering/")

### Load library
library(anyLib)
anyLib::anyLib(c("Biobase", 
                 "cluster",
                 "ComplexHeatmap",
                 "EMA",
                 "factoextra",
                 "genefilter",
                 "ggplot2",
                 "limma",
                 "NbClust",
                 "seriation",
                 "tibble"
))

### Load data
load(file='/Users/alain_lu/Dropbox/These_de_science/Data/Processed/Lima/exset2.RData')
exset2 ## Limma Matrix which contains several features (expression, sample annotation ...)
ex1 <- Biobase::exprs(exset2) ## extraction of expression matrix from exset2
dim(ex1)
class(ex1)
annot <- exset2@phenoData@data ## extraction of sample annotation from extset2 ## for merging an <-not$id == colnames(ex1) and colnames(ex1) == annot$id 
GS <- as.character(exset2@featureData@data$Gene.Symbol) ## gene list vector
fd <- exset2@featureData@data ## gene matrix with symbol, gene name, localisation
genemap <- read.csv2("/Users/alain_lu/Dropbox/These_de_science/Data/Unprocessed/genemap.csv", row.names = 1)
annot_plus_ad <- read.csv2(file = "/Users/alain_lu/Dropbox/These_de_science/Data/Unprocessed/MAARS_AD_full_20190131_12-34-49.csv", header = T, row.names = 1, sep = "\t")

### Clean feature data: only keep unique gene symbol
#fd <- cbind(fd, rownames(fd))
fd <- unique(fd)

### Keep only AD_LES data
## for clinical features
annot_ad_les <- annot[annot$clinical_group == "AD" & annot$lesional == "LES", ] # Clinical data for AD
dim(annot_ad_les)
id_sample_ad_les <- rownames(annot_ad_les)
## add more clinical feature by merging the annot to the annot_plus
dim(annot_plus_ad)
annot_plus_ad$involved.skin.biopsy.involved.skin.biopsy.MAARS.Sample.identifier..MAARS_Sample_identifier. %in% rownames(annot_ad_les) # the common column are ad_plus$involved.skin.biopsy.[...] and rownames(ad)
annot_plus_ad <- annot_plus_ad[annot_plus_ad$involved.skin.biopsy.involved.skin.biopsy.MAARS.Sample.identifier..MAARS_Sample_identifier. %in% rownames(annot_ad_les) == TRUE,]
rownames(annot_plus_ad) <- annot_plus_ad$involved.skin.biopsy.involved.skin.biopsy.MAARS.Sample.identifier..MAARS_Sample_identifier.
annot_ad_les <- merge(annot_ad_les, annot_plus_ad, by = "row.names")
colnames(annot_ad_les)
rownames(annot_ad_les) <- annot_ad_les$Row.names

## for genetic expression
ex1_ad_les <- ex1[ ,id_sample_ad_les]

### Compare groups with ttest

## Hight IgE vs normal IgE level
for(i in 1:nrow(ex1_ad_les)) {
  test <- t.test(as.numeric(ex1_ad_les[i, hight_IgE]), as.numeric(ex1_ad_les[i,-match(hight_IgE, colnames(ex1_ad_les))]))
  if(i == 1) res <- c(test$p.value, test$estimate)
  else res <- cbind(res, c(test$p.value, test$estimate))
}

### No p value under 0,9 -> this should be a problem !!

pvalue_ttest_hight_IgE <- res

padj_ttest_hight_IgE <- p.adjust(res[1,], method = "BH")

ttest_hight_IgE <- rbind(pvalue_ttest_hight_IgE, padj_ttest_hight_IgE)
colnames(ttest_hight_IgE) <- rownames(ex1_ad_les)
ttest_hight_IgE <- t(ttest_hight_IgE)
colnames(ttest_hight_IgE) <- c("pvalue_ttest", "mean_hight_IgE", "mean_normal_IgE","padj_ttest")
ttest_hight_IgE <- ttest_hight_IgE[ ,c(2,3,1,4)]
ttest_hight_IgE <- as.matrix(ttest_hight_IgE)

DEG_hight_igE_nl_IgE <- ttest_hight_IgE[ttest_hight_IgE[ ,4] < 0.05, ]
DEG_hight_igE_nl_IgE <- as.matrix(DEG_hight_igE_nl_IgE)
colnames(DEG_hight_igE_nl_IgE) <- c("GeneSymbol", "mean_male", "mean_female", "pvalue_ttest", "padj_ttest")

diff_means <- as.numeric(DEG_male_female_ctrl_MAARS[ ,2]) - as.numeric(DEG_male_female_ctrl_MAARS[ ,3]) 
DEG_male_female_ctrl_MAARS <- cbind(DEG_male_female_ctrl_MAARS, diff_means)

### Clustering whith complexe heatmap

col <- circlize::colorRamp2(c(0, 5.5, 11), c("green", "white", "red")) # Coloration function for Complexheatmap

## Clustering on most variable genes specifically in ad_les in comparision to healthy controls
var_ad_les_ctrl_top32 <- var_ad_les_ctrl_top100[1:32,] # Gene selection of variance twice superior to healthy donnor
ex1_ad_les_top32 <- ex1_ad_les_top100[rownames(var_ad_les_ctrl_top32), ]

hm <- Heatmap(as.matrix(ex1_ad_les_top32),
              clustering_distance_rows = "euclidean",
              clustering_distance_columns = "euclidean",
              column_split = 4,
              row_split = 4,
              name = "Intensity", 
              column_title = "Samples",
              row_title = "Genes",
              col=col,
              cluster_rows = T,
              cluster_columns = T,
              show_row_names = T,
              show_column_names = T,
              # top_annotation = ha
)
hm
