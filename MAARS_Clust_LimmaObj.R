############################################
## CLUSTERING USING THE LIMA OBJECT
############################################

### Define working directory
setwd("/Users/alain_lu/Dropbox/These_de_science/Data/Processed/Clustering/")


### Load data
load(file='/Users/alain_lu/Dropbox/These_de_science/Data/Processed/Lima/exset2.RData')
exset2 ## Limma Matrix which contains several features (expression, sample annotation ...)
ex1 <- exprs(exset2) ## extraction of expression matrix from exset2
class(ex1)
annot <- exset2@phenoData@data ## extraction of sample annotation from extset2 ## for merging an <-not$id == colnames(ex1) and colnames(ex1) == annot$id 
GS <- as.character(exset2@featureData@data$Gene.Symbol) ## gene list vector
fd <- exset2@featureData@data ## gene matrix with symbol, gene name, localisation


### Load library
anyLib::anyLib(c("ade4",
                 "circlize",
                 "ComplexHeatmap",
                 "dendextend",
                 "EMA",
                 "gdata",
                 "ggdendro",
                 "ggplot2",
                 "graphics",
                 "oligo",
                 "pd.hugene.2.1.st",
                 "pheatmap",
                 "ProjectTemplate",
                 "RColorBrewer",
                 "readxl", 
                 "reshape2"))

### Filter dataset
## keep only AD_NON_LES, AD_LES and CTRL

dim(annot)
annot <- annot[annot$status == "CTRL"|annot$status == "AD_NON_LES"|annot$status == "AD_LES", ]
dim(annot)

annot <- annot[order(rownames(annot)), ]
ex1 <- ex1[ ,order(colnames(ex1))]
ex1 <- ex1[ , colnames(ex1) %in% annot$sample_id]

## Discard probes with signal close to background level
ex1 <- EMA::expFilter(as.matrix(ex1), threshold = 5, p = 1/ncol(ex1)) # filter totally low expressed genes threshold < 5
density_ex1 <- density(as.matrix(ex1))
plot(density_ex1, main = "Distribution of normalized and filtered gene expression levels (> 5)", xlab = "log2 expression levels") 
+ abline(v = mean(as.matrix(ex1)), col = "red")
+ abline(v = median(as.matrix(ex1)), col = "blue") # Distribution of gene expression after filtering
dim(ex1)

##  Discard genes with low variability
non_const_ex1 <- EMA::genes.selection(ex1, thres.diff = 1.5, probs = 1/ncol(ex1))
density_ex1 <- density(as.matrix(ex1[non_const_ex1, ]))
plot(density_ex1, main = "Distribution of non-constant, normalized and filtered gene expression levels (> 5)", xlab = "log2 expression levels")
+ abline(v = mean(as.matrix(ex1[non_const_ex1, ])), col = "red")
+ abline(v = median(as.matrix(ex1[non_const_ex1, ])), col = "blue")
ex1 <- as.matrix(ex1[non_const_ex1,])
dim(ex1)

### Venn diagramm to determine the more stringent list of differentialy expressed genes between AD_LES, AD_NON_LES and CTRL



### Clustering usique selected genes by limma's comparision of ad_les and ad_non_les

## Using defaut function in R (by default  (disance : euclidean, hclust : complete)


col <- colorRampPalette(c("green", "black", "red"))(5)
heatmap(ex1_ad_clust, col=col, labRow = NA)

## Clustering whith complexe heatmap (https://jokergoo.github.io/ComplexHeatmap-reference/book/introduction.html)

col <- circlize::colorRamp2(c(0, 5.5, 11), c("green", "black", "red")) # Coloration function for Complexheatmap
essai <- ex1_ad_clust
essai <- merge(essai, fd, by = "row.names")

rownames(annot_ad)
colnames(essai) <- annot_ad$status
essai <- essai[, -c(166,167)]
essai <- essai[, -1]
# Clustering by raws and columns
hm <- Heatmap(essai,
        name = "Intensity", column_title = "Samples", 
        col=col,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = T,
        show_column_names = T,
        row_title = "Genes")
hm

#-------------------------------------------------------------------------------
# Classification par KMEANS
# pas de classification des colonnes mais seulement des lignes
#-------------------------------------------------------------------------------

# Kmeans + distance euclidiean

# 1) calcul de la matrice de distance
d = dist(t(ex1), method = 'euclidean') 
# 2) application de l'algorithme de regroupement
kgroups = kmeans(d, 3)

# 3) Représentation des résultats
Heatmap(ex1[names(sort(kgroups$cluster)), ],
        cluster_rows = F,
        cluster_columns = T,
        name = "Intensity", column_title = "Conditions", 
        col=col,
        show_row_names = F,
        show_column_names = F, 
        row_split = LETTERS[sort(kgroups$cluster)]
)

Heatmap(ex1,
        cluster_rows = F,
        cluster_columns = T,
        name = "Intensity", column_title = "Conditions", 
        col=col,
        show_row_names = F,
        show_column_names = F, 
)

# Kmeans + distance de corrélation

# 1) calcul de la matrice de distance
d = 1-cor(t(na.omit(expData)))
# 2) application de l'algorithme de regroupement
kgroups = kmeans(d, 10)

# 3) Représentation des résultats
Heatmap(na.omit(expData)[names(sort(kgroups$cluster)), ],
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        name = "Intensity", column_title = "Conditions", 
        col=col,
        show_row_names = FALSE,
        row_split = LETTERS[sort(kgroups$cluster)]
)


#Dendro
dt1 <- dist(t(ex1))
ct <- as.factor(cutree(hclust(dt1),k=3))

dd.row <- as.dendrogram(hclust(dt1))
ddata_x <- dendro_data(dd.row)
p1 <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
labs <- label(ddata_x)
m<-match(labs$label, rownames(annot))


ggd1 <- p1 + geom_text(data=label(ddata_x), size=2.5 ,
                       aes(label=label, x=x, y=-3,  angle = 0 ,hjust = 0 )) +
  coord_flip() +
  scale_y_reverse(expand=c(.1, 40)) + labs(color='Batches')
ggd1


c1 <- fcol(set = 'Set1' , annot = annot2 , ncol = 1)
c2 <- fcol(set = 'Set2' , annot = annot2 , ncol = 2)
c3 <- fcol(set = 'Set3' , annot = annot2 , ncol = 3)
ann_col <- as.list(c(c1,c2,c3))

ann_col
rn <- wh1

pheatmap( ex1,annotation_col = annot$status , annotation_colors = "red", cluster_cols = T, cluster_rows = T , gaps_col =  ,border_color = 0, filename = "/Users/alain_lu/Dropbox/These_de_science/Data/Processed/Lima/HM_1.pdf", fontsize = 10 , cutree_rows = 1, cutree_cols= 1 , scale = 'row'  ,cellwidth=1,cellheight=10,
          labels_row = as.character(paste(fd[rn,'Gene.Symbol'],' (',fd[rn,'Gene.Name'],')',sep='')) ,
          #labels_col = as.character( colnames(ex)) ,
          labels_col = rep('',ncol(ex1)) ,
          main = '')
browseURL("/Users/alain_lu/Dropbox/These_de_science/Data/Processed/Lima/HM_1.pdf")

pheatmap(mat = ex1, , cluster_rows = F, cluster_cols = T, scale = "col", border_color = 0, cutree_cols = 3,
         annotation_col = annot$sample_id, annotation_row = fd$GeneSymbol, 
         filename = "/Users/alain_lu/Dropbox/These_de_science/Data/Processed/Lima/HM_2.pdf")

browseURL("/Users/alain_lu/Dropbox/These_de_science/Data/Processed/Lima/HM_2.pdf")
