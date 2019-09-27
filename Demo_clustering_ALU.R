################################################################################
# Profils d’expression des gènes et méthodes de regroupement (clustering)
#-------------------------------------------------------------------------------
# Thomas DENECKER & Gaelle LELANDAIS
# Juin 2019
#
# OBJECTIFS 
# Cette démonstration a pour objectif de vous présenter les potentialité de R
# pour appliquer les méthodes de clustering et représenter les résultats.
#
# SOURCES DES DONNEES 
# Les données sont issus de la publication de Lelandais et al. (2016), cas d'étude
# O. tauri.
#
################################################################################

################################################################################
# Mise à jour R (si besoin) -> Ne marche pas sous mac
################################################################################

# installing/loading the package:
if(!require(installr)) {
  install.packages("installr"); require(installr)} #load / install+load installr

# using the package:
updateR()

################################################################################
# Installation des librairies nécessaires
################################################################################

# Cette librairie permet de vérifier la présence/absence des librairies et 
# de les installer si nécessaire.
install.packages("anyLib")
anyLib::anyLib(c("gridGraphics",
                 "grid",
                 "summarytools",
                 "ggplot2",
                 "reshape2",
                 "gplots", 
                 "RColorBrewer",
                 "ComplexHeatmap",
                 "igraph",
                 "RCy3", "RColorBrewer",
                 "ape",
                 "network",
                 "plotly", 
                 "circlize"
))


################################################################################
# Début de la démonstration !
################################################################################

#===============================================================================
# Partie 1 - Exploration des données
#===============================================================================

#-------------------------------------------------------------------------------
# Lecture des données
#-------------------------------------------------------------------------------

# Lecture des donnnées
expData = read.csv2("Dropbox/These_de_science/Data/Gaelle_script/clustering/GeneProfiles.txt", header = T, row.names = 1, sep = "\t",
                    stringsAsFactors = F)

# Création d'un tableau de chiffres
expData = data.matrix(expData)


# Visualisation des données (méthode summary tools) : dans RSTUDIO
view(dfSummary(expData))

# Enregistrement du fichier "bilan", au format HTML 
view(dfSummary(expData), file = "summaryFile.html") 

#-------------------------------------------------------------------------------
# QUESTION :: Combien de gènes différents sont présents dans le fichier ? 
#-------------------------------------------------------------------------------

# Méthode la plus prudente pour répondre à la question :
print(paste0("Il y a ",length(unique(rownames(expData)))," gènes différents."))

# Méthode la plus rapide pour répondre à la question :
print(paste0("Il y a ", nrow(expData)," gènes différents."))

#-------------------------------------------------------------------------------
# QUESTION :: Combien de points de mesures sont dans le fichier ? 
#-------------------------------------------------------------------------------

print(paste0("Il y a ", ncol(expData)," points de mesures dans le fichier."))

#-------------------------------------------------------------------------------
# REPRESENTATION GRAPHIQUE : FONCTIONS
# Nous allons représenter les profils d'expression de tous les gènes présents
# dans le fichier de données. 
#-------------------------------------------------------------------------------

plotGenesClassique <- function(expData, title = "", 
                               yMin = NULL, 
                               yMax = NULL, 
                               meanProfile = TRUE){
  
  # Vérification que les paramètres d'entrée de la fonction sont
  # sont bien complétés
  if(is.null(yMax) || is.null(yMin)){
    
    print("You must specify a min / max values for Y axis")
    
  }else{
    
    # Représentation du premier profil d'expression
    # Initialisation du graphique, spécification des axes et des titres
    plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
         ylim = c(floor(yMin), ceiling(yMax)),
         xlab = "Data point", ylab = "Gene expression level",
         main = title)
    
    # Ajout des autres profils
    for(i in 2:nrow(expData)){
      
      lines(1:ncol(expData), expData[i,], col = "grey")
      
      # fin du for()  
    }
    
    # Pour finir, Représentation du profil moyen
    if(meanProfile == TRUE){
      expMean = apply(expData, 2, mean, na.rm = T)
      lines(1:ncol(expData), expMean, col = "red", 
            lwd = 1.5, lty = "dashed")
    }
    
    # fin du else()   
  }
  
  # fin de la fonction plotGenes()  
}


plotGenesGgplot <- function(expData, title = "", 
                            yMin = NULL, 
                            yMax = NULL, 
                            meanProfile = TRUE){
  
  # Vérification que les paramètres d'entrée de la fonction sont
  # sont bien complétés
  if(is.null(yMax) || is.null(yMin)){
    
    print("You must specify a min / max values for Y axis")
    
  }else{
    
    # Le format est spécifique pour les graphes ggplot. La fonction melt permet
    # de la mettre dans le bon format
    
    expDataMelt <- melt(expData)  # convert to long format
    colnames(expDataMelt) =c("Genes", "Conditions", "Values")
    
    # Calcul de la moyenne
    meanData = apply(expData, 2, mean, na.rm=T)
    meanData = cbind.data.frame(Group = rep("Mean", 8), 
                                Conditions = names(meanData),
                                Values = meanData)
    
    # Création du graphe
    ggplot(data=expDataMelt) +
      geom_line(aes(x=Conditions, y=Values, group = Genes))  + 
      geom_line(data=meanData , mapping = aes(x= Conditions, y = Values, group = Group), color="red", size = 2) +
      ylab(label="Gene expression level") + 
      xlab("Data point") + 
      theme_classic() +
      ggtitle(title)
  }
  
}

plotGenesPlotly <- function(expData, title = "", 
                            yMin = NULL, 
                            yMax = NULL, 
                            meanProfile = TRUE){
  # Vérification que les paramètres d'entrée de la fonction sont
  # sont bien complétés
  if(is.null(yMax) || is.null(yMin)){
    
    print("You must specify a min / max values for Y axis")
    
  }else{
    
    # Calcul de la moyenne
    meanData = apply(expData, 2, mean, na.rm=T)
    meanData = cbind.data.frame(Group = rep("Mean", 8), 
                                Conditions = names(meanData),
                                Values = meanData)
    
    p <- plot_ly() %>%
      layout(title = title,
             xaxis = list(title = ""),
             yaxis = list (title = "Monthly Count of Products Sold"),
             showlegend = FALSE)
    
    ## Add the traces one at a time
    for(i in 1:100){
      p <- p %>% add_trace(x = paste("Condition", 1:8), y = expData[i,], name = rownames(expData)[i],
                           type = 'scatter',
                           mode = 'lines',
                           line = list(color = 'rgb(0, 0, 0)', width = 4))
    }
    
    p <- p %>% add_trace(x = paste("Condition", 1:8), y = meanData$Values, name = "Mean",
                         type = 'scatter',
                         mode = 'lines',
                         line = list(color = 'rgb(205, 12, 24)', width = 4))
    p }
}

#-------------------------------------------------------------------------------
# QUESTION :: En utilisant les fonctions précédentes, 
# représenter les profils d'expression des 100 premiers
# gènes présents dans le tableau de données "expData". 
#-------------------------------------------------------------------------------

# Utilisation de la fonction 1
plotGenesClassique(expData[1:100,], title = "100 genes - Classique", 
                   yMin = min(expData, na.rm = T), yMax = max(expData, na.rm = T))

# Utilisation de la fonction 2
plotGenesGgplot(expData[1:100,], title = "100 genes - ggplot2", 
                yMin = min(expData, na.rm = T), yMax = max(expData, na.rm = T))

# Méthode dynamique 
plotGenesPlotly(expData[1:100,], title = "100 genes", 
                yMin = min(expData, na.rm = T), yMax = max(expData, na.rm = T))


#-------------------------------------------------------------------------------
# QUESTION :: Représenter les profils d'expression tous les
# gènes présents dans le tableau de données "expData".
#-------------------------------------------------------------------------------

# Utilisation de la fonction 1
plotGenesClassique(expData, title = "all genes - Classique", 
                   yMin = min(expData, na.rm = T), yMax = max(expData, na.rm = T))

# Utilisation de la fonction 2
plotGenesGgplot(expData, title = "all genes - ggplot2", 
                yMin = min(expData, na.rm = T), yMax = max(expData, na.rm = T))

#===============================================================================
# Partie 2 - Visualiation des données par une heatmap
#===============================================================================

# Creation d'un vecteur de couleurs (pour la heatmap !)
# col <- colorRampPalette(c("green", "black", "red"))(5)
col <- c("green", "green", "green", "green", "forestgreen", "forestgreen",
         "black","black",
         "firebrick", "firebrick", "red", "red", "red", "red")

#-------------------------------------------------------------------------------
# Fonction disponible dans R
#-------------------------------------------------------------------------------

# Paramétrage par défaut (disance : euclidean, hclust : complete)
heatmap(na.omit(expData), col=col, labRow= NA)
# --> Le résultat est pas mal, il n'y a pas de légende, 
# et le paramétrage est difficile.

#-------------------------------------------------------------------------------
# AUTRE SOLUTION :
# ComplexHeatmap (2019)
# Documentation : https://jokergoo.github.io/ComplexHeatmap-reference/book/
#
# Avantage : une librairie très complète et très bien documentée
# Attention  : R doit être en 3.6 sinon des fonctionnalités ne marchent pas
# (--> MAJ de R au début du script !)
#-------------------------------------------------------------------------------

# Fonction de coloration pour la Heatmap du package ComplexHeatmap
col = colorRamp2(c(-3, 0, 3), c("green", "black", "red"))

# Première représentation, la plupart des paramètres sont conservés par défaut.
Heatmap(na.omit(expData),
        name = "Intensity", column_title = "Conditions", 
        col=col,
        show_row_names = FALSE,
        row_title = "Genes"
)


# Pour ce jeu de données la classification des colonnes n'a pas de sens,
# elle n'est plus faite avec la commande ci dessous :
Heatmap(na.omit(expData),
        name = "Intensity", column_title = "Conditions", 
        col=col,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        row_title = "Genes"
)

#-------------------------------------------------------------------------------
# Classification par KMEANS
# pas de classification des colonnes mais seulement des lignes
#-------------------------------------------------------------------------------

# Kmeans + distance euclidiean

# 1) calcul de la matrice de distance
d = dist(na.omit(expData), method = 'euclidean') 
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


#-------------------------------------------------------------------------------
# Classification par hclust (méthode HCL)
#-------------------------------------------------------------------------------

# Suppression des valeurs manquantes
expDataClean = na.omit(expData)

# hclust + euclidienne
Heatmap(expDataClean,
        cluster_rows = function(m) hclust(dist(m, method="euclidean"), method ="average"),
        cluster_columns = FALSE,
        name = "Intensity", column_title = "Conditions", 
        row_title = "Genes",
        col=col,
        show_row_names = FALSE)

# hclust + 1- cor (err)
Heatmap(expDataClean,
        cluster_rows = function(m) hclust(as.dist(1 - cor(t(m))), method ="average"),
        cluster_columns = FALSE,
        name = "Intensity", column_title = "Conditions", 
        row_title = "Genes",
        col=col,
        show_row_names = FALSE)

#-------------------------------------------------------------------------------
# Heatmap enrichie ...
#-------------------------------------------------------------------------------


# Ajout des fonctions des gènes
functionGene = read.csv2("Dropbox/These_de_science/Data/Gaelle_script/clustering/FunctionalAnnotation_BMCgenomics.csv", sep = "\t", row.names = 1)
completData = merge(na.omit(expData), functionGene, by = "row.names", all.x = T)
rownames(completData) = completData[,1]
completData = completData[,-1]

# Ajout des informations du fichier GFF
GFF = read.gff("Dropbox/These_de_science/Data/Gaelle_script/clustering/GCF_000214015.3_version_140606_genomic_DUO2.gff")
GFF$attributes = substring(GFF$attributes, 4)
rownames(GFF) = GFF$attributes
completData = merge(completData, GFF, by = "row.names", all.x = T)
rownames(completData) = completData[,1]
completData = completData[,-1]


# 1) calcul de la matrice de distance
d = dist(data.matrix(completData[,1:8]), method = 'euclidean') 
# 2) application de l'algorithme de regroupement
kgroups = kmeans(d, 10)

completDataSorted = completData[names(sort(kgroups$cluster)),]
completDataSorted = cbind.data.frame(completDataSorted,
                                     mean = apply(data.matrix(completDataSorted[,1:8]), 1, mean),
                                     cluster = sort(kgroups$cluster))



Heatmap(data.matrix(completDataSorted[,1:8]),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        name = "Intensity", column_title = "Conditions",row_title = "Genes",
        col=col,
        show_row_names = FALSE,
        row_split = LETTERS[sort(completDataSorted$cluster)],
        left_annotation = rowAnnotation(Groups = anno_block(gp = gpar(fill = brewer.pal(10, "Set3")),
                                                            labels = 1:10, 
                                                            labels_gp = gpar(col = "black", fontsize = 10)),
                                        Mean = anno_points(completDataSorted[,"mean", drop =F], 
                                                           pch = 16, size = unit(1, "mm"), 
                                                           axis_param = list(at = c(0, 2e5, 4e5, 6e5), 
                                                                             labels = c("0kb", "200kb", "400kb", "600kb")))
        ),
        bottom_annotation = HeatmapAnnotation(Values = anno_density(completDataSorted[,1:8],
                                                                    type = "violin", 
                                                                    height = unit(2, "cm"), 
                                                                    gp = gpar(fill = 1:10)))
)+
  Heatmap(completDataSorted$General.function..manual.inspection., 
          name = "Functions", col = rainbow(length(unique(completDataSorted$General.function..manual.inspection.))), 
          width = unit(5, "mm"))+
  Heatmap(completDataSorted$strand, 
          name = "Strand", col = rainbow(length(unique(completDataSorted$strand))), 
          width = unit(5, "mm"))+
  Heatmap(completDataSorted$seqid, 
          name = "Seq id", col = rainbow(length(unique(completDataSorted$seqid))), 
          width = unit(5, "mm"))



# Mais aussi avec plotly
plot_ly(
  x = colnames(completDataSorted)[1:8], y = rownames(completDataSorted),
  z = data.matrix(completDataSorted[,1:8]), type = "heatmap", 
  colors = colorRamp(c("green", "black", "red"))
)

#===============================================================================
# Partie 3
# Création d'un graphe de co-expression. Dans ce graphe, les gènes sont 
# représentés par des sommets et les arêtes relient deux gènes dont les profils 
# d'expression sont suffisament ressemblants.
#
#
# La ressemblance entre des profils d'expression peut se quantifier par le 
# calcul des distances entre les mesures d'expression de toutes les paires de 
# gènes (voir le cours). 
# Ces distances sont ensuite conservées dans une matice de distance.
# 
# Plusieurs types de distances existent : distance euclidienne, distance de 
# corrélation, etc. (voir la partie cours).
#===============================================================================

# La fonction pour calculer une matrice de distance avec R est "dist()". 
# Ci dessous un exemple de son utilisation :

distMat = as.matrix(dist(expData, method = "euclidean"))

#-------------------------------------------------------------------------------
# QUESTION :: Quels sont les deux gènes qui ont les profils d'expression sont les 
# plus différents ? Représenter les profils d'expression de ces deux gènes.
#-------------------------------------------------------------------------------


# Etape 1 : Rechercher la plus grande distance entre deux profils d'expression 
# de gènes 
# (REPONSE : 18.18849)

maxDist = max(distMat, na.rm = T)
maxDist
# 

# Etape 2 : Rechercher les gènes impliqué 
# (REPONSE : "ostta19g00150" et "ostta02g00080")

which(distMat == maxDist, arr.ind = TRUE)
distMat["ostta19g00150", "ostta02g00080"] == maxDist

# Etape 3 : Représentation des 2 profils avec la fonction plotGenesClassique
plotGenesClassique(expData[c("ostta19g00150", "ostta02g00080"),], 
                   title = "Genes with the highest distance between expression profiles", 
                   yMax = max(expData, na.rm = T), yMin = min(expData, na.rm = T), meanProfile = F)


#-------------------------------------------------------------------------------
# QUESTION :: Quels sont les deux gènes qui ont les profils d'expression sont les 
# plus proches ? Représenter les profils d'expression de ces deux gènes.
#-------------------------------------------------------------------------------

# Nous allons maintenant appliquer la même démarche pour la plus petite valeur 
# de la distance entre deux profils d'expression de gènes.

# Il est important de faire attention à ne pas tenir compte des distances écrites 
# sur la diagonale de la matrice de distance (distance un profil d'expression d'un gène 
# et lui même).

# Ainsi, les valeurs de la matrice sur la diagonale sont remplacées par NA
diag(distMat) = NA

# Pour voir le changement
distMat[1:5, 1:5]

# Recherche de la plus petite distance
minDist = min(distMat, na.rm = T)
minDist

# Les gènes dont les profils sont les plus ressemblants sont ensuite recherchés, et représentés
# graphiquement.
res = which(distMat == minDist, arr.ind = TRUE)
plotGenesClassique(expData[row.names(res),], 
                   title = "Genes with the lowest distance between expression profiles", 
                   yMax = max(expData, na.rm = T), yMin = min(expData, na.rm = T), meanProfile = F)

#-------------------------------------------------------------------------------
# QUESTION :: Comment les distances sont-elles réparties?
#-------------------------------------------------------------------------------

# Représentation histogramme de l'ensemble des distances présentes dans la
# matrice de distance. La fonction hist() eut être utilisée.
hist(as.dist(distMat), nclass = 50, col = "orange",
     main = "Distances between gene expression profiles",
     xlab = "distance value")
# Notez que la conversion as.dist() est importante, pour ne pas compter les distances deux
# fois (la distance en les gènes g1 et g2 est égale à la distance ente les gènes g2 et g1).


#-------------------------------------------------------------------------------
# QUESTION :: Déterminer la valeur de la distance pour laquelle 1% des distances
# sont en dessous de cette valeur. Cette valeur nous servira de seuil pour représenter
# un premier graphe de co-expression. La fonction "quantile()" peut être utile. 
#-------------------------------------------------------------------------------

#
# Ecrivez votre code ici
#
# (REPONSE : 1.019148)
Tval = quantile(as.dist(distMat), probs = 0.01, na.rm = T)
Tval


#-------------------------------------------------------------------------------
# QUESTION :: Selectionner les paires de gènes avec une distance sous ce seuil 
#-------------------------------------------------------------------------------

edgeList = NULL

for(i in 1:(nrow(distMat) - 1)){
  
  for(j in (i+1):ncol(distMat)){
    
    if(!is.na(distMat[i,j])){
      
      if(distMat[i,j] < Tval){
        
        edgeList = rbind(edgeList, c(row.names(distMat)[i],colnames(distMat)[j]))
        
      }
      
    }
    
    # fin du for()  
  }
  # fin du for()  
}

colnames(edgeList) = c("source", "target")

#-------------------------------------------------------------------------------
# QUESTION :: Combien de paires de gènes sont ainsi sélectionnées ?
#-------------------------------------------------------------------------------

# (REPONSE : 5487)
nrow(edgeList)

#-------------------------------------------------------------------------------
# Création du graphe de co-expression
#
# Nous allons réer le graphe de co-expression en utilisant la librairie R "igraph".
#-------------------------------------------------------------------------------


# Création du graphe (à partir de la liste des arêtes)
tauriGraph = graph_from_edgelist(edgeList, directed = F)

# Représentation du graphe co-expression, et des profils d'expression des gènes 
# correspondants
par(mfrow = c(1,2))
plot.igraph(tauriGraph, vertex.color = "red", vertex.size = 5, vertex.label = NA,
            arrow.size = 0.1, edge.color = "grey", main = "Co-expression gene network")
plotGenesClassique(expData[unique(c(edgeList[,1], edgeList[,2])),], 
                   title = "Gene expression profiles", 
                   yMax = max(expData, na.rm = T), yMin = min(expData, na.rm = T))

#-------------------------------------------------------------------------------
# Création du graphe de co-expression avec cytoscape
# Documentation : https://bioconductor.org/packages/release/bioc/vignettes/RCy3/inst/doc/Overview-of-RCy3.html
#-------------------------------------------------------------------------------

# /!\ /!\ Il faut lancer cytoscape /!\ /!\


# tester si la connection est bonne avec cytoscape 
cytoscapePing ()

# Création des informations pour le graphe
nodes <- cbind.data.frame( unique(c(edgeList[,1], edgeList[,2])))
row.names(nodes) = nodes[,1]
nodes <- merge(nodes, functionGene, by = "row.names")
nodes = nodes[,-1]
colnames(nodes)[1]= "id"

edges <- as.data.frame(edgeList)

createNetworkFromDataFrames(nodes,edges, title="Co-expression gene network", collection="DataFrame Example")
style.name = "MyStyle"
defaults.list <- list(NODE_SHAPE="ellipse",
                      EDGE_TRANSPARENCY=120)
createVisualStyle(style.name, defaults.list)
setNodeColorBypass(nodes$id,col2hex(as.color(nodes$General.function..manual.inspection.)))
setVisualStyle(style.name=style.name)

# Note 
# En cliquant sur les sommets dans cytoscape, nous avons de l'information qui 
# provient du fichier de fonctions