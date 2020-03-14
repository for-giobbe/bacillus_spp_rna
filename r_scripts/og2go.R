if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("topGO")
library("topGO")

geneID2GO <- readMappings(file = "../../Downloads/3sp_OG2GO.tab")
geneUniverse <- names(geneID2GO)
genesOfInterest <- read.table("../../Downloads/interestinggenes.txt",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata
resultFisher <- runTest(myGOdata, algorithm="classic", statistic="fisher")
allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 30)


# print out the genes that are annotated with the significantly enriched GO terms:
myterms <- allRes$GO.ID
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms))
{
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  myfactor <- mygenesforterm %in% genesOfInterest # find the genes that are in the list of genes of interest
  mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
  mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
  print(paste("Term",myterm,"genes:",mygenesforterm2))
}