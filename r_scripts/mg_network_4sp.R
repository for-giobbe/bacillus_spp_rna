library(WGCNA)
library(data.table)
library(tidyverse)
options(stringsAsFactors = FALSE);

############################################################################################################################################  expr. data ########################################

femData = read.csv("./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/4sp/BGM_RSEM_mf.TMM.EXPR.matrix.4sp.reformat.ogonly", sep =" ", header = TRUE)
dim(femData);
names(femData);

datExpr0 = as.data.frame(t(femData[, -c(1)]));
names(datExpr0) = femData$transcript;
rownames(datExpr0) = names(femData)[-c(1)];

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

############################################################################################################################################  trait data ########################################

traitData = read.csv("./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/mg_description.csv");
dim(traitData);
names(traitData);

traitData = traitData[, c(2,3:4) ];
dim(traitData);
names(traitData);

Samples = rownames(datExpr0);
traitRows = match(Samples, traitData$sample);
datTraits = traitData[traitRows, -1];
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage();

# Plot samples and traits tree
sampleTree2 = hclust(dist(datExpr0), method = "complete");
traitColors = numbers2colors(datTraits, signed = TRUE);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap");

############################################################################################################################################  network ########################################

allowWGCNAThreads() 
powers = c(c(1:10), seq(from = 5, to=25, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5,networkType = "signed")

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(datExpr0, power = 23,
                       TOMType = "signed", networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.01, corOptions = list(use = 'p', maxPOutliers = 0.1),
                       numericLabels = TRUE, pamRespectsDendro = FALSE, replaceMissingAdjacencies = TRUE,
                       saveTOMs = TRUE, saveTOMFileBase = "test",
                       verbose = 3, maxBlockSize = 6000, corType = "bicor")

# Plot the dendrogram and the module colors underneath
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

# Recalculate MEs with color labels
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "all.obs", method = c("pearson"));
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitPvalue  <- p.adjust(moduleTraitPvalue, method = "hochberg", n = length(moduleTraitPvalue))


############################################################################################################################################  heatmap! #####################

# Display the correlation values within a heatmap plot
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

### positive correlation
names(datExpr0)[moduleColors=="blue"]
names(datExpr0)[moduleColors=="black"]
names(datExpr0)[moduleColors=="turquoise"]
names(datExpr0)[moduleColors=="red"]

############################################################################################################################################  male gonad modules #############################

# Define variable containing the variable column of datTrait
mg = as.data.frame(datTraits$condition);
names(mg) = "mg"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "all.obs", method = c("spearman")));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, mg, use = "all.obs", method = c("spearman")));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(mg), sep="");
names(GSPvalue) = paste("p.GS.", names(mg), sep="");

module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for male gonads",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

############################################################################################################################################  Chooses the top hub gene in each module  #######

chooseTopHubInEachModule(
  datExpr0, 
  moduleColors, type = "signed",power =19)

############################################################################################################################################  module preservation BRO ######

BROData = read.csv("./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/4sp/BRO_RSEM.TMM.EXPR.matrix.4sp.reformat.og_gonad_only", sep =" ", header = TRUE)
datExprBRO= as.data.frame(t(BROData[, -c(1)]));
names(datExprBRO) = BROData$OG;
rownames(datExprBRO) = names(BROData)[-c(1)];
gsg_BRO = goodSamplesGenes(datExprBRO, verbose = 3);

if (!gsg_BRO$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg_BRO$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExprBRO)[!gsg_BRO$goodGenes], collapse = ", ")));
  if (sum(!gsg_BRO$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExprBRO)[!gsg_BRO$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExprBRO = datExprBRO[gsg_BRO$goodSamples, gsg_BRO$goodGenes]
}

colorsBGM = moduleColors

setLabels = c("BGM", "BRO");
multiExprBRO = list(BGM= list(data = datExpr0), BRO = list(data = datExprBRO));
multiColor = list(BGM = colorsBGM);

system.time( {
  mp_bro = modulePreservation(multiExprBRO, multiColor, networkType = "signed", maxModuleSize = 6000, #corOptions = "use = 'p', method = 'spearman'",
                              referenceNetworks = 1,
                              nPermutations = 10,
                              randomSeed = 1,
                              quickCor = 1,
                              verbose = 3)
} );

ref = 1
test = 2
statsObs = cbind(mp_bro$quality$observed[[ref]][[test]][, -1], mp_bro$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp_bro$quality$Z[[ref]][[test]][, -1], mp_bro$preservation$Z[[ref]][[test]][, -1]);

############################################################################################################################################  preservation network connectivity BRO TO DO ##################

preservationNetworkConnectivity_BRO <- preservationNetworkConnectivity(
  multiExprBAT,
  useSets = NULL, useGenes = NULL,
  corFnc = "cor", corOptions = "use='p'",
  networkType = "signed",
  power = 19,
  sampleLinks = NULL, nLinks = 5000,
  blockSize = 1000,
  setSeed = 12345,
  weightPower = 2,
  verbose = 2, indent = 0)



############################################################################################################################################  preservation stats BRO ##################

bro_preservation_stats <- print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

# Module labels and module sizes are also contained in the results
modColors = rownames(mp_bro$preservation$observed[[ref]][[test]])
moduleSizes = mp_bro$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp_bro$preservation$observed[[ref]][[test]][, 2], mp_bro$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Plots/BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}

############################################################################################################################################  module preservation BAT ######

BATData = read.csv("./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/4sp/BAT_RSEM.TMM.EXPR.matrix.4sp.reformat.og_gonad_only", sep =" ", header = TRUE)
datExprBAT= as.data.frame(t(BATData[, -c(1)]));
names(datExprBAT) = BATData$OG;
rownames(datExprBAT) = names(BATData)[-c(1)];
gsg_BAT = goodSamplesGenes(datExprBAT, verbose = 3);
gsg$allOK


if (!gsg_BAT$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg_BAT$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExprBAT)[!gsg_BAT$goodGenes], collapse = ", ")));
  if (sum(!gsg_BRO$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExprBAT)[!gsg_BAT$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExprBAT = datExprBRO[gsg_BRO$goodSamples, gsg_BRO$goodGenes]
}

colorsBGM = moduleColors

setLabels = c("BGM", "BAT");
multiExprBAT = list(BGM= list(data = datExpr0), BAT = list(data = datExprBAT));
multiColor = list(BGM = colorsBGM);


system.time( {
  mp_bat = modulePreservation(multiExprBAT, multiColor, networkType = "signed", maxModuleSize = 6000, #corOptions = "use = 'p', method = 'spearman'",
                              referenceNetworks = 1,
                              nPermutations = 10,
                              randomSeed = 1,
                              quickCor = 0,
                              verbose = 3)
} );


ref = 1
test = 2
statsObs = cbind(mp_bat$quality$observed[[ref]][[test]][, -1], mp_bat$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp_bat$quality$Z[[ref]][[test]][, -1], mp_bat$preservation$Z[[ref]][[test]][, -1]);
#We look at the main output: the preservation medianRank and Zsummary statistics.

############################################################################################################################################  preservation network connectivity BAT TO DO ##################

preservationNetworkConnectivity_BAT <- preservationNetworkConnectivity(
  multiExprBAT,
  useSets = NULL, useGenes = NULL,
  corFnc = "cor", corOptions = "use='p'",
  networkType = "signed",
  power = 19,
  sampleLinks = NULL, nLinks = 5000,
  blockSize = 1000,
  setSeed = 12345,
  weightPower = 2,
  verbose =TO  2, indent = 0)

############################################################################################################################################  preservation stats BAT ##################

bat_preservation_stats <- print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

# Module labels and module sizes are also contained in the results
modColors = rownames(mp_bat$preservation$observed[[ref]][[test]])
moduleSizes = mp_bat$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp_bat$preservation$observed[[ref]][[test]][, 2], mp_bat$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Plots/BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}

############################################################################################################################################  gene modules ##################

blue_genes <- names(datExpr0)[moduleColors=="blue"]
red_genes <- names(datExpr0)[moduleColors=="red"]
black_genes <- names(datExpr0)[moduleColors=="black"]
turquoise_genes <- names(datExpr0)[moduleColors=="turquoise"]

############################################################################################################################################  n connectivity calculation ##################

nncon_bgm <- nearestNeighborConnectivity(datExpr0, nNeighbors = 1000, power = 19, type = "signed", corFnc = "bicor", corOptions = "use = 'p'", blockSize = 1000)
names_bgm <- as.vector(names(datExpr0))
nncon_bgm_name <- "con"
names_bgm_name <- "og"
nncon_bgm_df <- data.frame(names_bgm,nncon_bgm)
colnames(nncon_bgm_df) <- c(names_bgm_name, nncon_bgm_name)

nncon_bro <- nearestNeighborConnectivity(datExprBRO, nNeighbors = 1000, power = 19, type = "signed", corFnc = "bicor", corOptions = "use = 'p'", blockSize = 1000)
names_bro <- as.vector(names(datExprBRO))
nncon_bro_name <- "con"
names_bro_name <- "og"
nncon_bro_df <- data.frame(names_bro,nncon_bro)
colnames(nncon_bro_df) <- c(names_bro_name, nncon_bro_name)

nncon_bat <- nearestNeighborConnectivity(datExprBAT, nNeighbors = 1000, power = 19, type = "signed", corFnc = "bicor", corOptions = "use = 'p'", blockSize = 1000)
names_bat <- as.vector(names(datExprBAT))
nncon_bat_name <- "con"
names_bat_name <- "og"
nncon_bat_df <- data.frame(names_bat,nncon_bat)
colnames(nncon_bat_df) <- c(names_bat_name, nncon_bat_name)

############################################################################################################################################  k connectivity calculation TO DO ############################

adj <- adjacency(datExpr0, selectCols = NULL,  type = "signed", 
                 power = 19,
                 corFnc = "cor", corOptions = "use = 'p', method = 'spearman'",
                 distFnc = "dist", distOptions = "method = 'euclidean'")

intramodularConnectivity <- intramodularConnectivity(adj, mergedColors, scaleByMax = FALSE)

fundamentalNetworkConcepts(adj, GS = NULL)

write.table(intramodularConnectivity, file = "./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/4sp/intramodularConnectivity_4sp.tab" , sep = " ", row.names = TRUE, col.names = TRUE)
pwd

############################################################################################################################################  n connectivity sex VS asex per-module ##################

nnconectivity_bgm_bro <- merge(nncon_bgm_df, nncon_bro_df[, c("og", "con")], by="og")

turquoise_nnconnectivity_bgm_bro <- subset(nnconectivity_bgm_bro, nnconectivity_bgm_bro$og %in% turquoise_genes)
turquoise_nnconnectivity_bgm_bro.spearman <- cor.test(turquoise_nnconnectivity_bgm_bro$con.x , turquoise_nnconnectivity_bgm_bro$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
turquoise_nnconnectivity_bgm_bro_plot <- ggplot(turquoise_nnconnectivity_bgm_bro, aes(x=con.x, y=con.y)) + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "nn connectivity  Bacillus rossius (asex)",  title = paste0("turquoise     pval: ", round(turquoise_nnconnectivity_bgm_bro.spearman$p.value,digits = 3), "     rho: ", round(turquoise_nnconnectivity_bgm_bro.spearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"), title=element_text(size=12,face="bold"))
blue_nnconnectivity_bgm_bro <- subset(nnconectivity_bgm_bro, nnconectivity_bgm_bro$og %in% blue_genes)
blue_nnconnectivity_bgm_bro.spearman <- cor.test(blue_nnconnectivity_bgm_bro$con.x , blue_nnconnectivity_bgm_bro$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
blue_nnconnectivity_bgm_bro_plot <- ggplot(blue_nnconnectivity_bgm_bro, aes(x=con.x, y=con.y)) + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "nn connectivity Bacillus rossius (asex)",  title = paste0("blue     pval: ", round(blue_nnconnectivity_bgm_bro.spearman$p.value,digits = 3), "     rho: ", round(blue_nnconnectivity_bgm_bro.spearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"), title=element_text(size=12,face="bold"))
black_nnconnectivity_bgm_bro <- subset(nnconectivity_bgm_bro, nnconectivity_bgm_bro$og %in% black_genes)
black_nnconnectivity_bgm_bro.spearman <- cor.test(black_nnconnectivity_bgm_bro$con.x , black_nnconnectivity_bgm_bro$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
black_nnconnectivity_bgm_bro_plot <- ggplot(black_nnconnectivity_bgm_bro, aes(x=con.x, y=con.y)) + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "nn connectivity Bacillus rossius (asex)",  title = paste0("black     pval: ", round(black_nnconnectivity_bgm_bro.spearman$p.value,digits = 3), "     rho: ", round(black_nnconnectivity_bgm_bro.spearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"), title=element_text(size=12,face="bold"))
red_nnconnectivity_bgm_bro <- subset(nnconectivity_bgm_bro, nnconectivity_bgm_bro$og %in% red_genes)
red_nnconnectivity_bgm_bro.spearman <- cor.test(red_nnconnectivity_bgm_bro$con.x , red_nnconnectivity_bgm_bro$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
red_nnconnectivity_bgm_bro_plot <- ggplot(red_nnconnectivity_bgm_bro, aes(x=con.x, y=con.y)) + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "nn connectivity Bacillus rossius (asex)",  title = paste0("red     pval: ", round(red_nnconnectivity_bgm_bro.spearman$p.value,digits = 3), "     rho: ", round(red_nnconnectivity_bgm_bro.spearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"), title=element_text(size=12,face="bold"))
nnconnectivity_bgm_bro.spearman <- cor.test(nnconectivity_bgm_bro$con.x , nnconectivity_bgm_bro$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
nnconnectivity_bgm_bro_plot <- ggplot(nnconectivity_bgm_bro, aes(x=con.x, y=con.y)) + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "nn connectivity Bacillus rossius (asex)",  title = paste0("total     pval: ", round(nnconnectivity_bgm_bro.spearman$p.value,digits = 3), "     rho: ", round(nnconnectivity_bgm_bro.spearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"), title=element_text(size=12,face="bold"))

nnconectivity_bgm_bat <- merge(nncon_bgm_df, nncon_bat_df[, c("og", "con")], by="og")

turquoise_nnconnectivity_bgm_bat <- subset(nnconectivity_bgm_bat, nnconectivity_bgm_bat$og %in% turquoise_genes)
turquoise_nnconnectivity_bgm_bat.spearman <- cor.test(turquoise_nnconnectivity_bgm_bat$con.x , turquoise_nnconnectivity_bgm_bat$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
turquoise_nnconnectivity_bgm_bat_plot <- ggplot(turquoise_nnconnectivity_bgm_bat, aes(x=con.x, y=con.y)) + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "nn connectivity  Bacillus rossius (asex)",  title = paste0("turquoise     pval: ", round(turquoise_nnconnectivity_bgm_bat.spearman$p.value,digits = 3), "     rho: ", round(turquoise_nnconnectivity_bgm_bat.spearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"), title=element_text(size=12,face="bold"))
blue_nnconnectivity_bgm_bat <- subset(nnconectivity_bgm_bat, nnconectivity_bgm_bat$og %in% blue_genes)
blue_nnconnectivity_bgm_bat.spearman <- cor.test(blue_nnconnectivity_bgm_bat$con.x , blue_nnconnectivity_bgm_bat$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
blue_nnconnectivity_bgm_bat_plot <- ggplot(blue_nnconnectivity_bgm_bat, aes(x=con.x, y=con.y)) + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "nn connectivity Bacillus rossius (asex)",  title = paste0("blue     pval: ", round(blue_nnconnectivity_bgm_bat.spearman$p.value,digits = 3), "     rho: ", round(blue_nnconnectivity_bgm_bat.spearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"), title=element_text(size=12,face="bold"))
black_nnconnectivity_bgm_bat <- subset(nnconectivity_bgm_bat, nnconectivity_bgm_bat$og %in% black_genes)
black_nnconnectivity_bgm_bat.spearman <- cor.test(black_nnconnectivity_bgm_bat$con.x , black_nnconnectivity_bgm_bat$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
black_nnconnectivity_bgm_bat_plot <- ggplot(black_nnconnectivity_bgm_bat, aes(x=con.x, y=con.y)) + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "nn connectivity Bacillus rossius (asex)",  title = paste0("black     pval: ", round(black_nnconnectivity_bgm_bat.spearman$p.value,digits = 3), "     rho: ", round(black_nnconnectivity_bgm_bat.spearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"), title=element_text(size=12,face="bold"))
red_nnconnectivity_bgm_bat <- subset(nnconectivity_bgm_bat, nnconectivity_bgm_bat$og %in% red_genes)
red_nnconnectivity_bgm_bat.spearman <- cor.test(red_nnconnectivity_bgm_bat$con.x , red_nnconnectivity_bgm_bat$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
red_nnconnectivity_bgm_bat_plot <- ggplot(red_nnconnectivity_bgm_bat, aes(x=con.x, y=con.y)) + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "nn connectivity Bacillus rossius (asex)",  title = paste0("red     pval: ", round(red_nnconnectivity_bgm_bat.spearman$p.value,digits = 3), "     rho: ", round(red_nnconnectivity_bgm_bat.spearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"), title=element_text(size=12,face="bold"))
nnconnectivity_bgm_bat.spearman <- cor.test(nnconectivity_bgm_bat$con.x , nnconectivity_bgm_bat$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
nnconnectivity_bgm_bat_plot <- ggplot(nnconectivity_bgm_bat, aes(x=con.x, y=con.y)) + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "nn connectivity Bacillus rossius (asex)",  title = paste0("total     pval: ", round(nnconnectivity_bgm_bat.spearman$p.value,digits = 3), "     rho: ", round(nnconnectivity_bgm_bat.spearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"), title=element_text(size=12,face="bold"))

grid.arrange(turquoise_nnconnectivity_bgm_bro_plot, blue_nnconnectivity_bgm_bro_plot, black_nnconnectivity_bgm_bro_plot, red_nnconnectivity_bgm_bro_plot, nnconnectivity_bgm_bro_plot, 
             turquoise_nnconnectivity_bgm_bat_plot, blue_nnconnectivity_bgm_bat_plot, black_nnconnectivity_bgm_bat_plot, red_nnconnectivity_bgm_bat_plot, nnconnectivity_bgm_bat_plot,
             ncol = 5 , nrow = 2)       

############################################################################################################################################  n connectivity VS LogFC per module ############################

formatted_con_dnds_BAT <- read.table(file = "./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/4sp/connectivicty_BATF_dnds_LogFC.4sp.tab", sep = " ", header= TRUE)
redtable_BAT <- subset(formatted_con_dnds_BAT, formatted_con_dnds_BAT$OG %in% red_genes)
bluetable_BAT <- subset(formatted_con_dnds_BAT, formatted_con_dnds_BAT$OG %in% blue_genes)
blacktable_BAT <- subset(formatted_con_dnds_BAT, formatted_con_dnds_BAT$OG %in% black_genes)
turquoisetable_BAT <- subset(formatted_con_dnds_BAT, formatted_con_dnds_BAT$OG %in% turquoise_genes)

BAT_redspearman <- cor.test(redtable_BAT$kWithin , redtable_BAT$BAT_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
BAT_red <- ggplot(redtable_BAT, aes(x=kWithin, y=BAT_logFC))    + geom_point(shape=16, color="black", size = 1)  + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "LogFC Bacillus atticus (asex)",  title = paste0("red     pval: ", round(BAT_redspearman$p.value,digits = 3), "     rho: ", round(BAT_redspearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))
BAT_bluespearman <- cor.test(bluetable_BAT$kWithin , bluetable_BAT$BAT_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
BAT_blue <- ggplot(bluetable_BAT, aes(x=kWithin, y=BAT_logFC))    + geom_point(shape=16, color="black", size = 1)  + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "LogFC Bacillus atticus (asex)",  title = paste0("blue     pval: ", round(BAT_bluespearman$p.value,digits = 3), "     rho: ", round(BAT_bluespearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))
BAT_blackspearman <- cor.test(blacktable_BAT$kWithin , blacktable_BAT$BAT_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
BAT_black <- ggplot(blacktable_BAT, aes(x=kWithin, y=BAT_logFC))    + geom_point(shape=16, color="black", size = 1)  + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "LogFC Bacillus atticus (asex)",  title = paste0("black     pval: ", round(BAT_blackspearman$p.value,digits = 3), "     rho: ", round(BAT_blackspearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))
BAT_turquoisespearman <- cor.test(turquoisetable_BAT$kWithin , turquoisetable_BAT$BAT_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
BAT_turquoise <- ggplot(turquoisetable_BAT, aes(x=kWithin, y=BAT_logFC))    + geom_point(shape=16, color="black", size = 1)  + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "LogFC Bacillus atticus (asex)",  title = paste0("turquoise     pval: ", round(BAT_turquoisespearman$p.value,digits = 3), "     rho: ", round(BAT_turquoisespearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))

formatted_con_dnds_BRO <- read.table(file = "./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/4sp/connectivicty_BROF_dnds_LogFC.4sp.tab", sep = " ", header= TRUE)
redtable_BRO <- subset(formatted_con_dnds_BRO, formatted_con_dnds_BRO$OG %in% red_genes)
bluetable_BRO <- subset(formatted_con_dnds_BRO, formatted_con_dnds_BRO$OG %in% blue_genes)
blacktable_BRO <- subset(formatted_con_dnds_BRO, formatted_con_dnds_BRO$OG %in% black_genes)
turquoisetable_BRO <- subset(formatted_con_dnds_BRO, formatted_con_dnds_BRO$OG %in% turquoise_genes)

BRO_redspearman <- cor.test(redtable_BRO$kWithin , redtable_BRO$BRO_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
BRO_red <- ggplot(redtable_BRO, aes(x=kWithin, y=BRO_logFC))    + geom_point(shape=16, color="black", size = 1)  + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "LogFC Bacillus atticus (asex)",  title = paste0("red     pval: ", round(BRO_redspearman$p.value,digits = 3), "     rho: ", round(BRO_redspearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))
BRO_bluespearman <- cor.test(bluetable_BRO$kWithin , bluetable_BRO$BRO_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
BRO_blue <- ggplot(bluetable_BRO, aes(x=kWithin, y=BRO_logFC))    + geom_point(shape=16, color="black", size = 1)  + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "LogFC Bacillus atticus (asex)",  title = paste0("blue     pval: ", round(BRO_bluespearman$p.value,digits = 3), "     rho: ", round(BRO_bluespearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))
BRO_blackspearman <- cor.test(blacktable_BRO$kWithin , blacktable_BRO$BRO_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
BRO_black <- ggplot(blacktable_BRO, aes(x=kWithin, y=BRO_logFC))    + geom_point(shape=16, color="black", size = 1)  + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "LogFC Bacillus atticus (asex)",  title = paste0("black     pval: ", round(BRO_blackspearman$p.value,digits = 3), "     rho: ", round(BRO_blackspearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))
BRO_turquoisespearman <- cor.test(turquoisetable_BRO$kWithin , turquoisetable_BRO$BRO_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
BRO_turquoise <- ggplot(turquoisetable_BRO, aes(x=kWithin, y=BRO_logFC))    + geom_point(shape=16, color="black", size = 1)  + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "LogFC Bacillus atticus (asex)",  title = paste0("turquoise     pval: ", round(BRO_turquoisespearman$p.value,digits = 3), "     rho: ", round(BRO_turquoisespearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))

grid.arrange(BAT_red, BAT_blue, BAT_black, BAT_turquoise, BRO_red, BRO_blue, BRO_black, BRO_turquoise, ncol = 4 , nrow = 2)

############################################################################################################################################  n connectivity sex VS asex - module VS non-module main fig ############################

# bgm VS bro

nnconectivity_bgm_bro <- merge(nncon_bgm_df, nncon_bro_df[, c("og", "con")], by="og")
nnconectivity_bgm_bro_network <- subset(nnconectivity_bgm_bro, nnconectivity_bgm_bro$og %in% network_genes)
nnconectivity_bgm_bro_nonnetwork <- subset(nnconectivity_bgm_bro, !(nnconectivity_bgm_bro$og %in% network_genes))

nnconectivity_bgm_bro_network_spearman <- cor.test(nnconectivity_bgm_bro_network$con.x , nnconectivity_bgm_bro_network$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
nnconectivity_bgm_bro_network_plot <- ggplot() +
  geom_point(data=nnconectivity_bgm_bro_network, aes(con.x, con.y), color="black", size = 2, alpha = 0.5, shape= 20) +
  labs(x="nn connectivity Bacillus grandii  (sex)", y = "nn connectivity Bacillus rossius (asex)", title = paste0(
    "network   pval: ", round(nnconectivity_bgm_bro_network_spearman$p.value,digits = 3), "     rho: ", round(nnconectivity_bgm_bro_network_spearman$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
  
nnconectivity_bgm_bro_nonnetwork_spearman <- cor.test(nnconectivity_bgm_bro_nonnetwork$con.x , nnconectivity_bgm_bro_nonnetwork$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
nnconectivity_bgm_bro_nonnetwork_plot <- ggplot() +
  geom_point(data=nnconectivity_bgm_bro_nonnetwork, aes(con.x, con.y), color="black", size = 2, alpha = 0.5, shape= 20) +
  labs(x="nn connectivity Bacillus grandii  (sex)", y = "nn connectivity Bacillus rossius (asex)", title = paste0(
    "non-network   pval: ", round(nnconectivity_bgm_bro_nonnetwork_spearman$p.value,digits = 3), "     rho: ", round(nnconectivity_bgm_bro_nonnetwork_spearman$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

# bgm VS bat

nnconectivity_bgm_bat <- merge(nncon_bgm_df, nncon_bat_df[, c("og", "con")], by="og")
nnconectivity_bgm_bat_network <- subset(nnconectivity_bgm_bat, nnconectivity_bgm_bat$og %in% network_genes)
nnconectivity_bgm_bat_nonnetwork <- subset(nnconectivity_bgm_bat, !(nnconectivity_bgm_bat$og %in% network_genes))

nnconectivity_bgm_bat_network_spearman <- cor.test(nnconectivity_bgm_bat_network$con.x , nnconectivity_bgm_bat_network$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
nnconectivity_bgm_bat_network_plot <- ggplot() +
  geom_point(data=nnconectivity_bgm_bat_network, aes(con.x, con.y), color="black", size = 2, alpha = 0.5, shape= 20) +
  labs(x="nn connectivity Bacillus grandii  (sex)", y = "nn connectivity Bacillus atticus (asex)", title = paste0(
    "network   pval: ", round(nnconectivity_bgm_bat_network_spearman$p.value,digits = 3), "     rho: ", round(nnconectivity_bgm_bat_network_spearman$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

nnconectivity_bgm_bat_nonnetwork_spearman <- cor.test(nnconectivity_bgm_bat_nonnetwork$con.x , nnconectivity_bgm_bat_nonnetwork$con.y, method = "spearman", continuity = FALSE, conf.level = 0.95)
nnconectivity_bgm_bat_nonnetwork_plot <- ggplot() +
  geom_point(data=nnconectivity_bgm_bat_nonnetwork, aes(con.x, con.y), color="black", size = 2, alpha = 0.5, shape= 20) +
  labs(x="nn connectivity Bacillus grandii  (sex)", y = "nn connectivity Bacillus atticus (asex)", title = paste0(
    "non-network   pval: ", round(nnconectivity_bgm_bat_nonnetwork_spearman$p.value,digits = 3), "     rho: ", round(nnconectivity_bgm_bat_nonnetwork_spearman$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

grid.arrange(nnconectivity_bgm_bro_network_plot, nnconectivity_bgm_bro_nonnetwork_plot, nnconectivity_bgm_bat_network_plot, nnconectivity_bgm_bat_nonnetwork_plot, ncol = 2 , nrow = 2)

############################################################################################################################################  k connectivity VS t ############################

redtable_t_BAT <- subset(redtable_BAT, t>0.01 & t<1)
bluetable_t_BAT <- subset(bluetable_BAT, t>0.01 & t<1)
blacktable_t_BAT <- subset(blacktable_BAT, t>0.01 & t<1)
turquoisetable_t_BAT <- subset(turquoisetable_BAT, t>0.01 & t<1)

BAT_redspearman_t <- cor.test(redtable_t_BAT$kWithin , redtable_t_BAT$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
BAT_red_t <- ggplot(redtable_t_BAT, aes(x=kWithin, y=t))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("red     pval: ", round(BAT_redspearman$p.value,digits = 3), "     rho: ", round(BAT_redspearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))
BAT_bluespearman_t <- cor.test(bluetable_t_BAT$kWithin , bluetable_t_BAT$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
BAT_blue_t <- ggplot(bluetable_t_BAT, aes(x=kWithin, y=t))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("blue     pval: ", round(BAT_bluespearman$p.value,digits = 3), "     rho: ", round(BAT_bluespearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))
BAT_blackspearman_t <- cor.test(blacktable_t_BAT$kWithin , blacktable_t_BAT$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
BAT_black_t <- ggplot(blacktable_t_BAT, aes(x=kWithin, y=t))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("black     pval: ", round(BAT_blackspearman$p.value,digits = 3), "     rho: ", round(BAT_blackspearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))
BAT_turquoisespearman_t <- cor.test(turquoisetable_t_BAT$kWithin , turquoisetable_t_BAT$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
BAT_turquoise_t <- ggplot(turquoisetable_t_BAT, aes(x=kWithin, y=t))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("turquoise     pval: ", round(BAT_turquoisespearman$p.value,digits = 3), "     rho: ", round(BAT_turquoisespearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))

redtable_t_BRO <- subset(redtable_BRO, t>0.01 & t<1)
bluetable_t_BRO <- subset(bluetable_BRO, t>0.01 & t<1)
blacktable_t_BRO <- subset(blacktable_BRO, t>0.01 & t<1)
turquoisetable_t_BRO <- subset(turquoisetable_BRO, t>0.01 & t<1)

BRO_redspearman_t <- cor.test(redtable_t_BRO$kWithin , redtable_t_BRO$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
BRO_red_t <- ggplot(redtable_t_BRO, aes(x=kWithin, y=t))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("red     pval: ", round(BRO_redspearman$p.value,digits = 3), "     rho: ", round(BRO_redspearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))
BRO_bluespearman_t <- cor.test(bluetable_t_BRO$kWithin , bluetable_t_BRO$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
BRO_blue_t <- ggplot(bluetable_t_BRO, aes(x=kWithin, y=t))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("blue     pval: ", round(BRO_bluespearman$p.value,digits = 3), "     rho: ", round(BRO_bluespearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))
BRO_blackspearman_t <- cor.test(blacktable_t_BRO$kWithin , blacktable_t_BRO$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
BRO_black_t <- ggplot(blacktable_t_BRO, aes(x=kWithin, y=t))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("black     pval: ", round(BRO_blackspearman$p.value,digits = 3), "     rho: ", round(BRO_blackspearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))
BRO_turquoisespearman_t <- cor.test(turquoisetable_t_BRO$kWithin , turquoisetable_t_BRO$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
BRO_turquoise_t <- ggplot(turquoisetable_t_BRO, aes(x=kWithin, y=t))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="kWithin Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("turquoise     pval: ", round(BRO_turquoisespearman$p.value,digits = 3), "     rho: ", round(BRO_turquoisespearman$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))

grid.arrange(BAT_red_t, BAT_blue_t, BAT_black_t, BAT_turquoise_t, BRO_red_t, BRO_blue_t, BRO_black_t, BRO_turquoise_t, ncol = 4 , nrow = 2)

############################################################################################################################################  n connectivity VS dNdS ############################

nncon_bgm <- nearestNeighborConnectivity(datExpr0, nNeighbors = 2875, power = 19, type = "signed", corFnc = "bicor", corOptions = "use = 'p'", blockSize = 6000)

names_bgm <- as.vector(names(datExpr0))
nncon_bgm_name <- "con"
names_bgm_name <- "OG"
nncon_bgm_df <- data.frame(names_bgm,nncon_bgm)
colnames(nncon_bgm_df) <- c(names_bgm_name, nncon_bgm_name)

nnconectivity_bgm_dNdS_bat <- merge(nncon_bgm_df, formatted_con_dnds_BAT[, c("OG", "dNdS","dN","dS","t")], by="OG")
redtable_nncon_bgm_dNdS_bat <- subset(nnconectivity_bgm_dNdS_bat, nnconectivity_bgm_dNdS_bat$OG %in% red_genes)
redtable_nncon_bgm_dNdS_bat <- subset(redtable_nncon_bgm_dNdS_bat, dS>0.0000 & dNdS<1)
bluetable_nncon_bgm_dNdS_bat <- subset(nnconectivity_bgm_dNdS_bat, nnconectivity_bgm_dNdS_bat$OG %in% blue_genes)
bluetable_nncon_bgm_dNdS_bat <- subset(bluetable_nncon_bgm_dNdS_bat, dS>0.0000 & dNdS<1)
blacktable_nncon_bgm_dNdS_bat <- subset(nnconectivity_bgm_dNdS_bat, nnconectivity_bgm_dNdS_bat$OG %in% black_genes)
blacktable_nncon_bgm_dNdS_bat <- subset(blacktable_nncon_bgm_dNdS_bat, dS>0.0000 & dNdS < 1)
turquoisetable_nncon_bgm_dNdS_bat <- subset(nnconectivity_bgm_dNdS_bat, nnconectivity_bgm_dNdS_bat$OG %in% turquoise_genes)
turquoisetable_nncon_bgm_dNdS_bat <- subset(turquoisetable_nncon_bgm_dNdS_bat, dS>0.0000 & dNdS < 1)
tottable_nncon_bgm_dNdS_bat <- subset(nnconectivity_bgm_dNdS_bat, dS>0.0000 & dNdS < 1)

networktable_nncon_bgm_dNdS_bat <- rbind(redtable_nncon_bgm_dNdS_bat,redtable_nncon_bgm_dNdS_bat, blacktable_nncon_bgm_dNdS_bat,turquoisetable_nncon_bgm_dNdS_bat)
network_genes <- networktable_nncon_bgm_dNdS_bat$OG
nonnetworktable_nncon_bgm_dNdS_bat <- subset(nnconectivity_bgm_dNdS_bat, !(nnconectivity_bgm_dNdS_bat$OG %in% network_genes))
nonnetworktable_nncon_bgm_dNdS_bat <- subset(nonnetworktable_nncon_bgm_dNdS_bat, dS>0.0000 & dNdS < 1)

redspearman_nncon_bgm_dNdS_bat <- cor.test(redtable_nncon_bgm_dNdS_bat$con , redtable_nncon_bgm_dNdS_bat$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
red_nncon_bgm_dNdS_bat <- ggplot(redtable_nncon_bgm_dNdS_bat, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 1) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus atticus (asex)",  title = paste0("red     pval: ", round(BAT_redspearman_nncon_bgm_dNdS_bat$p.value,digits = 3), "     rho: ", round(BAT_redspearman_nncon_bgm_dNdS_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
bluespearman_nncon_bgm_dNdS_bat <- cor.test(bluetable_nncon_bgm_dNdS_bat$con , bluetable_nncon_bgm_dNdS_bat$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
blue_nncon_bgm_dNdS_bat <- ggplot(bluetable_nncon_bgm_dNdS_bat, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 1) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus atticus (asex)",  title = paste0("blue     pval: ", round(BAT_bluespearman_nncon_bgm_dNdS_bat$p.value,digits = 3), "     rho: ", round(BAT_bluespearman_nncon_bgm_dNdS_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
blackspearman_nncon_bgm_dNdS_bat <- cor.test(blacktable_nncon_bgm_dNdS_bat$con , blacktable_nncon_bgm_dNdS_bat$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
black_nncon_bgm_dNdS_bat <- ggplot(blacktable_nncon_bgm_dNdS_bat, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 1) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus atticus (asex)",  title = paste0("black     pval: ", round(BAT_blackspearman_nncon_bgm_dNdS_bat$p.value,digits = 3), "     rho: ", round(BAT_blackspearman_nncon_bgm_dNdS_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
turquoisespearman_nncon_bgm_dNdS_bat <- cor.test(turquoisetable_nncon_bgm_dNdS_bat$con , turquoisetable_nncon_bgm_dNdS_bat$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
turquoise_nncon_bgm_dNdS_bat <- ggplot(turquoisetable_nncon_bgm_dNdS_bat, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 1) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus atticus (asex)",  title = paste0("turquoise     pval: ", round(BAT_turquoisespearman_nncon_bgm_dNdS_bat$p.value,digits = 3), "     rho: ", round(BAT_turquoisespearman_nncon_bgm_dNdS_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
totspearman_nncon_bgm_dNdS_bat <- cor.test(tottable_nncon_bgm_dNdS_bat$con , tottable_nncon_bgm_dNdS_bat$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
tot_nncon_bgm_dNdS_bat <- ggplot(tottable_nncon_bgm_dNdS_bat, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 1) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus atticus (asex)",  title = paste0("total     pval: ", round(BAT_totspearman_nncon_bgm_dNdS_bat$p.value,digits = 3), "     rho: ", round(BAT_totspearman_nncon_bgm_dNdS_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

networkspearman_nncon_bgm_dNdS_bat <- cor.test(networktable_nncon_bgm_dNdS_bat$con , networktable_nncon_bgm_dNdS_bat$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
network_nncon_bgm_dNdS_bat <- ggplot(networktable_nncon_bgm_dNdS_bat, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 1) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus atticus (asex)",  title = paste0("network     pval: ", round(networkspearman_nncon_bgm_dNdS_bat$p.value,digits = 3), "     rho: ", round(networkspearman_nncon_bgm_dNdS_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
nonnetworkspearman_nncon_bgm_dNdS_bat <- cor.test(nonnetworktable_nncon_bgm_dNdS_bat$con , nonnetworktable_nncon_bgm_dNdS_bat$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
nonnetwork_nncon_bgm_dNdS_bat <- ggplot(nonnetworktable_nncon_bgm_dNdS_bat, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 1) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus atticus (asex)",  title = paste0("non-network     pval: ", round(nonnetworkspearman_nncon_bgm_dNdS_bat$p.value,digits = 3), "     rho: ", round(nonnetworkspearman_nncon_bgm_dNdS_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

nnconectivity_bgm_dNdS_bro <- merge(nncon_bgm_df, formatted_con_dnds_BRO[, c("OG", "dNdS","dN","dS","t")], by="OG")
redtable_nncon_bgm_dNdS_bro <- subset(nnconectivity_bgm_dNdS_bro, nnconectivity_bgm_dNdS_bro$OG %in% red_genes)
redtable_nncon_bgm_dNdS_bro <- subset(redtable_nncon_bgm_dNdS_bro, dS>0.0000 & dNdS < 1)
bluetable_nncon_bgm_dNdS_bro <- subset(nnconectivity_bgm_dNdS_bro, nnconectivity_bgm_dNdS_bro$OG %in% blue_genes)
bluetable_nncon_bgm_dNdS_bro <- subset(bluetable_nncon_bgm_dNdS_bro, dS>0.0000 & dNdS < 1)
blacktable_nncon_bgm_dNdS_bro <- subset(nnconectivity_bgm_dNdS_bro, nnconectivity_bgm_dNdS_bro$OG %in% black_genes)
blacktable_nncon_bgm_dNdS_bro <- subset(blacktable_nncon_bgm_dNdS_bro, dS>0.0000 & dNdS < 1)
turquoisetable_nncon_bgm_dNdS_bro <- subset(nnconectivity_bgm_dNdS_bro, nnconectivity_bgm_dNdS_bro$OG %in% turquoise_genes)
turquoisetable_nncon_bgm_dNdS_bro <- subset(turquoisetable_nncon_bgm_dNdS_bro, dS>0.0000 & dNdS < 1)
tottable_nncon_bgm_dNdS_bro <- subset(nnconectivity_bgm_dNdS_bro, dS>0.0000 & dNdS < 1)

networktable_nncon_bgm_dNdS_bro <- rbind(redtable_nncon_bgm_dNdS_bro,redtable_nncon_bgm_dNdS_bro, blacktable_nncon_bgm_dNdS_bat,turquoisetable_nncon_bgm_dNdS_bat)
network_genes <- networktable_nncon_bgm_dNdS_bro$OG
nonnetworktable_nncon_bgm_dNdS_bro <- subset(networktable_nncon_bgm_dNdS_bro, !(nnconectivity_bgm_dNdS_bat$OG %in% network_genes))
nonnetworktable_nncon_bgm_dNdS_bro <- subset(nonnetworktable_nncon_bgm_dNdS_bro, dS>0.0000 & dNdS < 1)

redspearman_nncon_bgm_dNdS_bro <- cor.test(redtable_nncon_bgm_dNdS_bro$con , redtable_nncon_bgm_dNdS_bro$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
red_nncon_bgm_dNdS_bro <- ggplot(redtable_nncon_bgm_dNdS_bro, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 0.75) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus rossius (asex)",  title = paste0("red     pval: ", round(BRO_redspearman_nncon_bgm_dNdS_bro$p.value,digits = 3), "     rho: ", round(BRO_redspearman_nncon_bgm_dNdS_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
bluespearman_nncon_bgm_dNdS_bro <- cor.test(bluetable_nncon_bgm_dNdS_bro$con , bluetable_nncon_bgm_dNdS_bro$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
blue_nncon_bgm_dNdS_bro <- ggplot(bluetable_nncon_bgm_dNdS_bro, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 0.75) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus rossius (asex)",  title = paste0("blue     pval: ", round(BRO_bluespearman_nncon_bgm_dNdS_bro$p.value,digits = 3), "     rho: ", round(BRO_bluespearman_nncon_bgm_dNdS_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
blackspearman_nncon_bgm_dNdS_bro <- cor.test(blacktable_nncon_bgm_dNdS_bro$con , blacktable_nncon_bgm_dNdS_bro$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
black_nncon_bgm_dNdS_bro <- ggplot(blacktable_nncon_bgm_dNdS_bro, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 0.75) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus rossius (asex)",  title = paste0("black     pval: ", round(BRO_blackspearman_nncon_bgm_dNdS_bro$p.value,digits = 3), "     rho: ", round(BRO_blackspearman_nncon_bgm_dNdS_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
turquoisespearman_nncon_bgm_dNdS_bro <- cor.test(turquoisetable_nncon_bgm_dNdS_bro$con , turquoisetable_nncon_bgm_dNdS_bro$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
turquoise_nncon_bgm_dNdS_bro <- ggplot(turquoisetable_nncon_bgm_dNdS_bro, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 0.75) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus rossius (asex)",  title = paste0("turquoise     pval: ", round(BRO_turquoisespearman_nncon_bgm_dNdS_bro$p.value,digits = 3), "     rho: ", round(BRO_turquoisespearman_nncon_bgm_dNdS_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
totspearman_nncon_bgm_dNdS_bro <- cor.test(tottable_nncon_bgm_dNdS_bro$con , tottable_nncon_bgm_dNdS_bro$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
tot_nncon_bgm_dNdS_bro <- ggplot(tottable_nncon_bgm_dNdS_bro, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 0.75) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus rossius (asex)",  title = paste0("total     pval: ", round(BRO_totspearman_nncon_bgm_dNdS_bro$p.value,digits = 3), "     rho: ", round(BRO_totspearman_nncon_bgm_dNdS_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

networkspearman_nncon_bgm_dNdS_bro <- cor.test(networktable_nncon_bgm_dNdS_bro$con , networktable_nncon_bgm_dNdS_bro$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
network_nncon_bgm_dNdS_bro <- ggplot(networktable_nncon_bgm_dNdS_bro, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 1) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus rossius (asex)",  title = paste0("network     pval: ", round(networkspearman_nncon_bgm_dNdS_bro$p.value,digits = 3), "     rho: ", round(networkspearman_nncon_bgm_dNdS_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
nonnetworkspearman_nncon_bgm_dNdS_bro <- cor.test(nonnetworktable_nncon_bgm_dNdS_bro$con , nonnetworktable_nncon_bgm_dNdS_bro$dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
nonnetwork_nncon_bgm_dNdS_bro <- ggplot(nonnetworktable_nncon_bgm_dNdS_bro, aes(x=con, y=dNdS))    + geom_point(shape=16, color="black", size = 1) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus rossius (asex)",  title = paste0("non-network     pval: ", round(nonnetworkspearman_nncon_bgm_dNdS_bro$p.value,digits = 3), "     rho: ", round(nonnetworkspearman_nncon_bgm_dNdS_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

nnconectivity_bgm_dNdS_bgm <- merge(nncon_bgm_df, tot[, c("OG", "BGM_dNdS","BGM_dN","BGM_dS","BGM_t")], by="OG")
redtable_nncon_bgm_dNdS_bgm <- subset(nnconectivity_bgm_dNdS_bgm, nnconectivity_bgm_dNdS_bgm$OG %in% red_genes)
redtable_nncon_bgm_dNdS_bgm <- subset(redtable_nncon_bgm_dNdS_bgm, BGM_dS>0.0000 & BGM_dNdS < 1)
bluetable_nncon_bgm_dNdS_bgm <- subset(nnconectivity_bgm_dNdS_bgm, nnconectivity_bgm_dNdS_bgm$OG %in% blue_genes)
bluetable_nncon_bgm_dNdS_bgm <- subset(bluetable_nncon_bgm_dNdS_bgm, BGM_dS>0.0000 & BGM_dNdS < 1)
blacktable_nncon_bgm_dNdS_bgm <- subset(nnconectivity_bgm_dNdS_bgm, nnconectivity_bgm_dNdS_bgm$OG %in% black_genes)
blacktable_nncon_bgm_dNdS_bgm <- subset(blacktable_nncon_bgm_dNdS_bgm, BGM_dS>0.0000 & BGM_dNdS < 1)
turquoisetable_nncon_bgm_dNdS_bgm <- subset(nnconectivity_bgm_dNdS_bgm, nnconectivity_bgm_dNdS_bgm$OG %in% turquoise_genes)
turquoisetable_nncon_bgm_dNdS_bgm <- subset(turquoisetable_nncon_bgm_dNdS_bgm, BGM_dS>0.0000 & BGM_dNdS < 1)
tottable_nncon_bgm_dNdS_bgm <- subset(nnconectivity_bgm_dNdS_bgm, BGM_dS>0.0000 & BGM_dNdS < 1)

networktable_nncon_bgm_dNdS_bgm <- rbind(redtable_nncon_bgm_dNdS_bgm,bluetable_nncon_bgm_dNdS_bgm, blacktable_nncon_bgm_dNdS_bgm,turquoisetable_nncon_bgm_dNdS_bgm)
network_genes <- networktable_nncon_bgm_dNdS_bgm$OG
nonnetworktable_nncon_bgm_dNdS_bgm <- subset(networktable_nncon_bgm_dNdS_bgm, !(nnconectivity_bgm_dNdS_bat$OG %in% network_genes))
nonnetworktable_nncon_bgm_dNdS_bgm <- subset(nonnetworktable_nncon_bgm_dNdS_bgm, dS>0.0000 & dNdS < 1)

redspearman_nncon_bgm_dNdS_bgm <- cor.test(redtable_nncon_bgm_dNdS_bgm$con , redtable_nncon_bgm_dNdS_bgm$BGM_dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
red_nncon_bgm_dNdS_bgm <- ggplot(redtable_nncon_bgm_dNdS_bgm, aes(x=con, y=BGM_dNdS))    + geom_point(shape=16, color="black", size = 0.75) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus grandii (sex)",  title = paste0("red     pval: ", round(BRO_redspearman_nncon_bgm_dNdS_bgm$p.value,digits = 3), "     rho: ", round(BRO_redspearman_nncon_bgm_dNdS_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
bluespearman_nncon_bgm_dNdS_bgm <- cor.test(bluetable_nncon_bgm_dNdS_bgm$con , bluetable_nncon_bgm_dNdS_bgm$BGM_dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
blue_nncon_bgm_dNdS_bgm <- ggplot(bluetable_nncon_bgm_dNdS_bgm, aes(x=con, y=BGM_dNdS))    + geom_point(shape=16, color="black", size = 0.75) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus grandii (sex)",  title = paste0("blue     pval: ", round(BRO_bluespearman_nncon_bgm_dNdS_bgm$p.value,digits = 3), "     rho: ", round(BRO_bluespearman_nncon_bgm_dNdS_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
blackspearman_nncon_bgm_dNdS_bgm <- cor.test(blacktable_nncon_bgm_dNdS_bgm$con , blacktable_nncon_bgm_dNdS_bgm$BGM_dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
black_nncon_bgm_dNdS_bgm <- ggplot(blacktable_nncon_bgm_dNdS_bgm, aes(x=con, y=BGM_dNdS))    + geom_point(shape=16, color="black", size = 0.75) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus grandii (sex)",  title = paste0("black     pval: ", round(BRO_blackspearman_nncon_bgm_dNdS_bgm$p.value,digits = 3), "     rho: ", round(BRO_blackspearman_nncon_bgm_dNdS_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
turquoisespearman_nncon_bgm_dNdS_bgm <- cor.test(turquoisetable_nncon_bgm_dNdS_bgm$con , turquoisetable_nncon_bgm_dNdS_bgm$BGM_dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
turquoise_nncon_bgm_dNdS_bgm <- ggplot(turquoisetable_nncon_bgm_dNdS_bgm, aes(x=con, y=BGM_dNdS))    + geom_point(shape=16, color="black", size = 0.75) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus grandii (sex)",  title = paste0("turquoise     pval: ", round(BRO_turquoisespearman_nncon_bgm_dNdS_bgm$p.value,digits = 3), "     rho: ", round(BRO_turquoisespearman_nncon_bgm_dNdS_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
totspearman_nncon_bgm_dNdS_bgm <- cor.test(tottable_nncon_bgm_dNdS_bgm$con , tottable_nncon_bgm_dNdS_bgm$BGM_dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
tot_nncon_bgm_dNdS_bgm <- ggplot(tottable_nncon_bgm_dNdS_bgm, aes(x=con, y=BGM_dNdS))    + geom_point(shape=16, color="black", size = 0.75) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus grandii (sex)",  title = paste0("total     pval: ", round(BRO_totspearman_nncon_bgm_dNdS_bgm$p.value,digits = 3), "     rho: ", round(BRO_totspearman_nncon_bgm_dNdS_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

networkspearman_nncon_bgm_dNdS_bgm <- cor.test(networktable_nncon_bgm_dNdS_bgm$con , networktable_nncon_bgm_dNdS_bgm$BGM_dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
network_nncon_bgm_dNdS_bgm <- ggplot(networktable_nncon_bgm_dNdS_bgm, aes(x=con, y=BGM_dNdS))    + geom_point(shape=16, color="black", size = 1) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus grandii (sex)",  title = paste0("network     pval: ", round(networkspearman_nncon_bgm_dNdS_bgm$p.value,digits = 3), "     rho: ", round(networkspearman_nncon_bgm_dNdS_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
nonnetworkspearman_nncon_bgm_dNdS_bgm <- cor.test(nonnetworktable_nncon_bgm_dNdS_bgm$con , nonnetworktable_nncon_bgm_dNdS_bgm$BGM_dNdS, method = "spearman", continuity = FALSE, conf.level = 0.95)
nonnetwork_nncon_bgm_dNdS_bgm <- ggplot(nonnetworktable_nncon_bgm_dNdS_bgm, aes(x=con, y=BGM_dNdS))    + geom_point(shape=16, color="black", size = 1) + #geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "dNdS Bacillus grandii (sex)",  title = paste0("non-network     pval: ", round(nonnetworkspearman_nncon_bgm_dNdS_bgm$p.value,digits = 3), "     rho: ", round(nonnetworkspearman_nncon_bgm_dNdS_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))


grid.arrange(red_nncon_bgm_dNdS_bat, blue_nncon_bgm_dNdS_bat, black_nncon_bgm_dNdS_bat, turquoise_nncon_bgm_dNdS_bat, tot_nncon_bgm_dNdS_bat,
             red_nncon_bgm_dNdS_bro, blue_nncon_bgm_dNdS_bro, black_nncon_bgm_dNdS_bro, turquoise_nncon_bgm_dNdS_bro, tot_nncon_bgm_dNdS_bro,
             red_nncon_bgm_dNdS_bgm, blue_nncon_bgm_dNdS_bgm, black_nncon_bgm_dNdS_bgm, turquoise_nncon_bgm_dNdS_bgm, tot_nncon_bgm_dNdS_bgm,
             ncol = 5 , nrow = 3)


grid.arrange(network_nncon_bgm_dNdS_bat, nonnetwork_nncon_bgm_dNdS_bat,
             network_nncon_bgm_dNdS_bro, nonnetwork_nncon_bgm_dNdS_bro,
             network_nncon_bgm_dNdS_bgm, nonnetwork_nncon_bgm_dNdS_bgm,
             ncol = 2 , nrow = 3)
             

############################################################################################################################################  n connectivity VS CUB ############################

nnconectivity_bgm_CUB_bat <- merge(nncon_bgm_df, tot[, c("OG", "ENC_BAT","ENC_BRO","ENC_BGM")], by="OG")
redtable_nncon_bgm_CUB_bat <- subset(nnconectivity_bgm_CUB_bat, nnconectivity_bgm_CUB_bat$OG %in% red_genes)
bluetable_nncon_bgm_CUB_bat <- subset(nnconectivity_bgm_CUB_bat, nnconectivity_bgm_CUB_bat$OG %in% blue_genes)
blacktable_nncon_bgm_CUB_bat <- subset(nnconectivity_bgm_CUB_bat, nnconectivity_bgm_CUB_bat$OG %in% black_genes)
turquoisetable_nncon_bgm_CUB_bat <- subset(nnconectivity_bgm_CUB_bat, nnconectivity_bgm_CUB_bat$OG %in% turquoise_genes)

BAT_redspearman_nncon_bgm_CUB_bat <- cor.test(redtable_nncon_bgm_CUB_bat$con , redtable_nncon_bgm_CUB_bat$ENC_BAT, method = "spearman", continuity = FALSE, conf.level = 0.95)
BAT_red_nncon_bgm_CUB_bat <- ggplot(redtable_nncon_bgm_CUB_bat, aes(x=con, y=ENC_BAT))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus atticus (asex)",  title = paste0("red     pval: ", round(BAT_redspearman_nncon_bgm_CUB_bat$p.value,digits = 3), "     rho: ", round(BAT_redspearman_nncon_bgm_CUB_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

BAT_bluespearman_nncon_bgm_CUB_bat <- cor.test(bluetable_nncon_bgm_CUB_bat$con , bluetable_nncon_bgm_CUB_bat$ENC_BAT, method = "spearman", continuity = FALSE, conf.level = 0.95)
BAT_blue_nncon_bgm_CUB_bat <- ggplot(bluetable_nncon_bgm_CUB_bat, aes(x=con, y=ENC_BAT))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus atticus (asex)",  title = paste0("blue     pval: ", round(BAT_bluespearman_nncon_bgm_CUB_bat$p.value,digits = 3), "     rho: ", round(BAT_bluespearman_nncon_bgm_CUB_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

BAT_blackspearman_nncon_bgm_CUB_bat <- cor.test(blacktable_nncon_bgm_CUB_bat$con , blacktable_nncon_bgm_CUB_bat$ENC_BAT, method = "spearman", continuity = FALSE, conf.level = 0.95)
BAT_black_nncon_bgm_CUB_bat <- ggplot(blacktable_nncon_bgm_CUB_bat, aes(x=con, y=ENC_BAT))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus atticus (asex)",  title = paste0("black     pval: ", round(BAT_blackspearman_nncon_bgm_CUB_bat$p.value,digits = 3), "     rho: ", round(BAT_blackspearman_nncon_bgm_CUB_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

BAT_turquoisespearman_nncon_bgm_CUB_bat <- cor.test(turquoisetable_nncon_bgm_CUB_bat$con , turquoisetable_nncon_bgm_CUB_bat$ENC_BAT, method = "spearman", continuity = FALSE, conf.level = 0.95)
BAT_turquoise_nncon_bgm_CUB_bat <- ggplot(turquoisetable_nncon_bgm_CUB_bat, aes(x=con, y=ENC_BAT))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus atticus (asex)",  title = paste0("turquoise     pval: ", round(BAT_turquoisespearman_nncon_bgm_CUB_bat$p.value,digits = 3), "     rho: ", round(BAT_turquoisespearman_nncon_bgm_CUB_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

BAT_totspearman_nncon_bgm_CUB_bat <- cor.test(nnconectivity_bgm_CUB_bat$con , nnconectivity_bgm_CUB_bat$ENC_BAT, method = "spearman", continuity = FALSE, conf.level = 0.95)
BAT_tot_nncon_bgm_CUB_bat <- ggplot(nnconectivity_bgm_CUB_bat, aes(x=con, y=ENC_BAT))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus atticus (asex)",  title = paste0("total     pval: ", round(BAT_totspearman_nncon_bgm_CUB_bat$p.value,digits = 3), "     rho: ", round(BAT_totspearman_nncon_bgm_CUB_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

nnconectivity_bgm_CUB_bro <- merge(nncon_bgm_df, tot[, c("OG", "ENC_BRO","ENC_BRO","ENC_BGM")], by="OG")
redtable_nncon_bgm_CUB_bro <- subset(nnconectivity_bgm_CUB_bro, nnconectivity_bgm_CUB_bro$OG %in% red_genes)
bluetable_nncon_bgm_CUB_bro <- subset(nnconectivity_bgm_CUB_bro, nnconectivity_bgm_CUB_bro$OG %in% blue_genes)
blacktable_nncon_bgm_CUB_bro <- subset(nnconectivity_bgm_CUB_bro, nnconectivity_bgm_CUB_bro$OG %in% black_genes)
turquoisetable_nncon_bgm_CUB_bro <- subset(nnconectivity_bgm_CUB_bro, nnconectivity_bgm_dNdS_bro$OG %in% turquoise_genes)

BRO_redspearman_nncon_bgm_CUB_bro <- cor.test(redtable_nncon_bgm_CUB_bro$con , redtable_nncon_bgm_CUB_bro$ENC_BRO, method = "spearman", continuity = FALSE, conf.level = 0.95)
BRO_red_nncon_bgm_CUB_bro <- ggplot(redtable_nncon_bgm_CUB_bro, aes(x=con, y=ENC_BRO))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus atticus (asex)",  title = paste0("red     pval: ", round(BRO_redspearman_nncon_bgm_CUB_bro$p.value,digits = 3), "     rho: ", round(BRO_redspearman_nncon_bgm_CUB_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

BRO_bluespearman_nncon_bgm_CUB_bro <- cor.test(bluetable_nncon_bgm_CUB_bro$con , bluetable_nncon_bgm_CUB_bro$ENC_BRO, method = "spearman", continuity = FALSE, conf.level = 0.95)
BRO_blue_nncon_bgm_CUB_bro <- ggplot(bluetable_nncon_bgm_CUB_bro, aes(x=con, y=ENC_BRO))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus atticus (asex)",  title = paste0("blue     pval: ", round(BRO_bluespearman_nncon_bgm_CUB_bro$p.value,digits = 3), "     rho: ", round(BRO_bluespearman_nncon_bgm_CUB_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

BRO_blackspearman_nncon_bgm_CUB_bro <- cor.test(blacktable_nncon_bgm_CUB_bro$con , blacktable_nncon_bgm_CUB_bro$ENC_BRO, method = "spearman", continuity = FALSE, conf.level = 0.95)
BRO_black_nncon_bgm_CUB_bro <- ggplot(blacktable_nncon_bgm_CUB_bro, aes(x=con, y=ENC_BRO))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus atticus (asex)",  title = paste0("black     pval: ", round(BRO_blackspearman_nncon_bgm_CUB_bro$p.value,digits = 3), "     rho: ", round(BRO_blackspearman_nncon_bgm_CUB_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

BRO_turquoisespearman_nncon_bgm_CUB_bro <- cor.test(turquoisetable_nncon_bgm_CUB_bro$con , turquoisetable_nncon_bgm_CUB_bro$ENC_BRO, method = "spearman", continuity = FALSE, conf.level = 0.95)
BRO_turquoise_nncon_bgm_CUB_bro <- ggplot(turquoisetable_nncon_bgm_CUB_bro, aes(x=con, y=ENC_BRO))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus atticus (asex)",  title = paste0("turquoise     pval: ", round(BRO_turquoisespearman_nncon_bgm_CUB_bro$p.value,digits = 3), "     rho: ", round(BRO_turquoisespearman_nncon_bgm_CUB_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

BRO_totspearman_nncon_bgm_CUB_bro <- cor.test(nnconectivity_bgm_CUB_bro$con , nnconectivity_bgm_CUB_bro$ENC_BRO, method = "spearman", continuity = FALSE, conf.level = 0.95)
BRO_tot_nncon_bgm_CUB_bro <- ggplot(nnconectivity_bgm_CUB_bro, aes(x=con, y=ENC_BRO))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus atticus (asex)",  title = paste0("total     pval: ", round(BRO_totspearman_nncon_bgm_CUB_bro$p.value,digits = 3), "     rho: ", round(BRO_totspearman_nncon_bgm_CUB_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

nnconectivity_bgm_CUB_bgm <- merge(nncon_bgm_df, tot[, c("OG", "ENC_BGM","ENC_BRO","ENC_BGM")], by="OG")
redtable_nncon_bgm_CUB_bgm <- subset(nnconectivity_bgm_CUB_bgm, nnconectivity_bgm_CUB_bgm$OG %in% red_genes)
bluetable_nncon_bgm_CUB_bgm <- subset(nnconectivity_bgm_CUB_bgm, nnconectivity_bgm_CUB_bgm$OG %in% blue_genes)
blacktable_nncon_bgm_CUB_bgm <- subset(nnconectivity_bgm_CUB_bgm, nnconectivity_bgm_CUB_bgm$OG %in% black_genes)
turquoisetable_nncon_bgm_CUB_bgm <- subset(nnconectivity_bgm_CUB_bgm, nnconectivity_bgm_CUB_bgm$OG %in% turquoise_genes)

BGM_redspearman_nncon_bgm_CUB_bgm <- cor.test(redtable_nncon_bgm_CUB_bgm$con , redtable_nncon_bgm_CUB_bgm$ENC_BGM, method = "spearman", continuity = FALSE, conf.level = 0.95)
BGM_red_nncon_bgm_CUB_bgm <- ggplot(redtable_nncon_bgm_CUB_bgm, aes(x=con, y=ENC_BGM))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus grandii (sex)",  title = paste0("red     pval: ", round(BGM_redspearman_nncon_bgm_CUB_bgm$p.value,digits = 3), "     rho: ", round(BGM_redspearman_nncon_bgm_CUB_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

BGM_bluespearman_nncon_bgm_CUB_bgm <- cor.test(bluetable_nncon_bgm_CUB_bgm$con , bluetable_nncon_bgm_CUB_bgm$ENC_BGM, method = "spearman", continuity = FALSE, conf.level = 0.95)
BGM_blue_nncon_bgm_CUB_bgm <- ggplot(bluetable_nncon_bgm_CUB_bgm, aes(x=con, y=ENC_BGM))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus grandii (sex)",  title = paste0("blue     pval: ", round(BGM_bluespearman_nncon_bgm_CUB_bgm$p.value,digits = 3), "     rho: ", round(BGM_bluespearman_nncon_bgm_CUB_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

BGM_blackspearman_nncon_bgm_CUB_bgm <- cor.test(blacktable_nncon_bgm_CUB_bgm$con , blacktable_nncon_bgm_CUB_bgm$ENC_BGM, method = "spearman", continuity = FALSE, conf.level = 0.95)
BGM_black_nncon_bgm_CUB_bgm <- ggplot(blacktable_nncon_bgm_CUB_bgm, aes(x=con, y=ENC_BGM))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus grandii (sex)",  title = paste0("black     pval: ", round(BGM_blackspearman_nncon_bgm_CUB_bgm$p.value,digits = 3), "     rho: ", round(BGM_blackspearman_nncon_bgm_CUB_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

BGM_turquoisespearman_nncon_bgm_CUB_bgm <- cor.test(turquoisetable_nncon_bgm_CUB_bgm$con , turquoisetable_nncon_bgm_CUB_bgm$ENC_BGM, method = "spearman", continuity = FALSE, conf.level = 0.95)
BGM_turquoise_nncon_bgm_CUB_bgm <- ggplot(turquoisetable_nncon_bgm_CUB_bgm, aes(x=con, y=ENC_BGM))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus grandii (sex)",  title = paste0("turquoise     pval: ", round(BGM_turquoisespearman_nncon_bgm_CUB_bgm$p.value,digits = 3), "     rho: ", round(BGM_turquoisespearman_nncon_bgm_CUB_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

BGM_totspearman_nncon_bgm_CUB_bgm <- cor.test(nnconectivity_bgm_CUB_bgm$con , nnconectivity_bgm_CUB_bgm$ENC_BGM, method = "spearman", continuity = FALSE, conf.level = 0.95)
BGM_tot_nncon_bgm_CUB_bgm <- ggplot(nnconectivity_bgm_CUB_bgm, aes(x=con, y=ENC_BGM))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "CUB Bacillus grandii (sex)",  title = paste0("total     pval: ", round(BGM_totspearman_nncon_bgm_CUB_bgm$p.value,digits = 3), "     rho: ", round(BGM_totspearman_nncon_bgm_CUB_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

grid.arrange(BAT_red_nncon_bgm_CUB_bat, BAT_blue_nncon_bgm_CUB_bat, BAT_black_nncon_bgm_CUB_bat, BAT_turquoise_nncon_bgm_CUB_bat, BAT_tot_nncon_bgm_CUB_bat,
             BRO_red_nncon_bgm_CUB_bro, BRO_blue_nncon_bgm_CUB_bro, BRO_black_nncon_bgm_CUB_bro, BRO_turquoise_nncon_bgm_CUB_bro, BRO_tot_nncon_bgm_CUB_bro,
             BGM_red_nncon_bgm_CUB_bgm, BGM_blue_nncon_bgm_CUB_bgm, BGM_black_nncon_bgm_CUB_bgm, BGM_turquoise_nncon_bgm_CUB_bgm, BGM_tot_nncon_bgm_CUB_bgm,
             ncol = 5 , nrow = 3)

############################################################################################################################################  n connectivity VS t ############################

redtable_nncon_bgm_t_bat <- subset(nnconectivity_bgm_dNdS_bat, nnconectivity_bgm_dNdS_bat$OG %in% red_genes)
redtable_nncon_bgm_t_bat <- subset(redtable_nncon_bgm_t_bat, t<10)
redtable_nncon_bgm_t_bat <- subset(nnconectivity_bgm_dNdS_bat, nnconectivity_bgm_dNdS_bat$OG %in% blue_genes)
bluetable_nncon_bgm_t_bat <- subset(bluetable_nncon_bgm_t_bat, t<10)
blacktable_nncon_bgm_t_bat <- subset(nnconectivity_bgm_dNdS_bat, nnconectivity_bgm_dNdS_bat$OG %in% black_genes)
blacktable_nncon_bgm_t_bat <- subset(blacktable_nncon_bgm_t_bat, t<10)
turquoisetable_nncon_bgm_t_bat <- subset(nnconectivity_bgm_dNdS_bat, nnconectivity_bgm_dNdS_bat$OG %in% turquoise_genes)
turquoisetable_nncon_bgm_t_bat <- subset(turquoisetable_nncon_bgm_t_bat, t<10)
tottable_nncon_bgm_t_bat <- subset(nnconectivity_bgm_dNdS_bat, t<10)

networktable_nncon_bgm_t_bat <- rbind(redtable_nncon_bgm_t_bat,redtable_nncon_bgm_t_bat, blacktable_nncon_bgm_t_bat,turquoisetable_nncon_bgm_t_bat)
network_genes <- networktable_nncon_bgm_t_bat$OG
nonnetworktable_nncon_bgm_t_bat <- subset(nnconectivity_bgm_dNdS_bat, !(nnconectivity_bgm_dNdS_bat$OG %in% network_genes))

redspearman_nncon_bgm_t_bat <- cor.test(redtable_nncon_bgm_t_bat$con , redtable_nncon_bgm_t_bat$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
red_nncon_bgm_t_bat <- ggplot(redtable_nncon_bgm_t_bat, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("total     pval: ", round(redspearman_nncon_bgm_t_bat$p.value,digits = 3), "     rho: ", round(redspearman_nncon_bgm_t_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
bluespearman_nncon_bgm_t_bat <- cor.test(bluetable_nncon_bgm_t_bat$con , bluetable_nncon_bgm_t_bat$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
blue_nncon_bgm_t_bat <- ggplot(bluetable_nncon_bgm_t_bat, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("total     pval: ", round(bluespearman_nncon_bgm_t_bat$p.value,digits = 3), "     rho: ", round(bluespearman_nncon_bgm_t_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
blackspearman_nncon_bgm_t_bat <- cor.test(blacktable_nncon_bgm_t_bat$con , blacktable_nncon_bgm_t_bat$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
black_nncon_bgm_t_bat <- ggplot(blacktable_nncon_bgm_t_bat, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("total     pval: ", round(blackspearman_nncon_bgm_t_bat$p.value,digits = 3), "     rho: ", round(blackspearman_nncon_bgm_t_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
turquoisespearman_nncon_bgm_t_bat <- cor.test(turquoisetable_nncon_bgm_t_bat$con , turquoisetable_nncon_bgm_t_bat$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
turquoise_nncon_bgm_t_bat <- ggplot(turquoisetable_nncon_bgm_t_bat, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("total     pval: ", round(turquoisespearman_nncon_bgm_t_bat$p.value,digits = 3), "     rho: ", round(turquoisespearman_nncon_bgm_t_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
totspearman_nncon_bgm_t_bat <- cor.test(tottable_nncon_bgm_t_bat$con , tottable_nncon_bgm_t_bat$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
tot_nncon_bgm_t_bat <- ggplot(tottable_nncon_bgm_t_bat, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("total     pval: ", round(totspearman_nncon_bgm_t_bat$p.value,digits = 3), "     rho: ", round(totspearman_nncon_bgm_t_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

networkspearman_nncon_bgm_t_bat <- cor.test(networktable_nncon_bgm_t_bat$con , networktable_nncon_bgm_t_bat$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
network_nncon_bgm_t_bat <- ggplot(networktable_nncon_bgm_t_bat, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("network     pval: ", round(networkspearman_nncon_bgm_t_bat$p.value,digits = 3), "     rho: ", round(networkspearman_nncon_bgm_t_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
nonnetworkspearman_nncon_bgm_t_bat <- cor.test(nonnetworktable_nncon_bgm_t_bat$con , nonnetworktable_nncon_bgm_t_bat$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
nonnetwork_nncon_bgm_t_bat <- ggplot(nonnetworktable_nncon_bgm_t_bat, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus atticus (asex)",  title = paste0("non-network     pval: ", round(nonnetworkspearman_nncon_bgm_t_bat$p.value,digits = 3), "     rho: ", round(nonnetworkspearman_nncon_bgm_t_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

redtable_nncon_bgm_t_bro <- subset(nnconectivity_bgm_dNdS_bro, nnconectivity_bgm_dNdS_bro$OG %in% red_genes)
redtable_nncon_bgm_t_bro <- subset(redtable_nncon_bgm_t_bro, t<10)
bluetable_nncon_bgm_t_bro <- subset(nnconectivity_bgm_dNdS_bro, nnconectivity_bgm_dNdS_bro$OG %in% blue_genes)
bluetable_nncon_bgm_t_bro <- subset(bluetable_nncon_bgm_t_bro, t<10)
blacktable_nncon_bgm_t_bro <- subset(nnconectivity_bgm_dNdS_bro, nnconectivity_bgm_dNdS_bro$OG %in% black_genes)
blacktable_nncon_bgm_t_bro <- subset(blacktable_nncon_bgm_t_bro, t<10)
turquoisetable_nncon_bgm_t_bro <- subset(nnconectivity_bgm_dNdS_bro, nnconectivity_bgm_dNdS_bro$OG %in% turquoise_genes)
turquoisetable_nncon_bgm_t_bro <- subset(turquoisetable_nncon_bgm_t_bro, t<10)
tottable_nncon_bgm_t_bro <- subset(nnconectivity_bgm_dNdS_bro, t<10)

networktable_nncon_bgm_t_bro <- rbind(redtable_nncon_bgm_t_bro,redtable_nncon_bgm_t_bro, blacktable_nncon_bgm_t_bro,turquoisetable_nncon_bgm_t_bro)
network_genes <- networktable_nncon_bgm_t_bro$OG
nonnetworktable_nncon_bgm_t_bro <- subset(nnconectivity_bgm_dNdS_bro, !(nnconectivity_bgm_dNdS_bro$OG %in% network_genes))

redspearman_nncon_bgm_t_bro <- cor.test(redtable_nncon_bgm_t_bro$con , redtable_nncon_bgm_t_bro$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
red_nncon_bgm_t_bro <- ggplot(redtable_nncon_bgm_t_bro, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus rossius (asex)",  title = paste0("total     pval: ", round(redspearman_nncon_bgm_t_bro$p.value,digits = 3), "     rho: ", round(redspearman_nncon_bgm_t_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
bluespearman_nncon_bgm_t_bro <- cor.test(bluetable_nncon_bgm_t_bro$con , bluetable_nncon_bgm_t_bro$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
blue_nncon_bgm_t_bro <- ggplot(bluetable_nncon_bgm_t_bro, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus rossius (asex)",  title = paste0("total     pval: ", round(bluespearman_nncon_bgm_t_bro$p.value,digits = 3), "     rho: ", round(bluespearman_nncon_bgm_t_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
blackspearman_nncon_bgm_t_bro <- cor.test(blacktable_nncon_bgm_t_bro$con , blacktable_nncon_bgm_t_bro$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
black_nncon_bgm_t_bro <- ggplot(blacktable_nncon_bgm_t_bro, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus rossius (asex)",  title = paste0("total     pval: ", round(blackspearman_nncon_bgm_t_bro$p.value,digits = 3), "     rho: ", round(blackspearman_nncon_bgm_t_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
turquoisespearman_nncon_bgm_t_bro <- cor.test(turquoisetable_nncon_bgm_t_bro$con , turquoisetable_nncon_bgm_t_bro$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
turquoise_nncon_bgm_t_bro <- ggplot(turquoisetable_nncon_bgm_t_bro, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus rossius (asex)",  title = paste0("total     pval: ", round(turquoisespearman_nncon_bgm_t_bro$p.value,digits = 3), "     rho: ", round(turquoisespearman_nncon_bgm_t_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
totspearman_nncon_bgm_t_bro <- cor.test(tottable_nncon_bgm_t_bro$con , tottable_nncon_bgm_t_bro$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
tot_nncon_bgm_t_bro <- ggplot(tottable_nncon_bgm_t_bro, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus rossius (asex)",  title = paste0("total     pval: ", round(totspearman_nncon_bgm_t_bro$p.value,digits = 3), "     rho: ", round(totspearman_nncon_bgm_t_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

networkspearman_nncon_bgm_t_bro <- cor.test(networktable_nncon_bgm_t_bro$con , networktable_nncon_bgm_t_bro$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
network_nncon_bgm_t_bro <- ggplot(networktable_nncon_bgm_t_bro, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus rossius (asex)",  title = paste0("network     pval: ", round(networkspearman_nncon_bgm_t_bro$p.value,digits = 3), "     rho: ", round(networkspearman_nncon_bgm_t_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
nonnetworkspearman_nncon_bgm_t_bro <- cor.test(nonnetworktable_nncon_bgm_t_bro$con , nonnetworktable_nncon_bgm_t_bro$t, method = "spearman", continuity = FALSE, conf.level = 0.95)
nonnetwork_nncon_bgm_t_bro <- ggplot(nonnetworktable_nncon_bgm_t_bro, aes(x=con, y=sqrt(t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus rossius (asex)",  title = paste0("non-network     pval: ", round(nonnetworkspearman_nncon_bgm_t_bro$p.value,digits = 3), "     rho: ", round(nonnetworkspearman_nncon_bgm_t_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

redtable_nncon_bgm_t_bgm <- subset(redtable_nncon_bgm_dNdS_bgm, BGM_t < 6)
bluetable_nncon_bgm_t_bgm <- subset(bluetable_nncon_bgm_dNdS_bgm, BGM_t < 6)
blacktable_nncon_bgm_t_bgm <- subset(blacktable_nncon_bgm_dNdS_bgm, BGM_t < 6)
turquoisetable_nncon_bgm_t_bgm <- subset(turquoisetable_nncon_bgm_dNdS_bgm, BGM_t < 6)
tottable_nncon_bgm_t_bgm <- subset(nnconectivity_bgm_dNdS_bgm, BGM_t < 6)

networktable_nncon_bgm_t_bgm <- rbind(redtable_nncon_bgm_t_bgm,redtable_nncon_bgm_t_bgm, blacktable_nncon_bgm_t_bgm,turquoisetable_nncon_bgm_t_bgm)
network_genes <- networktable_nncon_bgm_t_bgm$OG
nonnetworktable_nncon_bgm_t_bgm <- subset(nnconectivity_bgm_dNdS_bgm, !(nnconectivity_bgm_dNdS_bgm$OG %in% network_genes))

redspearman_nncon_bgm_t_bgm <- cor.test(redtable_nncon_bgm_t_bgm$con , redtable_nncon_bgm_t_bgm$BGM_t, method = "spearman", continuity = FALSE, conf.level = 0.95)
red_nncon_bgm_t_bgm <- ggplot(redtable_nncon_bgm_t_bgm, aes(x=con, y=sqrt(BGM_t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus grandii (sex)",  title = paste0("total     pval: ", round(redspearman_nncon_bgm_t_bgm$p.value,digits = 3), "     rho: ", round(redspearman_nncon_bgm_t_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
bluespearman_nncon_bgm_t_bgm <- cor.test(bluetable_nncon_bgm_t_bgm$con , bluetable_nncon_bgm_t_bgm$BGM_t, method = "spearman", continuity = FALSE, conf.level = 0.95)
blue_nncon_bgm_t_bgm <- ggplot(bluetable_nncon_bgm_t_bgm, aes(x=con, y=sqrt(BGM_t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus grandii (sex)",  title = paste0("total     pval: ", round(bluespearman_nncon_bgm_t_bgm$p.value,digits = 3), "     rho: ", round(bluespearman_nncon_bgm_t_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
blackspearman_nncon_bgm_t_bgm <- cor.test(blacktable_nncon_bgm_t_bgm$con , blacktable_nncon_bgm_t_bgm$BGM_t, method = "spearman", continuity = FALSE, conf.level = 0.95)
black_nncon_bgm_t_bgm <- ggplot(blacktable_nncon_bgm_t_bgm, aes(x=con, y=sqrt(BGM_t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus grandii (sex)",  title = paste0("total     pval: ", round(blackspearman_nncon_bgm_t_bgm$p.value,digits = 3), "     rho: ", round(blackspearman_nncon_bgm_t_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
turquoisespearman_nncon_bgm_t_bgm <- cor.test(turquoisetable_nncon_bgm_t_bgm$con , turquoisetable_nncon_bgm_t_bgm$BGM_t, method = "spearman", continuity = FALSE, conf.level = 0.95)
turquoise_nncon_bgm_t_bgm <- ggplot(turquoisetable_nncon_bgm_t_bgm, aes(x=con, y=sqrt(BGM_t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus grandii (sex)",  title = paste0("total     pval: ", round(turquoisespearman_nncon_bgm_t_bgm$p.value,digits = 3), "     rho: ", round(turquoisespearman_nncon_bgm_t_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
totspearman_nncon_bgm_t_bgm <- cor.test(tottable_nncon_bgm_t_bgm$con , tottable_nncon_bgm_t_bgm$BGM_t, method = "spearman", continuity = FALSE, conf.level = 0.95)
tot_nncon_bgm_t_bgm <- ggplot(tottable_nncon_bgm_t_bgm, aes(x=con, y=sqrt(BGM_t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus grandii (sex)",  title = paste0("total     pval: ", round(totspearman_nncon_bgm_t_bgm$p.value,digits = 3), "     rho: ", round(totspearman_nncon_bgm_t_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

networkspearman_nncon_bgm_t_bgm <- cor.test(networktable_nncon_bgm_t_bgm$con , networktable_nncon_bgm_t_bgm$BGM_t, method = "spearman", continuity = FALSE, conf.level = 0.95)
network_nncon_bgm_t_bgm <- ggplot(networktable_nncon_bgm_t_bgm, aes(x=con, y=sqrt(BGM_t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus grandii (sex)",  title = paste0("network     pval: ", round(networkspearman_nncon_bgm_t_bgm$p.value,digits = 3), "     rho: ", round(networkspearman_nncon_bgm_t_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
nonnetworkspearman_nncon_bgm_t_bgm <- cor.test(nonnetworktable_nncon_bgm_t_bgm$con , nonnetworktable_nncon_bgm_t_bgm$BGM_t, method = "spearman", continuity = FALSE, conf.level = 0.95)
nonnetwork_nncon_bgm_t_bgm <- ggplot(nonnetworktable_nncon_bgm_t_bgm, aes(x=con, y=sqrt(BGM_t)))     + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="#f55a4f", se=F, size=1.5, alpha=0.8) +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "t Bacillus grandii (sex)",  title = paste0("non-network     pval: ", round(nonnetworkspearman_nncon_bgm_t_bgm$p.value,digits = 3), "     rho: ", round(nonnetworkspearman_nncon_bgm_t_bgm$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))

grid.arrange(red_nncon_bgm_t_bat, blue_nncon_bgm_t_bat, black_nncon_bgm_t_bat, turquoise_nncon_bgm_t_bat, tot_nncon_bgm_t_bat,
             red_nncon_bgm_t_bro, blue_nncon_bgm_t_bro, black_nncon_bgm_t_bro, turquoise_nncon_bgm_t_bro, tot_nncon_bgm_t_bro,
             red_nncon_bgm_t_bgm, blue_nncon_bgm_t_bgm, black_nncon_bgm_t_bgm, turquoise_nncon_bgm_t_bgm, tot_nncon_bgm_t_bgm,
             ncol = 5 , nrow = 3)

grid.arrange(network_nncon_bgm_t_bgm, nonnetwork_nncon_bgm_t_bgm,
             network_nncon_bgm_t_bro, nonnetwork_nncon_bgm_t_bro,
             network_nncon_bgm_t_bat, nonnetwork_nncon_bgm_t_bat,
             ncol = 2 , nrow = 3)
                 #defplot
############################################################################################################################################  n connectivity VS LogFC ############################

nncon_bgm <- nearestNeighborConnectivity(datExpr0, nNeighbors = 500, power = 19, type = "signed", corFnc = "bicor", corOptions = "use = 'p'", blockSize = 6000)
names_bgm <- as.vector(names(datExpr0))
nncon_bgm_name <- "con"
names_bgm_name <- "OG"
nncon_bgm_df <- data.frame(names_bgm,nncon_bgm)
colnames(nncon_bgm_df) <- c(names_bgm_name, nncon_bgm_name)

nnconectivity_bgm_LogFC_bat <- merge(nncon_bgm_df, formatted_con_dnds_BAT[, c("OG", "BAT_logFC")], by="OG")
redtable_nncon_bgm_LogFC_bat <- subset(nnconectivity_bgm_LogFC_bat, nnconectivity_bgm_dNdS_bat$OG %in% red_genes)
bluetable_nncon_bgm_LogFC_bat <- subset(nnconectivity_bgm_LogFC_bat, nnconectivity_bgm_dNdS_bat$OG %in% blue_genes)
blacktable_nncon_bgm_LogFC_bat <- subset(nnconectivity_bgm_LogFC_bat, nnconectivity_bgm_dNdS_bat$OG %in% black_genes)
turquoisetable_nncon_bgm_LogFC_bat <- subset(nnconectivity_bgm_LogFC_bat, nnconectivity_bgm_dNdS_bat$OG %in% turquoise_genes)

redspearman_nncon_bgm_LogFC_bat <- cor.test(redtable_nncon_bgm_LogFC_bat$BAT_logFC, redtable_nncon_bgm_LogFC_bat$con, method = "spearman", continuity = FALSE, conf.level = 0.95)
red_nncon_bgm_LogFC_bat <- ggplot(redtable_nncon_bgm_LogFC_bat, aes(x=con, y=BAT_logFC))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "LogFC Bacillus atticus (asex)",  title = paste0("red     pval: ", round(redspearman_nncon_bgm_LogFC_bat$p.value,digits = 3), "     rho: ", round(redspearman_nncon_bgm_LogFC_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))
bluespearman_nncon_bgm_LogFC_bat <- cor.test(bluetable_nncon_bgm_LogFC_bat$BAT_logFC, bluetable_nncon_bgm_LogFC_bat$con, method = "spearman", continuity = FALSE, conf.level = 0.95)
blue_nncon_bgm_LogFC_bat <- ggplot(bluetable_nncon_bgm_LogFC_bat, aes(x=con, y=BAT_logFC))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "LogFC Bacillus atticus (asex)",  title = paste0("blue     pval: ", round(bluespearman_nncon_bgm_LogFC_bat$p.value,digits = 3), "     rho: ", round(bluespearman_nncon_bgm_LogFC_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))
blackspearman_nncon_bgm_LogFC_bat <- cor.test(blacktable_nncon_bgm_LogFC_bat$BAT_logFC, blacktable_nncon_bgm_LogFC_bat$con, method = "spearman", continuity = FALSE, conf.level = 0.95)
black_nncon_bgm_LogFC_bat <- ggplot(blacktable_nncon_bgm_LogFC_bat, aes(x=con, y=BAT_logFC))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "LogFC Bacillus atticus (asex)",  title = paste0("black     pval: ", round(blackspearman_nncon_bgm_LogFC_bat$p.value,digits = 3), "     rho: ", round(blackspearman_nncon_bgm_LogFC_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))
turquoisespearman_nncon_bgm_LogFC_bat <- cor.test(turquoisetable_nncon_bgm_LogFC_bat$BAT_logFC, turquoisetable_nncon_bgm_LogFC_bat$con, method = "spearman", continuity = FALSE, conf.level = 0.95)
turquoise_nncon_bgm_LogFC_bat <- ggplot(turquoisetable_nncon_bgm_LogFC_bat, aes(x=con, y=BAT_logFC))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "LogFC Bacillus atticus (asex)",  title = paste0("turquoise     pval: ", round(turquoisespearman_nncon_bgm_LogFC_bat$p.value,digits = 3), "     rho: ", round(turquoisespearman_nncon_bgm_LogFC_bat$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))

nnconectivity_bgm_LogFC_bro <- merge(nncon_bgm_df, formatted_con_dnds_BRO[, c("OG", "BRO_logFC")], by="OG")
redtable_nncon_bgm_LogFC_bro <- subset(nnconectivity_bgm_LogFC_bro, nnconectivity_bgm_dNdS_bro$OG %in% red_genes)
bluetable_nncon_bgm_LogFC_bro <- subset(nnconectivity_bgm_LogFC_bro, nnconectivity_bgm_dNdS_bro$OG %in% blue_genes)
blacktable_nncon_bgm_LogFC_bro <- subset(nnconectivity_bgm_LogFC_bro, nnconectivity_bgm_dNdS_bro$OG %in% black_genes)
turquoisetable_nncon_bgm_LogFC_bro <- subset(nnconectivity_bgm_LogFC_bro, nnconectivity_bgm_dNdS_bro$OG %in% turquoise_genes)

redspearman_nncon_bgm_LogFC_bro <- cor.test(redtable_nncon_bgm_LogFC_bro$BRO_logFC, redtable_nncon_bgm_LogFC_bro$con, method = "spearman", continuity = FALSE, conf.level = 0.95)
red_nncon_bgm_LogFC_bro <- ggplot(redtable_nncon_bgm_LogFC_bro, aes(x=con, y=BRO_logFC))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "LogFC Bacillus rossius (asex)",  title = paste0("red     pval: ", round(redspearman_nncon_bgm_LogFC_bro$p.value,digits = 3), "     rho: ", round(redspearman_nncon_bgm_LogFC_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))
bluespearman_nncon_bgm_LogFC_bro <- cor.test(bluetable_nncon_bgm_LogFC_bro$BRO_logFC, bluetable_nncon_bgm_LogFC_bro$con, method = "spearman", continuity = FALSE, conf.level = 0.95)
blue_nncon_bgm_LogFC_bro <- ggplot(bluetable_nncon_bgm_LogFC_bro, aes(x=con, y=BRO_logFC))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "LogFC Bacillus rossius (asex)",  title = paste0("blue     pval: ", round(bluespearman_nncon_bgm_LogFC_bro$p.value,digits = 3), "     rho: ", round(bluespearman_nncon_bgm_LogFC_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))
blackspearman_nncon_bgm_LogFC_bro <- cor.test(blacktable_nncon_bgm_LogFC_bro$BRO_logFC, blacktable_nncon_bgm_LogFC_bro$con, method = "spearman", continuity = FALSE, conf.level = 0.95)
black_nncon_bgm_LogFC_bro <- ggplot(blacktable_nncon_bgm_LogFC_bro, aes(x=con, y=BRO_logFC))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "LogFC Bacillus rossius (asex)",  title = paste0("black     pval: ", round(blackspearman_nncon_bgm_LogFC_bro$p.value,digits = 3), "     rho: ", round(blackspearman_nncon_bgm_LogFC_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))
turquoisespearman_nncon_bgm_LogFC_bro <- cor.test(turquoisetable_nncon_bgm_LogFC_bro$BRO_logFC, turquoisetable_nncon_bgm_LogFC_bro$con, method = "spearman", continuity = FALSE, conf.level = 0.95)
turquoise_nncon_bgm_LogFC_bro <- ggplot(turquoisetable_nncon_bgm_LogFC_bro, aes(x=con, y=BRO_logFC))    + geom_point(shape=16, color="black", size = 1) + geom_smooth(method=lm, color="black") +
  labs(x="nn connectivity Bacillus grandii (sex)", y = "LogFC Bacillus rossius (asex)",  title = paste0("turquoise     pval: ", round(turquoisespearman_nncon_bgm_LogFC_bro$p.value,digits = 3), "     rho: ", round(turquoisespearman_nncon_bgm_LogFC_bro$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))

grid.arrange(red_nncon_bgm_LogFC_bat, blue_nncon_bgm_LogFC_bat, black_nncon_bgm_LogFC_bat, turquoise_nncon_bgm_LogFC_bat,
             red_nncon_bgm_LogFC_bro, blue_nncon_bgm_LogFC_bro, black_nncon_bgm_LogFC_bro, turquoise_nncon_bgm_LogFC_bro, ncol = 4 , nrow = 2)

############################################################################################################################################  sex-asex module-nonmodule dNdS ############################

dNdS_tot <- read.table(file = "./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/4sp/4sp_pur_sel.tab", sep = " ", header= TRUE)

module_dNdS <-  subset(dNdS_tot, nnconectivity_bgm_dNdS_bro$OG %in% network_genes)
module_dNdS_bgm <- subset(module_dNdS, BGM_dS>0.0001 & BGM_dNdS<1)
module_dNdS_bgm <- module_dNdS_bgm$BGM_dNdS
module_dNdS_bgm <- data.frame(group = "Bacillus grandii", value = module_dNdS_bgm)
module_dNdS_bat <- subset(module_dNdS, BAT_dS>0.0001 & BAT_dNdS<1)
module_dNdS_bat <- module_dNdS_bat$BAT_dNdS
module_dNdS_bat <- data.frame(group = "Bacillus atticus mocule", value = module_dNdS_bat)
module_dNdS_bro <- subset(module_dNdS, BRO_dS>0.0001 & BRO_dNdS<1)
module_dNdS_bro <- module_dNdS_bro$BRO_dNdS
module_dNdS_bro <- data.frame(group = "Bacillus rossius", value = module_dNdS_bro)
module_dNdS_plot.data <- rbind(module_dNdS_bgm,module_dNdS_bat,module_dNdS_bro)
module_dNdS_plot <- ggplot(module_dNdS_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "dNdS",  title = "genes in modules associated to male gonads")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"), title=element_text(size=16,face="bold"))  + ylim(0, 1)
module_dNdS_plot <- module_dNdS_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

nonmodule_dNdS <-  subset(dNdS_tot, !(nnconectivity_bgm_dNdS_bro$OG %in% network_genes))
nonmodule_dNdS_bgm <- subset(nonmodule_dNdS, BGM_dS>0.0001 & BGM_dNdS<1)
nonmodule_dNdS_bgm <- nonmodule_dNdS_bgm$BGM_dNdS
nonmodule_dNdS_bgm <- data.frame(group = "Bacillus grandii", value = nonmodule_dNdS_bgm)
nonmodule_dNdS_bat <- subset(nonmodule_dNdS, BAT_dS>0.0001 & BAT_dNdS<1)
nonmodule_dNdS_bat <- nonmodule_dNdS_bat$BAT_dNdS
nonmodule_dNdS_bat <- data.frame(group = "Bacillus atticus", value = nonmodule_dNdS_bat)
nonmodule_dNdS_bro <- subset(nonmodule_dNdS, BRO_dS>0.0001 & BRO_dNdS<1)
nonmodule_dNdS_bro <- nonmodule_dNdS_bro$BRO_dNdS
nonmodule_dNdS_bro <- data.frame(group = "Bacillus rossius", value = nonmodule_dNdS_bro)

nonmodule_dNdS_plot.data <- rbind(nonmodule_dNdS_bgm,nonmodule_dNdS_bat,nonmodule_dNdS_bro)
nonmodule_dNdS_plot <- ggplot(nonmodule_dNdS_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "dNdS",  title = "genes out of modules associated to male gonads")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"), title=element_text(size=16,face="bold"))  + ylim(0, 1)
nonmodule_dNdS_plot <- nonmodule_dNdS_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

ks.test(module_dNdS_bgm$value, module_dNdS_bat$value)
ks.test(module_dNdS_bgm$value, module_dNdS_bro$value)
ks.test(module_dNdS_bro$value, module_dNdS_bat$value)

ks.test(nonmodule_dNdS_bgm$value, nonmodule_dNdS_bat$value)
ks.test(nonmodule_dNdS_bgm$value, nonmodule_dNdS_bro$value)
ks.test(nonmodule_dNdS_bro$value, nonmodule_dNdS_bat$value)

ks.test(module_dNdS_bgm$value, nonmodule_dNdS_bgm$value)
ks.test(module_dNdS_bro$value, nonmodule_dNdS_bro$value)
ks.test(module_dNdS_bat$value, nonmodule_dNdS_bat$value)

grid.arrange(module_dNdS_plot, nonmodule_dNdS_plot, ncol = 2, nrow = 1)

############################################################################################################################################  sex-asex module-nonmodule t ############################

t_tot <- read.table(file = "./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/4sp/4sp_pur_sel.tab", sep = " ", header= TRUE)

module_t <-  subset(t_tot, t_tot$OG %in% network_genes)
module_t_bgm <- subset(module_t, BGM_t<10)
module_t_bgm <- module_t_bgm$BGM_t
module_t_bgm <- data.frame(group = "Bacillus grandii", value = module_t_bgm)
module_t_bat <- subset(module_t, BAT_t<10)
module_t_bat <- module_t_bat$BAT_t
module_t_bat <- data.frame(group = "Bacillus atticus", value = module_t_bat)
module_t_bro <- subset(module_t, BRO_t<10)
module_t_bro <- module_t_bro$BRO_t
module_t_bro <- data.frame(group = "Bacillus rossius", value = module_t_bro)
module_t_plot.data <- rbind(module_t_bgm,module_t_bat,module_t_bro)
module_t_plot <- ggplot(module_t_plot.data, aes(x=group, y=sqrt(value), fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "t",  title = "genes in modules associated to male gonads")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"), title=element_text(size=16,face="bold"))  + ylim(0, 1)
module_t_plot <- module_t_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

nonmodule_t <-  subset(t_tot, !(t_tot$OG %in% network_genes))
nonmodule_t_bgm <- subset(nonmodule_t, BGM_t<10)
nonmodule_t_bgm <- nonmodule_t_bgm$BGM_t
nonmodule_t_bgm <- data.frame(group = "Bacillus grandii", value = nonmodule_t_bgm)
nonmodule_t_bat <- subset(nonmodule_t, BAT_t<10)
nonmodule_t_bat <- nonmodule_t_bat$BAT_t
nonmodule_t_bat <- data.frame(group = "Bacillus atticus", value = nonmodule_t_bat)
nonmodule_t_bro <- subset(nonmodule_t, BRO_t<10)
nonmodule_t_bro <- nonmodule_t_bro$BRO_t
nonmodule_t_bro <- data.frame(group = "Bacillus rossius", value = nonmodule_t_bro)

nonmodule_t_plot.data <- rbind(nonmodule_t_bgm,nonmodule_t_bat,nonmodule_t_bro)
nonmodule_t_plot <- ggplot(nonmodule_t_plot.data, aes(x=group, y=sqrt(value), fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "t",  title = "genes out of modules associated to male gonads")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"), title=element_text(size=16,face="bold"))  + ylim(0, 1)
nonmodule_t_plot <- nonmodule_t_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

ks.test(module_t_bgm$value, module_t_bat$value)
ks.test(module_t_bgm$value, module_t_bro$value)
ks.test(module_t_bro$value, module_t_bat$value)

ks.test(nonmodule_t_bgm$value, nonmodule_t_bat$value)
ks.test(nonmodule_t_bgm$value, nonmodule_t_bro$value)
ks.test(nonmodule_t_bro$value, nonmodule_t_bat$value)

ks.test(module_t_bgm$value, nonmodule_t_bgm$value)
ks.test(module_t_bro$value, nonmodule_t_bro$value)
ks.test(module_t_bat$value, nonmodule_t_bat$value)


grid.arrange(module_t_plot, nonmodule_t_plot, ncol = 2, nrow = 1)

############################################################################################################################################  sex-asex module-nonmodule CUB ############################

CUB_tot <- read.table(file = "./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/4sp/def.tab", sep = " ", header= TRUE)

module_CUB <-  subset(CUB_tot, CUB_tot$OG %in% network_genes)
module_CUB_bgm <- module_CUB$ENC_BGM
module_CUB_bgm <- data.frame(group = "Bacillus grandii", value = module_CUB_bgm)
module_CUB_bat <- module_CUB$ENC_BAT
module_CUB_bat <- data.frame(group = "Bacillus atticus", value = module_CUB_bat)
module_CUB_bro <- module_CUB$ENC_BRO
module_CUB_bro <- data.frame(group = "Bacillus rossius", value = module_CUB_bro)
module_CUB_plot.data <- rbind(module_CUB_bgm,module_CUB_bat,module_CUB_bro)
module_CUB_plot <- ggplot(module_CUB_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "CUB",  title = "genes in modules associated to male gonads")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"), title=element_text(size=16,face="bold"))
module_CUB_plot <- module_CUB_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

nonmodule_CUB <-  subset(CUB_tot, !(CUB_tot$OG %in% network_genes))
nonmodule_CUB_bgm <- nonmodule_CUB$ENC_BGM
nonmodule_CUB_bgm <- data.frame(group = "Bacillus grandii", value = nonmodule_CUB_bgm)
nonmodule_CUB_bat <- nonmodule_CUB$ENC_BAT
nonmodule_CUB_bat <- data.frame(group = "Bacillus atticus", value = nonmodule_CUB_bat)
nonmodule_CUB_bro <- nonmodule_CUB$ENC_BRO
nonmodule_CUB_bro <- data.frame(group = "Bacillus rossius", value = nonmodule_CUB_bro)

nonmodule_CUB_plot.data <- rbind(nonmodule_CUB_bgm,nonmodule_CUB_bat,nonmodule_CUB_bro)
nonmodule_CUB_plot <- ggplot(nonmodule_CUB_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "CUB",  title = "genes out of modules associated to male gonads")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"), title=element_text(size=16,face="bold"))
nonmodule_CUB_plot <- nonmodule_CUB_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

ks.test(module_CUB_bgm$value, module_CUB_bat$value)
ks.test(module_CUB_bgm$value, module_CUB_bro$value)
ks.test(module_CUB_bro$value, module_CUB_bat$value)

ks.test(nonmodule_CUB_bgm$value, nonmodule_CUB_bat$value)
ks.test(nonmodule_CUB_bgm$value, nonmodule_CUB_bro$value)
ks.test(nonmodule_CUB_bro$value, nonmodule_CUB_bat$value)

ks.test(module_CUB_bgm$value, nonmodule_CUB_bgm$value)
ks.test(module_CUB_bro$value, nonmodule_CUB_bro$value)
ks.test(module_CUB_bat$value, nonmodule_CUB_bat$value)


grid.arrange(module_CUB_plot, nonmodule_CUB_plot, ncol = 2, nrow = 1)


############################################################################################################################################  sex-asex per-module CUB #############################

tot <- read.table(file = "./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/4sp/def.tab", sep = " ", header= TRUE)

red_tot <-  subset(tot, tot$OG %in% red_genes)
red_MILC_bgm <- red_tot$MILC_BGM
red_MILC_bgm <- data.frame(group = "Bacillus grandii", value = red_MILC_bgm)
red_MILC_bro <- red_tot$MILC_BRO
red_MILC_bro <- data.frame(group = "Bacillus rossius", value = red_MILC_bro)
red_MILC_bat <- red_tot$MILC_BAT
red_MILC_bat <- data.frame(group = "Bacillus atticus", value = red_MILC_bat)
red_MILC_plot.data <- rbind(red_MILC_bgm,red_MILC_bat,red_MILC_bro)
red_MILC_plot <- ggplot(red_MILC_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "CUB - MILC",  title = "red module")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) + ylim(0, 1)
red_MILC_plot <- red_MILC_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

blue_tot <-  subset(tot, tot$OG %in% blue_genes)
blue_MILC_bgm <- blue_tot$MILC_BGM
blue_MILC_bgm <- data.frame(group = "Bacillus grandii", value = blue_MILC_bgm)
blue_MILC_bro <- blue_tot$MILC_BRO
blue_MILC_bro <- data.frame(group = "Bacillus rossius", value = blue_MILC_bro)
blue_MILC_bat <- blue_tot$MILC_BAT
blue_MILC_bat <- data.frame(group = "Bacillus atticus", value = blue_MILC_bat)
blue_MILC_plot.data <- rbind(blue_MILC_bgm,blue_MILC_bat,blue_MILC_bro)
blue_MILC_plot <- ggplot(blue_MILC_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "CUB - MILC",  title = "blue module")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) + ylim(0, 1)
blue_MILC_plot <- blue_MILC_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

black_tot <-  subset(tot, tot$OG %in% black_genes)
black_MILC_bgm <- black_tot$MILC_BGM
black_MILC_bgm <- data.frame(group = "Bacillus grandii", value = black_MILC_bgm)
black_MILC_bro <- black_tot$MILC_BRO
black_MILC_bro <- data.frame(group = "Bacillus rossius", value = black_MILC_bro)
black_MILC_bat <- black_tot$MILC_BAT
black_MILC_bat <- data.frame(group = "Bacillus atticus", value = black_MILC_bat)
black_MILC_plot.data <- rbind(black_MILC_bgm,black_MILC_bat,black_MILC_bro)
black_MILC_plot <- ggplot(black_MILC_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "CUB - MILC",  title = "black module")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) + ylim(0, 1)
black_MILC_plot <- black_MILC_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

turquoise_tot <-  subset(tot, tot$OG %in% turquoise_genes)
turquoise_MILC_bgm <- turquoise_tot$MILC_BGM
turquoise_MILC_bgm <- data.frame(group = "Bacillus grandii", value = turquoise_MILC_bgm)
turquoise_MILC_bro <- turquoise_tot$MILC_BRO
turquoise_MILC_bro <- data.frame(group = "Bacillus rossius", value = turquoise_MILC_bro)
turquoise_MILC_bat <- turquoise_tot$MILC_BAT
turquoise_MILC_bat <- data.frame(group = "Bacillus atticus", value = turquoise_MILC_bat)
turquoise_MILC_plot.data <- rbind(turquoise_MILC_bgm,turquoise_MILC_bat,turquoise_MILC_bro)
turquoise_MILC_plot <- ggplot(turquoise_MILC_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "CUB - MILC",  title = "turquoise module")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) + ylim(0, 1)
turquoise_MILC_plot <- turquoise_MILC_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

MILC_bgm <- tot$MILC_BGM
MILC_bgm <- data.frame(group = "Bacillus grandii", value = MILC_bgm)
MILC_bro <- tot$MILC_BRO
MILC_bro <- data.frame(group = "Bacillus rossius", value = MILC_bro)
MILC_bat <- tot$MILC_BAT
MILC_bat <- data.frame(group = "Bacillus atticus", value = MILC_bat)
MILC_plot.data <- rbind(MILC_bgm,MILC_bat,MILC_bro)
MILC_plot <- ggplot(MILC_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "CUB - MILC",  title = "total")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) + ylim(0, 1)
MILC_plot <- MILC_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

grid.arrange(black_MILC_plot, blue_MILC_plot, turquoise_MILC_plot, red_MILC_plot, MILC_plot, ncol = 5 , nrow = 1)

############################################################################################################################################  sex-asex per-module t ############################

red_t <-  subset(dNdS_tot, dNdS_tot$OG %in% red_genes)
red_t_bgm <- subset(red_t, BGM_t<0.6)
red_t_bgm <- red_t_bgm$BGM_t
red_t_bgm <- data.frame(group = "Bacillus grandii", value = red_t_bgm)

red_t <-  subset(dNdS_tot, dNdS_tot$OG %in% red_genes)
red_t_bro <- subset(red_t, BRO_t<0.6)
red_t_bro <- red_t_bro$BRO_t
red_t_bro <- data.frame(group = "Bacillus rossius", value = red_t_bro)

red_t <-  subset(dNdS_tot, dNdS_tot$OG %in% red_genes)
red_t_bat <- subset(red_t, BAT_t<0.6)
red_t_bat <- red_t_bat$BAT_t
red_t_bat <- data.frame(group = "Bacillus atticus", value = red_t_bat)

red_t_plot.data <- rbind(red_t_bgm,red_t_bro,red_t_bat)
red_t_plot <- ggplot(red_t_plot.data, aes(x=group, y=sqrt(value), fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "subs. rate",  title = "red")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) + ylim(0, 1)
red_t_plot <- red_t_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

blue_t <-  subset(dNdS_tot, dNdS_tot$OG %in% blue_genes)
blue_t_bgm <- subset(blue_t, BGM_t<0.6)
blue_t_bgm <- blue_t_bgm$BGM_t
blue_t_bgm <- data.frame(group = "Bacillus grandii", value = blue_t_bgm)

blue_t <-  subset(dNdS_tot, dNdS_tot$OG %in% blue_genes)
blue_t_bro <- subset(blue_t, BRO_t<0.6)
blue_t_bro <- blue_t_bro$BRO_t
blue_t_bro <- data.frame(group = "Bacillus rossius", value = blue_t_bro)

blue_t <-  subset(dNdS_tot, dNdS_tot$OG %in% blue_genes)
blue_t_bat <- subset(blue_t, BAT_t<0.6)
blue_t_bat <- blue_t_bat$BAT_t
blue_t_bat <- data.frame(group = "Bacillus atticus", value = blue_t_bat)

blue_t_plot.data <- rbind(blue_t_bgm,blue_t_bro,blue_t_bat)
blue_t_plot <- ggplot(blue_t_plot.data, aes(x=group, y=sqrt(value), fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "subs. rate",  title = "blue")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) + ylim(0, 1)
blue_t_plot <- blue_t_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")


black_t <-  subset(dNdS_tot, dNdS_tot$OG %in% black_genes)
black_t_bgm <- subset(black_t, BGM_t<0.6)
black_t_bgm <- black_t_bgm$BGM_t
black_t_bgm <- data.frame(group = "Bacillus grandii", value = black_t_bgm)

black_t <-  subset(dNdS_tot, dNdS_tot$OG %in% black_genes)
black_t_bro <- subset(black_t, BRO_t<0.6)
black_t_bro <- black_t_bro$BRO_t
black_t_bro <- data.frame(group = "Bacillus rossius", value = black_t_bro)

black_t <-  subset(dNdS_tot, dNdS_tot$OG %in% black_genes)
black_t_bat <- subset(black_t, BAT_t<0.6)
black_t_bat <- black_t_bat$BAT_t
black_t_bat <- data.frame(group = "Bacillus atticus", value = black_t_bat)

black_t_plot.data <- rbind(black_t_bgm,black_t_bro,black_t_bat)
black_t_plot <- ggplot(black_t_plot.data, aes(x=group, y=sqrt(value), fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "subs. rate",  title = "black")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) + ylim(0, 1)
black_t_plot <- black_t_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

turquoise_t <-  subset(dNdS_tot, dNdS_tot$OG %in% turquoise_genes)
turquoise_t_bgm <- subset(turquoise_t, BGM_t<0.6)
turquoise_t_bgm <- turquoise_t_bgm$BGM_t
turquoise_t_bgm <- data.frame(group = "Bacillus grandii", value = turquoise_t_bgm)

turquoise_t <-  subset(dNdS_tot, dNdS_tot$OG %in% turquoise_genes)
turquoise_t_bro <- subset(turquoise_t, BRO_t<0.6)
turquoise_t_bro <- turquoise_t_bro$BRO_t
turquoise_t_bro <- data.frame(group = "Bacillus rossius", value = turquoise_t_bro)

turquoise_t <-  subset(dNdS_tot, dNdS_tot$OG %in% turquoise_genes)
turquoise_t_bat <- subset(turquoise_t, BAT_t<0.6)
turquoise_t_bat <- turquoise_t_bat$BAT_t
turquoise_t_bat <- data.frame(group = "Bacillus atticus", value = turquoise_t_bat)

turquoise_t_plot.data <- rbind(turquoise_t_bgm,turquoise_t_bro,turquoise_t_bat)
turquoise_t_plot <- ggplot(turquoise_t_plot.data, aes(x=group, y=sqrt(value), fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "subs. rate",  title = "turquoise")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) + ylim(0, 1)
turquoise_t_plot <- turquoise_t_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

t_bgm <- subset(tot, BGM_t<0.6)
t_bgm <- t_bgm$BGM_t
t_bgm <- data.frame(group = "Bacillus grandii", value = t_bgm)

t_bro <- subset(tot, BRO_t<0.6)
t_bro <- t_bro$BRO_t
t_bro <- data.frame(group = "Bacillus rossius", value = t_bro)

t_bat <- subset(tot, BAT_t<0.6)
t_bat <- t_bat$BAT_t
t_bat <- data.frame(group = "Bacillus atticus", value = t_bat)

t_plot.data <- rbind(t_bgm,t_bro,t_bat)
t_plot <- ggplot(t_plot.data, aes(x=group, y=sqrt(value), fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "subs. rate",  title = "total")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) + ylim(0, 1)
t_plot <- t_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

grid.arrange(black_dNdS_plot, blue_dNdS_plot, turquoise_dNdS_plot, red_dNdS_plot, ncol = 4 , nrow = 1)

############################################################################################################################################  sex-asex per-module dNdS ############################

dNdS_tot <- read.table(file = "./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/4sp/4sp_pur_sel.tab", sep = " ", header= TRUE)

red_dNdS <-  subset(dNdS_tot, nnconectivity_bgm_dNdS_bro$OG %in% red_genes)
red_dNdS_bgm <- subset(red_dNdS, BGM_dS>0.0001 & BGM_dNdS<1)
red_dNdS_bgm <- red_dNdS_bgm$BGM_dNdS
red_dNdS_bgm <- data.frame(group = "Bacillus grandii", value = red_dNdS_bgm)
red_dNdS_bat <- subset(red_dNdS, BAT_dS>0.0001 & BAT_dNdS<1)
red_dNdS_bat <- red_dNdS_bat$BAT_dNdS
red_dNdS_bat <- data.frame(group = "Bacillus atticus", value = red_dNdS_bat)
red_dNdS_bro <- subset(red_dNdS, BRO_dS>0.0001 & BRO_dNdS<1)
red_dNdS_bro <- red_dNdS_bro$BRO_dNdS
red_dNdS_bro <- data.frame(group = "Bacillus rossius", value = red_dNdS_bro)
red_dNdS_plot.data <- rbind(red_dNdS_bgm,red_dNdS_bat,red_dNdS_bro)
red_dNdS_plot <- ggplot(red_dNdS_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "dNdS",  title = "red module")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"), title=element_text(size=16,face="bold"))  + ylim(0, 1)
red_dNdS_plot <- red_dNdS_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

blue_dNdS <-  subset(dNdS_tot, nnconectivity_bgm_dNdS_bro$OG %in% blue_genes)
blue_dNdS_bgm <- subset(blue_dNdS, BGM_dS>0.0001 & BGM_dNdS<1)
blue_dNdS_bgm <- blue_dNdS_bgm$BGM_dNdS
blue_dNdS_bgm <- data.frame(group = "Bacillus grandii", value = blue_dNdS_bgm)
blue_dNdS_bat <- subset(blue_dNdS, BAT_dS>0.0001 & BAT_dNdS<1)
blue_dNdS_bat <- blue_dNdS_bat$BAT_dNdS
blue_dNdS_bat <- data.frame(group = "Bacillus atticus", value = blue_dNdS_bat)
blue_dNdS_bro <- subset(blue_dNdS, BRO_dS>0.0001 & BRO_dNdS<1)
blue_dNdS_bro <- blue_dNdS_bro$BRO_dNdS
blue_dNdS_bro <- data.frame(group = "Bacillus rossius", value = blue_dNdS_bro)
blue_dNdS_plot.data <- rbind(blue_dNdS_bgm,blue_dNdS_bat,blue_dNdS_bro)
blue_dNdS_plot <- ggplot(blue_dNdS_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "dNdS",  title = "blue module")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"), title=element_text(size=16,face="bold")) + ylim(0, 1)
blue_dNdS_plot <- blue_dNdS_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")


black_dNdS <-  subset(dNdS_tot, nnconectivity_bgm_dNdS_bro$OG %in% black_genes)
black_dNdS_bgm <- subset(black_dNdS, BGM_dS>0.0001 & BGM_dNdS<1)
black_dNdS_bgm <- black_dNdS_bgm$BGM_dNdS
black_dNdS_bgm <- data.frame(group = "Bacillus grandii", value = black_dNdS_bgm)
black_dNdS_bat <- subset(black_dNdS, BAT_dS>0.0001 & BAT_dNdS<1)
black_dNdS_bat <- black_dNdS_bat$BAT_dNdS
black_dNdS_bat <- data.frame(group = "Bacillus atticus", value = black_dNdS_bat)
black_dNdS_bro <- subset(black_dNdS, BRO_dS>0.0001 & BRO_dNdS<1)
black_dNdS_bro <- black_dNdS_bro$BRO_dNdS
black_dNdS_bro <- data.frame(group = "Bacillus rossius", value = black_dNdS_bro)
black_dNdS_plot.data <- rbind(black_dNdS_bgm,black_dNdS_bat,black_dNdS_bro)
black_dNdS_plot <- ggplot(black_dNdS_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "dNdS",  title = "black module")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"), title=element_text(size=16,face="bold")) + ylim(0, 1)
black_dNdS_plot <- black_dNdS_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")


turquoise_dNdS <-  subset(dNdS_tot, nnconectivity_bgm_dNdS_bro$OG %in% turquoise_genes)
turquoise_dNdS_bgm <- subset(turquoise_dNdS, BGM_dS>0.0001 & BGM_dNdS<1)
turquoise_dNdS_bgm <- turquoise_dNdS_bgm$BGM_dNdS
turquoise_dNdS_bgm <- data.frame(group = "Bacillus grandii", value = turquoise_dNdS_bgm)
turquoise_dNdS_bat <- subset(turquoise_dNdS, BAT_dS>0.0001 & BAT_dNdS<1)
turquoise_dNdS_bat <- turquoise_dNdS_bat$BAT_dNdS
turquoise_dNdS_bat <- data.frame(group = "Bacillus atticus", value = turquoise_dNdS_bat)
turquoise_dNdS_bro <- subset(turquoise_dNdS, BRO_dS>0.0001 & BRO_dNdS<1)
turquoise_dNdS_bro <- turquoise_dNdS_bro$BRO_dNdS
turquoise_dNdS_bro <- data.frame(group = "Bacillus rossius", value = turquoise_dNdS_bro)
turquoise_dNdS_plot.data <- rbind(turquoise_dNdS_bgm,turquoise_dNdS_bat,turquoise_dNdS_bro)
turquoise_dNdS_plot <- ggplot(turquoise_dNdS_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "dNdS",  title = "turquoise module")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"), title=element_text(size=16,face="bold")) + ylim(0, 1)
turquoise_dNdS_plot <- turquoise_dNdS_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

dNdS_bgm <- dNdS_tot$BGM_dNdS
dNdS_bgm <- data.frame(group = "Bacillus grandii", value = dNdS_bgm)
dNdS_bro <- dNdS_tot$BRO_dNdS
dNdS_bro <- data.frame(group = "Bacillus rossius", value = dNdS_bro)
dNdS_bat <- tot$BAT_dNdS
dNdS_bat <- data.frame(group = "Bacillus atticus", value = dNdS_bat)
dNdS_plot.data <- rbind(dNdS_bgm,dNdS_bat,dNdS_bro)
dNdS_plot <- ggplot(dNdS_plot.data, aes(x=group, y=value, fill=group)) + geom_boxplot(width=0.3, outlier.size=0) + labs(x=" ", y = "CUB - MILC",  title = "total")+ theme_classic() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"), title=element_text(size=16,face="bold")) + ylim(0, 1)
dNdS_plot <- dNdS_plot + scale_fill_manual(values=c("#04b4d6", "#ffc037", "#f55a4f")) + theme(legend.position = "none")

grid.arrange(black_dNdS_plot, blue_dNdS_plot, turquoise_dNdS_plot, red_dNdS_plot, dNdS_plot, ncol = 5 , nrow = 1)


############################################################################################################################################  molecular evolution main fig ############################

grid.arrange(module_dNdS_plot, nonmodule_dNdS_plot,
             ncol = 2, nrow = 1)

############################################################################################################################################  k connectivity VS LogFC main fig #########################################

formatted_con_dnds_BAT <- read.table(file = "./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/4sp/connectivicty_BATF_dnds_LogFC.4sp.tab", sep = " ", header= TRUE)

formatted_con_dnds_BRO <- read.table(file = "./Desktop/network analysis/expression_matrixes/exp_matrixes_orthologs_only/4sp/connectivicty_BROF_dnds_LogFC.4sp.tab", sep = " ", header= TRUE)

formatted_con_dnds_BAT_modules <- subset(formatted_con_dnds_BAT, formatted_con_dnds_BAT$OG %in% network_genes)
formatted_con_dnds_BAT_modules_pos <- subset(formatted_con_dnds_BAT_modules, BAT_logFC>0)
formatted_con_dnds_BAT_modules_neg <- subset(formatted_con_dnds_BAT_modules, BAT_logFC<0)

formatted_con_dnds_BAT_nonmodules <- subset(formatted_con_dnds_BAT, !(formatted_con_dnds_BAT$OG %in% network_genes))
formatted_con_dnds_BAT_nonmodules_pos <- subset(formatted_con_dnds_BAT_nonmodules, BAT_logFC>0)
formatted_con_dnds_BAT_nonmodules_neg <- subset(formatted_con_dnds_BAT_nonmodules, BAT_logFC<0)

formatted_con_dnds_BRO_modules <- subset(formatted_con_dnds_BRO, formatted_con_dnds_BRO$OG %in% network_genes)
formatted_con_dnds_BRO_modules_pos <- subset(formatted_con_dnds_BRO_modules, BRO_logFC>0)
formatted_con_dnds_BRO_modules_neg <- subset(formatted_con_dnds_BRO_modules, BRO_logFC<0)

formatted_con_dnds_BRO_nonmodules <- subset(formatted_con_dnds_BRO, !(formatted_con_dnds_BRO$OG %in% network_genes))
formatted_con_dnds_BRO_nonmodules_pos <- subset(formatted_con_dnds_BRO_nonmodules, BRO_logFC>0)
formatted_con_dnds_BRO_nonmodules_neg <- subset(formatted_con_dnds_BRO_nonmodules, BRO_logFC<0)

# kWithin sex VS LogFC asex

# bgm VS bat

bgm_kWithin_bat_LogFC_spearman_modules <- cor.test(formatted_con_dnds_BAT_modules$kWithin , formatted_con_dnds_BAT_modules$BAT_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
bgm_kWithin_bat_LogFC_plot_modules <- ggplot() +
  geom_point(data=formatted_con_dnds_BAT_modules_neg, aes(kWithin, BAT_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=formatted_con_dnds_BAT_modules_pos, aes(kWithin, BAT_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="kwithin B. grandii (sex)", y = "LogFC B. atticus (asex)", title = paste0(
    "module LogFC     pval: ", round(bgm_kWithin_bat_LogFC_spearman_modules$p.value,digits = 3), "     rho: ", round(bgm_kWithin_bat_LogFC_spearman_modules$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=13,face="bold"))

bgm_kWithin_bat_LogFC_spearman_nonmodules <- cor.test(formatted_con_dnds_BAT_nonmodules$kWithin , formatted_con_dnds_BAT_nonmodules$BAT_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
bgm_kWithin_bat_LogFC_plot_nonmodules <- ggplot() +
  geom_point(data=formatted_con_dnds_BAT_nonmodules_neg, aes(kWithin, BAT_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=formatted_con_dnds_BAT_nonmodules_pos, aes(kWithin, BAT_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="kwithin B. grandii (sex)", y = "LogFC B. atticus (asex)", title = paste0(
    "non-modules LogFC     pval: ", round(bgm_kWithin_bat_LogFC_spearman_nonmodules$p.value,digits = 3), "     rho: ", round(bgm_kWithin_bat_LogFC_spearman_nonmodules$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=13,face="bold"))

# bgm VS bro

bgm_kWithin_bro_LogFC_spearman_modules <- cor.test(formatted_con_dnds_BRO_modules$kWithin , formatted_con_dnds_BRO_modules$BRO_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
bgm_kWithin_bro_LogFC_plot_modules <- ggplot() +
  geom_point(data=formatted_con_dnds_BRO_modules_neg, aes(kWithin, BRO_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=formatted_con_dnds_BRO_modules_pos, aes(kWithin, BRO_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="kwithin B. grandii (sex)", y = "LogFC B. rossius (asex)", title = paste0(
    "module LogFC     pval: ", round(bgm_kWithin_bro_LogFC_spearman_modules$p.value,digits = 3), "     rho: ", round(bgm_kWithin_bro_LogFC_spearman_modules$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=13,face="bold"))

bgm_kWithin_bro_LogFC_spearman_nonmodules <- cor.test(formatted_con_dnds_BRO_nonmodules$kWithin , formatted_con_dnds_BRO_nonmodules$BRO_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
bgm_kWithin_bro_LogFC_plot_nonmodules <- ggplot() +
  geom_point(data=formatted_con_dnds_BRO_nonmodules_neg, aes(kWithin, BRO_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=formatted_con_dnds_BRO_nonmodules_pos, aes(kWithin, BRO_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="kwithin B. grandii (sex)", y = "LogFC B. rossius (asex)", title = paste0(
    "non-modules LogFC     pval: ", round(bgm_kWithin_bro_LogFC_spearman_nonmodules$p.value,digits = 3), "     rho: ", round(bgm_kWithin_bro_LogFC_spearman_nonmodules$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=13,face="bold"))

# kOut sex VS LogFC asex

# bgm VS bat

bgm_kOut_bat_LogFC_spearman_modules <- cor.test(formatted_con_dnds_BAT_modules$kOut , formatted_con_dnds_BAT_modules$BAT_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
bgm_kOut_bat_LogFC_plot_modules <- ggplot() +
  geom_point(data=formatted_con_dnds_BAT_modules_neg, aes(kOut, BAT_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=formatted_con_dnds_BAT_modules_pos, aes(kOut, BAT_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="kOut B. grandii (sex)", y = "LogFC B. atticus (asex)", title = paste0(
    "module LogFC     pval: ", round(bgm_kOut_bat_LogFC_spearman_modules$p.value,digits = 3), "     rho: ", round(bgm_kOut_bat_LogFC_spearman_modules$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=13,face="bold"))

bgm_kOut_bat_LogFC_spearman_nonmodules <- cor.test(formatted_con_dnds_BAT_nonmodules$kOut , formatted_con_dnds_BAT_nonmodules$BAT_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
bgm_kOut_bat_LogFC_plot_nonmodules <- ggplot() +
  geom_point(data=formatted_con_dnds_BAT_nonmodules_neg, aes(kOut, BAT_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=formatted_con_dnds_BAT_nonmodules_pos, aes(kOut, BAT_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="kOut B. grandii (sex)", y = "LogFC B. atticus (asex)", title = paste0(
    "non-modules LogFC     pval: ", round(bgm_kOut_bat_LogFC_spearman_nonmodules$p.value,digits = 3), "     rho: ", round(bgm_kOut_bat_LogFC_spearman_nonmodules$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=13,face="bold"))

# bgm VS bro

bgm_kOut_bro_LogFC_spearman_modules <- cor.test(formatted_con_dnds_BRO_modules$kOut , formatted_con_dnds_BRO_modules$BRO_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
bgm_kOut_bro_LogFC_plot_modules <- ggplot() +
  geom_point(data=formatted_con_dnds_BRO_modules_neg, aes(kOut, BRO_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=formatted_con_dnds_BRO_modules_pos, aes(kOut, BRO_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="kOut B. grandii (sex)", y = "LogFC B. rossius (asex)", title = paste0(
    "module LogFC     pval: ", round(bgm_kOut_bro_LogFC_spearman_modules$p.value,digits = 3), "     rho: ", round(bgm_kOut_bro_LogFC_spearman_modules$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=13,face="bold"))

bgm_kOut_bro_LogFC_spearman_nonmodules <- cor.test(formatted_con_dnds_BRO_nonmodules$kOut , formatted_con_dnds_BRO_nonmodules$BRO_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
bgm_kOut_bro_LogFC_plot_nonmodules <- ggplot() +
  geom_point(data=formatted_con_dnds_BRO_nonmodules_neg, aes(kOut, BRO_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=formatted_con_dnds_BRO_nonmodules_pos, aes(kOut, BRO_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="kOut B. grandii (sex)", y = "LogFC B. rossius (asex)", title = paste0(
    "non-modules LogFC     pval: ", round(bgm_kOut_bro_LogFC_spearman_nonmodules$p.value,digits = 3), "     rho: ", round(bgm_kOut_bro_LogFC_spearman_nonmodules$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=13,face="bold"))

# kTotal sex VS LogFC asex

# bgm VS bat

bgm_kTotal_bat_LogFC_spearman_modules <- cor.test(formatted_con_dnds_BAT_modules$kTotal , formatted_con_dnds_BAT_modules$BAT_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
bgm_kTotal_bat_LogFC_plot_modules <- ggplot() +
  geom_point(data=formatted_con_dnds_BAT_modules_neg, aes(kTotal, BAT_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=formatted_con_dnds_BAT_modules_pos, aes(kTotal, BAT_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="kTotal B. grandii (sex)", y = "LogFC B. atticus (asex)", title = paste0(
    "module LogFC pval: ", round(bgm_kTotal_bat_LogFC_spearman_modules$p.value,digits = 3), " rho: ", round(bgm_kTotal_bat_LogFC_spearman_modules$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=13,face="bold"))

bgm_kTotal_bat_LogFC_spearman_nonmodules <- cor.test(formatted_con_dnds_BAT_nonmodules$kTotal , formatted_con_dnds_BAT_nonmodules$BAT_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
bgm_kTotal_bat_LogFC_plot_nonmodules <- ggplot() +
  geom_point(data=formatted_con_dnds_BAT_nonmodules_neg, aes(kTotal, BAT_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=formatted_con_dnds_BAT_nonmodules_pos, aes(kTotal, BAT_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="kTotal B. grandii (sex)", y = "LogFC B. atticus (asex)", title = paste0(
    "non-modules LogFC pval: ", round(bgm_kTotal_bat_LogFC_spearman_nonmodules$p.value,digits = 3), " rho: ", round(bgm_kTotal_bat_LogFC_spearman_nonmodules$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=13,face="bold"))

# bgm VS bro

bgm_kTotal_bro_LogFC_spearman_modules <- cor.test(formatted_con_dnds_BRO_modules$kTotal , formatted_con_dnds_BRO_modules$BRO_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
bgm_kTotal_bro_LogFC_plot_modules <- ggplot() +
  geom_point(data=formatted_con_dnds_BRO_modules_neg, aes(kTotal, BRO_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=formatted_con_dnds_BRO_modules_pos, aes(kTotal, BRO_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="kTotal B. grandii (sex)", y = "LogFC B. rossius (asex)", title = paste0(
    "module LogFC pval: ", round(bgm_kTotal_bro_LogFC_spearman_modules$p.value,digits = 3), " rho: ", round(bgm_kTotal_bro_LogFC_spearman_modules$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=13,face="bold"))

bgm_kTotal_bro_LogFC_spearman_nonmodules <- cor.test(formatted_con_dnds_BRO_nonmodules$kTotal , formatted_con_dnds_BRO_nonmodules$BRO_logFC, method = "spearman", continuity = FALSE, conf.level = 0.95)
bgm_kTotal_bro_LogFC_plot_nonmodules <- ggplot() +
  geom_point(data=formatted_con_dnds_BRO_nonmodules_neg, aes(kTotal, BRO_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=formatted_con_dnds_BRO_nonmodules_pos, aes(kTotal, BRO_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="kTotal B. grandii (sex)", y = "LogFC B. rossius (asex)", title = paste0(
    "non-modules LogFC pval: ", round(bgm_kTotal_bro_LogFC_spearman_nonmodules$p.value,digits = 3), " rho: ", round(bgm_kTotal_bro_LogFC_spearman_nonmodules$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=13,face="bold"))

grid.arrange(bgm_kTotal_bro_LogFC_plot_modules, bgm_kTotal_bro_LogFC_plot_nonmodules, bgm_kTotal_bat_LogFC_plot_modules, bgm_kTotal_bat_LogFC_plot_nonmodules,
             bgm_kOut_bro_LogFC_plot_modules, bgm_kOut_bro_LogFC_plot_nonmodules, bgm_kOut_bat_LogFC_plot_modules, bgm_kOut_bat_LogFC_plot_nonmodules,
             bgm_kWithin_bro_LogFC_plot_modules, bgm_kWithin_bro_LogFC_plot_nonmodules, bgm_kWithin_bat_LogFC_plot_modules, bgm_kWithin_bat_LogFC_plot_nonmodules,
             ncol = 4 , nrow = 3)

############################################################################################################################################  n connectivity VS LogFC BGM VS BGM (sanity check) ####################################

# females BGMVSBGM

nnconectivity_LogFC_bgm_f_pos <- subset(nnconectivity_LogFC_bgm_f, BGM_F_logFC>0)
nnconectivity_LogFC_bgm_f_neg <- subset(nnconectivity_LogFC_bgm_f, BGM_F_logFC<0)

nnconectivity_LogFC_bgm_spearman_f_pos <- cor.test(nnconectivity_LogFC_bgm_f_pos$con , abs(nnconectivity_LogFC_bgm_f_pos$BGM_F_logFC), method = "spearman", continuity = FALSE, conf.level = 0.95)
nnconectivity_LogFC_bgm_spearman_f_neg <- cor.test(nnconectivity_LogFC_bgm_f_neg$con , abs(nnconectivity_LogFC_bgm_f_neg$BGM_F_logFC), method = "spearman", continuity = FALSE, conf.level = 0.95)

nnconectivity_LogFC_bgm_f_plot_abs <- ggplot() +
  geom_point(data=nnconectivity_LogFC_bgm_f_neg, aes(con, BGM_F_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=nnconectivity_LogFC_bgm_f_pos, aes(con, BGM_F_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="nn connectivity Bacillus grandii females (sex)", y = "absolute LogFC Bacillus grandii females (sex)", title = paste0(
    "positive LogFC     pval: ", round(nnconectivity_LogFC_bgm_spearman_f_pos$p.value,digits = 3), "     rho: ", round(nnconectivity_LogFC_bgm_spearman_f_pos$estimate,digits = 3), 
    "          negatve LogFC     pval: ", round(nnconectivity_LogFC_bgm_spearman_f_neg$p.value,digits = 3), "     rho: ", round(nnconectivity_LogFC_bgm_spearman_f_neg$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))

# males BGMVSBGM

nnconectivity_LogFC_bgm_m_pos <- subset(nnconectivity_LogFC_bgm_m, BGM_M_logFC>0)
nnconectivity_LogFC_bgm_m_neg <- subset(nnconectivity_LogFC_bgm_m, BGM_M_logFC<0)

nnconectivity_LogFC_bgm_spearman_m_pos <- cor.test(nnconectivity_LogFC_bgm_m_pos$con , abs(nnconectivity_LogFC_bgm_m_pos$BGM_M_logFC), method = "spearman", continuity = FALSE, conf.level = 0.95)
nnconectivity_LogFC_bgm_spearman_m_neg <- cor.test(nnconectivity_LogFC_bgm_m_neg$con , abs(nnconectivity_LogFC_bgm_m_neg$BGM_M_logFC), method = "spearman", continuity = FALSE, conf.level = 0.95)

nnconectivity_LogFC_bgm_m_plot_abs <- ggplot() +
  geom_point(data=nnconectivity_LogFC_bgm_m_neg, aes(con, BGM_M_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) +
  geom_point(data=nnconectivity_LogFC_bgm_m_pos, aes(con, BGM_M_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) +
  labs(x="nn connectivity Bacillus grandii males (sex)", y = "absolute LogFC Bacillus grandii males (sex)", title = paste0(
    "positive LogFC     pval: ", round(nnconectivity_LogFC_bgm_spearman_m_pos$p.value,digits = 3), "     rho: ", round(nnconectivity_LogFC_bgm_spearman_m_pos$estimate,digits = 3), 
    "          negatve LogFC     pval: ", round(nnconectivity_LogFC_bgm_spearman_m_neg$p.value,digits = 3), "     rho: ", round(nnconectivity_LogFC_bgm_spearman_m_neg$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))

grid.arrange(nnconectivity_LogFC_bgm_f_plot_abs, nnconectivity_LogFC_bgm_m_plot_abs, ncol = 2 , nrow = 1)

############################################################################################################################################  n connectivity VS LogFC main fig ####################################

# females BGMVSBGM

nnconectivity_LogFC_bgm_f_pos <- subset(nnconectivity_LogFC_bgm_f, BGM_F_logFC>0)
nnconectivity_LogFC_bgm_f_neg <- subset(nnconectivity_LogFC_bgm_f, BGM_F_logFC<0)

nnconectivity_LogFC_bgm_spearman_f_pos <- cor.test(nnconectivity_LogFC_bgm_f_pos$con , abs(nnconectivity_LogFC_bgm_f_pos$BGM_F_logFC), method = "spearman", continuity = FALSE, conf.level = 0.95)
nnconectivity_LogFC_bgm_spearman_f_neg <- cor.test(nnconectivity_LogFC_bgm_f_neg$con , abs(nnconectivity_LogFC_bgm_f_neg$BGM_F_logFC), method = "spearman", continuity = FALSE, conf.level = 0.95)

nnconectivity_LogFC_bgm_f_plot_abs <- ggplot() +
  geom_point(data=nnconectivity_LogFC_bgm_f_neg, aes(con, BGM_F_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) + ylim(-12,5) +
  geom_point(data=nnconectivity_LogFC_bgm_f_pos, aes(con, BGM_F_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) + ylim(-12,5) +
  labs(x="nn connectivity Bacillus grandii females (sex)", y = "absolute LogFC Bacillus grandii females (sex)", title = paste0(
    "positive LogFC     pval: ", round(nnconectivity_LogFC_bgm_spearman_f_pos$p.value,digits = 3), "     rho: ", round(nnconectivity_LogFC_bgm_spearman_f_pos$estimate,digits = 3), 
    "          negatve LogFC     pval: ", round(nnconectivity_LogFC_bgm_spearman_f_neg$p.value,digits = 3), "     rho: ", round(nnconectivity_LogFC_bgm_spearman_f_neg$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))

# females BGMVSBAT

nnconectivity_bgm_LogFC_bat_pos <- subset(nnconectivity_bgm_LogFC_bat, BAT_logFC>0)
nnconectivity_bgm_LogFC_bat_neg <- subset(nnconectivity_bgm_LogFC_bat, BAT_logFC<0)

nnconectivity_bgm_LogFC_batspearman_pos <- cor.test(nnconectivity_bgm_LogFC_bat_pos$con , abs(nnconectivity_bgm_LogFC_bat_pos$BAT_logFC), method = "spearman", continuity = FALSE, conf.level = 0.95)
nnconectivity_bgm_LogFC_batspearman_neg <- cor.test(nnconectivity_bgm_LogFC_bat_neg$con , abs(nnconectivity_bgm_LogFC_bat_neg$BAT_logFC), method = "spearman", continuity = FALSE, conf.level = 0.95)

nnconectivity_bgm_LogFC_bat_f_plot_abs <- ggplot() +
  geom_point(data=nnconectivity_bgm_LogFC_bat_neg, aes(con, BAT_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) + ylim(-12,5) +
  geom_point(data=nnconectivity_bgm_LogFC_bat_pos, aes(con, BAT_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) + ylim(-12,5) +
  labs(x="nn connectivity Bacillus grandii females (sex)", y = "absolute LogFC Bacillus atticus females (sex)", title = paste0(
    "positive LogFC     pval: ", round(nnconectivity_bgm_LogFC_batspearman_pos$p.value,digits = 3), "     rho: ", round(nnconectivity_bgm_LogFC_batspearman_pos$estimate,digits = 3), 
    "          negatve LogFC     pval: ", round(nnconectivity_bgm_LogFC_batspearman_neg$p.value,digits = 3), "     rho: ", round(nnconectivity_bgm_LogFC_batspearman_neg$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))

# females BGMVSBRO

nnconectivity_bgm_LogFC_bro_pos <- subset(nnconectivity_bgm_LogFC_bro, BRO_logFC>0)
nnconectivity_bgm_LogFC_bro_neg <- subset(nnconectivity_bgm_LogFC_bro, BRO_logFC<0)

nnconectivity_bgm_LogFC_brospearman_pos <- cor.test(nnconectivity_bgm_LogFC_bro_pos$con , abs(nnconectivity_bgm_LogFC_bro_pos$BRO_logFC), method = "spearman", continuity = FALSE, conf.level = 0.95)
nnconectivity_bgm_LogFC_brospearman_neg <- cor.test(nnconectivity_bgm_LogFC_bro_neg$con , abs(nnconectivity_bgm_LogFC_bro_neg$BRO_logFC), method = "spearman", continuity = FALSE, conf.level = 0.95)

nnconectivity_bgm_LogFC_bro_f_plot_abs <- ggplot() +
  geom_point(data=nnconectivity_bgm_LogFC_bro_neg, aes(con, BRO_logFC), color="#00A4CCFF", size = 2, alpha = 0.8, shape= 20) + ylim(-12,5) +
  geom_point(data=nnconectivity_bgm_LogFC_bro_pos, aes(con, BRO_logFC), color="#F95700FF", size = 2, alpha = 0.8, shape= 20) + ylim(-12,5) +
  labs(x="nn connectivity Bacillus grandii females (sex)", y = "absolute LogFC Bacillus rossius females (sex)", title = paste0(
    "positive LogFC     pval: ", round(nnconectivity_bgm_LogFC_brospearman_pos$p.value,digits = 3), "     rho: ", round(nnconectivity_bgm_LogFC_brospearman_pos$estimate,digits = 3), 
    "          negatve LogFC     pval: ", round(nnconectivity_bgm_LogFC_brospearman_neg$p.value,digits = 3), "     rho: ", round(nnconectivity_bgm_LogFC_brospearman_neg$estimate,digits = 3))) + 
  theme_classic() + theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))

grid.arrange(nnconectivity_LogFC_bgm_f_plot_abs, nnconectivity_bgm_LogFC_bat_f_plot_abs, nnconectivity_bgm_LogFC_bro_f_plot_abs, ncol = 1 , nrow = 3)

############################################################################################################################################  connectivity VS padj BGMVSBGM (sanity check) ####################################

nnconectivity_padj_bgm_spearman_f_abs <- cor.test(nnconectivity_padj_bgm_f$con , abs(nnconectivity_padj_bgm_f$BGM_F_padj), method = "spearman", continuity = FALSE, conf.level = 0.95)
nnconectivity_padj_bgm_f_plot_abs <- ggplot() +
  geom_point(data=nnconectivity_padj_bgm_f, aes(con, BGM_F_padj), color="black", size = 2, alpha = 0.8, shape= 20) +
  geom_smooth(x=nnconectivity_padj_bgm_f$con, y=nnconectivity_padj_bgm_f$BGM_F_padj) +
  labs(x="nn connectivity Bacillus grandii females (sex)", y = " padj Bacillus grandii males (sex)",  title = paste0("     pval: ", round(nnconectivity_padj_bgm_spearman_f_abs$p.value,digits = 3), "     rho: ", round(nnconectivity_padj_bgm_spearman_f_abs$estimate,digits = 3)))+ theme_classic() +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=18,face="bold"))

