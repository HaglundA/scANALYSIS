

source("/rds/general/user/ah3918/ephemeral/SCANALYSIS/single_cell_function.R")


#Create matrix from cellranger output

dataMatrix.1=Load.Data.Matrices(list(SAMPLE=paste("/rds/general/user/ah3918/ephemeral/CELLRANGER_DATA/21FEB20Test1630/outs/raw_feature_bc_matrix",sep="")))
dataMatrixDgC.1=dataMatrix.1[[1]]




#QC - barcodeRanks: Compute barcode rank statistics and identifry the knee and inflection points on the total count curve.

brRanks.1 <- barcodeRanks(dataMatrixDgC.1)

png("QC_plot1_27FEB19")
options(repr.plot.width=5, repr.plot.height=4)
plot(brRanks.1$rank, brRanks.1$total, log="xy", xlab="Barcodes Ranked", ylab="Total UMI per Barcode", pch=19, cex=0.3, main="HS")
ordered.1 <- order(brRanks.1$rank)
abline(h=metadata(brRanks.1)$knee, col="dodgerblue", lty=2)
abline(h=metadata(brRanks.1)$inflection, col="forestgreen", lty=2)
legend("topright", lty=2, col=c("dodgerblue", "forestgreen"),legend=c("knee", "inflection"), cex=1, pt.cex=0.1)

dev.off()

#Select subset

knee.1=brRanks.1@metadata$knee
inflection.1=brRanks.1@metadata$inflection
paste(inflection.1, knee.1)



#EmptyDrops-------------------

set.seed(1234)
emptyDropsOutput.1= emptyDrops(dataMatrixDgC.1, lower=inflection.1,retain=knee.1)


#This ensures emptyDrops contains vectors of p.values for each barcode, and
#filters out cells with FDR<0.01

emptyDropsKeep.1 = emptyDropsOutput.1$FDR <= 0.01
emptyDropsKeep.na.1 = emptyDropsKeep.1
emptyDropsKeep.na.1[is.na(emptyDropsKeep.1)]= FALSE


#check that emptydrops detects cells - number of barcodes containing cells
table(emptyDropsKeep.na.1)

#--

counts.1=dataMatrixDgC.1[,emptyDropsKeep.na.1]
countsPerCell.1 = Matrix::colSums(counts.1)
countsGerGene.1 = Matrix::rowSums(counts.1)
genesPerCell.1 = Matrix::colSums(counts.1>0)


paste("Median Counts per Cell = ",median(countsPerCell.1))
paste("Median Genes per Cell = ",median(genesPerCell.1))


png("countspercellhist_27FEB20.png")
plot(hist(countsPerCell.1,breaks=42, col="powderblue", border="black", main="HS", xlab="Counts Per Cell", xlim=c(0, 25000), ylim=c(0, 17000), cex.lab=1, cex.axis=0.75, las=1))
dev.off()

png("genesspercellhist_27FEB20.png")
plot(hist(genesPerCell.1,breaks=45,col="powderblue", border="black",main="HS", xlab="Genes Per Cell", xlim=c(0, 6000), ylim=c(0, 10000), cex.lab=1, cex.axis=0.75, las=1))
dev.off()


png("cell_metrics_27FEB2020.png")
plot(countsPerCell.1, genesPerCell.1, main="HS", xlab="Counts Per Cell", ylab="Genes Per Cell",pch=19, cex=0.01,  cex.lab=1, cex.axis=0.75, las=1)
dev.off()


#CREATING A SEURAT OBJECT--------------------

#Seurat is used for quality control, analysis and exploration of scRNA-seq data.
# The Seurat object represents the expression data (@Assay) aswell as the dimensionality reductions
# of the expression data. The Assay objects are transformed with various dimensional reduction
# techniques to the DimReduc object.




countsMx.1=as.matrix(counts.1)

#In the CreateSeuratObject function, min.cells is the minimum number of cells a gene needs to
# be expressed in to be included.

HS.Sample= CreateSeuratObject(countsMx.1, min.cells=5, min.features=5, project = "Control")

#We now need to identify mitochondrial genes and remove them from the cells - low Quality
# cells or dying cells often have extensive mitochondrial contamination.


# This function can be used to pull information from any of the slots in the Assay class.
# For example, pull one of the data matrices("counts", "data", or "scale.data"):

HS.Sample[["MT.PER"]] = PercentageFeatureSet(HS.Sample, pattern = "mt-")

png("MitochondrialGenes_hist_27FEB20.png")
plot(hist(HS.Sample@meta.data$MT.PER, breaks=45,col="powderblue", border="black", xlab="% Mitochondrial Genes", main="HS", cex.lab=1, cex.axis=0.75, las=1))
dev.off()


png("MitochondrialPercentage_FeatureScatter_27FEB20.png")
FeatureScatter(HS.Sample, feature1 = "nCount_RNA", feature2 = "MT.PER", cols="black", pt.size=0.1)
dev.off()

png("RNAfeatures_FeatureScatter_27FEB20.png")
FeatureScatter(HS.Sample, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols="black", pt.size=0.1)
dev.off()



#Mitochondrial cutoff of 10%

HS.MTfilt = subset(HS.Sample, subset = MT.PER < 10)

dim(HS.Sample)
dim(HS.MTfilt)
#(1.No of genes 2. Cells)

#SAVED R. DATA 27FEB20 12:05


#extract raw counts
HS.MTfilt.counts=GetAssayData(object = HS.MTfilt, assay = "RNA", slot = "counts")
HS.filtSeuratOb= CreateSeuratObject(counts=as.matrix(HS.MTfilt.counts), min.cells=5, min.features=5, project = "GRN")
HS.filtSeuratOb[["MT.PER"]] = PercentageFeatureSet(HS.filtSeuratOb, pattern = "mt-")

#Create variable with non-normalised filtered matrix
HS.FilteredCountsMatrix=GetAssayData(object = HS.filtSeuratOb, assay = "RNA", slot = "counts")

HS.MTfilt.data=GetAssayData(object = HS.filtSeuratOb, assay = "RNA", slot = "data")

png("Total_Expression_BeforeNorm_27FEB20.png")
options(repr.plot.width=5, repr.plot.height=4)
plot(hist(colSums(HS.MTfilt.data),breaks = 100, main = "HS - Total expression before normalisation", xlab = "Counts Per Cell"))
dev.off()

# NORMALISATION STEP
# "LogNormalise" is the default normalisation method that normalises gene expression for
# each cell by total expression and multiplies this by a scale factor (10,000 is the default),
# and log transforms (10) the result. (this is what is stored in the assay objeect "data").

HS.filtNormData=NormalizeData(HS.filtSeuratOb, normalization.method = "LogNormalize", scale.factor = 10000)

HS.MTfilt.data=GetAssayData(object = HS.filtNormData, assay = "RNA", slot = "data")

png("Total_Expression_AfterNorm_27FEB20.png")
options(repr.plot.width=5, repr.plot.height=4)
hist(colSums(HS.filtNormData),breaks = 100, main = "Total expression after normalisation", xlab = "Normalised Counts Per Cell")
dev.off()

#FINDING VARIABLE FEATURES --------
# After normalisation, the next step is to find genes that vary between single cells
# (genes that are constant among all cells are not as informative). FindVariableGenes calculates
# average expression and variance for each gene and places these into bins, then normalises each bin.
# It then identifies features that are outliers on a 'mean variability plot'.

HS.filtNormData = FindVariableFeatures(HS.filtNormData, selection.method = "mean.var.plot", verbose=TRUE)

# For the mean.var.plot method: Exact parameter settings may vary empirically from dataset to dataset,
# and based on visual inspection of the plot. Setting the y.cutoff parameter to 2 identifies features
# that are more than two standard deviations away from the average dispersion within a bin.
# The default X-axis function is the mean expression level, and for Y-axis it is the log(Variance/mean).
# All mean/variance calculations are not performed in log-space, but the results are reported
# in log-space.




# SCALE DATA-------------------------------------------------------------------------------
# ScaleData() scales and centers genes in the dataset, which standardizes the range of expression
# values for each gene. The function additionally regresses out unwanted sources of variation such as
# technical noise. Data needs to be scaled to remove confounding factors such as batch effect and
# cell cycle stage affect. These Should be adjusted to remove the unwanted sources of variation and to
# infer the "correct" gene expression pattern. A linear model is built using the number of
# detected molecules per cell and the percentage mitochondrial RNA.
# The scaled and normalised residuals are stored in 'scale.data' and are used for
# dimentionality reduction and clustering.


allGenes.HS = rownames(HS.filtNormData)
HS.filtNormScaled = ScaleData(HS.filtNormData, features = allGenes.HS)

#Run PCA

HS.filtNormScaledPCA = RunPCA(HS.filtNormScaled,npcs = 30, verbose=F, vars.to.regress="MT.PER")

png("elbowplot_NormscaledPCA_27FEB20.png")
options(repr.plot.width=5, repr.plot.height=4)
ElbowPlot(HS.filtNormScaledPCA)
dev.off()


png("dimheatmap_NormscaledPCA_27FEB20.png")
options(repr.plot.width=10, repr.plot.height=20)
DimHeatmap(HS.filtNormScaledPCA, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

png("dimplot_NormscaledPCA_27FEB20.png")
options(repr.plot.width=5, repr.plot.height=4)
DimPlot(object = HS.filtNormScaledPCA, dim.1 = 1, dim.2 = 2)
dev.off()


#Jackstraw to determine significant PCs

png("JackStrawPlot_27FEB20.png")
options(repr.plot.width=5, repr.plot.height=4)P
HS.filtNormScaledPCA = JackStraw(HS.filtNormScaledPCA, num.replicate = 300)
HS.filtNormScaledPCA = ScoreJackStraw(HS.filtNormScaledPCA, reduction = "pca", dims = 1:20)
JackStrawPlot(object = HS.filtNormScaledPCA, dims= 1:20)
dev.off()
#This takes a good 2 hours^^^^^^

#Save+readvariables
saveRDS(HS.filtNormScaledPCA, "HS_filtNormScaledPCA.rds")
HS.filtNormScaledPCA=readRDS("/rds/general/user/ah3918/ephemeral/SCANALYSIS/HS100bp/HS_filtNormScaledPCA.rds")


#findclusters

HS.filtNormScaledPCA = FindNeighbors(HS.filtNormScaledPCA, dims = 1:20)
HS.filtNormScaledPCAclust = FindClusters(HS.filtNormScaledPCA, resolution = 1.2, pc.use=1:20)
HS.filtNormScaledPCAclust = RunUMAP(HS.filtNormScaledPCAclust, dims = 1:20)

#Plot UMAP

png("UMAP_28FEB20.png")
png("UMAP_28FEB20.png")
options(repr.plot.width=5, repr.plot.height=4)
DimPlot(HS.filtNormScaledPCAclust, reduction = "umap", label = TRUE, pt.size = 1,label.size=5) + NoLegend()
dev.off()



#DOUBLET FINDER --------------------------------------------------------------
#
# DoubletFinder is a package that predicts doublets in scRNA-seq data. It can be broken up into 4 steps:
# 1) Generate artificial doublets from existing scRNA-seq data (this is the above from Seurat),
# 2) Pre-process merged real-artificial data
# 3) Perform PCA and use the PC distance matrix to find each cell's proportion of artificial k
#    nearest neighbors (pANN)
# 4) Rank order and threshold pANN values according to the expected number of doublets
#
#
# The data needs to be run through the seurat pre-processing steps first.
#
#
# pK Identification:



Con.res.list.HS = paramSweep_v3(HS.filtNormScaledPCAclust, PCs = 1:20)
Con.res.list.stats.HS = summarizeSweep(Con.res.list.HS, GT = FALSE)
bcmvn.HS= find.pK(Con.res.list.stats.HS)
pK.HS=as.numeric(as.vector(bcmvn.HS$pK[which.max(bcmvn.HS$BCmetric)]))


#Chooses pK based on highest BCmetric value

#Homotypic Doublet proportion estimate

annotations.HS = HS.filtNormScaledPCAclust@meta.data$RNA_snn_res.1.2
Prop.Homotypic.HS = modelHomotypic(annotations.HS)
nExp_poi.HS = round(0.07*length(HS.filtNormScaledPCAclust@active.ident))
nExp_poi.adj.HS= round(nExp_poi.HS*(1-Prop.Homotypic.HS))

#Run Doubletfinder with varying classification stringencies

HS.filtNormScaledPCAclust = doubletFinder_v3(HS.filtNormScaledPCAclust, PCs = 1:20, pN = 0.25, pK = pK.HS, nExp = nExp_poi.HS, reuse.pANN = FALSE,sct = FALSE)
HS.filtNormScaledPCAclust = doubletFinder_v3(HS.filtNormScaledPCAclust, PCs = 1:20, pN = 0.25, pK = pK.HS, nExp = nExp_poi.adj.HS, reuse.pANN = paste(colnames(HS.filtNormScaledPCAclust@meta.data[7])), sct = FALSE)

#remove doublets

HS.filtNormScaledPCAclust@meta.data[,"DF_Class"] = get('HS.filtNormScaledPCAclust')@meta.data[8]
HS.filtNormScaledPCAclust@meta.data$DF_Class[which(HS.filtNormScaledPCAclust@meta.data$DF_Class == "Doublet" & get('HS.filtNormScaledPCAclust')@meta.data[9] == "Singlet")] = "Doublet_LOW"
HS.filtNormScaledPCAclust@meta.data$DF_Class[which(HS.filtNormScaledPCAclust@meta.data$DF_Class == "Doublet")] = "Doublet_HIGH"
table(HS.filtNormScaledPCAclust@meta.data$DF_Class)

#plot
png("DoubletPlot_28FEB20.png")
options(repr.plot.width=5, repr.plot.height=4)
DimPlot(HS.filtNormScaledPCAclust, reduction = "umap", group.by = colnames(as.matrix(HS.filtNormScaledPCAclust@meta.data[1,]))[9], label = TRUE)
dev.off()

# drop doublets

HS.doublets=as.character(colnames(HS.filtNormScaledPCAclust)[which(HS.filtNormScaledPCAclust@meta.data$DF_Class == "Doublet_HIGH")])
HS.filtNormScaledPCAclustDF=HS.filtNormScaledPCAclust[, !(colnames(HS.filtNormScaledPCAclust) %in% HS.doublets), drop = FALSE]
saveRDS(HS.filtNormScaledPCAclustDF, "HS_filtNormScaledPCAclustDF.rds")
HS.filtNormScaledPCAclustDF=readRDS("/rds/general/user/ah3918/ephemeral/SCANALYSIS/HS100bp/HS_filtNormScaledPCAclustDF.rds")


#find markers for each cluster


HS.dataMarkers= FindAllMarkers(HS.filtNormScaledPCAclustDF, only.pos = TRUE, min.pct = 0.25, logfc.threHSold = 0.25)
HS.dataMarkers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
str(HS.dataMarkers)
HS.MarkersPerCluster.LIST=unstack(HS.dataMarkers, HS.dataMarkers$gene ~ HS.dataMarkers$cluster)
str(HS.MarkersPerCluster.LIST)
saveRDS(HS.MarkersPerCluster.LIST,"HS_markersPerCluster.rds")

#Load Zeisel Markers

zeisel.markers=readRDS("/rdsgpfs/general/ephemeral/user/ah3918/ephemeral/SCANALYSIS/zeisel_markers.rds")


#Fishers Exact test for zeisel markers/cluster

HS.FisherTest.zeisel=function(x)
{
  TMP=matrix(ncol=9,nrow=1)
  HS.uniqueGenesPerCluster=intersect(HS.uniqueGenesAllClusters,x)
  for(j in 1:9)
  {
    TMP.Genes=intersect(HS.uniqueGenesAllClusters,zeisel.markers[[j]])
    TMP.MAT=matrix(ncol=2,nrow=2)
    TMP.MAT[1,1]=length(intersect(TMP.Genes,HS.uniqueGenesPerCluster))
    TMP.MAT[1,2]=length(setdiff(TMP.Genes,HS.uniqueGenesPerCluster))
    TMP.MAT[2,1]=length(setdiff(HS.uniqueGenesPerCluster,TMP.Genes))
    TMP.MAT[2,2]=length(HS.uniqueGenesAllClusters)-TMP.MAT[1,1]-TMP.MAT[1,2]-TMP.MAT[2,1]
    TMP[1,j]=Fisher.test(TMP.MAT,alternative="greater")$p.value
  }
  TMP
}

HS.data.zeisel=HS.filtNormScaledPCAclustDF
HS.uniqueGenesAllClusters=unique(unlist(HS.MarkersPerCluster.LIST))
HS.FisherTest.zeisel.result=lapply(HS.MarkersPerCluster.LIST,HS.FisherTest.zeisel)
HS.matrixPval.CellTypeXCluster = matrix(unlist(HS.FisherTest.zeisel.result), ncol = 9, byrow = TRUE)
colnames(HS.matrixPval.CellTypeXCluster)=names(zeisel.markers)

HS.matrixPval.CellTypeXCluster


#Annotate Cell types

HS.CellTypeTest.zeisel=function(x)
{
  CellType=apply(x,1,function(y){colnames(x)[which(y == min(y) & min(y) < 0.05 )]})
  CellType.Names=c(rep("Unclassified",length(HS.FisherTest.zeisel.result)))
  CellType.Names[which(CellType=="Mural")]="Mural"
  CellType.Names[which(CellType=="Endothelial")]="Endothelial"
  CellType.Names[which(CellType=="Ependymal")]="Ependymal"
  CellType.Names[which(CellType=="Microglia")]="Microglia"
  CellType.Names[which(CellType=="Oligodendrocyte")]="Oligodendrocyte"
  CellType.Names[which(CellType=="Interneuron")]="Interneurone"
  CellType.Names[which(CellType=="Pyramidal")]="CA1.Pyramidal"
  CellType.Names[which(CellType=="S1.Pyramidal")]="S1.Pyramidal"
  CellType.Names[which(CellType=="Astrocyte")]="Astrocyte"
  CellType.Names
}


HS.LabelledCellTypeClusters.zeisel=HS.CellTypeTest.zeisel(HS.matrixPval.CellTypeXCluster)
names(HS.LabelledCellTypeClusters.zeisel)=names(HS.FisherTest.zeisel.result)
HS.data.zeisel = RenameIdents(HS.data.zeisel,HS.LabelledCellTypeClusters.zeisel)
levels(HS.data.zeisel)

#New UMAP with cell types

png("UMAP_CellTypes_28FEB20.png")
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(HS.data.zeisel, reduction = "umap", label = TRUE, pt.size = 1,label.size=4)
table(Idents(HS.data.zeisel))
dev.off()

#Save cluster identities

HS.data.zeisel[["ClusterIdent"]] = Idents(object = HS.data.zeisel)
HS.CellIdentities=as.character(HS.data.zeisel@meta.data$ClusterIdent)
names(HS.CellIdentities)=rownames(HS.data.zeisel@meta.data)

saveRDS(as.factor(HS.CellIdentities),"HS_CellID_Type.rds")

HS.FinalRawFilteredMatrix=GetAssayData(object = HS.data.zeisel, assay = "RNA", slot = "counts")
HS.FinalNormedFilteredMatrix=GetAssayData(object = HS.data.zeisel, assay = "RNA", slot = "data")

saveRDS(HS.FinalNormedFilteredMatrix,"HS.FinalNormedFilteredMatrix.rds")
saveRDS(HS.data.zeisel,"HS.data.zeisel.rds")

#Save non-normalised matrix to use in CONOS
saveRDS(HS.FinalRawFilteredMatrix,"HS_FinalRawFiltered_dgCMatrix.rds")





#========================================================================================
# Analysis, continued! extract cell types -
#========================================================================================

#Extract Astrocyte data



source("/rds/general/user/ah3918/ephemeral/SCANALYSIS/single_cell_function.R")

HS.data.zeisel<-
table(Idents(HS.data.zeisel)))

AstrocyteSeuratObj<-subset((HS.data.zeisel),idents = c("Astrocyte"))


Astrocyte.Expr.Matrix<-as.matrix(GetAssayData(HS.data.zeisel)[, WhichCells(HS.data.zeisel, ident = "Astrocyte")])

saveRDS(AstrocyteSeuratObj, "Astrocyte_Seurat_Object.rds")
saveRDS(Astrocyte.Expr.Matrix, "Astrocyte_Raw_Expression_Matrix.rds")


write.table(Astrocyte.Expr.Matrix, "Astrocyte_Raw_Expression_Matrix_TEST2.txt",quote=FALSE, sep="\t")



#Subsetting other cell types

EndothelialSeuratObj<-subset((HS.data.zeisel),idents = c("Endothelial"))
InterneuronSeuratObj<-subset((HS.data.zeisel),idents = c("Interneuron"))
PyramidalSeuratObj<-subset((HS.data.zeisel),idents = c("Pyramidal"))
MicrogliaSeuratObj<-subset((HS.data.zeisel),idents = c("Microglia"))
OligodendrocyteSeuratObj<-subset((HS.data.zeisel),idents = c("Oligodendrocyte"))

Endothelial.Expr.Matrix<-as.matrix(GetAssayData(HS.data.zeisel)[, WhichCells(HS.data.zeisel, ident = "Endothelial")])
Interneuron.Expr.Matrix<-as.matrix(GetAssayData(HS.data.zeisel)[, WhichCells(HS.data.zeisel, ident = "Interneuron")])
Pyramidal.Expr.Matrix<-as.matrix(GetAssayData(HS.data.zeisel)[, WhichCells(HS.data.zeisel, ident = "Pyramidal")])
Microglia.Expr.Matrix<-as.matrix(GetAssayData(HS.data.zeisel)[, WhichCells(HS.data.zeisel, ident = "Microglia")])
Oligodendrocyte.Expr.Matrix<-as.matrix(GetAssayData(HS.data.zeisel)[, WhichCells(HS.data.zeisel, ident = "Oligodendrocyte")])

write.table(Endothelial.Expr.Matrix, "18MARCH20_Endothelial_Raw_Expression_Matrix.txt",quote=FALSE, sep="\t")
write.table(Interneuron.Expr.Matrix, "18MARCH20_Interneuron_Raw_Expression_Matrix.txt",quote=FALSE, sep="\t")
write.table(Pyramidal.Expr.Matrix, "18MARCH20_Pyramidal_Raw_Expression_Matrix.txt",quote=FALSE, sep="\t")
write.table(Microglia.Expr.Matrix, "18MARCH20_Microglia_Raw_Expression_Matrix.txt",quote=FALSE, sep="\t")
write.table(Oligodendrocyte.Expr.Matrix, "18MARCH20_Oligodendrocyte_Raw_Expression_Matrix.txt",quote=FALSE, sep="\t")

saveRDS(InterneuronSeuratObj, "InterneuronSeuratObj.rds")
saveRDS(EndothelialSeuratObj, "EndothelialSeuratObj.rds")
saveRDS(PyramidalSeuratObj, "PyramidalSeuratObj.rds")
saveRDS(MicrogliaSeuratObj, "MicrogliaSeuratObj.rds")
saveRDS(OligodendrocyteSeuratObj, "OligodendrocyteSeuratObj.rds")



#========================================================================================
# removing CA1 and S1 Pyramidal labels
#========================================================================================


HS.data.zeisel<-readRDS("HS.data.zeisel.rds")

HS.FisherTest.zeisel=function(x)
{
  TMP=matrix(ncol=8,nrow=1)
  HS.uniqueGenesPerCluster=intersect(HS.uniqueGenesAllClusters,x)
  for(j in 1:8)
  {
    TMP.Genes=intersect(HS.uniqueGenesAllClusters,zeisel.markers[[j]])
    TMP.MAT=matrix(ncol=2,nrow=2)
    TMP.MAT[1,1]=length(intersect(TMP.Genes,HS.uniqueGenesPerCluster))
    TMP.MAT[1,2]=length(setdiff(TMP.Genes,HS.uniqueGenesPerCluster))
    TMP.MAT[2,1]=length(setdiff(HS.uniqueGenesPerCluster,TMP.Genes))
    TMP.MAT[2,2]=length(HS.uniqueGenesAllClusters)-TMP.MAT[1,1]-TMP.MAT[1,2]-TMP.MAT[2,1]
    TMP[1,j]=Fisher.test(TMP.MAT,alternative="greater")$p.value
  }
  TMP
}

HS.data.zeisel=HS.filtNormScaledPCAclustDF
HS.uniqueGenesAllClusters=unique(unlist(HS.MarkersPerCluster.LIST))
HS.FisherTest.zeisel.result=lapply(HS.MarkersPerCluster.LIST,HS.FisherTest.zeisel)
HS.matrixPval.CellTypeXCluster = matrix(unlist(HS.FisherTest.zeisel.result), ncol = 9, byrow = TRUE)
colnames(HS.matrixPval.CellTypeXCluster)=names(zeisel.markers)

HS.matrixPval.CellTypeXCluster


#Annotate Cell types

HS.CellTypeTest.zeisel=function(x)
{
  CellType=apply(x,1,function(y){colnames(x)[which(y == min(y) & min(y) < 0.05 )]})
  CellType.Names=c(rep("Unclassified",length(HS.FisherTest.zeisel.result)))
  CellType.Names[which(CellType=="Mural")]="Mural"
  CellType.Names[which(CellType=="Endothelial")]="Endothelial"
  CellType.Names[which(CellType=="Ependymal")]="Ependymal"
  CellType.Names[which(CellType=="Microglia")]="Microglia"
  CellType.Names[which(CellType=="Oligodendrocyte")]="Oligodendrocyte"
  CellType.Names[which(CellType=="Interneuron")]="Interneuron"
  CellType.Names[which(CellType=="Pyramidal")]="Pyramidal"
  CellType.Names[which(CellType=="Astrocyte")]="Astrocyte"
  CellType.Names
}


HS.LabelledCellTypeClusters.zeisel=HS.CellTypeTest.zeisel(HS.matrixPval.CellTypeXCluster)
names(HS.LabelledCellTypeClusters.zeisel)=names(HS.FisherTest.zeisel.result)
HS.data.zeisel = RenameIdents(HS.data.zeisel,HS.LabelledCellTypeClusters.zeisel)
levels(HS.data.zeisel)

HS.data.zeisel<-RunTSNE(HS.data.zeisel,dims=1:20)

#remove unclassifeds

NoUnclassifieds<-subset((HS.data.zeisel),idents=c("Unclassified"),invert=TRUE)

pdf("21MAR20_tSNE_CellTypes_FINAL_17FEBMAR.pdf")
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(NoUnclassifieds, reduction = "tsne", label = TRUE, pt.size = 1,label.size=4)
dev.off()

saveRDS(HS.data.zeisel,"HS.data.zeisel.rds")
