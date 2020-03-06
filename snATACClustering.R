# Portions of this R script were adapted from Satpathy*, Granja* et al. (2019) 
# https://github.com/GreenleafLab/10x-scATAC-2019

library("Seurat")
library("Matrix")


#Handle input arguments

args = commandArgs(trailingOnly = TRUE)
inputfile = args[1]

numcomponents = abs(as.integer(args[2])) #Typical values: Pass 1=25 or Pass2=50
mincellcount = abs(as.integer(args[3])) #200
initres = abs(as.double(args[4])) #0.8
resrate = abs(as.double(args[5])) #0.8

#Set rate to default value when 1 or greater so the program will finish and not get stuck in an infinite loop
if(resrate >= 1.0){
  resrate = 0.8
}

outdir = args[6]

topvariable = -1
if(length(args) > 6){
  topvariable = as.integer(args[7])
}


seedvalue = 1

###
###

data = read.table(inputfile, sep="\t", header = TRUE)

numcells = max(data[,2])+1
numregions = max(data[,3])

#Creates a sparse matrix where cells = columns and regions = rows
sparseM <- Matrix::sparseMatrix(
  i = data[,3], 
  j = data[,2]+1,
  x = data[,4], 
  dims = c(numregions, numcells))

#Binarize the counts
sparseM@x[sparseM@x > 0] <- 1

#Subset the matrix if using the top variable regions
if(topvariable > 0){
  library("edgeR")
  library("matrixStats")
  print(paste("Using top variable regions:", topvariable))
  cpmMatrix = cpm(sparseM, log=TRUE, prior.count = 3)
  variablecpmb = rowVars(cpmMatrix)
  
  #Subset the sparse matrix by the top variable
  sparseM = sparseM[order(-variablecpmb)[1:topvariable],]
}

#get tf-idf matrix
freqs <- t(t(sparseM)/Matrix::colSums(sparseM)) 
idf   <- as(log(1 + ncol(sparseM) / Matrix::rowSums(sparseM)), "sparseVector")
tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs

set.seed(seedvalue)
svd <- irlba::irlba(tfidf, numcomponents, numcomponents)


#Make SVD matrix
svdDiag <- matrix(0, nrow=numcomponents, ncol=numcomponents)
diag(svdDiag) <- svd$d
matSVD <- t(svdDiag %*% t(svd$v))
colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
rownames(matSVD) <- unique(data[order(data[,2]),1])


#Seurat 
colnames(sparseM) <- unique(data[order(data[,2]),1])
rownames(sparseM) <- paste0("RI:",seq_len(nrow(sparseM)))

obj <- CreateSeuratObject(sparseM, project='scATAC', min.cells=0, min.features=0) #????

obj[["svd"]] <- CreateDimReducObject(embeddings = matSVD, key = "SVD_", assay = DefaultAssay(obj))

pdf(paste(outdir,"/","svd_bysample.pdf",sep=""))
DimPlot(obj, reduction = "svd", pt.size = 0.75)
dev.off()

#suerat clustering
clusters <- FindNeighbors(obj, reduction = "svd", dims = 2:numcomponents)
clusters <- FindClusters(clusters, resolution = initres)

minsize = min(table(clusters@meta.data[[paste0("RNA_snn_res.",initres)]]))
curres = initres

while(minsize < mincellcount){
  curres = curres * resrate
  clusters <- FindClusters(clusters, resolution = curres)
  minsize = min(table(clusters@meta.data[[paste0("RNA_snn_res.",curres)]]))
  print(minsize)
}

pdf(paste(outdir,"/","svd.pdf",sep=""))
DimPlot(clusters, reduction = "svd", pt.size = 0.75)
dev.off()

clusters <- RunUMAP(clusters, reduction="svd", dims = 2:numcomponents)

pdf(paste(outdir,"/","umap.pdf",sep=""))
DimPlot(clusters, reduction = "umap")
dev.off()
#head(Idents(clusters), 5)

metadata = clusters@meta.data

clusterinfo = metadata[,'seurat_clusters', drop=FALSE]

clustersbysamples = data.frame(matrix(NA, nrow=nrow(clusterinfo), ncol=3 ))

#write clusters for each sample to split the BAM files
for (i in 1:nrow(clusterinfo)){
  strid = rownames(clusterinfo)[i]
  clusterid = as.numeric(clusterinfo[i,1])
  split = strsplit(strid, '__', fixed=TRUE)
  sampleid = split[[1]][1]
  cellid = paste("_",split[[1]][2],sep="")
  clustersbysamples[i,1] = sampleid
  clustersbysamples[i,2] = cellid
  clustersbysamples[i,3] = clusterid-1 #Subtracting by 1 to make it consistent with the plots
}


for(cursample in unique(clustersbysamples[,1])){
  curcluster = clustersbysamples[clustersbysamples[,1] == cursample, 2:3]
  outputfile =  paste(outdir,"/",cursample, "_clusters.txt",sep="")
  write.table(curcluster, outputfile, 
              quote=FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)
}


