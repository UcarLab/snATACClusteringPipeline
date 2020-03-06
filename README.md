# snATACClusteringPipeline

Pipeline for obtain clusters from snATAC-seq data.

For now, this pipeline is an implementation of the method describe by Satpathy, A.T, Granja, J.M., et al. (2019).

# Dependencies:
snATACClusteringTools: Obtain/Compile the jar from: https://github.com/UcarLab/snATACClusteringTools
Place the jar in the same directory at the shell script.

R (Tested with verion 3.6.2)
R Libraries: 

1. Seurat: https://cran.r-project.org/web/packages/Seurat/index.html
2. EdgeR: https://bioconductor.org/packages/release/bioc/html/edgeR.html
3. Matrix: https://cran.r-project.org/web/packages/Matrix/index.html
4. matrixStats: https://cran.r-project.org/web/packages/matrixStats/index.html
