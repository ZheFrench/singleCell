library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(velocyto.R)
library(glue)
print("Let's do it.")
#http://pklab.med.harvard.edu/velocyto/notebooks/R/DG1.nb.html
#http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html
#June 10, 2020
library(optparse)
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="Loom file ", metavar="Loom file")

)

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = TRUE);
opt  <- arguments$options
args <- arguments$args

print("> OPTS : ")
print(opt$file)
filename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(opt$file))

print(filename)
stop()
ldat <- ReadVelocity(file = opt$file)

bm <- as.Seurat(x = ldat)

bm[["RNA"]] <- bm[["spliced"]]

bm <- SCTransform(object = bm, assay = "spliced") #Use this function as an alternative to the NormalizeData, FindVariableFeatures, ScaleData workflow
bm <- RunPCA(object = bm, verbose = FALSE)

bm <- RunUMAP(object = bm,reduction = "pca", dims = 1:30)
bm <- RunTSNE(object = bm,reduction = "pca")

bm <- FindNeighbors(object = bm, reduction = "pca", dims = 1:30)
bm <- FindClusters(object = bm, resolution = 0.5)

bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02)

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)

cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)

png(glue("/data/villemin/code/singleCell/python/figures/{filename}.png"))
show.velocity.on.embedding.cor(emb = Embeddings(object = bm, reduction = "umap"), vel = Tool(object = bm, 
    slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
    cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
    do.par = FALSE, cell.border.alpha = 0.1,n.cores = 1)
dev.off()

DefaultAssay(bm) <- "RNA"

SaveH5Seurat(bm, filename = glue("/data/villemin/code/singleCell/python/{filename}.filtered.h5Seurat"))
Convert(glue("/data/villemin/code/singleCell/python/{filename}.h5Seurat"), dest = "h5ad")
