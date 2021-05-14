suppressMessages(require("scCB2"))
suppressMessages(require("Matrix"))
args <- commandArgs(trailingOnly = TRUE)

# command line arguments
samp <- args[1]
indir <- args[2]
outdir <- args[3]

if (!dir.exists(outdir)){
  dir.create(outdir)
} else {
  print("Dir already exists!")
}

mtx_file = paste0(samp,'_sparse_molecule_counts.mtx')
bc_file = paste0(samp,'_sparse_counts_barcodes.csv')
g_file = paste0(samp,'_sparse_counts_genes.csv')

ofile = paste0(samp, '_realdrops.csv')

barcodes = read.csv(paste0(indir,bc_file),stringsAsFactors = F,
                        header=F,sep=',',row.names=1,as.is=T)
barcodes = as.character(barcodes$V2)

genes = read.csv(paste0(indir,g_file),stringsAsFactors = F,
                      header=F,sep=',',row.names=1,as.is=T)
genes = as.character(genes$V2)

my.counts = readMM(paste0(indir,mtx_file))

rownames(my.counts) <- barcodes
colnames(my.counts) <- genes
my.counts = t(my.counts)

ind = !duplicated(genes)
my.counts = my.counts[ind,]

my.counts = as(my.counts,'dgCMatrix')

CBOut <- CB2FindCell(my.counts, FDR_threshold = 0.01, 
    lower = 100, Ncores = 2, verbose = TRUE)
summary(CBOut)

RealCell <- GetCellMat(CBOut)
str(RealCell)

real_cells = colnames(RealCell)

write.table( data.frame(real_cells), paste0(outdir,ofile), quote = F, row.names = F, col.names = F )

