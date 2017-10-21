library(ggplot2)
library(pheatmap)
#load("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/HepG2_hclust.RData")
#load("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/HepG2_hclust_5.RData")
#load("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/HepG2_single_cell_expr_swapped.RData")
load("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/HepG2_hclust_transcript_4.RData")
load("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/HepG2_single_cell_expr_swapped_transcript.RData")

##################################################################################################################
##################################################################################################################

library(vioplot)
my.vioplot <- function(yalist,cols,main,names){
  plot(0,0,type="n",xlim=c(0,(length(yalist)+1)), ylim=range(yalist,na.rm=T),  xaxt = 'n', xlab ="", ylab = "",main=main,bty="l")
  for (i in seq(length(yalist))) { vioplot(na.omit(yalist[[i]]), at = i, add = T, col = cols[i]) }
  axis(side=1,at=seq(1,length(yalist)),labels=names)
}
##################################################################################################################
##################################################################################################################

library(gridGraphics)
library(grid)
grab_grob <- function(){
    grid.echo()
  grid.grab()
}
#########################################################################
################### squeeze NFR into feature space ######################
squeezNFR <- function(plus,minus,nfr,bin.cnt){
  feat.cnt <- ncol(plus)/bin.cnt
  x.new <- NULL
  for(i in seq(feat.cnt)){
    half1.idx <- seq((i-1)*bin.cnt+1,((i-1)*bin.cnt + bin.cnt/2))
    half2.idx <- seq(((i-1)*bin.cnt + bin.cnt/2) + 1,i*bin.cnt)
    x.new <- cbind(x.new,minus[,half1.idx],nfr[,i],plus[,half2.idx])
  }
  return(x.new)
}
##################################################################################################################
######################################## convert FPKM to TPM #####################################################
FPKM2TPM <- function(fpkm){
  tpm <- fpkm / sum(fpkm) * 10^6
}
##################################################################################################################

############# Merge the smallest cluster #############

#clusts[which(clusts==5)] <- 3
#clusts[which(clusts==6)] <- 5
#k <- 5
k <- 6 ### with the new BPs clustering, I'm not merging the clusters anymore
k <- 5 ### we decided to stick to k=5 (cf. the trello board discussions)
k <- 4 ### With transcript data and 2000bp set to threshold for nearby TSSs

######################################################

cell <- "HepG2_newBPs_transcriptClustering"

plus <- read.table("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/transcript/HepG2_scell_plus_transcript_expression.txt",header=T,row.names = 1)
minus <- read.table("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/transcript/HepG2_scell_minus_transcript_expression.txt",header=T,row.names = 1)
plus$transcript_id <- NULL
minus$transcript_id <- NULL

bulk_plus_expr <- log2(1 + FPKM2TPM(as.numeric(readLines("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/HepG2/01_HepG2_newBPs_plus_mRNA.y"))))
bulk_minus_expr <- log2(1 + FPKM2TPM(as.numeric(readLines("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/HepG2/01_HepG2_newBPs_minus_mRNA.y"))))

hm.plus <- read.table("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/HepG2/HepG2_newBPs_plus_100bp.x")
hm.minus <- read.table("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/HepG2/HepG2_newBPs_minus_100bp.x")
hm.nfr <- read.table("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/HepG2/HepG2_HM_NFR_bin.x")

hm.k122ac.plus <- read.table("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/181/H3K122ac/HepG2_plus_100bp.x")
hm.k122ac.minus <- read.table("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/181/H3K122ac/HepG2_minus_100bp.x")
hm.k122ac.nfr <- as.numeric(readLines("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/181/H3K122ac/HepG2_HM_H3K122ac_NFR_bin.x"))

hm.plus <- cbind(hm.plus,hm.k122ac.plus)
hm.minus <- cbind(hm.minus,hm.k122ac.minus)
hm.nfr <- cbind(hm.nfr,hm.k122ac.nfr)
hm.bin.cnt <- 40
hm <- squeezNFR(as.matrix(hm.plus),as.matrix(hm.minus),as.matrix(hm.nfr),hm.bin.cnt)

#hm.plus <- read.table("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/HepG2/HepG2_plus_1bp.x")
#hm.minus <- read.table("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/HepG2/HepG2_minus_1bp.x")

#hm.plus <- read.table("../../DEEP_ChIP_HM/01/HepG2/HepG2_plus_10bp.x")
#hm.minus <- read.table("../../DEEP_ChIP_HM/01/HepG2/HepG2_minus_10bp.x")

dnase.plus <- as.matrix(read.table(paste("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/DNase/HepG2_newBPs_DNase_plus_100bp.x",sep="")))
dnase.minus <- as.matrix(read.table(paste("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/DNase/HepG2_newBPs_DNase_minus_100bp.x",sep="")))
dnase.nfr <- as.numeric(readLines("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/DNase/HepG2_newBPs_DNase_NFR_bin.x"))
dnase.bin.cnt <- 40

tf.plus <- read.table("/MMCI/MS/ExpRegulation/work/data/ENC_ChIP_HM/TFs/HepG2/HepG2_newBPs_TF_plus_100bp.x")
tf.minus <- read.table("/MMCI/MS/ExpRegulation/work/data/ENC_ChIP_HM/TFs/HepG2/minus/HepG2_newBPs_TF_minus_100bp.x")
tf.nfr <- read.table("/MMCI/MS/ExpRegulation/work/data/ENC_ChIP_HM/TFs/HepG2/HepG2_newBPs_TF_minus_NFRplus.x")
tf.names <- readLines("/MMCI/MS/ExpRegulation/work/data/ENC_ChIP_HM/TFs/HepG2/HepG2_TF_names.txt")

nfr.coord <- read.table("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/HepG2/di_NFR_mRNA_HepG2_newBPs_sorted.txt",header=F,col.names = c("chr","start","end"),stringsAsFactors = F)

##################################################################################


cage.plus <- log2(1+as.numeric(readLines("/MMCI/MS/ExpRegulation/work/data/ENC_ChIP_HM/plus/HepG2/HepG2_newBPs_plus_allBDP_CAGE.y")))
cage.minus <- log2(1+as.numeric(readLines("/MMCI/MS/ExpRegulation/work/data/ENC_ChIP_HM/plus/HepG2/HepG2_newBPs_minus_allBDP_CAGE.y")))
##########
gc.content <- read.table("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/GC_content/GC_content_40bins_NFR.txt")## It includes the NFR bin it as well
gc.content <- read.table("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/GC_content/GC_content_80bins_NFR_50bp.txt")## It includes the NFR bin it as well
GC_BINS <- 81
##########
gtf.bp.plus <- read.table("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/HepG2/di_plus_newBPs_mRNA_HepG2_fullGTF_sorted.txt",header=T,stringsAsFactors=F)
gtf.bp.minus <- read.table("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/HepG2/di_minus_newBPs_mRNA_HepG2_fullGTF_sorted.txt",header=T,stringsAsFactors=F)
########## group gene products into PC and NPC groups ##########
gtf.bp.plus.grouped <- gtf.bp.plus
gtf.bp.minus.grouped <- gtf.bp.minus

for(i in seq(nrow(gtf.bp.plus.grouped))){
  if(gtf.bp.plus.grouped$geneProduct[i] == "protein_coding"){
    gtf.bp.plus.grouped$geneProduct[i] <- "PC"
  }else{
    gtf.bp.plus.grouped$geneProduct[i] <- "NPC"
  }

  if(gtf.bp.minus.grouped$geneProduct[i] == "protein_coding"){
    gtf.bp.minus.grouped$geneProduct[i] <- "PC"
  }else{
    gtf.bp.minus.grouped$geneProduct[i] <- "NPC"
  }
}
gtf.bp.plus <- gtf.bp.plus.grouped
gtf.bp.minus <- gtf.bp.minus.grouped
######## Find and Remove crazy gene MRpL30, form all data ##############
ens.bp.plus <- read.table("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/HepG2/di_plus_newBPs_mRNA_HepG2_sorted.txt",col.names=c("ENS","gene","strand","chr","TSS","val"),stringsAsFactors=F)
ens.bp.minus <- read.table("/MMCI/MS/ExpRegulation/work/data/DEEP_ChIP_HM/01/HepG2/di_minus_newBPs_mRNA_HepG2_sorted.txt",col.names=c("ENS","gene","strand","chr","TSS","val"),stringsAsFactors=F)

crazy_gene_MRPL30 <- which(ens.bp.plus[,2] == "MRPL30")
################################################################################
#### TBP data for core promoter analysis #######
tbp.plus <- read.csv("/MMCI/MS/ExpRegulation/work/data/singleCell/TBP_affinity_newBPs/RegulatorTrail - TEPIC_newBPs_plus_input.tepic.zip.csv")
tbp.minus <- read.csv("/MMCI/MS/ExpRegulation/work/data/singleCell/TBP_affinity_newBPs/RegulatorTrail - TEPIC_newBPs_minus_input.tepic.zip.csv")
##################################################################################
matchingBPs <- NULL
for(i in seq(nrow(tbp.plus))){
  hit <- which(ens.bp.plus$ENS == tbp.plus[i,1]);
  if(length(hit) != 0){
    j <- which(tbp.minus[,1] == ens.bp.minus$ENS[hit]);
    if(length(j) != 0){
      matchingBPs <- rbind(matchingBPs,cbind(tbp.plus[i,],tbp.minus[j,]));
    }
  } 
}
#### Removing from bulk ####
bulk_plus_expr <- bulk_plus_expr[-crazy_gene_MRPL30]
bulk_minus_expr <- bulk_minus_expr[-crazy_gene_MRPL30]

#### Removing from NFR ####
nfr.coord <- nfr.coord[-crazy_gene_MRPL30,]
#### Removing from GTF ####
gtf.bp.plus <- gtf.bp.plus[-crazy_gene_MRPL30,]
gtf.bp.minus <- gtf.bp.minus[-crazy_gene_MRPL30,]

#### Removing from HM ####
hm.plus <- hm.plus[-crazy_gene_MRPL30,]
hm.minus <- hm.minus[-crazy_gene_MRPL30,]
hm.nfr <- hm.nfr[-crazy_gene_MRPL30,]
hm <- hm[-crazy_gene_MRPL30,]

#### Removing from TF ####
dnase.plus <- dnase.plus[-crazy_gene_MRPL30,]
dnase.minus <- dnase.minus[-crazy_gene_MRPL30,]
dnase.nfr <- dnase.nfr[-crazy_gene_MRPL30]
#### Removing from TF ####
tf.plus <- tf.plus[-crazy_gene_MRPL30,]
tf.minus <- tf.minus[-crazy_gene_MRPL30,]
tf.nfr <- tf.nfr[-crazy_gene_MRPL30,]

#### Removing from CAGE ####
cage.plus <- cage.plus[-crazy_gene_MRPL30]
cage.minus <- cage.minus[-crazy_gene_MRPL30]
#### Removing from GC content ####
gc.content <- gc.content[-crazy_gene_MRPL30,]

##################################################################################################################
####### swap the gtf files acording to higher expression obtained from single cell #######
gtf.bp.right.swapped <- gtf.bp.plus
gtf.bp.right.swapped[swapped$swap.idx,] <- gtf.bp.minus[swapped$swap.idx,]
gtf.bp.left.swapped <- gtf.bp.minus
gtf.bp.left.swapped[swapped$swap.idx,] <- gtf.bp.plus[swapped$swap.idx,]
##################################################################################################################
####### swap the HM data acording to higher expression obtained from single cell #######
hm.bin.cnt <- 40 + 1###Because NFR is included now
hm.cnt <- ncol(hm)/hm.bin.cnt
hm.mirrored.idx <- NULL
for(i in seq(hm.cnt)){
  hm.mirrored.idx <- c(hm.mirrored.idx,seq(i*hm.bin.cnt,(i-1)*hm.bin.cnt+1,-1))
}
hm[swapped$swap.idx,] <- hm[swapped$swap.idx,hm.mirrored.idx]
write.table(hm,"/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/HepG2_newBPs_HM_NFRbin_swapped.x",quote=F,row.names=F,col.names=F,sep="\t")
### Plot single cell RNA-seq expression (heatmap) ###
temp <- bulk_minus_expr
bulk_minus_expr[swapped$swap.idx] <- bulk_plus_expr[swapped$swap.idx]
bulk_plus_expr[swapped$swap.idx] <- temp[swapped$swap.idx]
##################
temp <- minus
minus[swapped$swap.idx,] <- plus[swapped$swap.idx,]
plus[swapped$swap.idx,] <- temp[swapped$swap.idx,]
##################
temp <- dnase.minus
dnase.minus[swapped$swap.idx,] <- dnase.plus[swapped$swap.idx,seq(dnase.bin.cnt,1,-1)]
dnase.plus[swapped$swap.idx,] <- temp[swapped$swap.idx,seq(dnase.bin.cnt,1,-1)]
dnase <- cbind(dnase.minus[,seq(1,dnase.bin.cnt/2)],dnase.nfr,dnase.plus[,seq(1+dnase.bin.cnt/2,dnase.bin.cnt)])
##################
temp <- cage.minus
cage.minus[swapped$swap.idx] <- cage.plus[swapped$swap.idx]
cage.plus[swapped$swap.idx] <- temp[swapped$swap.idx]
##################
pm.pairs <- swapped$res
scell.cnt <- ncol(pm.pairs)/2
clustering.idx <- NULL
cols <- NULL
cors <- list()
geneProduct.counts <- list()
downsample.idx <- sample(min(sapply(seq(k),function(i)length(which(clusts == i)))))
nfr.coord[,2] <- as.numeric(nfr.coord[,2])
nfr.coord[,3] <- as.numeric(nfr.coord[,3])
gc.content.nfr.ordered <- NULL
gc.content.ordered <- NULL
hm.nfr.sorted <- NULL
dnase.nfr.sorted <- NULL
print(c("min cluster size",length(downsample.idx)))
gc.content[swapped$swap.idx,] <- gc.content[swapped$swap.idx,seq(GC_BINS,1,-1)]
for(i in seq(length(unique(clusts)))){
  idx <- which(clusts==i)#[downsample.idx]
  geneProduct.counts[[i]] <- table(paste(gtf.bp.left.swapped$geneProduct[idx],gtf.bp.right.swapped$geneProduct[idx],sep="-->"))
  clustering.idx <- c(clustering.idx,idx)
  cors[[i]] <- sapply(seq(length(idx)),function(j)cor(as.numeric(pm.pairs[idx[j],seq(scell.cnt)]),as.numeric(pm.pairs[idx[j],seq((scell.cnt + 1),2 * scell.cnt)])))
  cols <- c(cols,rep(rainbow(k)[i],times=length(idx)))
  nfr.sorted.idx <- order(nfr.coord[idx,3] - nfr.coord[idx,2],decreasing=F)
  gc.sorted.idx <- order(rowSums(gc.content[idx,]),decreasing=T)
  print(i)
  print(nfr.coord[idx[nfr.sorted.idx],3] - nfr.coord[idx[nfr.sorted.idx],2])
  print(cor((nfr.coord[idx,3] - nfr.coord[idx,2]),(nfr.coord[idx[nfr.sorted.idx],3] - nfr.coord[idx[nfr.sorted.idx],2])))
  gc.content.nfr.ordered <- rbind(gc.content.nfr.ordered,gc.content[idx[nfr.sorted.idx],])
  gc.content.ordered <- rbind(gc.content.ordered,gc.content[idx[gc.sorted.idx],])
  hm.nfr.sorted <- rbind(hm.nfr.sorted,hm[idx,])
  dnase.nfr.sorted <- rbind(dnase.nfr.sorted,dnase[idx,])
}
print(geneProduct.counts)
my_palette <- colorRampPalette(c("white","steelblue"))(256)


##### Sort clusters by NFR #####
print(sum(sum(gc.content[clustering.idx,]-gc.content.ordered)))

###############################################################################################################
###############################################################################################################
################### correlation VS NFR ####################
binningFunc <- function(x,span_length,B=10){
  ord.idx <- order(span_length,decreasing=F)
  span_length <- span_length[ord.idx]
  x <- x[ord.idx]
  bucket_size <- ceiling(length(span_length) / B)

  bins <- list()
  binned.x <- list()
  for(i in seq(B)){
    idx <- seq((i-1)*bucket_size,i*bucket_size)
    exceed <- which(idx > length(span_length))
    if(length(exceed) > 1)
      idx <- idx[-exceed]
    bins[[i]] <- range(span_length[idx])
    binned.x[[i]] <- x[idx]
    print(length(idx))
  }
  x_labels <- sapply(seq(length(bins)),function(i)paste(bins[[i]][1],bins[[i]][2],sep="-"))
  return(list(x=binned.x,names=x_labels))
}
all.cors <- sapply(seq(nrow(plus)),function(i)cor(as.numeric(log2(1+plus[i,])),as.numeric(log2(1+minus[i,]))))
cors.binnd <- list()
nfr.size <- abs(nfr.coord$end - nfr.coord$start);
nfr.bins <- list()
B <- 5
bucket_size <- ceiling(length(nfr.size) / B)
ord.idx <- order(nfr.size,decreasing=F)
nfr.size <- nfr.size[ord.idx]
all.cors.ord <- all.cors[ord.idx]
for(i in seq(B))
  nfr.bins[[i]] <- c((i-1)*bucket_size,i*bucket_size)

for(i in seq(length(nfr.bins))){
  idx <- which(nfr.size > nfr.bins[[i]][1] & nfr.size <= nfr.bins[[i]][2]);
  cors.binnd[[i]] <- all.cors[idx]
}
cors.binned <- binningFunc(all.cors,nfr.size,B)
###############################################################################################################
###############################################################################################################
################### correlation VS transcript span ####################
trans_span_coords_plus <- read.table("transcript/HepG2_scell_plus_transcript_coordinates.txt")
trans_span_coords_minus <- read.table("transcript/HepG2_scell_minus_transcript_coordinates.txt")
span_length_plus <- trans_span_coords_plus[,4] - trans_span_coords_plus[,3]
span_length_minus <- trans_span_coords_minus[,4] - trans_span_coords_minus[,3]
p_m_span_diff <- abs(span_length_minus - span_length_plus)
B <- 10
end <- max(p_m_span_diff)
bin_length <- ceiling(end/B)
bins <- list()
for(i in seq(B))
    bins[[i]] <- c((i-1)*bin_length,i*bin_length)
cors.span.binned <- list();
for(i in seq(B)){
  idx <- which(p_m_span_diff > bins[[i]][1] & p_m_span_diff <= bins[[i]][2]);
  cors.span.binned[[i]] <- all.cors[idx]
}
cors.span.binned <- binningFunc(all.cors,p_m_span_diff,B)
boxplot_names <- sapply(seq(length(nfr.bins)),function(i)paste(nfr.bins[[i]][1],nfr.bins[[i]][2],sep="-"))
pdf("/MMCI/MS/ExpRegulation/work/data/singleCell/HepG2/NFR_transcriptSpan_plots_variableBins.pdf")
boxplot(cors.binned$x,na.rm=T,xlab="NFR Bins",names=cors.binned$names,ylab="L-R correlation")
boxplot(cors.span.binned$x,na.rm=T,names=cors.span.binned$names,las=2,main="transcript span",ylab="L-R correlation",cex.axis=.6)
dev.off()
###############################################################################################################
###############################################################################################################
################### avg. sc expr. vs bulk expr. ####################
sc <- rowMeans(apply(rbind(minus,plus),2,FUN=as.numeric))
sc <- log2(1 + sc)
bulk <- c(bulk_minus_expr,bulk_plus_expr)
cors_bulkVSsingle <- c(round(cor(bulk,sc,method='s'),2),round(cor(bulk,sc,method='p'),2))
pdf(paste(cell,"_avgSC_VS_bulk.pdf",sep=""))
plot(x=sc,y=bulk,pch=20,xlab="avg. single cell log2 expr. (TPM)",ylab="bulk log2 expr. (TPM)",main=paste("HepG2 (plus & minus)\n","Spearman_cor",cors_bulkVSsingle[1],"\nPearson_cor",cors_bulkVSsingle[2]))
abline(a=0,b=1)
dev.off()
####################################################################
par(cex.main=.5)
if(F){
heatmap(as.matrix(gc.content[clustering.idx,]),Rowv=NA,Colv=NA,labRow = "",labCol = "",RowSideColors = cols,main="HepG2 GC-content",col=my_palette)
gc.sc.ord <- grab_grob()
#gplots::heatmap.2(as.matrix(gc.content.nfr.ordered),Rowv=F,Colv=F,labRow = "",labCol = "",dendrogram="n",trace="n",RowSideColors = cols,col = my_palette,main="HepG2 GC-content (NFR sorted)",key=F)
heatmap(as.matrix(gc.content.nfr.ordered),Rowv=NA,Colv=NA,labRow = "",labCol = "",main="HepG2 GC-content (NFR sorted)",RowSideColors = cols,col=my_palette)
gc.nfr.ord <- grab_grob()
heatmap(as.matrix(gc.content.ordered),Rowv=NA,Colv=NA,labRow = "",labCol = "",main="HepG2 GC-content sorted",RowSideColors = cols,col=my_palette)
gc.gc.ord <- grab_grob()
}
gplots::heatmap.2(as.matrix(gc.content[clustering.idx,]),Rowv=F,Colv=F,dendrogram="n",trace="n",labRow = "",labCol = "",RowSideColors = cols,main="HepG2 GC-content",col=my_palette,key=F)
gc.sc.ord <- grab_grob()
gplots::heatmap.2(as.matrix(gc.content.nfr.ordered),Rowv=F,Colv=F,labRow = "",labCol = "",dendrogram="n",trace="n",RowSideColors = cols,col = my_palette,main="HepG2 GC-content (NFR sorted)",key=F)
gc.nfr.ord <- grab_grob()
gplots::heatmap.2(as.matrix(gc.content.ordered),Rowv=F,Colv=F,labRow = "",labCol = "",main="HepG2 GC-content sorted",dendrogram="n",trace="n",RowSideColors = cols,col=my_palette,key=F)
gc.gc.ord <- grab_grob()
pdf(paste(cell,"_",k,"_GC_content_NFRbin_ordered_50bp.pdf",sep=""))
gplots::heatmap.2(as.matrix(gc.content[clustering.idx,]),Rowv=F,Colv=F,dendrogram="n",trace="n",labRow = "",labCol = "",RowSideColors = cols,main="HepG2 GC-content",col=my_palette,key=T,density.info="n")
gplots::heatmap.2(as.matrix(gc.content.nfr.ordered),Rowv=F,Colv=F,labRow = "",labCol = "",dendrogram="n",trace="n",RowSideColors = cols,col = my_palette,main="HepG2 GC-content (NFR sorted)",key=T,density.info="n")
gplots::heatmap.2(as.matrix(gc.content.ordered),Rowv=F,Colv=F,labRow = "",labCol = "",main="HepG2 GC-content sorted",dendrogram="n",trace="n",RowSideColors = cols,col=my_palette,key=T,density.info="n")
gplots::heatmap.2(as.matrix(log2(1 + cbind(bulk_minus_expr[clustering.idx],bulk_plus_expr[clustering.idx]))),Rowv=F,Colv=F,labRow = "",labCol = "",main="HepG2 bulk RNA",dendrogram="n",trace="n",RowSideColors = cols,col=my_palette,key=T,density.info="n")
grid.newpage()
lay <- grid.layout(nrow = 1, ncol=3)
pushViewport(viewport(layout = lay))
grid.draw(editGrob(gc.sc.ord, vp=viewport(layout.pos.row = 1,layout.pos.col = 1, clip=TRUE)))
grid.draw(editGrob(gc.nfr.ord, vp=viewport(layout.pos.row = 1,layout.pos.col = 2, clip=TRUE)))
grid.draw(editGrob(gc.gc.ord, vp=viewport(layout.pos.row = 1,layout.pos.col = 3, clip=TRUE)))
dev.off()
###############################################################################################################
###############################################################################################################

#gplots::heatmap.2(as.matrix(pm.pairs[clustering.idx,]),Rowv=F,Colv=F,dendrogram = "n",trace="n",labRow = "",labCol = "",RowSideColors = cols,key=F,col = my_palette)#,lwid = c(.5,.5,4.5,.3),lmat=rbind(c(4,4,3,0),c(0,2,1,0)))
white_side_color <- rep("white",length(cage.minus)) ### For plotting purposes when using grid.draw
heatmap(as.matrix(pm.pairs[clustering.idx,]),Rowv=NA,Colv=NA,labRow = "",labCol = "",RowSideColors = cols,col = my_palette,revC=T,scale="none",main=cell)
sc.hamp <- grab_grob()
heatmap(cbind(cage.minus[clustering.idx],cage.plus[clustering.idx]),Rowv=NA,Colv=NA,labRow = "",labCol = "",col = my_palette,main="CAGE",scale="none",revC=T,RowSideColors = white_side_color)
cage.hmap <- grab_grob()
heatmap(log2(1 + cbind(bulk_minus_expr[clustering.idx],bulk_plus_expr[clustering.idx])),Rowv=NA,Colv=NA,labRow = "",labCol = "",col = my_palette,main="bulk RNA",scale="none",revC=T,RowSideColors = white_side_color)
bulk.hmap <- grab_grob()


cors.ttest_res <- matrix(NA,ncol=k,nrow=k)
for(i in seq(k))
  for(j in seq(k))
    cors.ttest_res[i,j] <- round(t.test(cors[[i]],cors[[j]],alternative="t")$p.value,3)
print("single cell BP gene corrlation p-values (t-test):")
print(cors.ttest_res)
write.table(cors.ttest_res,paste(cell,"_single_cell_hierarchical_clustering_",k,"_correlation_pvalues.txt"),quote=F,row.names=F,col.names=F)

pdf(paste(cell,"_single_cell_hierarchical_clustering_",k,"_correlation.pdf",sep=""))##_downsampling.pdf",sep=""))
grid.newpage()
lay <- grid.layout(nrow = 1, ncol=3,widths=c(4,.5,.5))
pushViewport(viewport(layout = lay))
grid.draw(editGrob(sc.hamp, vp=viewport(layout.pos.row = 1,layout.pos.col = 1, clip=TRUE)))
grid.draw(editGrob(bulk.hmap, vp=viewport(layout.pos.row = 1,layout.pos.col = 2, clip=TRUE)))
grid.draw(editGrob(cage.hmap, vp=viewport(layout.pos.row = 1,layout.pos.col = 3, clip=TRUE)))
upViewport(1)
my.vioplot(cors,rainbow(k),"HepG2 single cell Pearson Correlation",paste("clust.",seq(k)))
dev.off()
##################### Report gene product frequencies in each cluster ########################
names_geneProducts <- unique(c(names(geneProduct.counts[[1]]),names(geneProduct.counts[[2]]),names(geneProduct.counts[[3]]),names(geneProduct.counts[[4]])))
geneProducts_mat <- matrix(NA, nrow=length(names_geneProducts),ncol=4)
rownames(geneProducts_mat) <- names_geneProducts
colnames(geneProducts_mat) <- seq(4)
for(i in names_geneProducts){
  for(j in seq(k))
    geneProducts_mat[i,j] <- geneProduct.counts[[j]][i]
}
geneProducts_mat[which(is.na(geneProducts_mat)==T)] <- 0
total_n <- sum(geneProducts_mat)
probs <- colSums(geneProducts_mat) / total_n
n <- rowSums(geneProducts_mat)
get.expected <- function(p,n){
  return(p * n)
}
expected <- NULL;
for(i in seq(nrow(geneProducts_mat)))
    expected <- rbind(expected,get.expected(probs,n[i]))
observed_expected <- geneProducts_mat / expected
rownames(observed_expected) <- names_geneProducts
colnames(observed_expected) <- seq(4)
chi_sqrd <- (expected - geneProducts_mat)^2/expected
phyper.pvals <- matrix(NA,ncol=ncol(geneProducts_mat),nrow=nrow(geneProducts_mat));
for(i in seq(nrow(geneProducts_mat)))
  for(j in seq(ncol(geneProducts_mat)))
    phyper.pvals[i,j] <- phyper(geneProducts_mat[i,j],m=sum(geneProducts_mat[i,]),n=sum(geneProducts_mat)-sum(geneProducts_mat[i,]),k=sum(geneProducts_mat[,j]),log.p=F,lower.tail=F)
rownames(phyper.pvals) <- names_geneProducts
colnames(phyper.pvals) <- seq(4)
pdf(paste(cell,"_geneProduct_frequencies_",k,".pdf",sep=""))
par(mar=c(2,4,2,4),oma=c(1,4,2,8),main.cex=.8)
gplots::heatmap.2(geneProducts_mat,Rowv=F,Colv=F,dendrogram = "n",trace="n",col = my_palette,main="HepG2_geneProducts",notecol="black",notecex=1,cellnote=geneProducts_mat)
gplots::heatmap.2(phyper.pvals,Rowv=F,Colv=F,dendrogram = "n",trace="n",col = my_palette,main="HepG2_geneProducts\nhypergeo p-values",notecol="black",notecex=1,cellnote=round(phyper.pvals,2))
gplots::heatmap.2(observed_expected,Rowv=F,Colv=F,dendrogram = "n",trace="n",col = my_palette,main="HepG2_geneProducts\nobserved / expected",notecol="black",notecex=1,cellnote=round(observed_expected,2))
dev.off()
##################### make individual heatmaps #####################
library(pheatmap)
my_palette1 <- colorRampPalette(c("white","steelblue"))(256)
for(i in seq(length(unique(clusts)))){
  idx <- which(clusts==i)
  right <- pm.pairs[idx,seq((ncol(pm.pairs)/2+1),ncol(pm.pairs))]
  left <- pm.pairs[idx,seq(1,ncol(pm.pairs)/2)]

  left <- minus[idx,]
  right <- plus[idx,]
  colnames(left) <- NULL
  colnames(right) <- NULL
  pdf(paste(cell,"_individual_BP_heatmaps_clust_",i,".pdf",sep=""))
  for(j in seq(length(idx))){
    if(isTRUE(all.equal(as.numeric(left[j,]),as.numeric(right[j,]))))
      next;
    #pheatmap(rbind(as.numeric(left[j,]),as.numeric(right[j,])),main=i,cluster_rows=F)
    p_m_data <- rbind(as.numeric(left[j,]),as.numeric(right[j,]))
    rownames(p_m_data) <- c(ens.bp.minus$gene[which(ens.bp.minus$ENS == rownames(minus[j,]))],ens.bp.plus$gene[which(ens.bp.plus$ENS == rownames(plus[j,]))])
    gplots::heatmap.2(p_m_data,Rowv=F,Colv=T,dendrogram="n",trace="n",labCol = "",main=i,col=my_palette1,key=T,sepcolor="black",colsep=c(0,ncol(p_m_data)),rowsep=0:nrow(p_m_data),sepwidth=c(0.1,0.005),cexRow=.5)
  }
  dev.off()
}
#########################################################################
#########################################################################

tf <- tf.nfr
#hm.bin.cnt <- 4000
tf.bin.cnt <- ncol(tf)/length(tf.names)
#hm <- hm.plus
tf.mirrored.idx <- NULL
for(i in seq(length(tf.names))){
  tf.mirrored.idx <- c(tf.mirrored.idx,seq(i*tf.bin.cnt,(i-1)*tf.bin.cnt+1,-1))
}

tf[swapped$swap.idx,] <- tf[swapped$swap.idx,tf.mirrored.idx]

#########################################################################
#########################################################################
nfr.sizes <- list()
for(kk in seq(k))
  nfr.sizes[[kk]] <- as.numeric(nfr.coord$end[which(clusts == kk)]) - as.numeric(nfr.coord$start[which(clusts == kk)])
#########################################################################
#########################################################################

plot.stratified.tf <- function(tf,tf.names,clusts,k,func,func.name,tf.bin.cnt){
  cols <- rainbow(k)
  for(j in seq(length(tf.names))){
    tf.median <- matrix(NA,nrow=k,ncol=tf.bin.cnt)
    for(kk in seq(k))
      tf.median[kk,] <- apply(tf[which(clusts==kk),seq((j-1)*tf.bin.cnt+1,j*tf.bin.cnt)],2,FUN=func)
    for(i in seq(k)){
      idx <- which(clusts==i)
      if(i == 1)
        plot(tf.median[i,],type="l",col=cols[i],main=paste(tf.names[j],func.name,sep="_"),ylim=range(tf.median),lwd=3)
      else
        points(tf.median[i,],type="l",col=cols[i],lwd=3)
    }
  }
}


plot.stratified.hm <- function(hm,clusts,k,func,func.name,bin.cnt,whichHMs){
  #par(mfrow=c(3,2))
  cols <- rainbow(k)
  hm.heatmaps <- NULL
  hm.names <- c('H3K4me1','H3K4me3','H3K27me3','H3K36me3','H3K9me3','H3K27ac','H3K122ac')
  if(length(whichHMs) <=3)
    par(mfrow=c(3,1))
  for(j in whichHMs){
    hm.median <- matrix(NA,nrow=k,ncol=bin.cnt)
    for(kk in seq(k))
      hm.median[kk,] <- apply(hm[which(clusts==kk),seq((j-1)*bin.cnt+1,j*bin.cnt)],2,FUN=func)
    for(i in seq(k)){
      idx <- which(clusts==i)
      if(i == 1)
        plot(hm.median[i,],type="l",col=cols[i],ylim=range(hm.median),lwd=3,xaxt='n',bty="n",ann=F,yaxp=c(round(range(hm.median),1),1),main=hm.names[j])
        #plot(hm.median[i,],type="l",col=cols[i],main=paste(hm.names[j],func.name,sep="_"),ylim=range(hm.median),lwd=3,xaxt='n',bty="n",xlab="",ylab="")
      else
        points(hm.median[i,],type="l",col=cols[i],lwd=3)
    }
    print(paste("printing HM median ranges in HM",hm.names[j]))
    print(round(range(hm.median),1))
    hm.heatmaps[[j]] <- grab_grob()
  }
  return(hm.heatmaps)
}

plot.stratified.dnase.heatmap <- function(dnase,clusts,k,bin.cnt){
  my_palette1 <- colorRampPalette(c("white","steelblue"))(256)
  my_palette2 <- colorRampPalette(c("white","red"))(256)
  cols <- NULL
  par(cex.main=.6)
  SC_clustered <- NULL
  cols <- NULL
  medians <- list()
  for(kk in seq(k)){
      idx <- which(clusts==kk)
      cols <- c(cols,rep(rainbow(k)[kk],times=length(idx)))
      SC_clustered <- rbind(SC_clustered,dnase[idx,])
      medians[[kk]] <- apply(dnase[idx,],2,FUN=median)
  }
  plot(medians[[1]],lwd=3,xaxt='n',bty="n",ann=F,ylim=range(medians),col=rainbow(k)[1],type="l",yaxp=c(round(range(medians),1),1))#,main="DNase")
  print("printing DNase median ranges:")
  print(round(range(medians),1))
  sapply(seq(2,k),function(i)lines(medians[[i]],col=rainbow(k)[i],lwd=3))
  median.grob <- grab_grob()
  #gplots::heatmap.2(SC_clustered,Rowv=F,Colv=F,dendrogram = "n",trace="n",labRow = "",labCol = "",RowSideColors = cols,main="DNase",key=F)#default color palette, which is a yellow-red kinna scale
  heatmap(SC_clustered,Rowv=NA,Colv=NA,labRow = "",labCol = "",scale="none",revC=T,RowSideColors = white_side_color)#,main="DNase")
  #gplots::heatmap.2(hm.SC_clustered,Rowv=F,Colv=F,dendrogram = "n",trace="n",labRow = "",labCol = "",RowSideColors = cols,main=paste(hm.names[j],cell,sep="_"),col = my_palette1)#customized palette defined above
  dnase.grab <- grab_grob()
  return(list(heatmap=dnase.grab,median=median.grob))
}

plot.stratified.hm.heatmap <- function(hm,clusts,k,bin.cnt,whichHMs){
  hm.names <- c('H3K4me1','H3K4me3','H3K27me3','H3K36me3','H3K9me3','H3K27ac','H3K122ac')
  my_palette1 <- colorRampPalette(c("white","steelblue"))(256)
  my_palette2 <- colorRampPalette(c("white","red"))(256)
  cols <- NULL
  par(cex.main=.6)
  hm.heatmaps <- list()
  for(j in whichHMs){
    hm.SC_clustered <- NULL
    cols <- NULL
    for(kk in seq(k)){
      idx <- which(clusts==kk)
      cols <- c(cols,rep(rainbow(k)[kk],times=length(idx)))
      hm.SC_clustered <- rbind(hm.SC_clustered,hm[idx,seq((j-1)*bin.cnt+1,j*bin.cnt)])
    }
    if(j == 1){
      heatmap(hm.SC_clustered,Rowv=NA,Colv=NA,labRow = "",labCol = "",RowSideColors = cols,scale="none",revC=T)#,main=paste(hm.names[j]))#default color palette, which is a yellow-red kinna scale
    }else{
      heatmap(hm.SC_clustered,Rowv=NA,Colv=NA,labRow = "",labCol = "",RowSideColors = white_side_color,scale="none",revC=T)#,main=paste(hm.names[j]))#default color palette, which is a yellow-red kinna scale
    }
    #gplots::heatmap.2(hm.SC_clustered,Rowv=F,Colv=F,dendrogram = "n",trace="n",labRow = "",labCol = "",RowSideColors = cols,main=paste(hm.names[j]))#default color palette, which is a yellow-red kinna scale
    #gplots::heatmap.2(hm.SC_clustered,Rowv=F,Colv=F,dendrogram = "n",trace="n",labRow = "",labCol = "",RowSideColors = cols,main=paste(hm.names[j],cell,sep="_"),col = my_palette1)#customized palette defined above
    hm.heatmaps[[j]] <- grab_grob()
  }
  return(hm.heatmaps)
}
#######################################################################################
#######################################################################################
#######################################################################################
plot.dnase.heatmap <- function(dnase,clusts,k,bin.cnt){
  my_palette1 <- colorRampPalette(c("white","steelblue"))(256)
  my_palette2 <- colorRampPalette(c("white","red"))(256)
  par(cex.main=.6)
  heatmap(dnase,Rowv=NA,Colv=NA,labRow = "",labCol = "",scale="none",revC=T,RowSideColors = white_side_color)#,main="DNase")
  dnase.grab <- grab_grob()
  return(list(heatmap=dnase.grab))
}

plot.hm.heatmap <- function(hm,clusts,k,bin.cnt,whichHMs){
  hm.names <- c('H3K4me1','H3K4me3','H3K27me3','H3K36me3','H3K9me3','H3K27ac','H3K122ac')
  my_palette1 <- colorRampPalette(c("white","steelblue"))(256)
  my_palette2 <- colorRampPalette(c("white","red"))(256)
  cols <- NULL
  par(cex.main=.6)
  hm.heatmaps <- list()
  for(j in whichHMs){
    cols <- NULL
    for(kk in seq(k)){
      idx <- which(clusts==kk)
      cols <- c(cols,rep(rainbow(k)[kk],times=length(idx)))
      hm.SC_clustered <- hm[,seq((j-1)*bin.cnt+1,j*bin.cnt)]
    }
    if(j == 1){
      heatmap(hm.SC_clustered,Rowv=NA,Colv=NA,labRow = "",labCol = "",RowSideColors = cols,scale="none",revC=T)#,main=paste(hm.names[j]))#default color palette, which is a yellow-red kinna scale
    }else{
      heatmap(hm.SC_clustered,Rowv=NA,Colv=NA,labRow = "",labCol = "",RowSideColors = white_side_color,scale="none",revC=T)#,main=paste(hm.names[j]))#default color palette, which is a yellow-red kinna scale
    }
    hm.heatmaps[[j]] <- grab_grob()
  }
  return(hm.heatmaps)
}
#if(F){## I added this while I was in the liver meeting, because the TF data is not ready according to the newBPs (1244)
pdf(paste(cell,"_TF_stratified_swpapped_mirrored_NFRbin_100bp_",k,"clusters.pdf",sep=""))
par(mfrow=c(2,2))
plot.stratified.tf(tf,tf.names,clusts,k,median,func.name ="median",(tf.bin.cnt))
plot.stratified.tf(tf,tf.names,clusts,k,mean,func.name="mean",(tf.bin.cnt))
dev.off()
#}

pdf(paste(cell,"_HM_stratified_swpapped_mirrored_100bp_heatmapRawData_log2Scaled_NFRbin_yellowRedColor.pdf",sep=""))
hm.grob <- plot.stratified.hm.heatmap(log2(1+as.matrix(hm)),clusts,k,(hm.bin.cnt),seq(hm.cnt))
dev.off()
#####################
pdf(paste(cell,"_HM_stratified_swpapped_mirrored_NFRbin_100bp.pdf",sep=""))
#pdf(paste("HepG2_HM_stratified_swpapped_mirrored_1bp_",k,"clusters.pdf",sep=""))
plot.stratified.hm(hm,clusts,k,median,func.name ="median",(hm.bin.cnt),seq(hm.cnt))
plot.stratified.hm(hm,clusts,k,mean,func.name="mean",(hm.bin.cnt),seq(hm.cnt))
dev.off()
##################################################################################################################
if(F){
pdf(paste(cell,"_HM_stratified_swpapped_mirrored_NFRbin_100bp_enhancersOnly_",k,"clusters.pdf",sep=""))
plot.stratified.hm(hm,clusts,k,median,func.name ="median",(hm.bin.cnt),c(1,6,7))
plot.stratified.hm(hm,clusts,k,mean,func.name="mean",(hm.bin.cnt),c(1,6,7))
dev.off()
}
##################################################################################################################
##################################################################################################################
#hm_median_log.grobs <- plot.stratified.hm(log2(1+hm),clusts,k,median,func.name ="median",(hm.bin.cnt),seq(hm.cnt))
hm_median.grobs <- plot.stratified.hm(hm,clusts,k,median,func.name ="median",(hm.bin.cnt),seq(hm.cnt))
hm.grobs <- plot.stratified.hm.heatmap(log2(1+as.matrix(hm)),clusts,k,(hm.bin.cnt),seq(hm.cnt))
dnase.grob.res <- plot.stratified.dnase.heatmap(dnase,clusts,k,dnase.bin.cnt)
dnase_log.grob.res <- plot.stratified.dnase.heatmap(log2(1+dnase),clusts,k,dnase.bin.cnt)
dnase_median.grob <- dnase.grob.res$median
dnase.grob <- dnase_log.grob.res$heatmap

plot(seq(10),col="white",bty="n",xaxt="n",ann=F,yaxt="n")
blank_grob <- grab_grob()

pdf(paste(cell,"_HM_DNase_stratified_swpapped_",k,"clusters.pdf",sep=""))
### I discussed with Marcel, and we decided that it should be no problem to show the median on the non-log-transformed data
if(F){
  ##### Showing the log transformed data for both heatmaps and median signal
  grid.newpage()
  lay <- grid.layout(nrow = 2, ncol=8,heights=unit(c(.3,.8),rep("null",2)))#,widths=unit(rep(,8),"inches"))
  pushViewport(viewport(layout = lay))
  for(i in seq(7)){
    grid.draw(editGrob(hm_median_log.grobs[[i]], vp=viewport(layout.pos.row = 1,layout.pos.col = i, clip=TRUE)))
    grid.draw(editGrob(hm.grobs[[i]], vp=viewport(layout.pos.row = 2,layout.pos.col = i, clip=TRUE)))
  }
  grid.draw(editGrob(dnase_log.grob.res$median, vp=viewport(layout.pos.row = 1,layout.pos.col = 8, clip=TRUE)))
  grid.draw(editGrob(dnase_log.grob.res$heatmap, vp=viewport(layout.pos.row = 2,layout.pos.col = 8, clip=TRUE)))
  upViewport(1)
}
##### Showing the non-logged median results
grid.newpage()
lay <- grid.layout(nrow = 2, ncol=8,heights=unit(c(.2,.8),rep("null",2)))#,widths=unit(rep(,8),"inches"))
pushViewport(viewport(layout = lay))
for(i in seq(7)){
  grid.draw(editGrob(hm_median.grobs[[i]], vp=viewport(layout.pos.row = 1,layout.pos.col = i, clip=TRUE)))
  grid.draw(editGrob(hm.grobs[[i]], vp=viewport(layout.pos.row = 2,layout.pos.col = i, clip=TRUE)))
}
grid.draw(editGrob(dnase.grob.res$median, vp=viewport(layout.pos.row = 1,layout.pos.col = 8, clip=TRUE)))
grid.draw(editGrob(dnase_log.grob.res$heatmap, vp=viewport(layout.pos.row = 2,layout.pos.col = 8, clip=TRUE)))
upViewport(1)
dev.off()
##################################################################################################################
##################################################################################################################

  
cage.str.plus <- list()
cage.str.minus <- list()
for(i in seq(k)){
  idx <- which(clusts==i)
  cage.str.plus[[i]] <- cage.plus[idx]
  cage.str.minus[[i]] <- cage.minus[idx]
}
merged.strs <- c(rbind(cage.str.minus,cage.str.plus))
mann_withney_res <- NULL;
for(i in seq(k))
  mann_withney_res <- c(mann_withney_res,wilcox.test(cage.str.minus[[i]],cage.str.plus[[i]],alternative="t",paired=T)$p.value)
print("CAGE statistical test")
print(mann_withney_res)
writeLines(as.character(mann_withney_res),paste(cell,"_CAGE_pvalues_",k,"clusters.txt",sep=""))
df <- data.frame(cage=c(cage.plus,cage.minus),TSS=c(rep("plus",times=length(cage.plus)),rep("minus",times=length(cage.minus))),cluster=c(clusts,clusts))
df.plus <- data.frame(cage=c(cage.plus),TSS=rep("plus",times=length(cage.plus)),cluster=clusts)
df.minus <- data.frame(cage=c(cage.minus),TSS=rep("minus",times=length(cage.plus)),cluster=clusts)
pdf(paste(cell,"_CAGE_stratified_",k,"clusters.pdf",sep=""))
my.vioplot(merged.strs,rep(rainbow(k),each=2),"HepG2 CAGE (minus_plus)",names=paste(rep(c("minus","plus"),times=k),rep(seq(k),each=2),sep="_"))
ggplot(df,aes(TSS,cage,fill=as.factor(cluster)))+geom_violin()#draw_quantiles = c(0.25, 0.5, 0.75),colour = "#3366FF")
ggplot(df.plus,aes(cage,x=TSS,fill=as.factor(cluster)))+geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
ggplot(df.minus,aes(cage,x=TSS,fill=as.factor(cluster)))+geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
boxplot(cage.str.plus,col=rainbow(k),main="plus")
boxplot(cage.str.minus,col=rainbow(k),main="minus")
dev.off()

##################################################################################################################
##################################################################################################################

nfr.ttest_res <- matrix(NA,ncol=k,nrow=k)
for(i in seq(k))
  for(j in seq(k))
    nfr.ttest_res[i,j] <- round(t.test(nfr.sizes[[i]],nfr.sizes[[j]],alternative="t",paired=F)$p.value,3)
write.table(nfr.ttest_res,paste(cell,"_NFR_sizes_pvalues_",k,"clusters.txt",sep=""),quote=F,row.names=F,col.names=F)

pdf(paste(cell,"_clustering_NFR_sizes_",k,"clusters.pdf",sep=""))
my.vioplot(nfr.sizes,rainbow(k),"HepG2 NFR size",names=seq(k))
boxplot(nfr.sizes,col=rainbow(k),main="HepG2 NFR size",names=seq(k))
dev.off()
##################################################################################################################
##################################################################################################################
##################################################################################################################
if(F){
  for(j in seq(6))
    for(i in seq(k)){
      idx <- which(clusts==i)
      for(l in idx)
        if(l == 1)
          plot(as.numeric(hm[l,seq((j-1)*40+1,j*40)]),type="l",col=cols[i],main=hm.names[j],ylim = range(hm[,seq((j-1)*40+1,j*40)]))
        else
          points(as.numeric(hm[l,seq((j-1)*40+1,j*40)]),type="l",col=cols[i])
    }
}
