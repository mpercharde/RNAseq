
#####################################
#         PGC CNN/RDN SCRIPT        #
#    Michelle Percharde 2016        #
#     michelle.percharde@ucsf.edu   #
#####################################

source("http://bioconductor.org/biocLite.R")

library(limma)
library(gdata)
library(ggplot2)
library(gplots)
library("Rsamtools")
library("DESeq2")
library(GGally)

setwd("/Users/mpercharde/Documents/+BIOINF/PGC_RNAseq/")

#below is if using htseq-count on AWS. Better to use feature counts (Below) - much faster
# seqdata<-cbind((read.table("counts/trim.MPGC1_counts.txt", row.names=1, header=F, fill=T, sep="\t", check.names=F, quote="")),
#                (read.table("counts/trim.MPGC2_counts.txt", row.names=1, header=F, fill=T, sep="\t", check.names=F, quote="")),
#                (read.table("counts/trim.MPGC3_counts.txt", row.names=1, header=F, fill=T, sep="\t", check.names=F, quote="")),
#                (read.table("counts/trim.MSOMA1_counts.txt", row.names=1, header=F, fill=T, sep="\t", check.names=F, quote="")),
#                (read.table("counts/trim.MSOMA2_counts.txt", row.names=1, header=F, fill=T, sep="\t", check.names=F, quote="")),
#                (read.table("counts/trim.MSOMA3_counts.txt", row.names=1, header=F, fill=T, sep="\t", check.names=F, quote="")),
#                (read.table("counts/trim.FPGC1_counts.txt", row.names=1, header=F, fill=T, sep="\t", check.names=F, quote="")),
#                (read.table("counts/trim.FPGC2_counts.txt", row.names=1, header=F, fill=T, sep="\t", check.names=F, quote="")),
#                (read.table("counts/trim.FPGC3_counts.txt", row.names=1, header=F, fill=T, sep="\t", check.names=F, quote="")),
#                (read.table("counts/trim.FSOMA1_counts.txt", row.names=1, header=F, fill=T, sep="\t", check.names=F, quote="")),
#                (read.table("counts/trim.FSOMA2_counts.txt", row.names=1, header=F, fill=T, sep="\t", check.names=F, quote="")),
#                (read.table("counts/trim.FSOMA3_counts.txt", row.names=1, header=F, fill=T, sep="\t", check.names=F, quote="")))


#code below to use featurecounts, with downloaded bam files (recommended)
#check settings. Here, multimapping reads are discarded. Use TRUE and fraction=T to map multimappers fractionally

seqdata <- featureCounts(files=c("C1.sorted.bam","C2.sorted.bam","C3.sorted.bam","C4.sorted.bam","C5.sorted.bam",
                                "M1.sorted.bam","M2.sorted.bam","M3.sorted.bam","M4.sorted.bam","M5.sorted.bam"),
                        annot.ext="/Users/mpercharde/Documents/+BIOINF/path/to/genes_ercc.gtf",
                        isGTFAnnotationFile=T, GTF.featureType="exon", GTF.attrType="gene_id", useMetaFeatures=T, allowMultiOverlap=F,
                        isPairedEnd=F,nthreads=1, strandSpecific=1,countMultiMappingReads=F,chrAliases=NULL,reportReads=F)

#set colnames of samples
colnames(seqdata)<-c("MPGC1", "MPGC2", "MPGC3","MSOMA1","MSOMA2","MSOMA3","FPGC1","FPGC2","FPGC3","FSOMA1","FSOMA2","FSOMA3")  
# Xkr4	7	9	11	67	308	127	9	0	2	96	66	99

        write.table(seqdata, "htseq_counts_RAW.txt", sep="\t", quote=F, row.names=T)

####################   READ IN RAW COUNTS  ####################
seqdata<-read.table("htseq_counts_RAW.txt", header=T, row.names=1, quote="") #use if previously exported data table.
##############################################################

ercc<-seqdata[grep("^ERCC", rownames(seqdata)),] #get a table from subset of seqdata, only for ercc. WROTE TO FILE

#examine ERCC linearity between samples
pairs(~MPGC1+MPGC2+MPGC3+MSOMA1+MSOMA2+MSOMA3, data=ercc)
pairs(ercc[,1:6]) #male samples
pairs(ercc[,7:12]) #female samples


############################################################################################################
#Normalization and heatmap using Deseq - RAW DATA NORMALIZED TO READ DEPTH. #raw_expr is the table         #
############################################################################################################

#first filter the data to remove non-expressed rows, necessary. Generate tables for genes, and ERCCs.
is.expressed = apply(seqdata, 1, function(row) all(row !=0 ))   #remove any rows where there's a zero value
raw_expr <- seqdata[is.expressed,]
ercc_expr<-raw_expr[grep("^ERCC", rownames(raw_expr)),] #get a table from subset of seqdata, only for ercc. WROTE TO FILE
genes_expr <- raw_expr[!rownames(raw_expr) %in% rownames(ercc_expr),] #remove the ERCC data from the table
        write.table(raw_expr, "RAW_expressed.txt",sep="\t", quote=F, row.names=T)

#first read-depth normalize the data - using rlog method
library(DESeq2)
info<-read.csv("coldata.csv", header=T, row.names=1) #info table with sample name and condition. Example below:
        # SampleName Condition
# C1         C1    CreHet
# C2         C2    CreHet
# C3         C3    CreHet
# C4         C4    CreHet
# C5         C5    CreHet
# M1         M1      Null

matrix <- DESeqDataSetFromMatrix(countData = genes_expr,  colData = info, design = ~ Condition) #DESeq object made
matrixR <- rlog(matrix)   #this is a variance stabilization method. it also normalises for sequencing depth, and log transforms.
exprsR<-assay(matrixR)   #ONLY USRE THIS FOR LOOKING AT DATA, clustering, etc, USE RAW TABLE exprs FOR DE TESTING
head(exprsR)

        write.table(exprsR, "DEseq-READ-DEPTH-RLOG_norm.txt",sep="\t", quote=F, row.names=T)

#heatmap and UHC for all genes after read depth normalization
pcorr<-function(x) as.dist(1-cor(t(x),method="pearson")) #this is a function that allows pearson correlation to drive heatmap (Better than default euclidian distfun)                    
heatmap.2(as.matrix(exprsR), dendrogram="column", col=bluered, scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, Colv=T, distfun=pcorr)

########## TOPTABLE ANALYSIS USING THE READ DEPTH NORMALIZATION METHOD (NOT ERCC-see below) ##########

DEmatrix <- DESeq(matrix) #RAW data is processed to do Differential expression analysis NOT RLOG
toptable<-results(DEmatrix) #notes. P value is the Wald test p-value. 
mcols(toptable, use.names=T) #gives info on each column
# DataFrame with 6 rows and 2 columns
# type                                      description
# <character>                                      <character>
#         baseMean       intermediate        mean of normalized counts for all samples
# log2FoldChange      results log2 fold change (MAP): Condition Null vs CreHet
# lfcSE               results         standard error: Condition Null vs CreHet
# stat                results         Wald statistic: Condition Null vs CreHet
# pvalue              results      Wald test p-value: Condition Null vs CreHet
# padj                results                             BH adjusted p-values

summary(toptable) 
# out of 9973 with nonzero total read count
# blah blah blah etc

plotMA(toptable, main="DESeq2") #gives a basic MA plot. Use own code for prettier one

toptableF<-subset(toptable2, padj<0.1) #FDR <0.1 subbset
toptable2<-toptable[order((toptable$log2FoldChange)),] #order by increasing FC
        write.table(toptable2, "DESeq2-toptable-all.txt", sep="\t", quote=F, row.names=T) #all genes

##NB for other types of graphs etc, use code from the ERCC sections.

######################################################################################################
############################            voom version of ERCC norm                     ################  
#######################################################################################################

#NB expressed is data removing any zero values. # is.expressed = apply(seqdata, 1, function(row) all(row !=0 ))

#raw_expr - expressed raw genes
#ercc_expr - expressed ercc
names<-c("MPGC1", "MPGC2", "MPGC3","MSOMA1","MSOMA2","MSOMA3","FPGC1","FPGC2","FPGC3","FSOMA1","FSOMA2","FSOMA3")

library(limma)
library(edgeR)

condition <- factor(info[, "Condition"]) #Compare data by sample type
design <- model.matrix(~0+condition) #creates the matrix that tells you which condition it is. Try it - 1 = yes, 0 = no.
colnames(design) <- levels(condition)

#normalize the data using ERCC SPIKE-INs
N <- colSums(genes_expr)
nf <- calcNormFactors(ercc_expr, lib.size=N)
voom.norm <- voom(genes_expr, design, lib.size = N * nf, plot=T)

head(voom.norm$E) #- this is expression counts. Norm count + 0.5 log2. Check it looks ok on plot.
voom_genes<-voom.norm$E
        write.table(voom.norm$E, "VOOM_NORM-expr.genes_log2exprs.txt", sep="\t", quote=F, row.names=T)   #export ERCC norm count data 

#READ IN DATA ########################################################################################
voom_genes <-read.table("VOOM_NORM-expr.genes_log2exprs.txt", header=T, row.names=1, quote="")        
########################################################################################


############################### Density plots of raw vs ERCC-norm data ##############################

genes_expr <-raw_expr[!rownames(raw_expr) %in% rownames(ercc_expr),] #remove the ERCC data from the table
raw_av <- data.frame(MPGC=rowMeans(genes_expr[,1:3]), MSOMA=rowMeans(genes_expr[,4:6]), FPGC=rowMeans(genes_expr[,7:9]),FSOMA=rowMeans(genes_expr[,10:12]))
head(raw_av)

#################  now to with voom genes.   ##############

voom_av <- data.frame(MPGC=rowMeans(voom_genes[,1:3]), MSOMA=rowMeans(voom_genes[,4:6]), FPGC=rowMeans(voom_genes[,7:9]),FSOMA=rowMeans(voom_genes[,10:12]))
head(voom_av)

par( mfrow = c(1,2))

dat_female <- data.frame(dens = c(voom_av$FPGC, voom_av$FSOMA), key=rep(c("FPGCs", "FSOMA"), each = nrow(raw_av)))
ggplot(dat_female, aes(x = dens, fill = key)) + 
        geom_density(alpha = 0.3)
dat_male <- data.frame(dens = c(voom_av$MPGC, voom_av$MSOMA), key=rep(c("MPGCs", "MSOMA"), each = nrow(raw_av)))
ggplot(dat_male, aes(x = dens, fill = key)) + 
        geom_density(alpha = 0.3)

######################################### VOOM MDS ##################################################

plotMDS(voom.norm,xlim=c(-5,5), ylim=c(-5,5), col=key, pch=15)
legend(3.5,5, c("MPGCs", "MSOMA", "FPGCs", "FSOMA"), col=c("blue","grey","red","grey50"),  bty="n", pch=15) #clusters like other stuff - good.

key=c("blue","blue","blue","grey","grey","grey","red","red","red","grey50","grey50","grey50")

################# VOOM PCA - check clustering of data ##########################

pcaV<-prcomp(t(voom_genes), scale.=T, center=T)
head(pcaV$x)
head(pcaV$rotation)
scores=as.data.frame(pcaV$x)

######## create design info #######

colors=factor(c("blue","blue","blue","grey","grey","grey","red","red","red","grey50","grey50","grey50")) #these colors are by sample type

ggplot(data=scores, aes(x=PC1, y=PC2, label=rownames(scores))) +
        geom_point(color=colors, size=4) +
        geom_hline(yintercept = -80, color = "black", size=1) +
        geom_vline(xintercept = -120, color = "black", size=1) +
        # geom_text(colour = "black", alpha = 0.8, size = 4, hjust=-1) +
        theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12)) +
        theme(axis.title.x=element_text(vjust=-0.35, size=20), axis.title.y=element_text(vjust=0.35, size=20)) +
        theme(panel.background=element_rect(fill="white"), panel.grid.major=element_line(color="white", size=1)) +
        ggtitle("PC1 vs PC2 analysis -voom norm data") +
        theme(plot.title = element_text(size=20,lineheight=.8, vjust=2)) +
        ylim(c(-100,100))

#################################################################################
################################## VOOM CORRELATION HEATMAP #####################
#################################################################################

#use voom_genes (see read in above)
genes_corr<-cor(voom_genes, method="pearson")           
genes_corr
par(oma=c(3,4,3,4))
heatmap.2(genes_corr, scale="none", breaks=seq(0.5,1, length.out=41), col=colorpanel(40, "red3", "khaki1"), key=TRUE, symkey=F, 
          density.info="none", trace="none")

####### TOPTABLE ANALYSIS _ DE TESTING ##################

################################# MALE PGCs vs MALE SOMA #################################
fit <- lmFit(voom.norm,design) #fits a linear model for each gene. Takes the voom.norm object.
Mcontrasts <- makeContrasts(MPGC - MSOMA, levels = design) #compares up to three things, here is male PGCs vs SOMA
Mcontr.fit <- eBayes(contrasts.fit(fit, Mcontrasts)) #fits to the Bayesian model, comparing the chosen groups
Mdiff <- topTable(Mcontr.fit, coef=NULL, number=12974) #performs topTable analysis of genes for the defined genes. number is nrow(genes)
# > head(Mdiff)
# logFC  AveExpr         t       P.Value     adj.P.Val        B
# Dazl   11.005007 6.925225  27.30740 1.318678e-163 1.710853e-159 362.1509
# Ddx4    9.587726 6.125236  24.19770 5.370359e-129 3.483752e-125 282.9684
#blah etc

Mdiff_export <- Mdiff[order(-Mdiff$logFC),] #sort by descending logFC #12974 genes
Mdiff_export.padj05 = subset(Mdiff_export, adj.P.Val < 0.05) #changed to adj P <0.05 #10590 genes
Mdiff_export.padj05.logFC2 = subset(Mdiff_export.padj05, logFC > 2 | logFC < -2) #5785
        write.table(Mdiff_export.padj05.logFC2, "TopTable-Male-padj05,logFC2.txt", sep="\t", quote=F, row.names=T)
        write.table(Mdiff_export, "TopTable-Male-ALL.txt", sep="\t", quote=F, row.names=T)

nrow(subset(Mdiff_export.padj05, logFC < -2)) #1141 genes down
nrow(subset(Mdiff_export.padj05, logFC > 2)) #4644 genes up

################################# FEMALE PGCs vs FEMALE SOMA #################################
fit <- lmFit(voom.norm,design) #fits a linear model for each gene
Fcontrasts <- makeContrasts(FPGC - FSOMA, levels = design) #compares up to three things, here is male PGCs vs FED
Fcontr.fit <- eBayes(contrasts.fit(fit, Fcontrasts)) #fits to the Bayesian model, comparing the chosen groups
Fdiff <- topTable(Fcontr.fit, coef=NULL, number=12974) #performs topTable analysis of genes for the defined genes
head(Fdiff)
# logFC  AveExpr         t       P.Value     adj.P.Val        B
# Dazl    8.370203 6.925225  27.59554 5.072379e-167 6.580905e-163 371.3426
# Smc1b   8.499444 5.170298  25.74577 1.034526e-145 6.710970e-142 322.0247
# Hbb-y -10.319015 7.432861 -24.47225 6.918081e-132 2.991840e-128 289.2049
# Ddx4    8.449477 6.125236  24.03995 2.377051e-127 7.709965e-124 279.6534
# Stag3   7.075096 6.794783  23.94195 2.472693e-126 6.416144e-123 277.8116
# Stk31   8.314767 5.366119  22.83388 4.069228e-115 8.799028e-112 251.4958

Fdiff_export <- Fdiff[order(-Fdiff$logFC),] #sort by descending logFC 
Fdiff_export.padj05 = subset(Fdiff_export, adj.P.Val < 0.05) #changed to adj P <0.05 #11329 genes
Fdiff_export.padj05.logFC2 = subset(Fdiff_export.padj05, logFC > 2 | logFC < -2) #8453
write.table(Fdiff_export.padj05.logFC2, "TopTable-Female-padj05,logFC2.txt", sep="\t", quote=F, row.names=T)
write.table(Fdiff_export, "TopTable-Female-ALL.txt", sep="\t", quote=F, row.names=T)

nrow(subset(Fdiff_export.padj05, logFC < -2)) #449 genes down
nrow(subset(Fdiff_export.padj05, logFC > 2)) #7954 genes up


######################################################################################
################################# HEATMAPS YEAH ######################################
######################################################################################

#heatmaps of this stuff

#heatmaps of average data
voom_av <- data.frame(MPGC=rowMeans(voom_genes[,1:3]), MSOMA=rowMeans(voom_genes[,4:6]), FPGC=rowMeans(voom_genes[,7:9]),FSOMA=rowMeans(voom_genes[,10:12]))
# head(voom_av)

my_palette2 <- colorRampPalette(c("purple4", "white", "orange2"))(n=299)
heatmap.2(as.matrix(voom_genes), dendrogram="column", col=my_palette2, scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, Colv=T, labRow=F, distfun=pcorr)


ribo3 <-grep("Rpl|Rps", rownames(voom_genes)) #get row numbers which contain Rpl or Rps in the title.
voom_ribo <- voom_genes[ribo3,]

ribo_av<-grep("Rpl|Rps", rownames(voom_av)) #get row numbers which contain Rpl or Rps in the title.
voom_ribo_av<- voom_av[ribo_av,]
r6k_av <-voom_ribo_av[grep("Rps6k", rownames(voom_ribo_av)),]
voom_ribo_av2<-voom_ribo_av[!rownames(voom_ribo_av) %in% rownames(r6k_av),] #remove weird Rps6k

#heatmap for average male, average female.


################################################################################################################
#                                               MA Plots   = DE genes                                    #
################################################################################################################

library(gplots)

###### MA PLOT HERE #######

Mdiff$Color = "black"
Mdiff$Color[Mdiff$logFC > 2.0 ] = "orange2"
Mdiff$Color[Mdiff$logFC < -2 ] = "purple4"

Mdiff2<-as.data.frame(Mdiff)

ggplot(data=Mdiff2, aes(AveExpr,logFC)) +       
        geom_point(alpha=0.4, size=2, color=Mdiff2$Color) +
        geom_hline(yintercept = 0, colour = "black", size=1) +
        # geom_hline(yintercept = -2, colour = "black", size=1) +
        geom_hline(yintercept = 2.1, colour = "orange2", size=1.2) +
        geom_hline(yintercept = -2, colour = "purple4", size=1.2) +
        geom_vline(xintercept = -5, colour = "black", size=1) +
        xlim(c(-5, 15)) +
        xlab("log2 average expression") + ylab("log2 fold change") +
        theme(panel.background=element_rect(fill="white"), panel.grid.major=element_line(color="white", size=1)) +
        theme(axis.text.x=element_text(size=20, color="black"), axis.text.y=element_text(size=20, color="black")) +
        ggtitle("MA plot - MPGCs vs MSOMA") +
        theme(plot.title = element_text(size=18,lineheight=.8, vjust=2))

Fdiff$Color = "black"
Fdiff$Color[Fdiff$logFC > 2 ] = "orange2"
Fdiff$Color[Fdiff$logFC < -2 ] = "purple4"

Fdiff2<-as.data.frame(Fdiff)

# sel<-Fdiff2["Sp110",]

ggplot(data=Fdiff2, aes(AveExpr,logFC)) +       
        geom_point(alpha=0.4, size=2, color=Fdiff2$Color) +
        geom_hline(yintercept = 0, colour = "black", size=1) +
        # geom_hline(yintercept = -2, colour = "black", size=1) +
        geom_hline(yintercept = 2.1, colour = "orange2", size=1.2) +
        geom_hline(yintercept = -2, colour = "purple4", size=1.2) +
        geom_vline(xintercept = -5, colour = "black", size=1) +
        xlim(c(-5, 15)) +
        ylim(c(-10, 12)) +
        xlab("log2 average expression") + ylab("log2 fold change") +
        theme(panel.background=element_rect(fill="white"), panel.grid.major=element_line(color="white", size=1)) +
        theme(axis.text.x=element_text(size=20, color="black"), axis.text.y=element_text(size=20, color="black")) +
        ggtitle("MA plot - FPGCs vs FSOMA") +
        theme(plot.title = element_text(size=18,lineheight=.8, vjust=2))


################################################################################################################
#                                               Venn diagrams/overlaps                                      #
################################################################################################################


biocLite("graph")
biocLite("RBGL")
install.packages("Vennerable", repos="http://R-Forge.R-project.org", type="source")
library(Vennerable)
library(graph)
library(RBGL)
library(reshape)

Mdiff_up = subset(Mdiff_export.padj05.logFC2, logFC > 2) #4644
Fdiff_up = subset(Fdiff_export.padj05.logFC2, logFC > 2) #7954

vennData = list(rownames(Mdiff_up), rownames(Fdiff_up))
names(vennData) = c("MPGC genes", "FPGC genes")
str(vennData)
vennDiagram = Venn(vennData)
vennDiagram #see what the proportion of genes is

# plot a venn diagram - very basic.
plot(vennDiagram, doWeights=T) 

#replot the intersection using venneuler
biocLite("venneuler")
library(venneuler)
v <- venneuler(c(MPGCs=390, FPGCs=3700, "MPGCs&FPGCs"=4254))
plot(v)


#code to export the overlapping genes
mVenn<-rownames(Mdiff_up)
fVenn<-rownames(Fdiff_up)

is.element(mVenn, fVenn) #gives a logical vector of if they're the same (looks for F in M with this order)
intersect<-Mdiff_up[is.element(mVenn,fVenn),] #returns a list of the intersecting genes
write.table(intersect, "Venn-PGC-intersection.txt", sep="\t", quote=F, row.names=T)  

################################################################################################################
#                                               Venn diagrams/overlaps                                      #
################################################################################################################

#define expressed as count >=1 in all samples

Mexpressed = apply(seqdata[,1:3], 1, function(row) all(row >=1 ))
nrow(seqdata[Mexpressed,]) #14776
MSexpressed = apply(seqdata[,4:6], 1, function(row) all(row >=1 ))
nrow(seqdata[MSexpressed,]) #15703
Fexpressed = apply(seqdata[,7:9], 1, function(row) all(row >=1 ))
nrow(seqdata[Fexpressed,]) #14778
FSexpressed = apply(seqdata[,10:12], 1, function(row) all(row >=1 ))
nrow(seqdata[FSexpressed,]) #15847

#define expressed as count >=10 in all samples

Mexpressed = apply(seqdata[,1:3], 1, function(row) all(row >=10 ))
nrow(seqdata[Mexpressed,]) #12760
MSexpressed = apply(seqdata[,4:6], 1, function(row) all(row >=10 ))
nrow(seqdata[MSexpressed,]) #13717
Fexpressed = apply(seqdata[,7:9], 1, function(row) all(row >=10 ))
nrow(seqdata[Fexpressed,]) #12314
FSexpressed = apply(seqdata[,10:12], 1, function(row) all(row >=10 ))
nrow(seqdata[FSexpressed,]) #13403

######################################################################################################################
#Define universes of expressed genes
######################################################################################################################


####################   READ IN RAW COUNTS  ####################
seqdata<-read.table("htseq_counts_RAW.txt", header=T, row.names=1, quote="")
##############################################################

# is.expressed = apply(seqdata, 1, function(row) all(row !=0 ))   #remove any rows where there's a zero value

keep <- rowSums(cpm(seqdata)>1) >= 3
seqdata2<-seqdata[keep,]
nrow(seqdata2) #15342

ercc<-seqdata[grep("^ERCC", rownames(seqdata)),] #get a table from subset of seqdata, only for ercc. WROTE TO FILE
genes2 <-seqdata2[!rownames(seqdata2) %in% rownames(ercc),] #remove the ERCC data from the table


#normalize data
names<-c("MPGC1", "MPGC2", "MPGC3","MSOMA1","MSOMA2","MSOMA3","FPGC1","FPGC2","FPGC3","FSOMA1","FSOMA2","FSOMA3")
info<-read.csv("coldata.csv", header=T, row.names=1)
condition <- factor(info[, "Condition"]) #Compare data by sample type
design <- model.matrix(~0+condition) #creates the matrix that tells you which condition it is. Try it - 1 = yes, 0 = no.
colnames(design) <- levels(condition)

library(edgeR)
library(limma)
#normalize the data using ERCC SPIKE-INs
N <- colSums(genes2)
nf <- calcNormFactors(ercc, lib.size=N)
genes2.norm <- voom(genes2, design, lib.size = N * nf, plot=T)

head(genes2.norm$E) #- this is expression counts. Norm count + 0.5 log2
genes2.voom<-genes2.norm$E
# write.table(voom.norm$E, "VOOM_NORM-expr.genes_log2exprs.txt", sep="\t", quote=F, row.names=T)   #export count data 

#define expressed as count >=1 in all samples

Mexpressed = apply(genes2.voom[,1:3], 1, function(row) all(row >=1 ))
nrow(genes2.voom[Mexpressed,]) #11858 MPGC

# M1 <-lapply(genes2.voom[,1], length)>=5
# nrow(genes2.voom[M1,])

nrow(subset(genes2.voom, genes2.voom[,1] >=1)) #MPGC1 12644
nrow(subset(genes2.voom, genes2.voom[,2] >=1)) #MPGC2 12361
nrow(subset(genes2.voom, genes2.voom[,3] >=1)) #MPGC3 12427
nrow(subset(genes2.voom, genes2.voom[,4] >=1)) #MSOMA1 12080
nrow(subset(genes2.voom, genes2.voom[,5] >=1)) #MSOMA2 12049
nrow(subset(genes2.voom, genes2.voom[,6] >=1)) #MSOMA3 13172

nrow(subset(genes2.voom, genes2.voom[,7] >=1)) #FPGC1 12990
nrow(subset(genes2.voom, genes2.voom[,8] >=1)) #FPGC2 12883
nrow(subset(genes2.voom, genes2.voom[,9] >=1)) #FPGC3 12868
nrow(subset(genes2.voom, genes2.voom[,10] >=1)) #FSOMA1 11936
nrow(subset(genes2.voom, genes2.voom[,11] >=1)) #FSOMA2 11884
nrow(subset(genes2.voom, genes2.voom[,12] >=1)) #FSOMA3 11974

nrow(subset(genes2.voom, genes2.voom[,1] >=5)) #MPGC1 8040
nrow(subset(genes2.voom, genes2.voom[,2] >=5)) #MPGC2 7860
nrow(subset(genes2.voom, genes2.voom[,3] >=5)) #MPGC3 7829
nrow(subset(genes2.voom, genes2.voom[,4] >=5)) #MSOMA1 4500
nrow(subset(genes2.voom, genes2.voom[,5] >=5)) #MSOMA2 4452
nrow(subset(genes2.voom, genes2.voom[,6] >=5)) #MSOMA3 7157

nrow(subset(genes2.voom, genes2.voom[,7] >=5)) #FPGC1 8361
nrow(subset(genes2.voom, genes2.voom[,8] >=5)) #FPGC2 7824
nrow(subset(genes2.voom, genes2.voom[,9] >=5)) #FPGC3 8223
nrow(subset(genes2.voom, genes2.voom[,10] >=5)) #FSOMA1 2818
nrow(subset(genes2.voom, genes2.voom[,11] >=5)) #FSOMA2 2487
nrow(subset(genes2.voom, genes2.voom[,12] >=5)) #FSOMA3 3559



Mexpressed = apply(genes2.voom[,4:6], 1, function(row) all(row >=1 ))
nrow(genes2.voom[Mexpressed,]) #11775 MSOMA

Fexpressed = apply(genes2.voom[,7:9], 1, function(row) all(row >=1 ))
nrow(genes2.voom[Fexpressed,]) #12325 FPGC

Fexpressed = apply(genes2.voom[,10:12], 1, function(row) all(row >=1 ))
nrow(genes2.voom[Fexpressed,]) #11349

#define highly expressed as count >=5 in all samples

Mexpressed = apply(genes2.voom[,1:3], 1, function(row) all(row >=5 ))
nrow(genes2.voom[Mexpressed,]) #7519 MPGC

Mexpressed = apply(genes2.voom[,4:6], 1, function(row) all(row >=5 ))
nrow(genes2.voom[Mexpressed,]) #4260 MSOMA

Fexpressed = apply(genes2.voom[,7:9], 1, function(row) all(row >=5 ))
nrow(genes2.voom[Fexpressed,]) #7690 FPGC

Fexpressed = apply(genes2.voom[,10:12], 1, function(row) all(row >=5 ))
nrow(genes2.voom[Fexpressed,]) #2366









