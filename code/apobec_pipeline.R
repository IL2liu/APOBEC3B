#################################################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## October 31, 2014
#################################################

## APOBEC and Genomic instability in breast cancer
## genomic instability = chromosomal instability (CIN70)
## Carter, S. L., Eklund, A. C., Kohane, I. S., Harris, L. N., & Szallasi, Z. (2006). A signature of chromosomal instability inferred from gene expression profiles predicts clinical outcome in multiple human cancers. Nature genetics, 38(9), 1043–1048. doi:10.1038/ng1861.

## install devel version of MetaGx and PharmacGx
# library(devtools)
# devtools::install_github("bhklab/survcomp", ref="master")
# devtools::install_github("bhklab/genefu", ref="master")
# devtools::install_github("bhklab/MetaGx", ref="master")

library(parallel)
library(Biobase)
library(biomaRt)
library(genefu)
library(inSilicoDb)
library(MetaGx)
library(org.Hs.eg.db)
library(plotrix)
library(cgdsr)
library(Hmisc)
library(xtable)
library(gdata)
library(nnet)
library(lmtest)
library(SuppDists)

## prevent strings to be converted into factors
options(stringsAsFactors=FALSE)

## number of cpu cores available for the analysis pipeline
## set to 'NULL' if all the available cores should be used
nbcore <- 16
availcore <- parallel::detectCores()
if (is.null(nbcore) || nbcore > availcore) { nbcore <- availcore }
options("mc.cores"=nbcore)
  
## GSEA parameters
## path to the GSEA java executable
gsea.exec <- file.path("gsea", "gsea2-2.1.0.jar")
## number of random geneset permutations
gsea.nperm <- 1000
## minimum size for a geneset to be analyzed
min.geneset.size <- 15
## maximum size for a geneset to be analyzed
max.geneset.size <- 250
## path to the geneset gmt file
genesets.filen <- file.path("gsea", "c5.all.v4.0.entrez.gmt")

## functions

## recursive application of functions
fold <- function(f, x, y, ...){
    if (missing(...)) { f(x, y) } else { f(x, fold(f, y, ...)) }
}

## subtyping
## Curtis et al
# sbtm <- "intClust"
# ll <- as.list(paste("iC", 1:10, sep=""))
# names(ll) <- unlist(ll)
# sbtoi <- c(list("Global"=paste("iC", 1:10, sep="")), ll)
## Haibe-Kains et al
sbtm <- "scmod2"
sbtoi <- list("Global"=c("Basal", "Her2", "LumB", "LumA"), "Basal"="Basal", "Her2"="Her2", "Luminals"=c("LumB", "LumA"), "LuminalB"="LumB", "LuminalA"="LumA")
# sbtoi <- list("Global"=c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif"), "Basal"="ER-/HER2-", "Her2"="HER2+", "Luminals"=c("ER+/HER2- High Prolif", "ER+/HER2- Low Prolif"), "LuminalB"="ER+/HER2- High Prolif", "LuminalA"="ER+/HER2- Low Prolif")
sbtcol <- rainbow(length(sbtoi), alpha=0.6)
names(sbtcol) <- names(sbtoi)

## create directory to save the results
saveres <- sprintf("saveres_%s", sbtm)
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE) }
  
## genes of interest
## apobec family
goi <- c("APOBEC3B"="geneid.9582", "APOBEC3G"="geneid.60489", "APOBEC3F"="geneid.200316", "APOBEC3C"="geneid.27350", "APOBEC3A"="geneid.200315", "APOBEC3H"="geneid.164668", "APOBEC3D"="geneid.140564", "APOBEC3A_B"="geneid.100913187")
## proliferation-related genes
goi.prolif <- c("AURKA"="geneid.6790", "MKI67"="geneid.4288", "CCNB1"="geneid.891", "CDKN2A"="geneid.1029")

#################################################
## create CIN70 signature
#################################################
## extracted from GeneSigDB: http://compbio.dfci.harvard.edu/genesigdb/signaturedetail.jsp?signatureId=16921376-GeneList

cin70 <- read.csv(file.path("data", "CIN70.csv"), quote="'", stringsAsFactors=FALSE)
## manual curation
curation <- rbind(c("CDC2", "CDK1"),
  c("C20orf24/TGIF2", "TGIF2-C20orf24"),
  c("CNAP1", "NCAPD2"),
  c("CDC45L", "CDC45"),
  c("ch-TOG", "CKAP5"),
  c("KS2", "CKS2"),
  c("CTPS", "CTPS1"),
  c("BRRN1", "NCAPH"),
  c("MTB", "NCAPG2"),
  c("TOPK", "PBK"),
  c("FLJ10036", "ZWILCH"),
  c("GPIandMGC13096", "GPI"),
  c("SFRS2", "SRSF2"),
  c("STK6", "AURKA"),
  c("KIAA0286", "TMEM194A"))
cin70[match(curation[ , 1], cin70[ , "Gene.Symbol"]), "Gene.Symbol"] <- curation[ , 2]

## annotations
## get the gene symbols from entrez gene id using org.Hs.eg.db
gs <- toTable(org.Hs.egSYMBOL)
gs <- gs[!is.na(gs[ , "symbol"]) & !duplicated(gs[ , "symbol"]), , drop=FALSE]
rownames(gs) <- gs[ , "symbol"]
gs <- gs[cin70[ , "Gene.Symbol"], "gene_id"]

## create signature
cin70.sig <- list("CIN70"=cbind("probe"=cin70[ , "Gene.Symbol"], "EntrezGene.ID"=gs, "coefficient"=1))
iix <- !is.na(gs)
cin70t <- cbind("probe"=paste("geneid", gs[iix], sep="."), "EntrezGene.ID"=gs[iix], "coefficient"=1)
write.csv(cin70.sig, file=file.path(saveres, "CIN70_new.csv"))


#################################################
## get breast cancer data from TCGA
#################################################

myfna <- file.path(saveres, "tcga_bc_data.RData")
if (!file.exists(myfna)) {
  source(file.path("code", "apobec_data_TCGA.R"))
  save(list="tcga", compress=TRUE, file=myfna)
}

#################################################
## get breast cancer data from METABRIC
#################################################

myfna <- file.path(saveres, "metabric_bc_data.RData")
if (!file.exists(myfna)) {
  source(file.path("code", "apobec_data_METABRIC.R"))
  save(list="metabric", compress=TRUE, file=myfna)
}


#################################################
## tcga analyses
#################################################

dataset.name <- "TCGA"
load(file.path(saveres, "tcga_bc_data.RData"))
dataset <- tcga
rm(list="tcga")
gc()

## mutations and CNV calls
nsheets <- c("PIK3CA.E542K", "PIK3CA.E545K", "PIK3CA.H1047R", "P53", "Other")
tt <- ids <- NULL
for (i in 1:length(nsheets)) {
  tt <- c(tt, list(gdata::read.xls(xls=file.path("data", "TCGA TP53 and PIK3CA mutation calls (cbio).xlsx"), sheet=i, header=TRUE, stringsAsFactors=FALSE)))
  if (i == 1) {
    ids <- tt[[i]][ , 1]
  } else {
    ids <- union(ids, tt[[1]][ , 1])
  }
}
tt <- sapply(tt, function (x, y) { return (x[match(y, x[ , 1]), -1, drop=FALSE]) }, y=ids)
tt <- do.call(cbind, tt)
rownames(tt) <- ids
tt2 <- data.frame(matrix(NA, ncol=ncol(tt), nrow=nrow(dataset$clin), dimnames=list(rownames(dataset$clin), colnames(tt))))
myx <- intersect(rownames(dataset$clin), rownames(tt))
tt2[myx, ] <- tt[myx, ]
mutations <- tt2
mutations[!is.na(mutations) & (mutations == "" | mutations == "#N/A" | (mutations != "0" & mutations != "1"))] <- NA

## add the hypermutation phenotype
hypermut <- gdata::read.xls(xls=file.path("data", "TCGA breast hypermutator status from nik-zainal.xlsx"), sheet=1)
nn <- substr(as.character(hypermut[ , 1]), 1, 12)
hypermut <- hypermut[ , 2, drop=TRUE]
names(hypermut) <- nn
iix <- intersect(rownames(mutations), names(hypermut))
mutations <- cbind(mutations, "HYPERMUTATION"=NA)
mutations[iix, "HYPERMUTATION"] <- hypermut[iix]
## APOBEC3B CNV calls
tt <- gdata::read.xls(xls=file.path("data", "TCGA CNV calls 6-4-14.xlsx"), sheet=1, header=FALSE, stringsAsFactors=FALSE)
nn <- as.character(tt[ , 1])
tt <- tt[ , 2]
names(tt) <- nn
cnvcalls <- array(NA, dim=nrow(dataset$clin), dimnames=list(rownames(dataset$clin)))
nn <- intersect(names(cnvcalls), nn)
cnvcalls[nn] <- tt[nn]
cnvcalls <- factor(cnvcalls)
levels(cnvcalls) <- c("A3B del/del", "A3B del/wt", "A3B wt/wt")

source(file.path("code", "apobec_analysis.R"))

rm(list=c("dataset", "dataset.name", "cnvcalls", "mutations"))

#################################################
## metabric analyses
#################################################

dataset.name <- "METABRIC"
load(file.path(saveres, "metabric_bc_data.RData"))
dataset <- metabric
rm(list="metabric")
gc()

tt <- gdata::read.xls(xls=file.path("data", "METABRIC CNV and p53 calls.xlsx"), sheet=1, stringsAsFactors=FALSE)
rownames(tt) <- tt[ , "ID"]
tt2 <- data.frame(matrix(NA, ncol=ncol(tt), nrow=nrow(dataset$clin), dimnames=list(rownames(dataset$clin), colnames(tt))))
myx <- intersect(rownames(dataset$clin), rownames(tt))
tt2[myx, ] <- tt[myx, ]
## APOBEC3B CNV calls
cnvcalls <- tt2[ , "CNV.Call"]
cnvcalls <- factor(cnvcalls)
levels(cnvcalls) <- c("A3B del/del", "A3B del/wt", "A3B wt/wt")
names(cnvcalls) <- rownames(tt2)
# levels(cnvcalls) <- c("deletion", "heterozygous", "amplified")
## TP53 mutations
mutations <- tt2[ , "X53.MUT", drop=FALSE]
mutations[!is.na(mutations) & (mutations == "" | mutations == "#N/A")] <- NA
# mutations <- apply(mutations, 2, as.factor)
## 0=WT/slient, 1=MUT
colnames(mutations) <- c("P53")
## add the hypermutation phenotype
mutations <- cbind(mutations, "HYPERMUTATION"=NA)

source(file.path("code", "apobec_analysis.R"))

rm(list=c("dataset", "dataset.name", "cnvcalls", "mutations"))

## save session info
write(toLatex(sessionInfo(), locale = FALSE), file="sessionInfoR.tex", append=FALSE)

## end





