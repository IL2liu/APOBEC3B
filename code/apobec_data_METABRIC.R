#################################################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## October 31, 2014
#################################################


#################################################
## Download breast cancer data from METABRIC
#################################################

myfn <- file.path(saveres, "metabric_bc.RData")
message("Retrieve gene expression data for METABRIC")
if (!file.exists(myfn)) {
  inSilicoDb::InSilicoLogin(login="bhaibeka@gmail.com", password="747779bec8a754b91076d6cc1f700831")
  metabric <- inSilicoDb::getDatasets(dataset="ISDB10278", norm="ORIGINAL", curation="25010", features="PROBE")[[1]]
  metabric <- MetaGx::probeGeneMapping(eset=metabric)
  inSilicoDb::InSilicoLogout()
  data.ge <- t(Biobase::exprs(metabric))
  annot <-  Biobase::fData(metabric)
  data.clin <- Biobase::pData(metabric)
  colnames(data.clin) <- gsub(" ", ".", colnames(data.clin))
  ## select tumor tissues
  tix <- !is.na(data.clin[ , "tissue"]) & data.clin[ , "tissue"] == "TUMOR"
  data.ge.normal <- data.ge[!tix, , drop=FALSE]
  data.clin.normal <- data.clin[!tix, , drop=FALSE]
  data.ge <- data.ge[tix, , drop=FALSE]
  data.clin <- data.clin[tix, , drop=FALSE]
  ## revert sample names to original metabric identifiers
  rownames(data.ge) <- rownames(data.clin) <- gsub("MB", "MB-", as.character(data.clin[ , "id"]))
  ## save data
  save(list=c("data.ge", "data.ge.normal", "annot", "data.clin", "data.clin.normal"), compress=TRUE, file=myfn)
} else { load(myfn) }

annot <- data.frame("probe"=rownames(annot), "EntrezGene.ID"=MetaGx::stripWhiteSpace(as.character(annot[ , "ENTREZID"])), annot)
rownames(annot) <- as.character(annot$probe)

## molecular subtyping using MetaGx
message("molecular subtyping for METABRIC")
## build exprs object
esetge <- ExpressionSet(assayData=t(data.ge), phenoData=AnnotatedDataFrame(data=data.clin), featureData=AnnotatedDataFrame(data=annot))
esetge <- subtypeClassification(eset=esetge, model=sbtm)
ss <- MetaGx::getSubtype(eset=esetge, method="class")
ss.proba <- MetaGx::getSubtype(eset=esetge, method="fuzzy")

# pdf(file.path(saveres, "tcga_bc_subtyping.pdf"), width=7, height=7)
# ss <- genefu::subtype.cluster.predict(sbt.model=scmod1.robust, data=data.ge, annot=annot, do.mapping=TRUE, plot=TRUE, verbose=FALSE)
# dev.off()

## double check subtyping
print(table("SCM"=ss, "ER"=data.clin[ , "er"]))
print(table("SCM"=ss, "HER2"=data.clin[ , "her2"]))
print(table("SCM"=ss, "GRADE"=data.clin[ , "grade"]))
tt <- read.csv("data/table_S2_revised.txt", sep="\t")
rownames(tt) <- as.character(tt[ , "METABRIC_ID"])
tt <- tt[rownames(data.clin), , drop=FALSE]
print(table(ss, tt[ , "IntClustMemb"]))

# print(table("SCM"=ss$subtype2, "GRADE"=data.clin[ , "grade"]))
survd <- survcomp::censor.time(surv.time=as.numeric(data.clin[ , "t.rfs"]) / 365, surv.event=as.numeric(data.clin[ , "e.rfs"]), time.cens=0)
dd <- data.frame("time"=survd[[1]], "event"=survd[[2]], "subtype"=ss)
pdf(file.path(saveres, "metabric_bc_subtyping_km.pdf"), width=7, height=7)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ subtype), data.s=dd, x.label="Time (years)", y.label="Probability of DFS", main.title=sprintf("Molecular subtypes (METABRIC)"), sub.title=NULL, leg.text=levels(dd$subtype), leg.pos="bottomright", leg.inset=0.0, v.line=NULL, h.line=NULL, .col=rainbow(length(levels(dd$subtype))), .lty=1, show.n.risk=TRUE, n.risk.step=3, n.risk.cex=0.85, bty="n", leg.bty="n", verbose=FALSE)
dev.off()


metabric <- list("annot"=annot, "ge"=data.ge, "clin"=data.clin, "subtype"=ss, "subtype.proba"=ss.proba, "ge.normal"=data.ge.normal, "clin.normal"=data.clin.normal)
rm(list=c("data.ge", "annot", "data.clin", "ss", "data.ge.normal", "data.clin.normal"))
gc(TRUE)