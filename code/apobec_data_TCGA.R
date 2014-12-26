#################################################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## October 31, 2014
#################################################

#################################################
## Download breast cancer data from TCGA
#################################################

## all human gene symbols
gs <- toTable(org.Hs.egSYMBOL)
## get all gene symbols
gs2 <- as.character(unique(gs[ , "symbol"]))

myfn <- file.path(saveres, "tcga_bc.RData")
if (!file.exists(myfn)) {
  ## retrieve the breast cancer data
  # Create CGDS object
  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/public-portal/")
  ## test the connection
  # message("Test connection to cBioPortal")
  # test(mycgds)
  # message("")

  # Get list of cancer studies at server
  ## get the breast cancer data
  cancerstudy <- cgdsr::getCancerStudies(mycgds)
  mycancerstudy <- cancerstudy[which(cancerstudy[ , 1] == "brca_tcga"), 1] 
  ## get the collection of cases
  caselist <- cgdsr::getCaseLists(mycgds, mycancerstudy)
  # mycaselist <- caselist[1, 1]
  mycaselist <- "brca_tcga_3way_complete"
  
  ## get the gene expression with RNA-seq
  message("Retrieve gene expression data for TCGA")
  myfn2 <- file.path(saveres, "tcga_bc_ge.RData")
  if (!file.exists(myfn2)) {
    geneticprofile <- cgdsr::getGeneticProfiles(mycgds, mycancerstudy)
    # mygeneticprofile <- geneticprofile[3, 1]
    mygeneticprofile <- "brca_tcga_rna_seq_v2_mrna"
    ## parallel version
    splitix <- parallel::splitIndices(nx=length(gs2), ncl=ceiling(length(gs2) / 100))
    mcres <- parallel::mclapply(splitix, function(x, mycgds, gs2, mygeneticprofile, mycaselist) {
    	try(dd <- cgdsr::getProfileData(x=mycgds, genes=gs2[x], geneticProfiles=mygeneticprofile, caseList=mycaselist), silent=FALSE)
      if (class(try) == "try-error") {
        dd <- NULL
      } else {
        cix <- apply(dd, 2, function(x) { return(all(is.na(x))) })
      	dd <- dd[ , !cix, drop=FALSE]
        if(nrow(dd) == 0 || ncol(dd) == 0) {
          dd <- NULL
        } else {
          dd <- data.matrix(dd)
          dd[!is.na(dd) & dd == "NaN"] <- NA
        }
      }
    	return (dd)
    }, mycgds=mycgds, gs2=gs2, mygeneticprofile=mygeneticprofile, mycaselist=mycaselist, mc.cores=nbcore)
    mcres <- mcres[!sapply(mcres, is.null)]
    data.ge <- do.call(cbind, mcres)
    # rix <- apply(data.ge, 1, function(x) { return(all(is.na(x))) })
    # data.ge <- data.ge[!rix, , drop=FALSE]
    data.ge <- data.matrix(data.ge)
    data.ge <- log2(data.ge + 1)
    save(list="data.ge", compress=TRUE, file=myfn2)
  } else { load(myfn2) }

  ## retrieval of copy number variations
  message("Retrieve CNV data for TCGA")
  myfn2 <- file.path(saveres, "tcga_bc_cnv.RData")
  if (!file.exists(myfn2)) {
  # mygeneticprofile <- geneticprofile[5, 1]
    mygeneticprofile <- "brca_tcga_log2CNA"
    ## get APOBEC3B copy number alteration
    apobec3b.cnv <- cgdsr::getProfileData(x=mycgds, genes="APOBEC3B", geneticProfiles=mygeneticprofile, caseList=mycaselist)
    nn <- rownames(apobec3b.cnv)
    apobec3b.cnv <- apobec3b.cnv[ , 1, drop=TRUE]
    names(apobec3b.cnv) <- nn
    ## parallel version tpo get all the copy number alterations
    splitix <- parallel::splitIndices(nx=length(gs2), ncl=ceiling(length(gs2) / 100))
    mcres <- parallel::mclapply(splitix, function(x, mycgds, gs2, mygeneticprofile, mycaselist) {
    	try(dd <- cgdsr::getProfileData(x=mycgds, genes=gs2[x], geneticProfiles=mygeneticprofile, caseList=mycaselist))
      if (class(try) == "try-error") {
        dd <- NULL
      } else {
        cix <- apply(dd, 2, function(x) { return(all(is.na(x))) })
      	dd <- dd[ , !cix, drop=FALSE]
        if(nrow(dd) == 0 || ncol(dd) == 0) {
          dd <- NULL
        } else {
          dd <- data.matrix(dd)
          dd[!is.na(dd) & dd == "NaN"] <- NA
        }
      }
    	return (dd)
    }, mycgds=mycgds, gs2=gs2, mygeneticprofile=mygeneticprofile, mycaselist=mycaselist, mc.cores=nbcore)
    mcres <- mcres[!sapply(mcres, is.null)]
    data.cnv <- do.call(cbind, mcres)
    # rix <- apply(data.cnv, 1, function(x) { return(all(is.na(x))) })
    # data.cnv <- data.cnv[!rix, , drop=FALSE]
    data.cnv <- data.matrix(data.cnv)
    save(list="data.cnv", compress=TRUE, file=myfn2)
  } else { load(myfn2) }
  
  ## retrieval of copy number alterations
  message("Retrieve CNA data for TCGA")
  myfn2 <- file.path(saveres, "tcga_bc_cna.RData")
  if (!file.exists(myfn2)) {
    # mygeneticprofile <- geneticprofile[1, 1]
    mygeneticprofile <- "brca_tcga_gistic"
    ## get APOBEC3B copy number alteration
    apobec3b.cna <- cgdsr::getProfileData(x=mycgds, genes="APOBEC3B", geneticProfiles=mygeneticprofile, caseList=mycaselist)
    nn <- rownames(apobec3b.cna)
    apobec3b.cna <- apobec3b.cna[ , 1, drop=TRUE]
    names(apobec3b.cna) <- nn
    ## parallel version tpo get all the copy number alterations
    splitix <- parallel::splitIndices(nx=length(gs2), ncl=ceiling(length(gs2) / 100))
    mcres <- parallel::mclapply(splitix, function(x, mycgds, gs2, mygeneticprofile, mycaselist) {
    	try(dd <- cgdsr::getProfileData(x=mycgds, genes=gs2[x], geneticProfiles=mygeneticprofile, caseList=mycaselist))
      if (class(try) == "try-error") {
        dd <- NULL
      } else {
        cix <- apply(dd, 2, function(x) { return(all(is.na(x))) })
      	dd <- dd[ , !cix, drop=FALSE]
        if(nrow(dd) == 0 || ncol(dd) == 0) {
          dd <- NULL
        } else {
          dd <- data.matrix(dd)
          dd[!is.na(dd) & dd == "NaN"] <- NA
        }
      }
    	return (dd)
    }, mycgds=mycgds, gs2=gs2, mygeneticprofile=mygeneticprofile, mycaselist=mycaselist, mc.cores=nbcore)
    mcres <- mcres[!sapply(mcres, is.null)]
    data.cna <- do.call(cbind, mcres)
    # rix <- apply(data.cna, 1, function(x) { return(all(is.na(x))) })
    # data.cna <- data.cna[!rix, , drop=FALSE]
    data.cna <- data.matrix(data.cna)
    save(list="data.cna", compress=TRUE, file=myfn2)
  } else { load(myfn2) }

  ## retrieval of mutations
  message("Retrieve mutation data for TCGA")
  myfn2 <- file.path(saveres, "tcga_bc_mut.RData")
  if (!file.exists(myfn2)) {
    # mygeneticprofile <- getGeneticProfiles(mycgds,mycancerstudy)[8, 1]
    mygeneticprofile <- "brca_tcga_mutations"
    ## get APOBEC3B mutations
    apobec3b.mut <- cgdsr::getProfileData(x=mycgds, genes="APOBEC3B", geneticProfiles=mygeneticprofile, caseList=mycaselist)
    nn <- rownames(apobec3b.mut)
    apobec3b.mut <- apobec3b.mut[ , 1, drop=TRUE]
    names(apobec3b.mut) <- nn
    ## parallel version
    splitix <- parallel::splitIndices(nx=length(gs2), ncl=ceiling(length(gs2) / 100))
    mcres <- parallel::mclapply(splitix, function(x, mycgds, gs2, mygeneticprofile, mycaselist) {
    	try(dd <- cgdsr::getProfileData(x=mycgds, genes=gs2[x], geneticProfiles=mygeneticprofile, caseList=mycaselist))
      if (class(try) == "try-error") {
        dd <- NULL
      } else {
        cix <- apply(dd, 2, function(x) { return(all(is.na(x))) })
      	dd <- dd[ , !cix, drop=FALSE]
        if(nrow(dd) == 0 || ncol(dd) == 0) {
          dd <- NULL
        } else {
          nn <- dimnames(dd)
          dd <- apply(dd, 2, as.character)
          dimnames(dd) <- nn
          dd[!is.na(dd) & dd == "NaN"] <- NA
        }
      }
    	return(dd)
    }, mycgds=mycgds, gs2=gs2, mygeneticprofile=mygeneticprofile, mycaselist=mycaselist, mc.cores=nbcore)
    mcres <- mcres[!sapply(mcres, is.null)]
    data.mut <- do.call(cbind, mcres)
    # rix <- apply(data.mut, 1, function(x) { return(all(is.na(x))) })
    # data.mut <- data.mut[!rix, , drop=FALSE]
    ## concatenate mutations for each gene
    dupln <- sort(unique(colnames(data.mut)[duplicated(colnames(data.mut))]))
    uniqn <- colnames(data.mut)[!is.element(colnames(data.mut), dupln)]
    dd <- matrix(NA, ncol=length(uniqn) + length(dupln), nrow=nrow(data.mut), dimnames=list(rownames(data.mut), c(uniqn, dupln)))
    dd[rownames(data.mut), uniqn] <- data.mut[ , uniqn, drop=FALSE]
    ## merge genes with potentially multiple mutations
    splitix <- parallel::splitIndices(nx=length(dupln), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, dupln, data.mut) {
      rr <- sapply(dupln[x], function (x, y) {
        yy <- y[ , which(colnames(y) == x), drop=FALSE]
        yy <- apply(yy, 1, function (x) {
          if (all(is.na(x))) {
            return (NA)
          } else {
            return (paste(sort(unique(as.character(x[!is.na(x)]))), collapse="///"))
          }
        })
        return (yy)
      }, y=data.mut)
      return (rr)
    }, dupln=dupln, data.mut=data.mut)
    ddupl <- do.call(cbind, mcres)
    dd[rownames(ddupl), colnames(ddupl)] <- ddupl
    data.mut <- dd
    save(list="data.mut", compress=TRUE, file=myfn2)
  } else { load(myfn2) }

  ## get the clinical data
  message("Retrieve clinical data for TCGA")
  cn <- c("age", "er", "her2", "pgr", "node", "grade", "size", "t.rfs", "e.rfs", "t.os", "e.os", "treatment", "tissue", "dataset", "id")
  # data.clin <- cgdsr::getClinicalData(mycgds, mycaselist)
  ## load curated TCGA clinical data form Mauro Delorenzi's lab
  load(file.path("data", "tcga_breast_clin_curated.RData"))
  data.clin <- C1
  ## problem with her2 status
  data.clin[ , "her2"] <- NA
  colnames(data.clin) <- gsub("_", ".", colnames(C1))
  rownames(data.clin) <- gsub("-", ".", rownames(data.clin))
  cn <- cn[!is.element(cn, colnames(data.clin))]
  tt <- matrix(NA, nrow=nrow(data.clin), ncol=length(cn), dimnames=list(rownames(data.clin), cn))
  data.clin <- data.frame(data.clin, tt)
  data.clin[ , "t.rfs"] <- (as.numeric(data.clin[ , "t.rfs"]) / 12) * 365
  data.clin[ , "t.os"] <- (as.numeric(data.clin[ , "t.os"]) / 12) * 365

  ## merge all data
  pid <- sort(fold(union, rownames(data.ge), rownames(data.cnv), rownames(data.cna), rownames(data.mut), rownames(data.clin)))
  ## keep only gene symbols for which there exist an entrez gene id
  gid <- sort(fold(union, colnames(data.ge), colnames(data.cnv), colnames(data.cna), colnames(data.mut)))
  gs <- toTable(org.Hs.egSYMBOL)
  gs <- gs[!is.na(gs[ , "symbol"]) & !duplicated(gs[ , "symbol"]) & !is.na(gs[ , "gene_id"]) & !duplicated(gs[ , "gene_id"]), , drop=FALSE]
  rownames(gs) <- gs[ , "symbol"]
  gid <- gid[is.element(gsub("[.]", "-", gid), rownames(gs))]
  ## gene expression dtaa
  dd <- matrix(NA, ncol=length(gid), nrow=length(pid), dimnames=list(pid, gid))
  dd[rownames(data.ge), intersect(gid, colnames(data.ge))] <- data.ge[ , intersect(gid, colnames(data.ge))]
  data.ge <- dd
  ## copy number variation data
  dd <- matrix(NA, ncol=length(gid), nrow=length(pid), dimnames=list(pid, gid))
  dd[rownames(data.cnv), intersect(gid, colnames(data.cnv))] <- data.cnv[ , intersect(gid, colnames(data.cnv))]
  data.cnv <- dd
  ## copy number alteration data
  dd <- matrix(NA, ncol=length(gid), nrow=length(pid), dimnames=list(pid, gid))
  dd[rownames(data.cna), intersect(gid, colnames(data.cna))] <- data.cna[ , intersect(gid, colnames(data.cna))]
  data.cna <- dd
  ## mutation data
  dd <- matrix(NA, ncol=length(gid), nrow=length(pid), dimnames=list(pid, gid))
  dd[rownames(data.mut), intersect(gid, colnames(data.mut))] <- data.mut[ , intersect(gid, colnames(data.mut))]
  data.mut <- dd
  ## clinical data
  dd <- data.frame(matrix(NA, ncol=ncol(data.clin), nrow=length(pid), dimnames=list(pid, colnames(data.clin))))
  newlev <- sapply(data.clin, levels)
  dd <- genefu::setcolclass.df(df=dd, colclass=sapply(data.clin, class), factor.levels=newlev)
  dd[rownames(data.clin), ] <- data.clin
  data.clin <- dd
  data.clin[ , "tissue"] <- "TUMOR"
  data.clin[ , "dataset"] <- "TCGA"
  ## replace "." by "-" in patient and gene ids
  colnames(data.ge) <- colnames(data.cnv) <- colnames(data.cna) <- colnames(data.mut) <- paste("geneid", gs[gsub("[.]", "-", gid), "gene_id"], sep=".")
  rownames(data.ge) <- rownames(data.cnv) <- rownames(data.cna) <- rownames(data.mut) <- rownames(data.clin) <- gsub("[.]", "-", pid)
  annot <- gs[gsub("[.]", "-", gid), , drop=FALSE]
  colnames(annot) <- c("ENTREZID", "SYMBOL")
  rownames(annot) <- paste("geneid", as.character(annot[ , "ENTREZID"]), sep=".")
  ## save data
  save(list=c("data.ge", "data.cnv", "data.cna", "data.mut", "annot", "data.clin"), compress=TRUE, file=myfn)
} else { load(myfn) }

annot <- data.frame("probe"=rownames(annot), "EntrezGene.ID"=MetaGx::stripWhiteSpace(as.character(annot[ , "ENTREZID"])), annot)
rownames(annot) <- as.character(annot$probe)

## molecular subtyping using MetaGx
message("molecular subtyping for TCGA")
## build exprs object
esetge <- ExpressionSet(assayData=t(data.ge), phenoData=AnnotatedDataFrame(data=data.clin), featureData=AnnotatedDataFrame(data=annot))
esetge <- subtypeClassification(eset=esetge, model=sbtm)
ss <- MetaGx::getSubtype(eset=esetge, method="class")
ss.proba <- MetaGx::getSubtype(eset=esetge, method="fuzzy")

# colnames(aa)[colnames(aa) == "ENTREZID"] <- "EntrezGene.ID"
# pdf(file.path(saveres, "tcga_bc_subtyping.pdf"), width=7, height=7)
# ss <- genefu::subtype.cluster.predict(sbt.model=scmod1.robust, data=data.ge, annot=annot, do.mapping=TRUE, plot=TRUE, verbose=FALSE)
# dev.off()

## double check subtyping
print(table("SCM"=ss, "ER"=data.clin[ , "er"]))
print(table("SCM"=ss, "HER2"=data.clin[ , "her2"]))
# print(table("SCM"=ss$subtype2, "GRADE"=data.clin[ , "grade"]))
survd <- survcomp::censor.time(surv.time=as.numeric(data.clin[ , "t.rfs"]) / 365, surv.event=as.numeric(data.clin[ , "e.rfs"]), time.cens=0)
dd <- data.frame("time"=survd[[1]], "event"=survd[[2]], "subtype"=ss)
pdf(file.path(saveres, "tcga_bc_subtyping_km.pdf"), width=7, height=7)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ subtype), data.s=dd, x.label="Time (years)", y.label="Probability of DFS", main.title=sprintf("Molecular subtypes (TCGA)"), sub.title=NULL, leg.text=levels(dd$subtype), leg.pos="bottomright", leg.inset=0.0, v.line=NULL, h.line=NULL, .col=rainbow(length(levels(dd$subtype))), .lty=1, show.n.risk=TRUE, n.risk.step=3, n.risk.cex=0.85, bty="n", leg.bty="n", verbose=FALSE)
dev.off()

tcga <- list("annot"=annot, "ge"=data.ge, "cnv"=data.cnv, "cna"=data.cna, "mut"=data.mut, "clin"=data.clin, "subtype"=ss, "subtype.proba"=ss.proba)
rm(list=c("data.ge", "data.cnv", "data.cna", "data.mut", "annot", "data.clin", "ss"))
gc(TRUE)


## end

