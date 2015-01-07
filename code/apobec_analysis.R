#################################################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## October 31, 2014
#################################################


## subtypes
sbt2 <- dataset$subtype

#################################################
## computation of published signature genes
#################################################

## CIN70
cin70.score <- list("score"=genefu::sig.score(x=cin70t, data=dataset$ge, annot=dataset$annot, do.mapping=FALSE)$score, "risk"=NULL)
cin70.score$score <- (genefu::rescale(cin70.score$score, q=0.05, na.rm=TRUE) - 0.5) * 2
cin70.score$risk <- as.numeric(cut2(x=cin70.score$score, cuts=median(cin70.score$score, na.rm=TRUE))) - 1
names(cin70.score$risk) <- names(cin70.score$score)

## GGI
grade <- dataset$clin[ , "grade"]
if (sum(grade == 1, na.rm=TRUE) > 3 && sum(grade == 3, na.rm=TRUE) > 3) {
  ggi.score <- genefu::ggi(data=dataset$ge, annot=dataset$annot, do.mapping=TRUE, hg=grade)
} else {
  ggi.score <- list("score"=genefu::ggi(data=dataset$ge, annot=dataset$annot, do.mapping=TRUE)$score, "risk"=NULL)
  # ggi.score$score <- (ggi.score$score - median(ggi.score$score, na.rm=TRUE)) / (quantile(ggi.score$score, probs=0.75) - quantile(ggi.score$score, probs=0.255))
 #  ggi.score$risk <- as.numeric(ggi.score$score > 0)
 #  names(ggi.score$risk) <- names(ggi.score$score)
  ggi.score$score <- (genefu::rescale(ggi.score$score, q=0.05, na.rm=TRUE) - 0.5) * 2
  ggi.score$risk <- as.numeric(cut2(x=ggi.score$score, cuts=median(ggi.score$score, na.rm=TRUE))) - 1
  names(ggi.score$risk) <- names(ggi.score$score)
}

## PIK3CAGS
pik3cags.score <- list("score"=genefu::pik3cags(data=dataset$ge, annot=dataset$annot, do.mapping=TRUE), "risk"=NULL)
pik3cags.score$score <- (genefu::rescale(pik3cags.score$score, q=0.05, na.rm=TRUE) - 0.5) * 2
pik3cags.score$risk <- as.numeric(cut2(x=pik3cags.score$score, cuts=median(pik3cags.score$score, na.rm=TRUE))) - 1
names(pik3cags.score$risk) <- names(pik3cags.score$score)

## gene70
gene70.score <- genefu::gene70(data=dataset$ge, annot=dataset$annot, do.mapping=TRUE, std="robust")

## oncotype DX
sig.oncotypedx.save <- sig.oncotypedx
sig.oncotypedx["GSTM1", "EntrezGene.ID"] <- "2948"
## the affymetrix probes for GSTM1 are ambiguious with GSTM2, GSTM3, and GSTM4; GSTM3 is present in the data
## http://xavierlab2.mgh.harvard.edu/EnrichmentProfiler/primary/Enrichment/204550_x_at.html
oncotype.score <- genefu::oncotypedx(data=dataset$ge, annot=dataset$annot, do.mapping=TRUE)
sig.oncotypedx <- sig.oncotypedx.save

## RORS
rors.score <- genefu::rorS(data=dataset$ge, annot=dataset$annot, do.mapping=TRUE)

## endopredict
endopredict.score <- genefu::endoPredict(data=dataset$ge, annot=dataset$annot, do.mapping=TRUE)

## gene modules (desmedt et al, CCR, 2008)
module.score <- NULL
for (i in 1:length(genefu::mod1)) {
  mod.score <- list("score"=genefu::sig.score(x=genefu::mod1[[i]], data=dataset$ge, annot=dataset$annot, do.mapping=TRUE)$score, "risk"=NULL)
  mod.score$score <- (genefu::rescale(mod.score$score, q=0.05, na.rm=TRUE) - 0.5) * 2
  mod.score$risk <- as.numeric(cut2(x=mod.score$score, cuts=median(mod.score$score, na.rm=TRUE))) - 1
  names(mod.score$risk) <- names(mod.score$score)
  module.score <- c(module.score, list(mod.score))
}
names(module.score) <- paste(toupper(names(genefu::mod1)), "MOD", sep="")

## combine signatures
bc.signatures <- c(list("CIN70"=cin70.score, "GGI"=ggi.score, "PIK3CAGS"=pik3cags.score, "GENE70"=gene70.score, "ONCOTYPE"=oncotype.score, "RORS"=rors.score, "ENDOPREDICT"=endopredict.score), module.score)

## association with node, grade, er, her2
cn <- c("node", "grade", "er", "her2")
for (i in 1:length(cn)) {
  x <- dataset$clin[ , cn[i]]
  x[!is.na(x) & x == "NA"] <- NA
  x <- factor(x)
  names(x) <- rownames(dataset$clin)
  if (sum(complete.cases(x)) >= 3) {
    pdf(file.path(saveres, dataset.name, sprintf("%s_subtype_%s.pdf", cn[i], dataset.name)), height=4, width=6)
      tt <- table(x, sbt2)
      names(dimnames(tt)) <- c(cn[i], "subtype")
      if(length(tt) > 1 & !all(tt == 0)) {
        tts <- vcd::assocstats(tt)
        plot.new()
        plotrix::addtable2plot(x=par("usr")[1] + 0.15, y=par("usr")[4] - 0.5, xjust=0, yjust=0, table=as.matrix(tt), bty="o", display.rownames=TRUE, display.colnames=TRUE, hlines=TRUE, vlines=TRUE, title=sprintf("%s vs. subtype\n(%i patients, %s)", cn[i], sum(tt), dataset.name), cex=1.25, xpad=1, ypad=1)
        title(sub=sprintf("Cramer's V=%.2g, p=%.1E", tts$cramer, tts$chisq_tests[1,3]), cex=1)
      }
    dev.off()
  }
}

#################################################
## association with APOBEC3B expression
#################################################

goidir <- file.path(saveres, dataset.name, "exprs")
if (!file.exists(goidir)) { dir.create(goidir, showWarnings=FALSE, recursive=TRUE) }

yylim <- floor(range(dataset$ge[ , intersect(c(goi, goi.prolif), colnames(dataset$ge))], na.rm=TRUE))

#################################################
## subtype-specific co-expression  analysis
#################################################

## compute spearman correlation for co-expresison
myfn <- file.path(goidir, sprintf("coexpression.res.%s", dataset.name))
if (!file.exists(myfn)) {
	for (ii in 1:length(goi)) {
	    if (goi[ii] %in% colnames(dataset$ge)) {
			 a.exprs <- dataset$ge[ , goi[ii]]
			if (sum(complete.cases(a.exprs)) >= 3) {
			  res.all <- NULL
			  for (j in 1:length(sbtoi)) {
			    myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
			    res <- t(apply(dataset$ge[myx, , drop=FALSE], 2, function (x, y) {
					 rr <- c("rho"=NA, "p"=NA)
					 if (sum(complete.cases(x, y)) >= 3) {
						 rr <- cor.test(x, y, use="pairwise.complete.obs", method="spearman")
						 rr <- c("rho"=rr$estimate, "p"=rr$p.value)
					 }
					 return (rr)
				 }, y=a.exprs[myx]))
				 colnames(res) <- c("rho", "p")
				 res <- data.frame("Symbol"=dataset$annot[rownames(res), "SYMBOL"], res)
				 res <- res[order(res[ , "rho"], decreasing=TRUE, na.last=TRUE), , drop=FALSE]
				 res.all <- c(res.all, list(res))
			  }
			  names(res.all) <- names(sbtoi)
			  if (length(res.all) > 0) {
			    WriteXLS::WriteXLS(x="res.all", ExcelFileName=file.path(goidir, sprintf("coexpression_%s_%s.xls", names(goi)[ii], dataset.name)), row.names=TRUE, BoldHeaderRow=TRUE)
			  }
			}
		}
	}
	file.create(myfn)
}	

## association with subtypes
pdf(file.path(goidir, sprintf("subtype_apobec_exprs_%s.pdf", dataset.name)), height=6, width=8)
for (i in 1:length(goi)) {
  x <- dataset$subtype
  if (goi[i] %in% colnames(dataset$ge)) {
    a.exprs <- dataset$ge[ , goi[i]]
    iix <- intersect(names(x), names(a.exprs))
    kwt <- kruskal.test(a.exprs[iix] ~ x[iix])
    boxplot(a.exprs[iix] ~ x[iix], ylab=sprintf("%s expression", names(goi)[i]), sub=sprintf("Kruskal-Walis test p-value = %.1E", kwt$p.value), main=sprintf("%s vs subtypes in %s\n%i patients", names(goi)[i], dataset.name, length(iix)), ylim=yylim)
  }
}
dev.off()

## association with proliferation-related genes
pdf(file.path(goidir, sprintf("prolif_apobec_exprs_%s.pdf", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
for (ii in 1:length(goi.prolif)) {
  if (goi.prolif[ii] %in% colnames(dataset$ge)) {
    x <- dataset$ge[ , goi.prolif[[ii]]]
    if (sum(complete.cases(x)) > 10) {
      for (i in 1:length(goi)) {
        if (goi[i] %in% colnames(dataset$ge)) {
          a.exprs <- dataset$ge[ , goi[i]]
          par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
          iix <- intersect(names(x), names(a.exprs))
          for (j in 1:length(sbtoi)) {
            myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
            myx <- intersect(myx, iix)
            if(length(myx) > 10) {
              cc <- cor.test(x=a.exprs[myx], y=x[myx], methpd="spearman", conf.level=0.95)
              cc <- list("rho"=cc$estimate, "lower"=cc$conf.in[1], "upper"=cc$conf.in[2], "p"=cc$p.value)
            } else { cc <- list("rho"=NA, "lower"=NA, "upper"=NA, "p"=NA) }
            plot(x=a.exprs[myx], y=x[myx], ylab=sprintf("%s expression", names(goi.prolif)[ii]), xlab=sprintf("%s expression", names(goi)[i]), cex.main=0.9, main=sprintf("%s vs %s in %s (%s)\n%i patients", names(goi.prolif)[ii], names(goi)[i], names(sbtoi)[j], dataset.name, length(myx)), sub=sprintf("Spearman = %.3g 95%%CI [%.3g, %.3g], p = %.1E", cc$rho, cc$lower, cc$upper, cc$p), pch=20, col="darkgrey", xlim=yylim, ylim=yylim)
          }
        }
      }
    }
  }
}
dev.off()
  
## breast cancer published signatures
pdf(file.path(goidir, sprintf("signatures_apobec_exprs_%s.pdf", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
for (ii in 1:length(bc.signatures)) {
  x <- bc.signatures[[ii]]$score
  if (sum(complete.cases(x)) > 10) {
    for (i in 1:length(goi)) {
      if (goi[i] %in% colnames(dataset$ge)) {
        a.exprs <- dataset$ge[ , goi[i]]
        par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
        iix <- intersect(names(x), names(a.exprs))
        for (j in 1:length(sbtoi)) {
          myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
          myx <- intersect(myx, iix)
          if(length(myx) > 10) {
            cc <- cor.test(x=a.exprs[myx], y=x[myx], methpd="spearman", conf.level=0.95)
            cc <- list("rho"=cc$estimate, "lower"=cc$conf.in[1], "upper"=cc$conf.in[2], "p"=cc$p.value)
          } else { cc <- list("rho"=NA, "lower"=NA, "upper"=NA, "p"=NA) }
          plot(x=a.exprs[myx], y=x[myx], ylab=sprintf("%s score", names(bc.signatures)[ii]), xlab=sprintf("%s expression", names(goi)[i]), cex.main=0.9, main=sprintf("%s vs %s in %s (%s)\n%i patients", names(bc.signatures)[ii], names(goi)[i], names(sbtoi)[j], dataset.name, length(myx)), sub=sprintf("Spearman = %.3g 95%%CI [%.3g, %.3g], p = %.1E", cc$rho, cc$lower, cc$upper, cc$p), pch=20, col="darkgrey", xlim=yylim, ylim=yylim)
        }
      }
    }
  }
}
dev.off()


## association with age
x <- as.numeric(dataset$clin[ , "age"])
names(x) <- rownames(dataset$clin)
if (sum(complete.cases(x)) >= 3) {
  pdf(file.path(goidir, sprintf("%s_apobec_exprs_%s.pdf", "age", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
  for (i in 1:length(goi)) {
    if (goi[i] %in% colnames(dataset$ge)) {
      a.exprs <- dataset$ge[ , goi[i]]
      par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
      iix <- intersect(rownames(dataset$clin), names(a.exprs))
      for (j in 1:length(sbtoi)) {
        myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
        myx <- intersect(myx, iix)
        if(length(myx) > 10) {
          cc <- cor.test(x=a.exprs[myx], y=x[myx], methpd="spearman", conf.level=0.95)
          cc <- list("rho"=cc$estimate, "lower"=cc$conf.in[1], "upper"=cc$conf.in[2], "p"=cc$p.value)
        } else { cc <- list("rho"=NA, "lower"=NA, "upper"=NA, "p"=NA) }
        plot(x=a.exprs[myx], y=x[myx], ylab="Age at diagnosis", xlab=sprintf("%s expression", names(goi)[i]), cex.main=0.9, main=sprintf("%s vs %s in %s (%s)\n%i patients", "Age", names(goi)[i], names(sbtoi)[j], dataset.name, length(myx)), sub=sprintf("Spearman = %.3g 95%%CI [%.3g, %.3g], p = %.1E", cc$rho, cc$lower, cc$upper, cc$p), pch=20, col="darkgrey", xlim=yylim)
      }
    }
  }
  dev.off()
}

## association with tumor size
x <- as.numeric(dataset$clin[ , "size"])
names(x) <- rownames(dataset$clin)
if (sum(complete.cases(x)) >= 3) {
  pdf(file.path(goidir, sprintf("%s_apobec_exprs.pdf", "size", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
  for (i in 1:length(goi)) {
    if (goi[i] %in% colnames(dataset$ge)) {
      a.exprs <- dataset$ge[ , goi[i]]
      par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
      iix <- intersect(rownames(dataset$clin), names(a.exprs))
      for (j in 1:length(sbtoi)) {
        myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
        myx <- intersect(myx, iix)
        if(sum(complete.cases(a.exprs[myx], x[myx])) > 3) {
          cc <- cor.test(x=a.exprs[myx], y=x[myx], methpd="spearman", conf.level=0.95)
          cc <- list("rho"=cc$estimate, "lower"=cc$conf.in[1], "upper"=cc$conf.in[2], "p"=cc$p.value)
          plot(x=a.exprs[myx], y=x[myx], ylab="Tumor size", xlab=sprintf("%s expression", names(goi)[i]), cex.main=0.9, main=sprintf("%s vs %s in %s (%s)\n%i patients", "Tumor size", names(goi)[i], names(sbtoi)[j], dataset.name, length(myx)), sub=sprintf("Spearman = %.3g 95%%CI [%.3g, %.3g], p = %.1E", cc$rho, cc$lower, cc$upper, cc$p), pch=20, col="darkgrey", xlim=yylim)
        } else { cc <- list("rho"=NA, "lower"=NA, "upper"=NA, "p"=NA) }
      }
    }
  }
  dev.off()
}

## association with node, grade, er, her2
iix <- intersect(rownames(dataset$clin), names(a.exprs))
cn <- c("node", "grade", "er", "her2", "age.bin", "size.bin")
for (i in 1:length(cn)) {
  x <- dataset$clin[ , cn[i]]
  x[!is.na(x) & x == "NA"] <- NA
  x <- factor(x, levels=c("", sort(as.character(unique(x)))))
  names(x) <- rownames(dataset$clin)
  if (sum(complete.cases(x)) >= 3) {
    pdf(file.path(goidir, sprintf("%s_apobec_exprs_%s.pdf", cn[i], dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
    for (ii in 1:length(goi)) {
      if (goi[[ii]] %in% colnames(dataset$ge)) {
        a.exprs <- dataset$ge[ , goi[[ii]]]
        doplot <- FALSE
        for (j in 1:length(sbtoi)) {
          myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
          myx <- intersect(myx, iix)
          tt <- table(x[myx])
          tt <- tt[names(tt) != ""]
          if(length(tt) > 1 & all(tt > 3)) {
            if (!doplot) { par(mfrow=c(3, ceiling(length(sbtoi) / 3))) }
            doplot <- TRUE
            wt <- kruskal.test(a.exprs[myx] ~ x[myx])
            boxplot(a.exprs[myx] ~ x[myx], outline=FALSE, ylab=sprintf("%s expression", names(goi)[ii]), xlab=cn[i], cex.main=0.9, main=sprintf("%s vs %s in %s (%s)\n%i patients", cn[i], names(goi)[ii], names(sbtoi)[j], dataset.name, length(myx)), sub=sprintf("Kruskal-Wallis test p-value = %.1E", wt$p.value), ylim=yylim)
            legend("topleft", title=cn[i], legend=sprintf("%s (n=%i)", names(tt), tt), bty="n")
          } else { wt <- list("p.value"=NA) }
        }
      }
    }
    dev.off()
  }
}
## combined boxplot
pdf(file.path(goidir, sprintf("clin_apobec_exprs_%s.pdf", dataset.name)), height=3, width=2.5)
for (ii in 1:length(goi)) {
	if (goi[[ii]] %in% colnames(dataset$ge)) {
	  a.exprs <- dataset$ge[ , goi[[ii]]]
	  for (j in 1:length(sbtoi)) {
		myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
		myx <- intersect(myx, iix)
		## build a list of data points for each category
		rr.all <- NULL
		for (i in 1:length(cn)) {
			x <- dataset$clin[ , cn[i]]
			x[!is.na(x) & (x == "NA" | x == "")] <- NA
			if (!is.factor(x)) {
				x <- factor(x, levels=c(sort(as.character(unique(x)))))
			}
			names(x) <- rownames(dataset$clin)
			if (sum(complete.cases(x[myx])) >= 3) {
				ll <- levels(x)
				ixx <- table(x)[ll] > 0
				ll <- ll[!is.na(ixx) & ixx > 0]
				rr <- lapply(ll, function (x, y, z) {
					return (y[!is.na(z) & z == x])
				}, y=a.exprs[myx], z=x[myx])
				names(rr) <- paste(cn[i], ll, sep=".")
				rr.all <- c(rr.all, list(NA), rr)
			}
		}
		rr.all <- rr.all[-1]
      # wt <- kruskal.test(a.exprs[myx] ~ x[myx])
		par(xaxt="n", las=1, mar=c(7, 4, 3, 1) + 0.1, cex=0.7)
      mp <- boxplot(rr.all, outline=FALSE, ylab=sprintf("%s expression", names(goi)[ii]), xlab="Clinical parameters", cex.main=0.9, main=sprintf("%s in %s", names(goi)[ii], names(sbtoi)[j]), ylim=yylim, col=sbtcol[j])
	   axis(1, at=seq(1, length(rr.all), by=1), labels=FALSE)
	   text(x=(1:length(mp$names)) + 0.25, y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=toupper(mp$names), srt=45, xpd=NA, cex=0.7, font=1)
	}
 }
}
dev.off()






## combined boxplot per subtypes
pdf(file.path(goidir, sprintf("clin3_apobec_exprs_%s.pdf", dataset.name)), height=3, width=9)
for (ii in 1:length(goi)) {
	if (goi[[ii]] %in% colnames(dataset$ge)) {
		a.exprs <- dataset$ge[ , goi[[ii]]]
		rr.all3 <- mycol3 <- NULL
		for (i in 1:length(cn)) {
			x <- dataset$clin[ , cn[i]]
			x[!is.na(x) & (x == "NA" | x == "")] <- NA
			if (!is.factor(x)) {
				x <- factor(x, levels=c(sort(as.character(unique(x)))))
			}
			names(x) <- rownames(dataset$clin)
			rr.all2 <- mycol2 <- NULL
			for (j in 1:length(sbtoi)) {
				rr.all <- NULL
				for (k in 1:length(levels(x))) {
					iix2 <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
					iix3 <- names(x)[!is.na(x) & x == levels(x)[k]]
					myx <- fold(intersect, iix3, iix2, iix)
					## build a list of data points for each category
					if (sum(complete.cases(x[myx])) >= 3) {
						rr <- list(a.exprs[myx])
						names(rr) <- paste(cn[i], levels(x)[k], names(sbtoi)[j], sep=".")
					} else {
						rr <- -1
					}
					rr.all <- c(rr.all, rr)
				}
				rr.all2 <- c(rr.all2, list(NA), rr.all)
				mycol2 <- c(mycol2, "white", rep(sbtcol[j], each=length(levels(x))))
			}
			rr.all3 <- c(rr.all3, list(NA), rr.all2)
			mycol3 <- c(mycol3, "white", mycol2)
		}
		rr.all3 <- rr.all3[-c(1:2)]
		mycol3 <- mycol3[-c(1:2)]
		myx <- sapply(rr.all3, function (x) { return (!(!is.na(x[1]) & x[1] == -1))})
		rr.all3 <- rr.all3[myx]
		mycol3 <- mycol3[myx]
		par(xaxt="n", las=1, mar=c(7, 4, 3, 1) + 0.1, cex=0.7)
		mp <- boxplot(rr.all3, outline=FALSE, ylab=sprintf("%s expression", names(goi)[ii]), xlab="", cex.main=0.9, main=sprintf("%s", names(goi)[ii]), col=mycol3)
		axis(1, at=seq(1, length(rr.all3), by=1), labels=FALSE)
		text(x=(1:length(mp$names)) + 0.25, y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=toupper(mp$names), srt=45, xpd=NA, cex=0.7, font=1)
	}
}
dev.off()








## association with mutations
iix <- intersect(rownames(dataset$ge), rownames(mutations))
pdf(file.path(goidir, sprintf("mutations_apobec_exprs_%s.pdf", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
for (ii in 1:ncol(mutations)) {
  x <- mutations[ , ii]
  names(x) <- rownames(mutations)
  xx <- as.numeric(x)
  xx <- factor(x, levels=c("", sort(unique(as.character(xx)))))
  names(xx) <- names(x)
  if (sum(complete.cases(x)) > 10) {
    for (i in 1:length(goi)) {
      if (goi[i] %in% colnames(dataset$ge)) {
        a.exprs <- dataset$ge[ , goi[i]]
        par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
        for (j in 1:length(sbtoi)) {
          myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
          myx <- intersect(myx, iix)
          tt <- table(xx[myx])
          tt <- tt[names(tt) != ""]
          if(length(tt) > 1 & all(tt > 3)) {
            wt <- kruskal.test(a.exprs[myx] ~ xx[myx])
            boxplot(a.exprs[myx] ~ xx[myx], outline=FALSE, ylab=sprintf("%s expression", names(goi)[i]), xlab=sprintf("%s mutation status", colnames(mutations)[ii]), cex.main=0.9, main=sprintf("%s mut vs %s in %s\n(%i patients, %s)", colnames(mutations)[ii], names(goi)[i], names(sbtoi)[j], length(myx), dataset.name), sub=sprintf("Kruskal-Wallis test p-value = %.1E", wt$p.value), ylim=yylim)
            legend("topleft", title=sprintf("%s mut", colnames(mutations)[ii]), legend=sprintf("%s (n=%i)", names(tt), tt), bty="n")
          } else { wt <- list("p.value"=NA) }
        }
      }
    }
  }
}
dev.off()
## combined boxplot per mutation type
mutn <- c("TP53", "PIK3CA", "PIK3CA.E542K", "PIK3CA.E545K")
mutn <- mutn[mutn %in% colnames(mutations)]
if (length(mutn) > 0) {
	
	pdf(file.path(goidir, sprintf("mutations2_apobec_exprs_%s.pdf", dataset.name)), height=3, width=ceiling(length(mutn) * 2.5))
	for (ii in 1:length(goi)) {
		if (goi[ii] %in% colnames(dataset$ge)) {
			a.exprs <- dataset$ge[ , goi[ii]]
			# par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
			rr.all3 <- mycol3 <- NULL
			for (i in 1:length(mutn)) {
				x <- mutations[ , mutn[i]]
				x[!is.na(x) & (x == "NA" | x == "")] <- NA
				if (!is.factor(x)) {
					x <- factor(x, levels=c(sort(as.character(unique(x)))))
				}
				levels(x) <- c("WT", "MUT")
				names(x) <- rownames(mutations)
				rr.all2 <- mycol2 <- NULL
				for (k in 1:length(levels(x))) {
					rr.all <- NULL
					for (j in 1:length(sbtoi)) {
						iix2 <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
						iix3 <- names(x)[!is.na(levels(x)) & levels(x) == levels(x)[k]]
						myx <- fold(intersect, iix3, iix2, iix)
						## build a list of data points for each category
						if (sum(complete.cases(x[myx])) >= 3) {
							rr <- list(a.exprs[myx])
							names(rr) <- paste(mutn[i], levels(x)[k], names(sbtoi)[j], sep=".")
						} else {
							rr <- NA
						}
						rr.all <- c(rr.all, rr)
					}
					rr.all2 <- c(rr.all2, list(NA), rr.all)
					mycol2 <- c(mycol2, "white", sbtcol)
				}
				rr.all3 <- c(rr.all3, list(NA), rr.all2)
				mycol3 <- c(mycol3, "white", mycol2)	
			}
			rr.all3 <- rr.all3[-c(1:2)]
			mycol3 <- mycol3[-c(1:2)]
			# wt <- kruskal.test(a.exprs[myx] ~ x[myx])
			par(xaxt="n", las=1, mar=c(7, 4, 3, 1) + 0.1, cex=0.7)
			mp <- boxplot(rr.all3, outline=FALSE, ylab=sprintf("%s expression", names(goi)[ii]), xlab="", cex.main=0.9, main=sprintf("%s", names(goi)[ii]), col=mycol3, ylim=yylim)
			axis(1, at=seq(1, length(rr.all3), by=1), labels=FALSE)
			text(x=(1:length(mp$names)) + 0.25, y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=toupper(mp$names), srt=45, xpd=NA, cex=0.7, font=1)
		}
	}
	dev.off()

	## combined boxplot per mutation
	pdf(file.path(goidir, sprintf("mutations3_apobec_exprs_%s.pdf", dataset.name)), height=3, width=ceiling(length(mutn) * 2.5))
	for (ii in 1:length(goi)) {
		if (goi[ii] %in% colnames(dataset$ge)) {
			a.exprs <- dataset$ge[ , goi[ii]]
			# par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
			rr.all3 <- mycol3 <- NULL
			for (i in 1:length(mutn)) {
				x <- mutations[ , mutn[i]]
				x[!is.na(x) & (x == "NA" | x == "")] <- NA
				if (!is.factor(x)) {
					x <- factor(x, levels=c(sort(as.character(unique(x)))))
				}
				levels(x) <- c("WT", "MUT")
				names(x) <- rownames(mutations)
				rr.all2 <- mycol2 <- NULL
				for (j in 1:length(sbtoi)) {
					rr.all <- NULL
					for (k in 1:length(levels(x))) {
						iix2 <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
						iix3 <- names(x)[!is.na(x) & x == levels(x)[k]]
						myx <- fold(intersect, iix3, iix2, iix)
						## build a list of data points for each category
						if (sum(complete.cases(x[myx])) >= 3) {
							rr <- list(a.exprs[myx])
							names(rr) <- paste(mutn[i], levels(x)[k], names(sbtoi)[j], sep=".")
						} else {
							rr <- -1
						}
						rr.all <- c(rr.all, rr)
					}
					rr.all2 <- c(rr.all2, list(NA), rr.all)
					mycol2 <- c(mycol2, "white", rep(sbtcol[j], each=length(levels(x))))
				}
				rr.all3 <- c(rr.all3, list(NA), rr.all2)
				mycol3 <- c(mycol3, "white", mycol2)	
			}
			rr.all3 <- rr.all3[-c(1:2)]
			mycol3 <- mycol3[-c(1:2)]
			myx <- sapply(rr.all3, function (x) { return (!(!is.na(x[1]) & x[1] == -1))})
			rr.all3 <- rr.all3[myx]
			mycol3 <- mycol3[myx]
			# wt <- kruskal.test(a.exprs[myx] ~ x[myx])
			par(xaxt="n", las=1, mar=c(7, 4, 3, 1) + 0.1, cex=0.7)
			mp <- boxplot(rr.all3, outline=FALSE, ylab=sprintf("%s expression", names(goi)[ii]), xlab="", cex.main=0.9, main=sprintf("%s", names(goi)[ii]), col=mycol3, ylim=yylim)
			axis(1, at=seq(1, length(rr.all3), by=1), labels=FALSE)
			text(x=(1:length(mp$names)) + 0.25, y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=toupper(mp$names), srt=45, xpd=NA, cex=0.7, font=1)
		}
	}
	dev.off()
}


## Association with grade independently of subtypes (multivariate linear regression model with A3B expression + grade + subtypes)
if (!all(is.na(dataset$clin[ , "grade"]))) {
  sink(file.path(goidir, "apobec_exprs_grade_subtype_lm.txt"))
  for (i in 1:length(goi)) {
    if (goi[i] %in% colnames(dataset$ge)) {
      cat(sprintf("%s\n| %s |\n%s\n", paste(rep("-", nchar(names(goi)[i]) + 4), collapse=""), names(goi)[i], paste(rep("-", nchar(names(goi)[i]) + 4), collapse="")))
      a.exprs <- dataset$ge[ , goi[i]]
      dd <- data.frame("Gene"=a.exprs, "grade"=factor(dataset$clin[ , "grade"], levels=c("1", "2", "3"), ordered=TRUE), "subtype"=sbt2)
      dd <- dd[ , c(TRUE, apply(dd[ , -1], 2, function (x) { return (!all(is.na(x))) })), drop=FALSE]
      if (sum(complete.cases(dd)) > 10) {
        print(summary(glm(formula=Gene ~ ., data=dd, family="gaussian")))
      }
      cat("\n")
    }
  }
  sink()
}

## association with survival
## survival data
survd <- survcomp::censor.time(surv.time=as.numeric(dataset$clin[ , "t.rfs"]) / 365, surv.event=as.numeric(dataset$clin[ , "e.rfs"] == "1"), time.cens=0)
names(survd[[1]]) <- names(survd[[2]]) <- rownames(dataset$clin)
## cox with interaction term for subtypes
rr.all <- NULL
for (i in 1:length(goi)) {
	if (goi[i] %in% colnames(dataset$ge)) {
		a.exprs <- dataset$ge[ , goi[i]]
		cc.ix <- complete.cases(survd[[1]], survd[[2]], a.exprs)
		if (sum(cc.ix) > 10) {
			## force normal ditribution to compute D index
			sx <- sort(a.exprs[cc.ix], decreasing=FALSE)
			kap <- sqrt(8/pi)
			z <- kap^-1 * SuppDists::normOrder(N=sum(cc.ix))
			#ties?
			dup <- duplicated(sx)
			if(any(dup)) {
				udup <- unique(sx[dup])
				for(j in 1:length(udup)) { z[sx == udup[j]] <- mean(z[sx == udup[j]]) }
			}
			oo <- order(a.exprs[cc.ix], decreasing=FALSE)
			z2 <- array(NA, dim=length(a.exprs), dimnames=list(names(a.exprs)))
			z2[cc.ix][oo] <- z
			s <- factor(!is.na(sbt2) & is.element(sbt2, sbtoi[["Luminals"]]))
			levels(s) <- c("nonLuminals", "Luminals")
			dd <- data.frame("time"=survd[[1]], "event"=survd[[2]], "x"=z2, "s"=s)
			rr <- summary(coxph(Surv(time, event) ~ x + s + x * s, data=dd))
			rr <- cbind(rownames(rr$coefficients), round(rr$conf.int[ , 1], digits=2), round(rr$conf.int[ , 3:4], digits=2), rr$n, sprintf("%.1E", rr$coefficients[ , 5]))
			rr <- cbind(c(names(goi)[i], rep("", nrow(rr) - 1)), rr)
			colnames(rr) <- c("Gene", "variable", "hr", "lower", "upper", "n", "p")
			rr.all <- rbind(rr.all, "", data.frame(rr))
		}
	}
}
## save to HTML
xt <- xtable(rr.all)
print.xtable(xt, include.rownames=FALSE, type="html", file=file.path(goidir, sprintf("cox_apobec_exprs_interaclums.html")))

## survival curves
## ternary
pdf(file.path(goidir, sprintf("km3_apobec_exprs_%s.pdf", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
for (i in 1:length(goi)) {
  if (goi[i] %in% colnames(dataset$ge)) {
    a.exprs <- dataset$ge[ , goi[i]]
    par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
    for (j in 1:length(sbtoi)) {
       myx <- !is.na(sbt2) & is.element(sbt2, sbtoi[[j]])
       x <- Hmisc::cut2(a.exprs[myx], cuts=quantile(x=a.exprs[myx], probs=c(0.33, 0.66), na.rm=TRUE))
       levels(x) <- c("Low", "Intermediate", "High")
       dd <- data.frame("time"=survd[[1]][myx], "event"=survd[[2]][myx], "x"=x)
       km.coxph.plot(formula.s=formula(Surv(time, event) ~ x), data.s=dd, x.label="Time (years)", y.label="Probability of DFS", main.title=sprintf("%s expression (%s)\n%s", names(goi)[i], dataset.name, names(sbtoi)[j]), sub.title=NULL, leg.text=levels(x), leg.pos="bottomright", leg.inset=0.0, v.line=NULL, h.line=NULL, .col=c("darkblue", "darkgreen", "darkred"), .lty=c(1, 1), show.n.risk=TRUE, n.risk.step=3, n.risk.cex=0.85, bty="n", leg.bty="n", verbose=FALSE)
    }
  }
}
dev.off()
## binary
pdf(file.path(goidir, sprintf("km2_apobec_exprs_%s.pdf", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
for (i in 1:length(goi)) {
  if (goi[i] %in% colnames(dataset$ge)) {
    a.exprs <- dataset$ge[ , goi[i]]
    par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
    rr.all <- rr2.all <- NULL
    for (j in 1:length(sbtoi)) {
       myx <- !is.na(sbt2) & is.element(sbt2, sbtoi[[j]])
       x <- Hmisc::cut2(a.exprs[myx], cuts=quantile(x=a.exprs[myx], probs=c(0.5), na.rm=TRUE))
       levels(x) <- c("Low", "High")
       dd <- data.frame("time"=survd[[1]][myx], "event"=survd[[2]][myx], "x"=x, "xx"=a.exprs[myx])
       km.coxph.plot(formula.s=formula(Surv(time, event) ~ x), data.s=dd, x.label="Time (years)", y.label="Probability of DFS", main.title=sprintf("%s expression (%s)\n%s", names(goi)[i], dataset.name, names(sbtoi)[j]), sub.title=NULL, leg.text=levels(x), leg.pos="bottomright", leg.inset=0.0, v.line=NULL, h.line=NULL, .col=c("darkblue", "darkred"), .lty=c(1, 1), show.n.risk=TRUE, n.risk.step=3, n.risk.cex=0.85, bty="n", leg.bty="n", verbose=FALSE)
       ## cox models
       rr <- summary(coxph(formula=Surv(time, event) ~ x, data=dd))
       rr <- c(rr$conf.int[1], rr$conf.int[3:4], rr$n, rr$logtest[3])
       names(rr) <- c("hr", "lower", "upper", "n", "p")
       rr.all <- rbind(rr.all, rr)
       rr2 <- survcomp::D.index(x=dd$xx, surv.time=dd$time, surv.event=dd$event, na.rm=TRUE)
       rr2 <- c(rr2$d.index, rr2$lower, rr2$upper, rr2$n, rr2$p.value)
       names(rr2) <- c("hr", "lower", "upper", "n", "p")
       rr2.all <- rbind(rr2.all, rr2)
    }
    ## save to HTML
    rownames(rr.all) <- rownames(rr2.all) <- names(sbtoi)
    tt <- rr2.all
    rownames(tt) <- paste(rownames(tt), " ", sep="")
    rr <- rbind("Binary"=NA, rr.all, " "=NA, "Continuous"=NA, tt)
    xt <- xtable(rr, digits=c(0, 3, 3, 3, 0, -1))
    print.xtable(xt, type="html", file=file.path(goidir, sprintf("cox_%s_exprs.html", names(goi)[i])))
  }
}
dev.off()


## Association with survival independently of clinical variables, subtypes, and proliferation-related genes
survd <- survcomp::censor.time(surv.time=as.numeric(dataset$clin[ , "t.rfs"]) / 365, surv.event=as.numeric(dataset$clin[ , "e.rfs"]), time.cens=0)
names(survd[[1]]) <- names(survd[[2]]) <- rownames(dataset$clin)
sink(file.path(goidir, "cox_multiv_apobec_exprs_lm.txt"))
for (i in 1:length(goi)) {
  if (goi[i] %in% colnames(dataset$ge)) {
    cat(sprintf("%s\n| %s |\n%s\n\n", paste(rep("-", nchar(names(goi)[i]) + 4), collapse=""), names(goi)[i], paste(rep("-", nchar(names(goi)[i]) + 4), collapse="")))
    a.exprs <- dataset$ge[ , goi[i]]
    dd <- data.frame("time"=survd[[1]], "event"=survd[[2]], "gene"=a.exprs, "grade"=factor(dataset$clin[ , "grade"], levels=c("1", "2", "3"), ordered=TRUE), "age"=as.numeric(dataset$clin[ , "age"]), "er"=factor(dataset$clin[ , "er"], levels=c("0", "1"), ordered=TRUE), "node"=factor(dataset$clin[ , "node"], levels=c("0", "1"), ordered=TRUE), "size"=as.numeric(dataset$clin[ , "size"]))
    dd <- dd[ , c(TRUE, TRUE, TRUE, apply(dd[ , -(1:3)], 2, function (x) { return (!all(is.na(x))) })), drop=FALSE]
    for (j in 1:length(sbtoi)) {
      cat(sprintf("%s:\n%s\n", names(sbtoi)[j], paste(rep("-", nchar(names(sbtoi)[j]) + 1), collapse="")))
      myx <- !is.na(sbt2) & is.element(sbt2, sbtoi[[j]])
      if (sum(complete.cases(dd[myx, , drop=FALSE])) > 10) {
        print(summary(coxph(formula=Surv(time, event) ~ ., data=dd[myx, , drop=FALSE])))
      }
      cat("\n")
    }
  }
}
sink()

## Association with survival independently of clinical variables, subtypes, and proliferation-related genes
survd <- survcomp::censor.time(surv.time=as.numeric(dataset$clin[ , "t.rfs"]) / 365, surv.event=as.numeric(dataset$clin[ , "e.rfs"]), time.cens=0)
names(survd[[1]]) <- names(survd[[2]]) <- rownames(dataset$clin)
for (ii in 1:length(goi.prolif)) {
  if (goi.prolif[ii] %in% colnames(dataset$ge)) {
    x <- dataset$ge[ , goi.prolif[ii]]
    sink(file.path(goidir, sprintf("cox_multiv_%s_apobec_exprs_lm_%s.txt", names(goi.prolif)[ii], dataset.name)))
    for (i in 1:length(goi)) {
      if (goi[i] %in% colnames(dataset$ge)) {
        cat(sprintf("%s\n| %s |\n%s\n\n", paste(rep("-", nchar(names(goi)[i]) + 4), collapse=""), names(goi)[i], paste(rep("-", nchar(names(goi)[i]) + 4), collapse="")))
        a.exprs <- dataset$ge[ , goi[i]]
        dd <- data.frame("time"=survd[[1]], "event"=survd[[2]], "gene"=a.exprs, "prolif"=x, "grade"=factor(dataset$clin[ , "grade"], levels=c("1", "2", "3"), ordered=TRUE), "age"=as.numeric(dataset$clin[ , "age"]), "er"=factor(dataset$clin[ , "er"], levels=c("0", "1"), ordered=TRUE), "node"=factor(dataset$clin[ , "node"], levels=c("0", "1"), ordered=TRUE), "size"=as.numeric(dataset$clin[ , "size"]))
        dd <- dd[ , c(TRUE, TRUE, TRUE, TRUE, apply(dd[ , -(1:4)], 2, function (x) { return (!all(is.na(x))) })), drop=FALSE]
        for (j in 1:length(sbtoi)) {
          cat(sprintf("%s:\n%s\n", names(sbtoi)[j], paste(rep("-", nchar(names(sbtoi)[j]) + 1), collapse="")))
          myx <- !is.na(sbt2) & is.element(sbt2, sbtoi[[j]])
          if (sum(complete.cases(dd[myx, , drop=FALSE])) > 10) {
            ## model without the apobec gene
            ff <- sprintf("Surv(time, event) ~ %s", paste(colnames(dd)[!is.element(colnames(dd), c("time", "event", "gene"))], collapse=" + "))
            rr0 <- coxph(formula=formula(ff), data=dd[myx, , drop=FALSE])
            ## model with the apobec gene
            ff <- sprintf("Surv(time, event) ~ %s", paste(colnames(dd)[!is.element(colnames(dd), c("time", "event"))], collapse=" + "))
            rr1 <- coxph(formula=formula(ff), data=dd[myx, , drop=FALSE])
            print(anova(rr0, rr1))
            cat("\n")
            print(summary(coxph(formula=formula(ff), data=dd[myx, , drop=FALSE])))
          }
          cat("\n")
        }
      }
    }
    sink()
  }
}


#################################################
## association with APOBEC3B CNV calls
#################################################

goidir <- file.path(saveres, dataset.name, "cnvcalls")
if (!file.exists(goidir)) { dir.create(goidir, showWarnings=FALSE, recursive=TRUE) }


yylim <- floor(range(dataset$ge[ , intersect(c(goi, goi.prolif), colnames(dataset$ge))], na.rm=TRUE))

## correlation between apobec cnv calls and expression of apobec genes
pdf(file.path(goidir, sprintf("apobec_cnvcalls_%s.pdf", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
xx <- as.character(cnvcalls)
xx[!is.na(xx) & xx == "NA"] <- NA
xx <- factor(xx, levels=c("", sort(unique(as.character(xx)))))
names(xx) <- names(cnvcalls)
iix <- intersect(rownames(dataset$ge), names(cnvcalls))
for (i in 1:length(goi)) {
  ## expression levels of APOBEC genes
  if (goi[i] %in% colnames(dataset$ge)) {
    x <- dataset$ge[ , goi[i]]
    ## expression wrt CNV calls
    par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
    for (j in 1:length(sbtoi)) {
      myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
      myx <- intersect(myx, iix)
      tt <- table(xx[myx])
      tt <- tt[names(tt) != ""]
      if(length(tt) > 1 & all(tt > 3)) {
        wt <- kruskal.test(x[myx] ~ cnvcalls[myx])
        boxplot(x[myx] ~ xx[myx], outline=FALSE, ylab=sprintf("%s expression", names(goi)[i]), xlab="CNV calls", cex.main=0.9, main=sprintf("%s in %s (%s)\n%i patients", names(goi)[i], names(sbtoi)[j], dataset.name, length(myx)), sub=sprintf("Kruskal-Wallis test p-value = %.1E", wt$p.value))
        legend("topleft", title="CNV calls", legend=sprintf("%s (n=%i)", names(tt), tt), bty="n")
      } else { wt <- list("p.value"=NA) }
    }
  }
}
dev.off()
## combined boxplot per subtypes

pdf(file.path(goidir, sprintf("apobec3_cnvcalls_%s.pdf", dataset.name)), height=3, width=4)
for (ii in 1:length(goi)) {
	if (goi[ii] %in% colnames(dataset$ge)) {
		a.exprs <- dataset$ge[ , goi[ii]]
		x <- cnvcalls
		x[!is.na(x) & (x == "NA" | x == "")] <- NA
		rr.all2 <- mycol2 <- NULL
		for (j in 1:length(sbtoi)) {
			rr.all <- NULL
			for (k in 1:length(levels(x))) {
				iix2 <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
				iix3 <- names(x)[!is.na(x) & x == levels(x)[k]]
				myx <- fold(intersect, iix3, iix2, iix)
				## build a list of data points for each category
				if (sum(complete.cases(x[myx])) >= 3) {
					rr <- list(a.exprs[myx])
					names(rr) <- paste(levels(x)[k], names(sbtoi)[j], sep=".")
				} else {
					rr <- -1
				}
				rr.all <- c(rr.all, rr)
			}
			rr.all2 <- c(rr.all2, list(NA), rr.all)
			mycol2 <- c(mycol2, "white", rep(sbtcol[j], each=length(levels(x))))
		}
		rr.all2 <- rr.all2[-c(1)]
		mycol2 <- mycol2[-c(1)]
		myx <- sapply(rr.all2, function (x) { return (!(!is.na(x[1]) & x[1] == -1))})
		rr.all2 <- rr.all2[myx]
		mycol2 <- mycol2[myx]
		par(xaxt="n", las=1, mar=c(7, 4, 3, 1) + 0.1, cex=0.7)
		mp <- boxplot(rr.all2, outline=FALSE, ylab=sprintf("%s expression", names(goi)[ii]), xlab="", cex.main=0.9, main=sprintf("%s", names(goi)[ii]), col=mycol2, ylim=yylim)
		axis(1, at=seq(1, length(rr.all2), by=1), labels=FALSE)
		text(x=(1:length(mp$names)) + 0.25, y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=toupper(mp$names), srt=45, xpd=NA, cex=0.7, font=1)
	}
}
dev.off()

## correlation between apobec cnv calls and expression of proliferation geness
pdf(file.path(goidir, sprintf("prolif_cnvcalls_%s.pdf", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
xx <- as.character(cnvcalls)
xx[!is.na(xx) & xx == "NA"] <- NA
xx <- factor(xx, levels=c("", sort(unique(as.character(xx)))))
names(xx) <- names(cnvcalls)
iix <- intersect(rownames(dataset$ge), names(xx))
for (i in 1:length(goi.prolif)) {
  ## expression levels of APOBEC genes
  if (goi.prolif[i] %in% colnames(dataset$ge)) {
    x <- dataset$ge[ , goi.prolif[i]]
    ## expression wrt CNV calls
    par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
    for (j in 1:length(sbtoi)) {
      myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
      myx <- intersect(myx, iix)
      tt <- table(xx[myx])
      tt <- tt[names(tt) != ""]
      if(length(tt) > 1 & all(tt > 3)) {
        wt <- kruskal.test(x[myx] ~ cnvcalls[myx])
        boxplot(x[myx] ~ xx[myx], outline=FALSE, ylab=sprintf("%s expression", names(goi.prolif)[i]), xlab="CNV calls", cex.main=0.9, main=sprintf("%s in %s (%s)\n%i patients", names(goi.prolif)[i], names(sbtoi)[j], dataset.name, length(myx)), sub=sprintf("Kruskal-Wallis test p-value = %.1E", wt$p.value))
        legend("topleft", title="CNV calls", legend=sprintf("%s (n=%i)", names(tt), tt), bty="n")
      } else { wt <- list("p.value"=NA) }
    }
  }
}
dev.off()
## combined boxplot per subtypes
pdf(file.path(goidir, sprintf("prolif3_cnvcalls_%s.pdf", dataset.name)), height=3, width=4)
for (ii in 1:length(goi.prolif)) {
	if (goi.prolif[ii] %in% colnames(dataset$ge)) {
		a.exprs <- dataset$ge[ , goi.prolif[ii]]
		x <- cnvcalls
		x[!is.na(x) & (x == "NA" | x == "")] <- NA
		rr.all2 <- mycol2 <- NULL
		for (j in 1:length(sbtoi)) {
			rr.all <- NULL
			for (k in 1:length(levels(x))) {
				iix2 <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
				iix3 <- names(x)[!is.na(x) & x == levels(x)[k]]
				myx <- fold(intersect, iix3, iix2, iix)
				## build a list of data points for each category
				if (sum(complete.cases(x[myx])) >= 3) {
					rr <- list(a.exprs[myx])
					names(rr) <- paste(levels(x)[k], names(sbtoi)[j], sep=".")
				} else {
					rr <- -1
				}
				rr.all <- c(rr.all, rr)
			}
			rr.all2 <- c(rr.all2, list(NA), rr.all)
			mycol2 <- c(mycol2, "white", rep(sbtcol[j], each=length(levels(x))))
		}
		rr.all2 <- rr.all2[-c(1)]
		mycol2 <- mycol2[-c(1)]
		myx <- sapply(rr.all2, function (x) { return (!(!is.na(x[1]) & x[1] == -1))})
		rr.all2 <- rr.all2[myx]
		mycol2 <- mycol2[myx]
		par(xaxt="n", las=1, mar=c(7, 4, 3, 1) + 0.1, cex=0.7)
		mp <- boxplot(rr.all2, outline=FALSE, ylab=sprintf("%s expression", names(goi.prolif)[ii]), xlab="", cex.main=0.9, main=sprintf("%s", names(goi.prolif)[ii]), col=mycol2, ylim=yylim)
		axis(1, at=seq(1, length(rr.all2), by=1), labels=FALSE)
		text(x=(1:length(mp$names)) + 0.25, y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=toupper(mp$names), srt=45, xpd=NA, cex=0.7, font=1)
	}
}
dev.off()

## signature scores
pdf(file.path(goidir, sprintf("signatures_cnvcalls_%s.pdf", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
xx <- as.character(cnvcalls)
xx[!is.na(xx) & xx == "NA"] <- NA
xx <- factor(xx, levels=c("", sort(unique(as.character(xx)))))
names(xx) <- names(cnvcalls)
for (ii in 1:length(bc.signatures)) {
  x <- bc.signatures[[ii]]$score
  if (sum(complete.cases(x)) > 10) {
    par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
    iix <- intersect(names(x), names(xx))
    for (j in 1:length(sbtoi)) {
      myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
      myx <- intersect(myx, iix)
      tt <- table(xx[myx])
      tt <- tt[names(tt) != ""]
      if(length(tt) > 1 & all(tt > 1)) {
        wt <- kruskal.test(x[myx] ~ xx[myx])
        boxplot(x[myx] ~ xx[myx], outline=FALSE, ylab=sprintf("%s score", names(bc.signatures)[ii]), xlab="CNV calls", cex.main=0.9, main=sprintf("%s vs %sin %s (%s)\n%i patients", names(bc.signatures)[ii], "A3B CNV call", names(sbtoi)[j], dataset.name, length(myx)), sub=sprintf("Kruskal-Wallis test p-value = %.1E", wt$p.value))
        tt <- table(cnvcalls[myx])
        legend("topleft", title="CNV calls", legend=sprintf("%s (n=%i)", names(tt), tt), bty="n")
      } else { wt <- list("p.value"=NA) }
    }
  }
}
dev.off()
## combined boxplot per subtypes
pdf(file.path(goidir, sprintf("signatures3_cnvcalls_%s.pdf", dataset.name)), height=3, width=4)
for (ii in 1:length(bc.signatures)) {
	a.exprs <- bc.signatures[[ii]]$score
	x <- cnvcalls
	x[!is.na(x) & (x == "NA" | x == "")] <- NA
	rr.all2 <- mycol2 <- NULL
	for (j in 1:length(sbtoi)) {
		rr.all <- NULL
		for (k in 1:length(levels(x))) {
			iix2 <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
			iix3 <- names(x)[!is.na(x) & x == levels(x)[k]]
			myx <- fold(intersect, iix3, iix2, iix)
			## build a list of data points for each category
			if (sum(complete.cases(x[myx])) >= 3) {
				rr <- list(a.exprs[myx])
				names(rr) <- paste(levels(x)[k], names(sbtoi)[j], sep=".")
			} else {
				rr <- -1
			}
			rr.all <- c(rr.all, rr)
		}
		rr.all2 <- c(rr.all2, list(NA), rr.all)
		mycol2 <- c(mycol2, "white", rep(sbtcol[j], each=length(levels(x))))
	}
	rr.all2 <- rr.all2[-c(1)]
	mycol2 <- mycol2[-c(1)]
	myx <- sapply(rr.all2, function (x) { return (!(!is.na(x[1]) & x[1] == -1))})
	rr.all2 <- rr.all2[myx]
	mycol2 <- mycol2[myx]
	par(xaxt="n", las=1, mar=c(7, 4, 3, 1) + 0.1, cex=0.7)
	mp <- boxplot(rr.all2, outline=FALSE, ylab=sprintf("%s expression", names(bc.signatures)[ii]), xlab="", cex.main=0.9, main=sprintf("%s", names(bc.signatures)[ii]), col=mycol2)
	axis(1, at=seq(1, length(rr.all2), by=1), labels=FALSE)
	text(x=(1:length(mp$names)) + 0.25, y=par("usr")[3] - (par("usr")[4] * 0.025), pos=2, labels=toupper(mp$names), srt=45, xpd=NA, cex=0.7, font=1)
}
dev.off()

## association with survival
pdf(file.path(goidir, sprintf("km_cnvcalls_%s.pdf", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
survd <- survcomp::censor.time(surv.time=as.numeric(dataset$clin[ , "t.rfs"]) / 365, surv.event=as.numeric(dataset$clin[ , "e.rfs"]), time.cens=0)
names(survd[[1]]) <- names(survd[[2]]) <- rownames(dataset$clin)
par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
for (j in 1:length(sbtoi)) {
   myx <- !is.na(sbt2) & is.element(sbt2, sbtoi[[j]])
   dd <- data.frame("time"=survd[[1]][myx], "event"=survd[[2]][myx], "CNV"=cnvcalls[myx])
   try(km.coxph.plot(formula.s=formula(Surv(time, event) ~ CNV), data.s=dd, x.label="Time (years)", y.label="Probability of DFS", main.title=sprintf("APOBEC3B CNV calls (%s)\n%s", dataset.name, names(sbtoi)[j]), sub.title=NULL, leg.text=levels(cnvcalls), leg.pos="bottomright", leg.inset=0.0, v.line=NULL, h.line=NULL, .col=c("darkblue", "darkgreen", "darkred"), .lty=c(1, 1), show.n.risk=TRUE, n.risk.step=3, n.risk.cex=0.85, bty="n", leg.bty="n", verbose=FALSE))
}
dev.off()

## association with node, grade, er, her2
iix <- intersect(rownames(dataset$clin), names(cnvcalls))
cn <- c("node", "grade", "er", "her2")
names(x) <- rownames(dataset$clin)
for (i in 1:length(cn)) {
  x <- dataset$clin[ , cn[i]]
  x[!is.na(x) & x == "NA"] <- NA
  x <- factor(x)
  names(x) <- rownames(dataset$clin)
  if (sum(complete.cases(x)) >= 3) {
    pdf(file.path(goidir, sprintf("%s_cnvcalls_%s.pdf", cn[i], dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
    par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
    for (j in 1:length(sbtoi)) {
      myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
      myx <- intersect(myx, iix)
      tt <- table(x[myx], cnvcalls[myx])
      names(dimnames(tt)) <- c(cn[i], "cnvcalls")
      if(length(tt) > 1 & !all(tt == 0)) {
        tts <- vcd::assocstats(tt)
        plot.new()
        plotrix::addtable2plot(x=par("usr")[1] + 0.15, y=par("usr")[4] - 0.5, xjust=0, yjust=0, table=as.matrix(tt), bty="o", display.rownames=TRUE, display.colnames=TRUE, hlines=TRUE, vlines=TRUE, title=sprintf("%s vs. %s in %s\n(%i patients, %s)", cn[i], "CNV calls", names(sbtoi)[j], length(myx), dataset.name), cex=1.25, xpad=1, ypad=1)
        title(sub=sprintf("Cramer's V=%.2g, p=%.1E", tts$cramer, tts$chisq_tests[1,3]), cex=1)
      }
    }
    dev.off()
  }
}


## association with age
pdf(file.path(goidir, sprintf("%s_cnvcalls_%s.pdf", "age", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
xx <- as.character(cnvcalls)
xx[!is.na(xx) & xx == "NA"] <- NA
xx <- factor(xx, levels=c("", sort(unique(as.character(xx)))))
names(xx) <- names(cnvcalls)
x <- as.numeric(dataset$clin[ , "age"])
names(x) <- rownames(dataset$clin)
if (sum(complete.cases(xx, x)) >= 3) {
  iix <- intersect(rownames(dataset$clin), names(cnvcalls))
  for (j in 1:length(sbtoi)) {
    myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
    myx <- intersect(myx, iix)
    tt <- table(xx)
    if(sum(complete.cases(x, xx)) >= 3) {
      wt <- kruskal.test(x[myx] ~ cnvcalls[myx])
      boxplot(x[myx] ~ xx[myx], outline=FALSE, ylab=sprintf("Age at diagnosis"), xlab="CNV calls", cex.main=0.9, main=sprintf("%s in %s (%s)\n%i patients", "Age", names(sbtoi)[j], dataset.name, length(myx)), sub=sprintf("Kruskal-Wallis test p-value = %.1E", wt$p.value))
      legend("topleft", legend=sprintf("%s (n=%i)", names(tt), tt), bty="n")
    } else { wt <- list("p.value"=NA) }
  }
  dev.off()
}
## create tables with summary stats
if (sum(complete.cases(xx, x)) >= 3) {
  iix <- intersect(rownames(dataset$clin), names(cnvcalls))
  ss.all <- NULL
  for (j in 1:length(sbtoi)) {
    myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
    myx <- intersect(myx, iix)
    tt <- table(xx)
    if(sum(complete.cases(x, xx)) >= 3) {
      # wt <- kruskal.test(x[myx] ~ cnvcalls[myx])
      # boxplot(x[myx] ~ xx[myx], outline=FALSE, ylab=sprintf("Age at diagnosis"), xlab="CNV calls", cex.main=0.9, main=sprintf("%s in %s (%s)\n%i patients", "Age", names(sbtoi)[j], dataset.name, length(myx)), sub=sprintf("Kruskal-Wallis test p-value = %.1E", wt$p.value))
      # legend("topleft", legend=sprintf("%s (n=%i)", names(tt), tt), bty="n")
      ss <- NULL
      for (ii in levels(cnvcalls)) {
        xxx <- x[!is.na(xx[myx]) & xx[myx] == ii]
        ss <- cbind(ss, c(summary(xxx), "SD"=sd(xxx, na.rm=TRUE), "MAD"=mad(xxx, na.rm=TRUE)))
      }
      colnames(ss) <- paste("APOBEC3B.CNV", levels(cnvcalls), sep=".")
      ss.all <- c(ss.all, list(data.frame(ss)))
      names(ss.all)[length(ss.all)] <- names(sbtoi)[j]
    }
  }
  if (length(ss.all) > 0) {
    WriteXLS::WriteXLS(x="ss.all", ExcelFileName=file.path(goidir, sprintf("%s_cnvcalls_%s.xls", "age", dataset.name)), row.names=TRUE, BoldHeaderRow=TRUE)
  }
}

## association with tumor size
xx <- as.character(cnvcalls)
xx[!is.na(xx) & xx == "NA"] <- NA
xx <- factor(xx, levels=c("", sort(unique(as.character(xx)))))
names(xx) <- names(cnvcalls)
x <- as.numeric(dataset$clin[ , "size"])
names(x) <- rownames(dataset$clin)
if (sum(complete.cases(xx, x)) >= 3) {
  pdf(file.path(goidir, sprintf("%s_cnvcalls_%s.pdf", "size", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
  par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
  iix <- intersect(rownames(dataset$clin), names(cnvcalls))
  for (j in 1:length(sbtoi)) {
    myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
    myx <- intersect(myx, iix)
    tt <- table(xx)
    if(sum(complete.cases(x, xx)) >= 3) {
      wt <- kruskal.test(x[myx] ~ cnvcalls[myx])
      boxplot(x[myx] ~ xx[myx], outline=FALSE, ylab=sprintf("Tumor size"), xlab="CNV calls", cex.main=0.9, main=sprintf("%s in %s (%s)\n%i patients", "Tumor size", names(sbtoi)[j], dataset.name, length(myx)), sub=sprintf("Kruskal-Wallis test p-value = %.1E", wt$p.value))
      legend("topleft", title="CNV calls", legend=sprintf("%s (n=%i)", names(tt), tt), bty="n")
    } else { wt <- list("p.value"=NA) }
  }
  dev.off()
}
## create tables with summary stats
if (sum(complete.cases(xx)) >= 3) {
  iix <- intersect(rownames(dataset$clin), names(cnvcalls))
  ss.all <- NULL
  for (j in 1:length(sbtoi)) {
    myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
    myx <- intersect(myx, iix)
    tt <- table(xx)
    if(sum(complete.cases(x, xx)) >= 3) {
      ss <- NULL
      for (ii in levels(cnvcalls)) {
        xxx <- x[!is.na(xx[myx]) & xx[myx] == ii]
        ss <- cbind(ss, c(summary(xxx), "SD"=sd(xxx, na.rm=TRUE), "MAD"=mad(xxx, na.rm=TRUE)))
      }
      colnames(ss) <- paste("APOBEC3B.CNV", levels(cnvcalls), sep=".")
      ss.all <- c(ss.all, list(data.frame(ss)))
      names(ss.all)[length(ss.all)] <- names(sbtoi)[j]
    }
  }
  if (length(ss.all) > 0) {
    WriteXLS::WriteXLS(x="ss.all", ExcelFileName=file.path(goidir, sprintf("%s_cnvcalls_%s.xls", "size", dataset.name)), row.names=TRUE, BoldHeaderRow=TRUE)
  }
}

## association with mutation status
pdf(file.path(goidir, sprintf("mutation_cnvcalls_%s.pdf", dataset.name)), height=15, width=ceiling(length(sbtoi) / 3) * 5)
par(mfrow=c(3, ceiling(length(sbtoi) / 3)))
for (ii in 1:ncol(mutations)) {
  x <- mutations[ , ii]
  names(x) <- rownames(mutations)
  if (sum(complete.cases(x)) >= 3) {
    for (j in 1:length(sbtoi)) {
      myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
      myx <- intersect(myx, iix)
      tt <- table(x[myx], cnvcalls[myx])
      names(dimnames(tt)) <- c(colnames(mutations)[ii], "cnvcalls")
      if(length(tt) > 1 & !all(tt == 0)) {
        tts <- vcd::assocstats(tt)
        plot.new()
        plotrix::addtable2plot(x=par("usr")[1] + 0.15, y=par("usr")[4] - 0.5, xjust=0, yjust=0, table=as.matrix(tt), bty="o", display.rownames=TRUE, display.colnames=TRUE, hlines=TRUE, vlines=TRUE, title=sprintf("%s vs. %s in %s\n(%i patients, %s)", colnames(mutations)[ii], "CNV calls", names(sbtoi)[j], length(myx), dataset.name), cex=1.25, xpad=1, ypad=1)
        title(sub=sprintf("Cramer's V=%.2g, p=%.1E", tts$cramer, tts$chisq_tests[1,3]), cex=1)
      }
    }
  }
}
dev.off()

#################################################
## subtype-specific GSEA for association with APOBEC3B CNV calls
#################################################

## subtype-specific GSEA
gsea.out <- file.path(goidir, "GSEA")
if(!file.exists(gsea.out)) { dir.create(gsea.out, showWarnings=FALSE, recursive=TRUE) }
gsea.out.rank <- file.path(goidir, "GSEA", "rankings")
if(!file.exists(gsea.out.rank)) { dir.create(gsea.out.rank, showWarnings=FALSE, recursive=TRUE) }
gsea.out.res <- file.path(goidir, "GSEA", "reports")
if(!file.exists(gsea.out.res)) { dir.create(gsea.out.res, showWarnings=FALSE, recursive=TRUE) }

## create gmt, file containing the geneset definitions
myfn <- file.path(gsea.out, "genesets.RData")
genesets.filen2 <- unlist(strsplit(x=basename(genesets.filen), split="[.]"))
genesets.filen2 <- paste(c(genesets.filen2[1:(length(genesets.filen2) - 1)], "format", genesets.filen2[length(genesets.filen2)]), collapse=".")
genesets.filen2 <- file.path(gsea.out, genesets.filen2)
if(!file.exists(myfn)) {
  genesets <- MetaGx::formatGMT(infile=genesets.filen, outfile=genesets.filen2, replace=TRUE)
  save(list="genesets", compress=TRUE, file=myfn)
} else { load(myfn) }


## create rankings
myfn <- file.path(gsea.out, "rankings.RData")
if (!file.exists(myfn)) {
  iix <- intersect(rownames(dataset$ge), names(cnvcalls))
  rank.files <- rank.files.0vs2 <- rank.files.1vs2 <- rank.files.01vs2 <- NULL
  for (j in 1:length(sbtoi)) {
    myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
    myx <- intersect(myx, iix)
    yy <- cnvcalls
	 levels(yy) <- c(0, 1, 2)
    ## parallel execution
    splitix <- parallel::splitIndices(nx=ncol(dataset$ge), ncl=nbcore)
    splitix <- splitix[sapply(splitix, length) > 0]
    mcres <- parallel::mclapply(splitix, function(x, data, y) {
      res <- apply(data[ , x, drop=FALSE], 2, function (x, y) {
        ccix <- complete.cases(x, y)
        if(sum(ccix) >= 10) {
          ## multinomial model (for outcome as unordered factor)
          tempff <- file.path(gsea.out.rank, basename(tempfile(pattern="res_", tmpdir="", fileext=".tmp")))
          sink(file=tempff, type="output")
          ## reference is wild type
          y <- relevel(y, ref="2")
          res1 <- nnet::multinom(y ~ x)
          rr <- summary(res1)
          if ("x" %in% colnames(rr$coefficients)) {
            res <- - (rr$coefficients / rr$standard.errors)[ , "x"]
          } else { res <- c(NA, NA) }
          sink()
          unlink(x=tempff, force=TRUE)
        } else { res <- c(NA, NA) }
        return(res)
      }, y=y)
      return(res)
    }, data=dataset$ge[myx, , drop=FALSE], y=yy[myx])
    sst1 <- t(do.call(cbind, mcres))
    
    ## 0/1 vs 2
    levels(yy) <- c(0, 0, 1)
    mcres <- parallel::mclapply(splitix, function(x, data, y) {
      res <- apply(data[ , x, drop=FALSE], 2, function (x, y) {
        ccix <- complete.cases(x, y)
        if(sum(ccix) >= 10) {
          ## logistic regression 0/1 vs 2
          res1 <- glm(y ~ x, family = "binomial")
          rr <- summary(res1)
          if ("x" %in% rownames(rr$coefficients)) {
            res <- rr$coefficients["x", "z value"]
          } else { res <- NA }
        } else { res <- NA }
        return(res)
      }, y=y)
      return(res)
    }, data=dataset$ge[myx, , drop=FALSE], y=yy[myx])
    sst2 <- do.call(c, mcres)
    
    ss.all <- cbind(sst1, sst2)
    colnames(ss.all) <- c("CNV.0.vs.2", "CNV.1.vs.2", "CNV.01.vs.2")
    ss.all <- ss.all[colnames(dataset$ge), , drop=FALSE]
    
    ## 0vs2
    ss <- ss.all[ , "CNV.0.vs.2"]
    ssn <- "0vs2"
    ss <- sort(ss, decreasing=TRUE, na.last=NA)
    ss[ss == Inf] <- .Machine$double.xmax
    ss[ss == -Inf] <- -.Machine$double.xmax
    rankg <- data.frame("entrezid"=dataset$annot[names(ss), "EntrezGene.ID"], "statistic"=ss)
    fff <- file.path(gsea.out.rank, sprintf("%s_%s.rnk", names(sbtoi)[j], ssn))
    write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    rank.files.0vs2 <- c(rank.files.0vs2, fff)
    
    ## 1vs2
    ss <- ss.all[ , "CNV.1.vs.2"]
    ssn <- "1vs2"
    ss <- sort(ss, decreasing=TRUE, na.last=NA)
    ss[ss == Inf] <- .Machine$double.xmax
    ss[ss == -Inf] <- -.Machine$double.xmax
    rankg <- data.frame("entrezid"=dataset$annot[names(ss), "EntrezGene.ID"], "statistic"=ss)
    fff <- file.path(gsea.out.rank, sprintf("%s_%s.rnk", names(sbtoi)[j], ssn))
    write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    rank.files.1vs2 <- c(rank.files.1vs2, fff)
    
    ## 01vs2
    ss <- ss.all[ , "CNV.01.vs.2"]
    ssn <- "01vs2"
    ss <- sort(ss, decreasing=TRUE, na.last=NA)
    ss[ss == Inf] <- .Machine$double.xmax
    ss[ss == -Inf] <- -.Machine$double.xmax
    rankg <- data.frame("entrezid"=dataset$annot[names(ss), "EntrezGene.ID"], "statistic"=ss)
    fff <- file.path(gsea.out.rank, sprintf("%s_%s.rnk", names(sbtoi)[j], ssn))
    write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
    rank.files.01vs2 <- c(rank.files.01vs2, fff)
    
  }
  names(rank.files.0vs2) <- names(rank.files.1vs2) <- names(rank.files.01vs2) <- names(sbtoi)
  rank.files <- list("0vs2"=rank.files.0vs2, "1vs2"=rank.files.1vs2, "01vs2"=rank.files.01vs2)
  save(list=c("rank.files"), compress=TRUE, file=myfn)
} else { load(myfn) }


## prerank GSEA
for (i in 1:length(rank.files)) {
  myfn <- file.path(gsea.out, sprintf("gseas_%s.RData", names(rank.files)[i]))
  if (!file.exists(myfn)) {
    gsea.out.res2 <- file.path(gsea.out.res, names(rank.files)[i])
    if(!file.exists(gsea.out.res2)) { dir.create(gsea.out.res2, showWarnings=FALSE, recursive=TRUE) }
    splitix <- as.list(1:length(rank.files[[i]]))
    splitix <- splitix[sapply(splitix, length) > 0]
    gsea.res <- parallel::mclapply(splitix, function(x, rank.files, genesets.filen, gsea.nperm, gsea.out, max.geneset.size, min.geneset.size) {
      res <- MetaGx::prerankGSEA(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=rank.files[x], nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.report=names(rank.files)[x], gsea.out=gsea.out, replace.res=FALSE, gsea.seed=987654321)
    }, rank.files=rank.files[[i]], genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, gsea.out=gsea.out.res2, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size)
    save(list="gsea.res", compress=TRUE, file=myfn)
  }
}

#################################################
## association with hypermutation
#################################################
if (!all(is.na(mutations[ "HYPERMUTATION"]))) {
  goidir <- file.path(saveres, dataset.name, "hypermutation")
  if (!file.exists(goidir)) { dir.create(goidir, showWarnings=FALSE, recursive=TRUE) }
}

if (!all(is.na(mutations[ "HYPERMUTATION"]))) {
	
	## subtype-specific GSEA
	gsea.out <- file.path(goidir, "GSEA")
	if(!file.exists(gsea.out)) { dir.create(gsea.out, showWarnings=FALSE, recursive=TRUE) }
	gsea.out.rank <- file.path(goidir, "GSEA", "rankings")
	if(!file.exists(gsea.out.rank)) { dir.create(gsea.out.rank, showWarnings=FALSE, recursive=TRUE) }
	gsea.out.res <- file.path(goidir, "GSEA", "reports")
	if(!file.exists(gsea.out.res)) { dir.create(gsea.out.res, showWarnings=FALSE, recursive=TRUE) }

	## create gmt, file containing the geneset definitions
	myfn <- file.path(gsea.out, "genesets.RData")
	genesets.filen2 <- unlist(strsplit(x=basename(genesets.filen), split="[.]"))
	genesets.filen2 <- paste(c(genesets.filen2[1:(length(genesets.filen2) - 1)], "format", genesets.filen2[length(genesets.filen2)]), collapse=".")
	genesets.filen2 <- file.path(gsea.out, genesets.filen2)
	if(!file.exists(myfn)) {
	  genesets <- MetaGx::formatGMT(infile=genesets.filen, outfile=genesets.filen2, replace=TRUE)
	  save(list="genesets", compress=TRUE, file=myfn)
	} else { load(myfn) }

	## create rankings
	myfn <- file.path(gsea.out, "rankings.RData")
	if (!file.exists(myfn)) {
	  iix <- intersect(rownames(dataset$ge), names(cnvcalls))
	  rank.files <- NULL
	  for (j in 1:length(sbtoi)) {
	    myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
	    myx <- intersect(myx, iix)
		 yy <- mutations[ , "HYPERMUTATION"]
		 yy <- as.numeric(yy)
		 names(yy) <- rownames(mutations)
	    ## parallel execution
	    splitix <- parallel::splitIndices(nx=ncol(dataset$ge), ncl=nbcore)
	    splitix <- splitix[sapply(splitix, length) > 0]
	    mcres <- parallel::mclapply(splitix, function(x, data, y) {
	      res <- apply(data[ , x, drop=FALSE], 2, function (x, y) {
	        ccix <- complete.cases(x, y)
	        if(sum(ccix) >= 10) {
	          ## logistic regression
	          res1 <- glm(y ~ x, family = "binomial")
	          rr <- summary(res1)
	          if ("x" %in% rownames(rr$coefficients)) {
	            res <- rr$coefficients["x", "z value"]
	          } else { res <- NA }
	        } else { res <- NA }
	        return(res)
	      }, y=y)
	      return(res)
	    }, data=dataset$ge[myx, , drop=FALSE], y=yy[myx])
		 ss <- do.call(c, mcres)
		 ss <- ss[colnames(dataset$ge)]
       ss <- sort(ss, decreasing=TRUE, na.last=NA)
       ss[ss == Inf] <- .Machine$double.xmax
       ss[ss == -Inf] <- -.Machine$double.xmax
       rankg <- data.frame("entrezid"=dataset$annot[names(ss), "EntrezGene.ID"], "statistic"=ss)
       fff <- file.path(gsea.out.rank, sprintf("%s.rnk", names(sbtoi)[j]))
       write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
       rank.files <- c(rank.files, fff)
	  }
	  names(rank.files) <- names(sbtoi)
	  save(list=c("rank.files"), compress=TRUE, file=myfn)
	} else { load(myfn) }

	## prerank GSEA
	myfn <- file.path(gsea.out, "gseas.RData")
	if (!file.exists(myfn)) {
		splitix <- as.list(1:length(rank.files))
		splitix <- splitix[sapply(splitix, length) > 0]
		gsea.res <- parallel::mclapply(splitix, function(x, rank.files, genesets.filen, gsea.nperm, gsea.out, max.geneset.size, min.geneset.size) {
		  res <- MetaGx::prerankGSEA(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=rank.files[x], nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.report=names(rank.files)[x], gsea.out=gsea.out, replace.res=FALSE, gsea.seed=987654321)
		}, rank.files=rank.files, genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, gsea.out=gsea.out.res, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size)
		save(list="gsea.res", compress=TRUE, file=myfn)
	}
}


#################################################
## association with APOBEC3B calls controlling for hypermutation
#################################################
if (!all(is.na(mutations[ "HYPERMUTATION"]))) {
  goidir <- file.path(saveres, dataset.name, "a3bcnv_control_hypermutation")
  if (!file.exists(goidir)) { dir.create(goidir, showWarnings=FALSE, recursive=TRUE) }
}

if (!all(is.na(mutations[ "HYPERMUTATION"]))) {
  ## subtype-specific GSEA
  gsea.out <- file.path(goidir, "GSEA")
  if(!file.exists(gsea.out)) { dir.create(gsea.out, showWarnings=FALSE, recursive=TRUE) }
  gsea.out.rank <- file.path(goidir, "GSEA", "rankings")
  if(!file.exists(gsea.out.rank)) { dir.create(gsea.out.rank, showWarnings=FALSE, recursive=TRUE) }
  gsea.out.res <- file.path(goidir, "GSEA", "reports")
  if(!file.exists(gsea.out.res)) { dir.create(gsea.out.res, showWarnings=FALSE, recursive=TRUE) }

  ## create gmt, file containing the geneset definitions
  myfn <- file.path(gsea.out, "genesets.RData")
  genesets.filen2 <- unlist(strsplit(x=basename(genesets.filen), split="[.]"))
  genesets.filen2 <- paste(c(genesets.filen2[1:(length(genesets.filen2) - 1)], "format", genesets.filen2[length(genesets.filen2)]), collapse=".")
  genesets.filen2 <- file.path(gsea.out, genesets.filen2)
  if(!file.exists(myfn)) {
    genesets <- MetaGx::formatGMT(infile=genesets.filen, outfile=genesets.filen2, replace=TRUE)
    save(list="genesets", compress=TRUE, file=myfn)
  } else { load(myfn) }

  ## create rankings
	myfn <- file.path(gsea.out, "rankings.RData")
	if (!file.exists(myfn)) {
		iix <- intersect(rownames(dataset$ge), names(cnvcalls))
		rank.files <- rank.files.0vs2 <- rank.files.1vs2 <- rank.files.01vs2 <- NULL
		for (j in 1:length(sbtoi)) {
		  myx <- names(sbt2)[!is.na(sbt2) & is.element(sbt2, sbtoi[[j]])]
		  myx <- intersect(myx, iix)
		  yy <- cnvcalls
		  levels(yy) <- c(0, 1, 2)
		  zz <- mutations[ , "HYPERMUTATION"]
		  zz <- as.numeric(zz)
		  names(zz) <- rownames(mutations)
		  ## parallel execution
		  splitix <- parallel::splitIndices(nx=ncol(dataset$ge), ncl=nbcore)
		  splitix <- splitix[sapply(splitix, length) > 0]
		  mcres <- parallel::mclapply(splitix, function(x, data, y, z) {
		    res <- apply(data[ , x, drop=FALSE], 2, function (x, y, z) {
		      ccix <- complete.cases(x, y)
		      if(sum(ccix) >= 10) {
		        ## multinomial model (for outcome as unordered factor)
		        tempff <- file.path(gsea.out.rank, basename(tempfile(pattern="res_", tmpdir="", fileext=".tmp")))
		        sink(file=tempff, type="output")
		        ## reference is wild type
		        y <- relevel(y, ref="2")
		        res1 <- nnet::multinom(y ~ x + z)
		        rr <- summary(res1)
		        if ("x" %in% colnames(rr$coefficients)) {
		          res <- - (rr$coefficients / rr$standard.errors)[ , "x"]
		        } else { res <- c(NA, NA) }
		        sink()
		        unlink(x=tempff, force=TRUE)
		      } else { res <- c(NA, NA) }
		      return(res)
		    }, y=y, z=z)
		    return(res)
		  }, data=dataset$ge[myx, , drop=FALSE], y=yy[myx], z=zz[myx])
		  sst1 <- t(do.call(cbind, mcres))

		  ## 0/1 vs 2
		  levels(yy) <- c(0, 0, 1)
		  mcres <- parallel::mclapply(splitix, function(x, data, y, z) {
		    res <- apply(data[ , x, drop=FALSE], 2, function (x, y, z) {
		      ccix <- complete.cases(x, y)
		      if(sum(ccix) >= 10) {
		        ## logistic regression 0/1 vs 2
		        res1 <- glm(y ~ x + z, family = "binomial")
		        rr <- summary(res1)
		        if ("x" %in% rownames(rr$coefficients)) {
		          res <- rr$coefficients["x", "z value"]
		        } else { res <- NA }
		      } else { res <- NA }
		      return(res)
		    }, y=y, z=z)
		    return(res)
		  }, data=dataset$ge[myx, , drop=FALSE], y=yy[myx], z=zz[myx])
		  sst2 <- do.call(c, mcres)

		  ss.all <- cbind(sst1, sst2)
		  colnames(ss.all) <- c("CNV.0.vs.2", "CNV.1.vs.2", "CNV.01.vs.2")
		  ss.all <- ss.all[colnames(dataset$ge), , drop=FALSE]

		  ## 0vs2
		  ss <- ss.all[ , "CNV.0.vs.2"]
		  ssn <- "0vs2"
		  ss <- sort(ss, decreasing=TRUE, na.last=NA)
		  ss[ss == Inf] <- .Machine$double.xmax
		  ss[ss == -Inf] <- -.Machine$double.xmax
		  rankg <- data.frame("entrezid"=dataset$annot[names(ss), "EntrezGene.ID"], "statistic"=ss)
		  fff <- file.path(gsea.out.rank, sprintf("%s_%s.rnk", names(sbtoi)[j], ssn))
		  write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
		  rank.files.0vs2 <- c(rank.files.0vs2, fff)

		  ## 1vs2
		  ss <- ss.all[ , "CNV.1.vs.2"]
		  ssn <- "1vs2"
		  ss <- sort(ss, decreasing=TRUE, na.last=NA)
		  ss[ss == Inf] <- .Machine$double.xmax
		  ss[ss == -Inf] <- -.Machine$double.xmax
		  rankg <- data.frame("entrezid"=dataset$annot[names(ss), "EntrezGene.ID"], "statistic"=ss)
		  fff <- file.path(gsea.out.rank, sprintf("%s_%s.rnk", names(sbtoi)[j], ssn))
		  write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
		  rank.files.1vs2 <- c(rank.files.1vs2, fff)

		  ## 01vs2
		  ss <- ss.all[ , "CNV.01.vs.2"]
		  ssn <- "01vs2"
		  ss <- sort(ss, decreasing=TRUE, na.last=NA)
		  ss[ss == Inf] <- .Machine$double.xmax
		  ss[ss == -Inf] <- -.Machine$double.xmax
		  rankg <- data.frame("entrezid"=dataset$annot[names(ss), "EntrezGene.ID"], "statistic"=ss)
		  fff <- file.path(gsea.out.rank, sprintf("%s_%s.rnk", names(sbtoi)[j], ssn))
		  write.table(rankg, file=fff, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
		  rank.files.01vs2 <- c(rank.files.01vs2, fff)
	  }
	  names(rank.files.0vs2) <- names(rank.files.1vs2) <- names(rank.files.01vs2) <- names(sbtoi)
	  rank.files <- list("0vs2"=rank.files.0vs2, "1vs2"=rank.files.1vs2, "01vs2"=rank.files.01vs2)
	  save(list=c("rank.files"), compress=TRUE, file=myfn)
	} else { load(myfn) }

  ## prerank GSEA
  for (i in 1:length(rank.files)) {
    myfn <- file.path(gsea.out, sprintf("gseas_%s.RData", names(rank.files)[i]))
    if (!file.exists(myfn)) {
      gsea.out.res2 <- file.path(gsea.out.res, names(rank.files)[i])
      if(!file.exists(gsea.out.res2)) { dir.create(gsea.out.res2, showWarnings=FALSE, recursive=TRUE) }
      splitix <- as.list(1:length(rank.files[[i]]))
      splitix <- splitix[sapply(splitix, length) > 0]
      gsea.res <- parallel::mclapply(splitix, function(x, rank.files, genesets.filen, gsea.nperm, gsea.out, max.geneset.size, min.geneset.size) {
        res <- MetaGx::prerankGSEA(exe.path=gsea.exec, gmt.path=genesets.filen, rank.path=rank.files[x], nperm=gsea.nperm, scoring.scheme="weighted", make.sets=TRUE, include.only.symbols=FALSE, plot.top.x=20, set.max=max.geneset.size, set.min=min.geneset.size, zip.report=FALSE, gsea.report=names(rank.files)[x], gsea.out=gsea.out, replace.res=FALSE, gsea.seed=987654321)
      }, rank.files=rank.files[[i]], genesets.filen=genesets.filen2, gsea.nperm=gsea.nperm, gsea.out=gsea.out.res2, max.geneset.size=max.geneset.size, min.geneset.size=min.geneset.size)
      save(list="gsea.res", compress=TRUE, file=myfn)
    }
  }
}

## additional analyses
library(vcd)
table(sbt2, cnvcalls)
vcd::assocstats(table(sbt2, cnvcalls))



## end

