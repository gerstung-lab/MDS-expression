#+ Preliminaries, echo=FALSE
options(width=120)
pdf.options(pointsize=8)
knit_hooks$set(smallMar = function(before, options, envir) {
			if (before) par(mar=c(3,3,1,1), bty="n", mgp=c(2,.5,0)) 
		})
opts_chunk$set(dev=c('my_png','pdf'), fig.ext=c('png','pdf'), fig.width=3, fig.height=3, smallMar=TRUE)
my_png <-  function(file, width, height, pointsize=12, ...) {
	png(file, width = 1.5*width, height = 1.5*height, units="in", res=72*1.5, pointsize=pointsize, ...)
}
#opts_knit$set(root.dir = file.path(getwd(),".."))

#' TCGA AML expression data analysis
#' ===========================

#' #### Libraries
library(limma)
library(org.Hs.eg.db)
library(hgu133plus2.db )
library(RColorBrewer)
set1 <- brewer.pal(9,"Set1")
library(cgdsr)
library(CoxHD)
library(mg14)
library(xtable)
library(Hmisc)

#' 1. Differential expression analysis
#' ----------------------
#' 
#' ### 1.1 Load data 
#' #### Get entrez ids for genes of interest
entrez <- mappedRkeys(hgu133plus2ENTREZID)

#' #### Load expression data from cBio portal
#+ cBioData, cache=TRUE
mycgds <-  CGDS("http://www.cbioportal.org/public-portal/")
tcgaAML <- getCancerStudies(mycgds)[1,1]
cases <- getCaseLists(mycgds,tcgaAML)[8,1]
g <-  lapply(split(as.numeric(entrez), seq_along(entrez)%/%500), function(genes) getProfileData(mycgds,genes,getGeneticProfiles(mycgds,tcgaAML)[2,1],cases)) ## load in batches of 500
g <- do.call("cbind", g)
tcgaExpr <- log(t(g)+.5)
rownames(tcgaExpr) <- sub("\\.","-",sub("^[0-9]+\\.","",rownames(tcgaExpr))) ## Fix rownames
#rm(e,s,g)

#' The query will have dropped a few genes. Match back to entrezids...
s <- AnnotationDbi::select(org.Hs.eg.db, rownames(tcgaExpr), "ENTREZID", "SYMBOL")
m <- match(rownames(tcgaExpr), s$SYMBOL)
sum(is.na(m))
rownames(tcgaExpr) <- s$ENTREZID[m]
colnames(tcgaExpr) <- gsub("\\.","-",colnames(tcgaExpr))
tcgaExpr <- tcgaExpr[rowSums(is.na(tcgaExpr))==0 & !is.na(rownames(tcgaExpr)),]

#' Check the log2 FPKM densities across samples
#+ TCGA-logFPKM, fig.width=2, fig.height=1.8
plot(density(tcgaExpr[,1]))
for(j in 1:ncol(tcgaExpr)) lines(density(tcgaExpr[,j]))

#' #### Load curated clinical and mutation data
library(xlsx)
tmp <- tempfile()
download.file("http://www.nature.com/ncomms/2015/150109/ncomms6901/extref/ncomms6901-s6.xlsx", tmp)
numericalize <- function(x){if(class(x)!='factor') return(x) else if(all(is.na(as.numeric(levels(x))))) return(x) else return(as.numeric(as.character(x)))}
tcgaClinical <- as.data.frame(sapply(read.xlsx(tmp, sheetName="Clinical", startRow=3), numericalize, simplify=FALSE))
colnames(tcgaClinical)[3:12] <- capitalize(colnames(tcgaClinical)[3:12])
rownames(tcgaClinical) <- as.character(tcgaClinical$TCGA_ID), numericalize, simplify=FALSE))
colnames(tcgaClinical)[3:12] <- capitalize(colnames(tcgaClinical)[3:12])
tcgaGenetic <- as.data.frame(sapply(read.xlsx(tmp, sheetName="Genetic", startRow=3)
tcgaGenetic$TCGA_ID <- factor(as.character(tcgaGenetic$TCGA_ID), levels = levels(tcgaClinical$TCGA_ID))
g <- as.character(tcgaGenetic$Hugo_Symbol)
g[tcgaGenetic$Hugo_Symbol=="FLT3" & tcgaGenetic$Variant_Type == 'INS'] <- "FLT3_ITD"
g[tcgaGenetic$Hugo_Symbol=="FLT3" & tcgaGenetic$Variant_Type == 'SNP'] <- "FLT3_TKD"
tcgaMutation <- (table(tcgaGenetic$TCGA_ID,g)) > 0

#' ### 1.2 Compute PCA
tcgaPca <- prcomp(t(tcgaExpr))

#' #### Covariates
t <- cbind(tcgaMutation +0, tcgaClinical[,14:24])
t <- t[,colSums(t[colnames(tcgaExpr),],na.rm=TRUE)>=5]
tcgaCovariates <- as.matrix(cbind(Offset=1,t, Gender=tcgaClinical$Gender, Age=tcgaClinical$AOD/10))[colnames(tcgaExpr),]
groups <- factor(c("Offset", rep("Genetics",22), "Translocations",rep("CNA", 5), rep("Translocations",2), rep("Demographics",2)), levels=c("Offset","Genetics","CNA","Translocations","Demographics"))
tcgaCovariates <- tcgaCovariates[,order(groups)]
groups <- groups[order(groups)]
col1 <- c("grey",set1[c(3,5,2,7)])
names(col1) <- levels(groups)


#' #### PCA overview
#' Now plot a large panel with small multiples of the first two PCs overlaid with the mutation status of each gene.
#+ TCGA-PCA-panel,  fig.width=5, fig.height=4
par(bty="n", mgp = c(0,0.5,0), mar=c(1,1,1.5,0)+.1, las=1, tcl=-.25, font.main=3, mfrow=c(6,6), xpd=NA)
i<-0
for (geneId in colnames(tcgaCovariates)[-1]){
	i<-i+1
	plot(tcgaPca$x[rownames(tcgaCovariates),], cex=0.5, 
			pch=NA, 
			xlab=ifelse(i==1,"PC1",""), ylab=ifelse(i==1,"PC2",""),main=geneId, font.main=ifelse(grepl("[[:lower:]]",geneId),1,3), cex.main=1.33, cex.lab=1.2,xaxt="n", yaxt="n", ylim=c(-65,65))
	if(geneId != "Age")
		w <- rownames(tcgaCovariates)[which(tcgaCovariates[,geneId] == 1)]
	else
		w <- rownames(tcgaCovariates)[which(tcgaCovariates[,geneId] > median(tcgaCovariates[,geneId], na.rm=TRUE))]
	points(tcgaPca$x[!rownames(tcgaPca$x) %in% w,], pch=ifelse(is.na(tcgaCovariates[!rownames(tcgaPca$x) %in% w,geneId]),1,19), cex=0.5, col="grey", lwd=0.5)
	points(tcgaPca$x[w,], pch=16, cex=0.85, col=col1[groups[i+1]], lwd=0.05)
	u <- matrix(par("usr"), ncol=2)
	if(i==1){
		arrows(u[1,1],u[1,2], u[2,1],u[1,2],length=0.02)
		arrows(u[1,1],u[1,2], u[1,1],u[2,2],length=0.02)
	}
	text(u[2],u[3] + (u[4]-u[3])*.9, labels=paste("n=",length(w), sep=""), bty="n", cex=1.2, pos=2)
}
plot.new()
legend("center",c("Missing","Wildtype","Mutant"), pch=c(1,19,19), col=c("grey","grey","black"), bty="n", pt.cex=c(0.5,0.5,rep(0.85,4)),cex=1.2, pt.lwd=0.5)


#' ### 1.2 Explained variance R2/F-statistic
#' Here we want to determine all genes which are associated with any covariate. This will be based on an F-statistic.
#' lmFit also tests wether the offset is different from zero (trivially true).
#' Create design matrix of covariates; impute missing values by mean. 
poorMansImpute <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}
tcgaDesign <- apply(tcgaCovariates,2,poorMansImpute)

#' #### Fit the linear model
glm = lmFit(as.matrix(tcgaExpr), design = tcgaDesign) 
glm = eBayes(glm)

#' #### Random model
#' Compare to a model where all values of the covariates are randomly permuted. This helps reassure FDR estimates.
set.seed(42)
rlm <- lmFit(tcgaExpr[,rownames(tcgaDesign)], apply(tcgaDesign, 2, sample))
rlm <- eBayes(rlm)
F.stat <- classifyTestsF(rlm[,-1],fstat.only=TRUE)
rlm$F <- as.vector(F.stat)
df1 <- attr(F.stat,"df1")
df2 <- attr(F.stat,"df2")
if(df2[1] > 1e6){ # Work around bug in R 2.1
	rlm$F.p.value <- pchisq(df1*rlm$F,df1,lower.tail=FALSE)
}else
	rlm$F.p.value <- pf(rlm$F,df1,df2,lower.tail=FALSE)

#' #### Explained variance by different categories
#' The F-statistic is directly related to the R2.
F.stat <- classifyTestsF(glm[,2:22],fstat.only=TRUE) ## All genetics & cytogenetics
df1 <- attr(F.stat,"df1")
df2 <- attr(F.stat,"df2")
F.p.value <- pchisq(df1*F.stat,df1,lower.tail=FALSE)

R.stat <- classifyTestsF(rlm[,2:22],fstat.only=TRUE) ## Random

Rall = 1 - 1/(1 + glm$F * (ncol(tcgaDesign)-1)/(nrow(tcgaDesign)-ncol(tcgaDesign)))
Rgenetics = 1 - 1/(1 + F.stat * 21/(nrow(tcgaDesign)-ncol(tcgaDesign)))
Pgenetics = 1 - 1/(1 + R.stat * 21/(nrow(tcgaDesign)-ncol(tcgaDesign)))
names(Rgenetics) <- names(Pgenetics) <- names(Rall) <-  rownames(tcgaExpr)

#' #### Plot the variance explained by genetics
#+ TCGA-GE-R2, fig.width=2, fig.height=1.8
par(bty="n", mgp = c(2,.33,0), mar=c(3,2.5,1,1)+.1, las=1, tcl=-.25, xpd=NA)
d <- density(Pgenetics,bw=1e-3)
f <- 1#nrow(gexpr)/512
plot(d$x, d$y * f, col='grey', xlab=expression(paste("Explained variance per gene ", R^2)), main="", lwd=2, type="l", ylab="", xlim=c(0,0.7))
title(ylab="Density", line=1.5)
d <- density(Rgenetics, bw=1e-3)
r <- min(Rgenetics[p.adjust(F.p.value,"BH")<0.05])
x0 <- which(d$x>r)
polygon(d$x[c(x0[1],x0)], c(0,d$y[x0])* f, col=paste(set1[1],"44",sep=""), border=NA)
lines(d$x, d$y* f, col=set1[1], lwd=2)
#points(d$x[x0[1]], d$y[x0[1]]*f, col=set1[1], pch=16)
text(d$x[x0[1]], d$y[x0[1]]*f, pos=4, paste(sum(Rgenetics > r), "genes q < 0.05"))
legend("topright", bty="n", col=c(set1[1], "grey"), lty=1, c("Observed","Random"), lwd=2)

#' #### Predictions
glmPrediction <- glm$coefficients %*% t(tcgaDesign)
rlmPrediction <- rlm$coefficients %*% t(tcgaDesign)
glmProjection <-  t(glm$coefficients[,-1] ) %*% tcgaPca$rotation[,1:20]

#' ### 1.4. Test results
#' Prepare the test results using a hierarchical procedure: 
#' 1. Adjust down transcripts using F-stat
#' 2. Adjust along covariates
#+ TCGA-testResults
testResults <- decideTests(glm, method="hierarchical",adjust.method="BH", p.value=0.05)[,-1]
significantGenes <- sapply(1:ncol(testResults), function(j){
			c <- glm$coefficients[testResults[,j]!=0,j+1]
			table(cut(c, breaks=c(-5,seq(-1.5,1.5,l=7),5)))
		})
colnames(significantGenes) <- colnames(testResults)


#' #### Table of top 50 deregulated genes contains many HOX genes
#+ TCGA-top50, results='asis'
t <- head(sort(Rgenetics, d=TRUE), 50)
g <- apply(testResults[names(t),], 1, function(x) paste(colnames(testResults)[x!=0], collapse=", "))
print(xtable(cbind(AnnotationDbi::select(org.Hs.eg.db, names(t), c("SYMBOL","GENENAME"), "ENTREZID"),`R^2`=t, `Mutations`=g)), type="html", )

#' #### Top genes
#+ TCGA-GE-predict, fig.width=2, fig.height=1.8
par(bty="n", mgp = c(1.5,.33,0), mar=c(2.5,2.5,1,1)+.1, las=1, tcl=-.25)
for(w in names(head(sort(Rgenetics, decreasing=TRUE),4))){
gene <- AnnotationDbi::select(org.Hs.eg.db, w, "SYMBOL", "ENTREZID")$SYMBOL
plot(glmPrediction[w,], tcgaExpr[w,rownames(tcgaDesign)], ylab=parse(text=paste("Observed  ~ italic(",gene,") ~ expression")), xlab=parse(text=paste("Predicted  ~ italic(",gene,") ~ expression")), pch=16, cex=.8)
par(xpd=FALSE)
abline(0,1)
u <- par("usr")
par(xpd=NA)
y <- glm$coefficients[w,-1]+glm$coefficients[w,1]
u <- par("usr")
x0 <- rep(u[4]-(u[4]-u[3])/8,ncol(tcgaDesign)-1)
y0 <- u[4] + 0.05*(u[4]-u[3]) - rank(-y)/length(y) * (u[4]-u[3])/1.2
d <- density(y)
lines(d$x, d$y/5+u[4]-(u[4]-u[3])/8, col="grey")
lines(d$x, -d$y/5+u[4]-(u[4]-u[3])/8, col="grey")
points(x=y, y=x0+violinJitter(y, magnitude=0.25)$y,pch=19, col=col1[groups[-1]])
text(x=glm$coefficients[w,1], y= u[4], "Model coefficients (logFC)", cex=0.8)
v <- glm$p.value[w,-1] < 0.01
rotatedLabel(y[v], x0[v]+0.1, labels=colnames(tcgaDesign)[-1][v], font=ifelse(grepl("[[:lower:]]", colnames(tcgaDesign)[-1]),1,3)[v], cex=.66, pos=1)
axis(at=-1:1 + glm$coefficients[w,1], labels=-1:1, side=3, cex.axis=.8, line=-1, mgp = c(1.5,.05,0), tcl=-.15)
text(u[2],u[3] + (u[4]-u[3])/10, substitute(paste(R^2==r),list(r=round(Rgenetics[w],2))), pos=2)
}


#' #### Barplot
#+ TCGA-GE-quantiles,  fig.width=4, fig.height=2
par(bty="n", mgp = c(2.5,.33,0), mar=c(3.5,3.3,2,0)+.1, las=2, tcl=-.25)
b <- barplot(significantGenes, las=2, ylab = "Differentially expressed genes", col=brewer.pal(8,"RdYlBu"), legend.text=FALSE , border=0, xaxt="n")#, col = set1[simple.annot[names(n)]], border=NA)
rotatedLabel(x0=b, y0=rep(10, ncol(significantGenes)), labels=colnames(significantGenes), cex=.7, srt=45, font=ifelse(grepl("[[:lower:]]", colnames(tcgaDesign))[-1], 1,3), col=col1[groups[-1]])
clip(0,30,0,1000)
#text(b+0.2, colSums(n)+50, colSums(n), pos=3, cex=.7, srt=90)
x0 <- max(b) - 1.5
image(x=x0+c(0,0.8), y=par("usr")[4]+seq(-100,100,l=9), z=matrix(1:8, ncol=8), col=brewer.pal(8,"RdYlBu"), add=TRUE)
text(x=x0+1.5, y=par("usr")[4]+seq(-50,50,l=3), format(seq(-1,1,l=3),2), cex=0.66)
lines(x=rep(x0+.8,2), y=par("usr")[4]+c(-75,75))
segments(x0+.8,par("usr")[4]+seq(-75,75,l=7),x0+.9,par("usr")[4]+seq(-75,75,l=7))
text(x0+.8, par("usr")[4]+125, "log2 FC", cex=.66)
rotatedLabel(b-0.1, colSums(significantGenes), colSums(significantGenes), pos=3, cex=, srt=45)

#' #### Associated mutations per transcript:
#+ TCGA-GE-transcripts,  fig.width=2, fig.height=2
par(bty="n", mgp = c(2.5,.33,0), mar=c(3,3.3,3,0)+.1, las=1, tcl=-.25)
t <- table(rowSums(abs(testResults[,1:16])))
b <- barplot(t[-1],ylab="Differentially expressed genes", col=rev(brewer.pal(7, "Spectral")[-(4:5)]), border=NA)
rotatedLabel(b-0.1, t[-1], t[-1], pos=3, cex=1, srt=45)
title(xlab="Associated drivers", line=2)


#' #### Chromosomal distribution
chr = factor(sapply(AnnotationDbi::mget(rownames(tcgaExpr), org.Hs.egCHR, ifnotfound=NA), `[`,1), levels=c(1:22, "X","Y","MT"))
chromTable <- apply(testResults,2, function(x) table(chr[x!=0]))

#+ TCGA-GE-CHR,  fig.width=6, fig.height=6
par(bty="n", mgp = c(0.5,0.5,0), las=1, tcl=-.25, font.main=3, mfrow=c(6,6), xpd=NA, mar=c(0,0,1.5,0))
for(j in 1:ncol(testResults)){
	n <- sum(testResults[,j]!=0)
	pie(chromTable[,j], col=colorRampPalette(brewer.pal(11,'Spectral'))(24), border="white",  radius=0.8, init.angle=90, labels=ifelse(chromTable[,j]/sum(chromTable[,j]) > 0.02, paste("",rownames(chromTable), "(" ,chromTable[,j], ")",sep=""),""))
	title(main = colnames(chromTable)[j], font.main = ifelse(grepl("[[:lower:]]", colnames(chromTable)[j]), 1,3), cex.main=1.33)
	symbols(0,0,circles=.3, inches=FALSE, col="white", bg="white", lty=0, add=TRUE)
	#symbols(0,0,circles=.8*(1-sqrt(n/max(colSums(testResults!=0)))), col="white", add=TRUE, lty=0, bg="white", inches=FALSE)
	#cat(n,"\n")
}
t <- table(chr)
pie(t, col=colorRampPalette(brewer.pal(11,'Spectral'))(24), border="white",  radius=0.8, cex.main=1.33, labels=ifelse(t/sum(t) > 0.02, names(t),""), init.angle=90)
symbols(0,0,circles=.3, inches=FALSE, col="white", bg="white", lty=0, add=TRUE)
title(main = "# Genes", font.main = ifelse(grepl("[[:lower:]]", colnames(chromTable)[j]), 1,3))




#' 2. Survival analyses
#' -----------------
#' ### 2.1 Prepare data
library(CoxHD)
tcgaData <- data.frame(tcgaMutation, tcgaClinical[14:24])
tcgaData <- tcgaData[colSums(tcgaData, na.rm=TRUE)>=5]
tcgaData <- ImputeXMissing(data.frame(tcgaData, Gender=tcgaClinical$Gender, scale(tcgaClinical[,c(7:12)])))
rownames(tcgaData) <- rownames(tcgaMutation)
tcgaData <- tcgaData[rowSums(is.na(tcgaData))==0 & rownames(tcgaData) %in% colnames(tcgaExpr),]
dataFrame <- data.frame(tcgaData, scale(tcgaPca$x[rownames(tcgaData),1:20]))
survivalGroups <- rep("Blood", ncol(dataFrame))
survivalGroups[colnames(tcgaData)%in%colnames(tcgaMutation)] <- "Genetics"
survivalGroups[grep("^PC", colnames(dataFrame))] <- "Expression"
survivalGroups[grep("^[a-z]", colnames(dataFrame))] <- "CNA"
survivalGroups[grep("(t_)|_t", colnames(dataFrame))] <- "Translocations"
survivalGroups[grep("(Gender)|(AOD)", colnames(dataFrame))] <- "Demographics"
survivalGroups <- factor(survivalGroups, levels=c("Genetics","CNA","Translocations","Expression","Demographics","Blood"))
survivalCol=set1[c(3,2,5,4,7,1)]#  c(brewer.pal(8, "Dark2")[1:3], brewer.pal(8, "Set1")[2:1])
names(survivalCol) <- levels(survivalGroups)

library(survival)
tcgaSurvival <- Surv(tcgaClinical$OS + .5, tcgaClinical$Status)[match(rownames(tcgaData),tcgaClinical$TCGA_ID)]

#' ### 2.2 Fit model and compute variance components of the predicted log-hazard
#+ TCGA-varianceComponents,  fig.width=3.5, fig.height=2
coxRFX <- CoxRFX(dataFrame, tcgaSurvival, which.mu = NULL)
VarianceComponents(coxRFX, groups=survivalGroups)
PlotVarianceComponents(coxRFX, col=survivalCol, groups=survivalGroups)
points(0,0,pch=16, col="white", cex=15)
title(main="Variance components AML")

#' Using a jacknife-type subsampling approach
#+ TCGA-coxKnife, cache=TRUE, fig.width=3.5, fig.height=2, warn=FALSE
coxRFXKnife <- sapply(1:50, function(i){
			w <- sample(1:nrow(dataFrame), floor(nrow(dataFrame)*2/3))
			coef(CoxRFX(dataFrame[w,], tcgaSurvival[w],  which.mu = NULL))
		})
survConcordance(tcgaSurvival ~ as.matrix(dataFrame) %*% rowMeans(coxRFXKnife, na.rm=TRUE))
c <- coxRFX
c$coefficients <- rowMeans(coxRFXKnife, na.rm=TRUE)
c$var <- cov(t(coxRFXKnife), use='p')
VarianceComponents(c)
PlotVarianceComponents(c, groups=survivalGroups,col=survivalCol)

#' #### Stability selection
#+ TCGA-CoxCPSS, cache=TRUE, warning=FALSE
w <- !is.na(tcgaSurvival) & tcgaSurvival[,1]>0
CoxCPSS(dataFrame[w,], tcgaSurvival[w], control="BH")

#' #### Stepwise model selection
#+ TCGA-CoxBIC
c <- coxph(tcgaSurvival ~ 1, data=dataFrame)
scopeStep <- as.formula(paste("tcgaSurvival ~", paste(colnames(dataFrame), collapse="+")))
coxBIC <- step(c, scope=scopeStep, k = log(sum(!is.na(tcgaSurvival))), trace=0)
summary(coxBIC)


#' #### Random forests
#+ TCGA-RSF, cache=TRUE
library(randomForestSRC)
m <- match(rownames(tcgaPca$x),tcgaClinical$TCGA_ID)
rsf <- rfsrc(Surv(time, status) ~ ., data=data.frame(time=tcgaClinical$OS[m], status = tcgaClinical$Status[m], dataFrame ), ntree=100)

#' Variable importance
#+ TCGA-RSF-boxplot,  fig.width=1.5, fig.height=2
par(bty="n", mgp = c(2,.33,0), mar=c(4,3,1,0.5)+.1, las=2, tcl=-.25, las=3, xpd=NA)
boxplot(rsf$importance ~ survivalGroups, border= survivalCol, staplewex=0, pch=16, cex=0.75, ylab="Variable importance", lty=1, xaxt="n")
rotatedLabel(x0=1:nlevels(survivalGroups), y0=rep(-0.002,nlevels(survivalGroups)), labels=levels(survivalGroups), srt=45)
#' The plot confirms the result that expression, blood counts, and clinical variables are most influental.

#' #### Subset
set.seed(42)
cvIdx <- sample(1:5, nrow(dataFrame), replace=TRUE)

#+ TCGA-concordance, cache=TRUE, warning=FALSE, fig.height=2, fig.width=2.5
subsets <- list(Genetics="Genetics", Cytogenetics=c("Translocations","CNA"), Blood="Blood", Demographics="Demographics", Expression="Expression", `Gen+Cyt` = c("Genetics","Translocations","CNA"), `Gen+Cyt+Blo+Exp`=c("Genetics","Translocations","CNA","Blood","Expression"), All=unique(survivalGroups))
colSubsets <- set1
names(colSubsets) <- c("Blood","Cytogenetics", "Genetics","Expression","Gen+Cyt","Gen+Cyt+Blo+Exp","Demographics","All","Std. Risk")
concordance <- sapply(1:5, function(i){
			v <- cvIdx == i
			c(sapply(subsets, function(s,v){
						w <- survivalGroups %in% s
						fit <- CoxRFX(dataFrame[!v, w]+0, tcgaSurvival[!v], which.mu=NULL)
						p <- as.matrix(dataFrame[v,w]) %*% coef(fit)
						survConcordance(tcgaSurvival[v]~p)$concordance
					}, v=v), 
			Std.Risk =  survConcordance(tcgaSurvival[v]~ c(3,1,2)[tcgaClinical$C_Risk[match(rownames(tcgaData),tcgaClinical$TCGA_ID)][v]])$concordance
			)
		})
rownames(concordance) <- c(names(subsets), "Std. Risk")
par(mar=c(5,4,1,1), mgp=c(3,0.5,0))
m <- rowMeans(concordance, na.rm=TRUE)
e <- apply(concordance,1,var)/ncol(concordance)
o <- order(m)
barplot(m[o], col=colSubsets[names(m[o])], names.arg=rep("", nrow(concordance)), ylim=c(0.5,0.75), xpd=FALSE, ylab="Concordance (5x CV)") -> b
points(jitter(rep(b,5)), concordance[o,], col="grey", pch=16, cex=.5)
rotatedLabel(b,rep(0.49, nrow(concordance)), rownames(concordance)[o])
segments(b, m[o]+sqrt(e)[o], b , m[o]-sqrt(e)[o])



#' Session
#' -----
#+ sessionInfo, eval=TRUE
sessionInfo()