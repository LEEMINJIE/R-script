# R-script
TNBC R script

BiocManager::install(version="3.12") 	# BiocManager v1.30.10 # do once
library(BiocManager)
BiocManager::version() 			# check the version

#BiocManager::install(c("affyio", "affy", "limma", "oligo"), ask=F) # run every time
# This may be slow because of the large number of dependencies
library(affy)              # activate library, run every time
library(affyio)            # activate library, run every time
library(oligo)             # activate library, run every time
library(limma)             # activate library, run every time

#BiocManager::install(c("GEOquery"), ask=F)
library(GEOquery)
#BiocManager::install(c("Biobase"), ask=F)
library(Biobase)
#BiocManager::install(c("knitr"), ask=F)
library(knitr)
#BiocManager::install(c("simpleaffy"), ask=F)
library(simpleaffy)
#BiocManager::install(c("GenomicRanges"), ask=F)
library(GenomicRanges)
#BiocManager::install(c("Organism.dplyr"), ask=F)
library(Organism.dplyr)
#BiocManager::install(c("illuminaio"), ask=F)
library(illuminaio)
#install.packages("XML", ask=F)
library(XML)
#install.packages("statmod", ask=F)
library(statmod)
#install.packages("cgdsr", ask=F)
library(cgdsr)
#install.packages("R2HTML", ask=F)
library(R2HTML)
#install.packages("hash", ask=F)
library(hash)
#install.packages("ggplot2", ask=F)
library(ggplot2)
#install.packages("ggrepel", ask=F)
library(ggrepel)
#install.packages("gplots", ask=F)
library(gplots)
#install.packages("heatmap.plus", ask=F)
library(heatmap.plus)
#install.packages("beeswarm", ask=F)
library(beeswarm)
#install.packages("xtable", ask=F) # let this re-start, hit "yes"
library(xtable)
#install.packages("rmarkdown")  # let this re-start, hit "yes"
library(rmarkdown)
#install.packages("Hmisc")
library(Hmisc)
#install_github("kassambara/factoextra")
library("factoextra")
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
library("ggpubr")

sessionInfo() # to check versions

library(biomaRt) # biomaRt v2.46.2 # 
citation("biomaRt")
#BiocManager::install("hgu133a.db", ask=F) # hgu133a.db v3.2.3
library(hgu133a.db) # https://bioconductor.org/packages/release/data/annotation/html/hgu133a.db.html
citation("hgu133a.db")
#BiocManager::install("huex10sttranscriptcluster.db", ask=F)  
library(huex10sttranscriptcluster.db)
library(data.table)
#install.packages("gridExtra")
library(gridExtra)

huex <- contents(hgu133aACCNUM)
length(huex) # list of 22,283 of linked probes and genes
mapped_probes <- mappedkeys(huex) # Convert to a list
xxhuex <- as.list(huex[mapped_probes]) # compare mapped probes # 22,283 elements

#----------------------------------------------- write csv file for GSE25066 dataset ----------------------------------------------- 
geo <- getGEO("GSE25066",GSEMatrix=TRUE,getGPL=FALSE)
if (length(geo) > 1) idx <- grep("GPL96", attr(geo, "names")) else idx <- 1
expression.set <- geo[[idx]]
str(expression.set)

expression.set2 <- exprs(expression.set) #show as matrix
str(expression.set2)
head(expression.set2)

data <- pData(expression.set)[,] #dataframe with 508 rows and 80 columns
write.csv(x=data, file="GSE25066.csv") # dataframe is saved as csv file in the CELFILES directory


#----------------------------------------------- read data ----------------------------------------------- 
##Part 2 - QC
names <- 1:508
data <- read.csv("GSE25066.csv")
colnames(data)
data = subset(data, select = -c(X)) # remove the first column
celfiles <- list.files(getwd(), pattern="CEL") # get all celfiles name 
data$celfiles= celfiles # add a column called “celfiles” to the dataframe
 
#rename metadata column name
sample.no. <- data[,2] 				# accession number
pam50_class<- data[,71] 			# pam50_class 
nodal_status <- data[,59] 	   		# clinical_nodal_status.ch1
pcr_rd <- data[,72]				 #pathologic_response_pcr_rd:ch1
rcb <- data[,73] 				# pathologic_response_rcb_class:ch1 (residual cancer burden)
her2.status <- data[,70] 			# her2_status:ch1
er.status <- data[,65]				 # er_status_ihc:ch1
pr.status  <- data[,74] 			# pr_status_ihc:ch1
er_status_ihc <- data[,64] 			#er_status_ihc_esr1_for.indeterminate.ch1
erbb2_status <-  data[,66] 			#erbb2_status.ch1
clinical_t_stage  <- data[,60] 			#"clinical_t_stage.ch1"     
ajcc <- data[,58]				#"clinical_ajcc_stage.ch1"  
grade   <- data[,69] 				# grade.ch1
ggi_class <- data[,68]				# ggi_class.ch1
set_class <- data[,77]				#set_class.ch1
rcb_prediction <- data[,75]			# rcb_0_i_prediction.ch1
chemosensitivity.prediction <- data[,57] 	#chemosensitivity_prediction:ch1
dlda30  <- data[,61] 				#dlda30_prediction:ch1
drfs <- data[,62]				#drfs_1_event_0_censored.ch1 (distance relapse free survival)
						#censoring status; 0=censored, 1=event.
drfs_years <- data[,63]				#drfs_even_time_years.ch1
celfiles  <- data[,81]                          #celfies

# make a smaller dataframe (sdata) with only some metadata are included
sdata <- data.frame(sample.no., her2.status, er.status, pr.status,er_status_ihc, erbb2_status, ajcc, nodal_status,set_class, rcb_prediction, chemosensitivity.prediction, pcr_rd, rcb, clinical_t_stage, grade,drfs_years ,celfiles) 

#filter out TNBC subtype 
tnbc <- sdata %>% filter(her2.status == "N", er.status == "N", pr.status == "N") # filter out TNBC subtypes
dim(tnbc) # check the total number of tnbc (should be 178)
tnbc <- tnbc[order(tnbc$pcr_rd),] # sort groups in response to pcr / rd
tnbc$subtypes = 'TNBC'

# remove tnbc from the dataframe
remove <- tnbc$sample.no.
sdatawithouttnbc <- sdata[!sdata$sample.no.%in%remove,] #remove tnbc from the dataframe
data2 <- sdatawithouttnbc
data2$subtypes = 'non-TNBC'
nontnbc <- data2

#merge tnbc with the sorted dataframe
mydata_tnbc_nontnbc <- rbind(tnbc,nontnbc)

#rawdata
rawdata_tnbc_nontnbc  <- read.celfiles(mydata_tnbc_nontnbc$celfiles)
eset <-exprs(rawdata_tnbc_nontnbc)
length(rownames(eset)) #506944 probes

# draw heatmap of raw data
png(file="rawheatmap-all.png", width=2000, height=2000, units="px")
heatmap(eset[1:100,], Colv=as.dendrogram(hclust(dist(t(eset)))), cexRow=2, cexCol=2, labCol = names, labRow = c(), margins=c(42,8))
dev.off() 
# draw boxplot of raw data expression levels
png(file="rawboxplot-all.png", width=1200, height=1200, units="px")
boxplot(log2(eset), las=2, ylab='log2 of expression level', names=names, col="red", margins = c(10, 10), cex=1)
dev.off() # log2 here because no transformation completed

# draw PCA of raw data expression levels
pca1<- prcomp(t(eset))
pcvalue1 = round(100*pca1$sdev[1]^2/(sum(pca1$sdev^2)),1)
pcvalue2 = round(100*pca1$sdev[2]^2/(sum(pca1$sdev^2)),1)
pcvalue1width  = 40*pcvalue1 # pixel width PC1				
pcvalue1height = 40*pcvalue2 # pixel height PC2				

png(file="rawPCA-all.png", units="px")
plot(pca1$x, col=mycol, pch=mypch,xlab=pcvalue1, ylab=pcvalue2, main='PCA of all raw data', cex=1, cex.axis=1, cex.lab=1)
text(pca1$x, names, pos=3.9, cex =0.5)
#legend("right", legend= names, col=mycol, pch=mypch, cex = 0.55, inset=c(-0.1,-0.1),xpd = TRUE) 
dev.off()  

#----------------------------------------------- normalisation ----------------------------------------------- 
#normalise the rawdata
normdata_tnbc_nontnbc <-rma(rawdata_tnbc_nontnbc)
eset2 <-exprs(normdata_tnbc_nontnbc)
length(rownames(eset2)) #22283 probes

# draw heatmap of normalised data # just 100 probes to get an idea
png(file="finalheatmap-all.png", width=2000, height=2000, units="px")
heatmap(eset2[1:100,], Colv=as.dendrogram(hclust(dist(data2))), cexRow=3, cexCol=3, col=colorRampPalette(c("red", "yellow", "green"))(n = 299), labCol = names, labRow = c(), margins=c(42,8))
dev.off()

# draw boxplot of normalised data expression levels
png(file="finalboxplot-all.png", units="px")
boxplot(log2(eset2), las=2, ylab='Normalised expression level', cex.lab= 1.5, cex.axis= 1,names=names, col="red", cex=1)
dev.off() # log2-transformation completed by rma() above already


#----------------------------------------------- allocating breast cancer subtypes ----------------------------------------------- 
# filter out lumA
lumA <- sdata %>% filter(er.status == "P", pr.status == "P", her2.status == "N") # filter out lumA subtypes
dim(lumA) # check the total number of lumA (216)

PNN <- sdata %>% filter(er.status == "P", pr.status == "N", her2.status == "N") # filter out subtypes
dim(PNN) # check the total number of PNN (72)
NPN <- sdata %>% filter(er.status == "N", pr.status == "P", her2.status == "N") # filter out subtypes
dim(NPN) # check the total number of NPN (17)
int1 <- sdata %>% filter(er.status == "I", er_status_ihc == "P", erbb2_status == "N") 
int2 <- sdata %>% filter(er.status == "I", er_status_ihc == "N", erbb2_status == "N")
lumA <-rbind(lumA,PNN,NPN,int1,int2)
lumA$subtypes = 'lumA'

lumB <- sdata %>% filter(er.status == "P", pr.status == "P", her2.status == "P") # filter out lumB subtypes
dim(lumB) # check the total number of lumB (2)
lumB <- lumB[order(lumB$pcr_rd),] # sort groups in response to pcr / rd
lumB$subtypes = 'lumB'

her2 <- sdata %>% filter(er.status == "N", pr.status == "N", her2.status == "P") # filter out her2 subtypes
dim(her2) # check the total number of her2 (3)
NPP <- sdata %>% filter(er.status == "N", pr.status == "P", her2.status == "P") # filter out subtypes
dim(NPP) # check the total number of NPP (1)
her2 <-rbind(her2,NPP)
her2$subtypes = 'her2'
her2 <- her2[order(her2$pcr_rd),] # sort groups in response to pcr / rd

#----------------------------------------------- PCA plots (after normalisation) ------------------------------------------
newdata <- rbind(tnbc,lumA,lumB,her2)
na <- sdata[!sdata$sample.no.%in%newdata$sample.no.,]  
na$subtypes <- 'N/A'
newdata <- rbind(newdata,na)
newdata <- read.celfiles(newdata$celfiles)
normnewdata<-rma(newdata)  #normalise the newdata
eset3<- exprs(normnewdata)
pca_proc1<- prcomp(t(eset3))

pcvalue1 = round(100*pca_proc1$sdev[1]^2/(sum(pca_proc1$sdev^2)),1)
pcvalue2 = round(100*pca_proc1$sdev[2]^2/(sum(pca_proc1$sdev^2)),1)
pcvalue3 = round(100*pca_proc1$sdev[3]^2/(sum(pca_proc1$sdev^2)),1)
pcvalue4 = round(100*pca_proc1$sdev[4]^2/(sum(pca_proc1$sdev^2)),1)
 
pcvalue1hw = 40*pcvalue1 # pixel PC1
pcvalue2hw = 40*pcvalue2 # pixel PC2
pcvalue3hw = 40*pcvalue3 # pixel PC3
pcvalue4hw = 40*pcvalue4 # pixel PC4
 
mycol <- c(rep("blue",178),rep("orange",309),rep("green",2),rep("red",4),rep("black",15)) # group samples
mypch = rep(20,508) # pattern
names <- 1:508

#PC1 vs PC2
png(file="PCA-1v2-IHCsubtypes.png", width=pcvalue1hw, height=pcvalue2hw, units="px")
plot(pca_proc1$x[,1], pca_proc1$x[,2], col=mycol, pch=mypch, xlab=paste0("PC1: ",pcvalue1,"% variance"), ylab=paste0("PC2: ",pcvalue2,"% variance"),xlim=c(-65,65), ylim=c(-60,60),main='PC1 vs PC2', cex=1, cex.axis=1, cex.lab=1)
legend("topleft", col=unique(mycol), legend = c('TNBC','LumA','LumB','Her2','N/A'),pch = 20, bty='n', cex=.75)
dev.off()

#PC2 vs PC3
png(file="PCA-2v3-IHCsubtypes.png", width=pcvalue1hw, height=pcvalue2hw, units="px")
plot(pca_proc1$x[,2], pca_proc1$x[,3], col=mycol, pch=mypch,  xlab=paste0("PC2: ",pcvalue2,"% variance"), ylab=paste0("PC3: ",pcvalue3,"% variance"),xlim=c(-65,65), ylim=c(-60,60), main='PC2 vs PC3', cex=1, cex.axis=1, cex.lab=1) 
legend("topleft", col=unique(mycol), legend =  c('TNBC','LumA','LumB','Her2','N/A'),pch = 20, bty='n', cex=.75)
dev.off()

#PC3 vs PC4
png(file="PCA-3v4-IHCsubtypes.png", width=pcvalue1hw, height=pcvalue2hw, units="px")
plot(pca_proc1$x[,3], pca_proc1$x[,4], col=mycol, pch=mypch,  xlab=paste0("PC3: ",pcvalue3,"% variance"), ylab=paste0("PC4: ",pcvalue4,"% variance"),xlim=c(-65,65), ylim=c(-60,60), main='PC3 vs PC4', cex=1, cex.axis=1, cex.lab=1) 
legend("topleft", col=unique(mycol), legend =  c('TNBC','LumA','LumB','Her2','N/A'),pch = 20, bty='n', cex=.75)
dev.off()

#PC1 vs PC3
png(file="PCA-1v3-IHCsubtypes.png", width=pcvalue1hw, height=pcvalue2hw, units="px")
plot(pca_proc1$x[,1], pca_proc1$x[,3], col=mycol, pch=mypch,  xlab=paste0("PC1: ",pcvalue1,"% variance"), ylab=paste0("PC3: ",pcvalue3,"% variance"),xlim=c(-65,65), ylim=c(-60,60), main='PC1 vs PC3', cex=1, cex.axis=1, cex.lab=1) 
legend("topleft", col=unique(mycol), legend =  c('TNBC','LumA','LumB','Her2','N/A'),pch = 20, bty='n', cex=.75)
dev.off()

#PC1 vs PC4
png(file="PCA-1v4-IHCsubtypes.png", width=pcvalue1hw, height=pcvalue2hw, units="px")
plot(pca_proc1$x[,1], pca_proc1$x[,4], col=mycol, pch=mypch,  xlab=paste0("PC1: ",pcvalue1,"% variance"), ylab=paste0("PC4: ",pcvalue4,"% variance"),xlim=c(-65,65), ylim=c(-60,60), main='PC1 vs PC4', cex=1, cex.axis=1, cex.lab=1) 
legend("topleft", col=unique(mycol), legend =  c('TNBC','LumA','LumB','Her2','N/A'),pch = 20, bty='n', cex=.75)
dev.off()

#PC2 vs PC4
png(file="PCA-2v4-IHCsubtypes.png", width=pcvalue1hw, height=pcvalue2hw, units="px")
plot(pca_proc1$x[,2], pca_proc1$x[,4], col=mycol, pch=mypch,  xlab=paste0("PC2: ",pcvalue2,"% variance"), ylab=paste0("PC4: ",pcvalue4,"% variance"),xlim=c(-65,65), ylim=c(-60,60), main='PC2 vs PC4', cex=1, cex.axis=1, cex.lab=1) 
legend("topleft", col=unique(mycol), legend =  c('TNBC','LumA','LumB','Her2','N/A'),pch = 20, bty='n', cex=.75)
dev.off()

#----------------------------------------------- DGE analysis (based on patient response) ------------------------------------------

# DGE analysis of TNBC vs non-TNBC
# volcano plots for each pairwise comparison
library(limma) # check limma is working, use install.packages(“limma”) if it is not 
mypch2_tnbc_nontnbc <- c(rep(1,178),rep(2,330)) # 1:330 in mydata_tnbc_nontnbc= non-TNBC , 331:508 = TNBC 
design_tnbc_nontnbc <- model.matrix(~-1+factor(mypch2_tnbc_nontnbc))
colnames(design_tnbc_nontnbc )<- c('T', 'N') #  group labels
top0_tnbc_nontnbc <-  eBayes(contrasts.fit(lmFit(normdata_tnbc_nontnbc, design_tnbc_nontnbc), makeContrasts(T - N, levels=design_tnbc_nontnbc)))
length(rownames(eset2))
top1_tnbc_nontnbc <- topTable(top0_tnbc_nontnbc, coef='T - N',adjust='BH',number=22283,sort.by='P')

#upregulated probes in non-TNBC
up.reg1_tnbc_nontnbc  <- top1_tnbc_nontnbc%>%filter(top1_tnbc_nontnbc$logFC < -1.5)
up.reg1_tnbc_nontnbc  <- up.reg1_tnbc_nontnbc %>%filter(up.reg1_tnbc_nontnbc$ adj.P.Val < 0.05)
names_upreg_nontnbc <- rownames(up.reg1_tnbc_nontnbc) #31 probes
#upregulated probes in TNBC
up.reg2_tnbc_nontnbc <- top1_tnbc_nontnbc%>%filter(top1_tnbc_nontnbc$logFC > 1.5)
up.reg2_tnbc_nontnbc <- up.reg2_tnbc_nontnbc%>%filter(up.reg2_tnbc_nontnbc$ adj.P.Val  <0.05)
names_upreg_tnbc <- rownames(up.reg2_tnbc_nontnbc) #9 probes

# ** get the symbol of the differentially expressed probes **
name1 <- select(hgu133a.db,names_upreg_nontnbc,c('SYMBOL','GENENAME'))
symbol_upreg_nontnbc <- name1$SYMBOL
name2 <- select(hgu133a.db,names_upreg_tnbc,c('SYMBOL','GENENAME'))
symbol_upreg_tnbc <- name2$SYMBOL

png(file="volcano-TNBC_vs_non-TNBC.png",units="px")
plot(top1_tnbc_nontnbc$logFC,-log10(top1_tnbc_nontnbc$adj.P.Val), main='Differential Gene Expression between TNBC and Non-TNBC Subtypes', xlab='log2FC(TNBC/Non-TNBC)',ylab='-log10(BH P value)', cex=0.5, cex.main=1,cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.5,70,label='Upregulated in \nTNBC',cex=1)
text(-1.5,70,label='Upregulated in \nNon-TNBC',cex=1)
abline(v=-1.5, col="red")
abline(v=1.5, col="red")
abline(h=1.3, col="green")
dev.off()


#filter out each population ----- TNBC(RD), TNBC(pCR), non-TNBC(RD), non-TNBC(pCR)
rd_tnbc <- tnbc%>%filter(pcr_rd == 'RD') #113
pcr_tnbc <- tnbc%>%filter(pcr_rd == 'pCR') #57
rd_data2 <- data2%>%filter(pcr_rd == 'RD') #276
pcr_data2 <- data2%>%filter(pcr_rd == 'pCR') #42
mydata_rd_pcr_tnbc_nontnbc <- rbind(rd_tnbc,pcr_tnbc,rd_data2,pcr_data2)
other <- sdata[!sdata$sample.no.%in%mydata_rd_pcr_tnbc_nontnbc$sample.no.,] #20
other$subtypes <- 'other'

mydata_rd_pcr_tnbc_nontnbc <- rbind(rd_tnbc,pcr_tnbc,rd_data2,pcr_data2,other)
rawdata_rd_pcr_tnbc_nontnbc <- read.celfiles(mydata_rd_pcr_tnbc_nontnbc$celfiles)
normdata_rd_pcr_tnbc_nontnbc <-rma(rawdata_rd_pcr_tnbc_nontnbc) 

# make volcano plots for each pairwise comparison
# DGE analysis of TNBC(RD) vs non-TNBC(RD)
mypch2_rd_pcr_tnbc_nontnbc <- c(rep(1,113),rep(2,57),rep(3,276),rep(4,42),rep(5,20))	
design_rd_pcr_tnbc_nontnbc <- model.matrix(~-1+factor(mypch2_rd_pcr_tnbc_nontnbc))
colnames(design_rd_pcr_tnbc_nontnbc)<- c('T.rd', 'T.pcr','N.rd', 'N.pcr','other') 

top0_rd_tnbc_nontnbc <-  eBayes(contrasts.fit(lmFit(normdata_rd_pcr_tnbc_nontnbc, design_rd_pcr_tnbc_nontnbc), makeContrasts(T.rd - N.rd, levels=design_rd_pcr_tnbc_nontnbc)))
top1_rd_tnbc_nontnbc <- topTable(top0_rd_tnbc_nontnbc, coef='T.rd - N.rd',adjust='BH',number=22283,sort.by='P') 

#upregulated probe_id in non-TNBC(RD)
up.reg1_rd_tnbc_nontnbc  <- top1_rd_tnbc_nontnbc%>%filter(top1_rd_tnbc_nontnbc$logFC < -1.5)
up.reg1_rd_tnbc_nontnbc <- up.reg1_rd_tnbc_nontnbc%>%filter(up.reg1_rd_tnbc_nontnbc$ adj.P.Val <0.05)
names_upreg_rd_nontnbc<- rownames(up.reg1_rd_tnbc_nontnbc)
length(names_upreg_rd_nontnbc) #34
#upregulated probe_id in TNBC(RD)
up.reg2_rd_tnbc_nontnbc  <- top1_rd_tnbc_nontnbc%>%filter(top1_rd_tnbc_nontnbc$logFC > 1.5)
up.reg2_rd_tnbc_nontnbc  <- up.reg2_rd_tnbc_nontnbc%>%filter(up.reg2_rd_tnbc_nontnbc$ adj.P.Val <0.05)
names_upreg_rd_tnbc <- rownames(up.reg2_rd_tnbc_nontnbc)
length(names_upreg_rd_tnbc)  #9

png(file="volcano-TNBC(RD)_vs_non-TNBC(RD).png",units="px")
plot(top1_rd_tnbc_nontnbc$logFC,-log10(top1_rd_tnbc_nontnbc$adj.P.Val), main= 'Differential Gene Expression \n between TNBC (RD) and Non-TNBC Subtypes (RD)', xlab='log2FC(TNBC(RD)/Non-TNBC(RD))', ylab='-log10(BH P value)', cex=0.5, cex.main=1,cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.5,50,label='Upregulated in \n TNBC(RD)',cex=1)
text(-1.3,50,label='Upregulated in \n Non-TNBC(RD)',cex=1)
abline(v=1.5, col="red")
abline(v=-1.5, col="red")
abline(h=1.3, col="green")
dev.off()

# DGE analysis of TNBC(pCR) vs non-TNBC(pCR)
# make volcano plots for each pairwise comparison
mypch2_rd_pcr_tnbc_nontnbc <- c(rep(1,113),rep(2,57),rep(3,276),rep(4,42),rep(5,20))	
design_rd_pcr_tnbc_nontnbc <- model.matrix(~-1+factor(mypch2_rd_pcr_tnbc_nontnbc))
colnames(design_rd_pcr_tnbc_nontnbc)<- c('T.rd', 'T.pcr','N.rd', 'N.pcr','other') 

top0_pcr_tnbc_nontnbc <-  eBayes(contrasts.fit(lmFit(normdata_rd_pcr_tnbc_nontnbc, design_rd_pcr_tnbc_nontnbc), makeContrasts(T.pcr - N.pcr, levels=design_rd_pcr_tnbc_nontnbc)))
top1_pcr_tnbc_nontnbc <- topTable(top0_pcr_tnbc_nontnbc, coef='T.pcr - N.pcr',adjust='BH',number=22283,sort.by='P') 

#upregulated probe_id in non-TNBC(pCR)
up.reg1_pcr_tnbc_nontnbc <- top1_pcr_tnbc_nontnbc%>%filter(top1_pcr_tnbc_nontnbc$logFC < -1.5)
up.reg1_pcr_tnbc_nontnbc <- up.reg1_pcr_tnbc_nontnbc%>%filter(up.reg1_pcr_tnbc_nontnbc$ adj.P.Val <0.05)
names_pcr_nontnbc  <- rownames(up.reg1_pcr_tnbc_nontnbc)
length(names_pcr_nontnbc) #12
#upregulated probe_id inTNBC(pCR)
up.reg2_pcr_tnbc_nontnbc <- top1_pcr_tnbc_nontnbc%>%filter(top1_pcr_tnbc_nontnbc$logFC > 1.5)
up.reg2_pcr_tnbc_nontnbc <- up.reg2_pcr_tnbc_nontnbc%>%filter(up.reg2_pcr_tnbc_nontnbc$ adj.P.Val<0.05)
names_pcr_tnbc <- rownames(up.reg2_pcr_tnbc_nontnbc)
length(names_pcr_tnbc) #3

png(file="volcano-TNBC(pCR)_vs_non-TNBC(pCR).png",units="px")
plot(top1_pcr_tnbc_nontnbc$logFC,-log10(top1_pcr_tnbc_nontnbc$adj.P.Val), main= 'Differential Gene Expression \n between TNBC (pCR) and Non-TNBC Subtypes (pCR)', xlab='log2FC(TNBC(pCR)/Non-TNBC(pCR))', ylab='-log10(BH P value)', cex=0.5, cex.main=1.2,cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.5,10,label='Upregulated in \n TNBC(pCR)',cex=1)
text(-1.3,10,label='Upregulated in \n non-TNBC(pCR)',cex=1)
abline(v=1.5, col="red")
abline(v=-1.5, col="red")
abline(h=1.3, col="green")
dev.off()

# DGE analysis of TNBC(RD) vs TNBC(pCR)
# make volcano plots for each pairwise comparison
mypch2_rd_pcr_tnbc_nontnbc <- c(rep(1,113),rep(2,57),rep(3,276),rep(4,42),rep(5,20))	
design_rd_pcr_tnbc_nontnbc <- model.matrix(~-1+factor(mypch2_rd_pcr_tnbc_nontnbc))
colnames(design_rd_pcr_tnbc_nontnbc)<- c('T.rd', 'T.pcr','N.rd', 'N.pcr','other') 

top0_tnbc_response <-  eBayes(contrasts.fit(lmFit(normdata_rd_pcr_tnbc_nontnbc, design_rd_pcr_tnbc_nontnbc), makeContrasts(T.rd - T.pcr, levels=design_rd_pcr_tnbc_nontnbc)))
top1_tnbc_response <- topTable(top0_tnbc_response, coef='T.rd - T.pcr',adjust='BH',number=22283,sort.by='P') 

#upregulated probes in TNBC(RD)
up.reg1_tnbc_response <- top1_tnbc_response%>%filter(top1_tnbc_response$ adj.P.Val <0.05)
up.reg1_tnbc_response<- up.reg1_tnbc_response%>%filter(up.reg1_tnbc_response$logFC > 1.5)
names_upreg_tnbc_rd <- rownames(up.reg1_tnbc_response)
#upregulated probes in TNBC(pCR)
up.reg2_tnbc_response <- top1_tnbc_response%>%filter(top1_tnbc_response$ adj.P.Val<0.05)
up.reg2_tnbc_response<- up.reg2_tnbc_response%>%filter(up.reg2_tnbc_response$logFC < -1.5)
names_upreg_tnbc_pcr <- rownames(up.reg2_tnbc_response)

png(file="volcano-TNBC(RD)_vs_TNBC(pCR).png",units="px")
plot(top1_tnbc_response$logFC,-log10(top1_tnbc_response$adj.P.Val),main= 'Differential Gene Expression between TNBC (RD) and TNBC (pCR)', xlab='log2FC(TNBC-rd/TNBC-pCR)', ylab='-log10(BH P value)', cex=0.5, cex.main=1, cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.5,3,label='Upregulated in \n TNBC(RD)',cex=1)
text(-1.5,3,label='Upregulated in \n TNBC(pCR)',cex=1)
abline(v=1.5, col="red")
abline(v=-1.5, col="red")
abline(h=1.3, col="green")
dev.off()

# DGE analysis of non-TNBC(RD) vs non-TNBC(pCR)
# make volcano plots for each pairwise comparison
mypch2_rd_pcr_tnbc_nontnbc <- c(rep(1,113),rep(2,57),rep(3,276),rep(4,42),rep(5,20))	
design_rd_pcr_tnbc_nontnbc <- model.matrix(~-1+factor(mypch2_rd_pcr_tnbc_nontnbc))
colnames(design_rd_pcr_tnbc_nontnbc)<- c('T.rd', 'T.pcr','N.rd', 'N.pcr','other') 

top0_nontnbc_response <-  eBayes(contrasts.fit(lmFit(normdata_rd_pcr_tnbc_nontnbc, design_rd_pcr_tnbc_nontnbc), makeContrasts(N.rd - N.pcr, levels=design_rd_pcr_tnbc_nontnbc)))
top1_nontnbc_response <- topTable(top0_nontnbc_response, coef='N.rd - N.pcr',adjust='BH',number=22283,sort.by='P') 

#upregulated probe_id in non-TNBC(RD)
up.reg1_nontnbc_response  <- top1_nontnbc_response %>%filter(top1_nontnbc_response$logFC < -1.5)
up.reg1_nontnbc_response  <-  up.reg1_nontnbc_response%>%filter(up.reg1_nontnbc_response$adj.P.Val <0.05)
names_upreg_rd_nonTNBC <- rownames(up.reg1_nontnbc_response)
length(names_upreg_rd_nonTNBC) #0
#upregulated probe_id in non-TNBC(pCR)
up.reg2_nontnbc_response  <- top1_nontnbc_response %>%filter(top1_nontnbc_response$logFC > 1.5)
up.reg2_nontnbc_response  <- up.reg2_nontnbc_response%>%filter(up.reg2_nontnbc_response $ adj.P.Val) <0.05)
names_upreg_pcr_nonTNBC <- rownames(up.reg2_nontnbc_response)
length(names_upreg_pcr_nonTNBC) #1

png(file="volcano-nonTNBC(RD)_vs_nonTNBC(pCR).png",units="px")
plot(top1_nontnbc_response$logFC,-log10(top1_nontnbc_response$adj.P.Val), main= 'Differential Gene Expression \n between Non-TNBC (RD) and Non-TNBC (pCR)', xlab='log2FC(non-TNBC(RD)/non-TNBC(pCR))', ylab='-log10(BH P value)', cex=0.5,cex.main=1.2, cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.5,6.6,label='Upregulated in \n Non-TNBC(RD)',cex=0.9)
text(-1.5,6.6,label='Upregulated in \n Non-TNBC(pCR)',cex=0.9)
abline(v=1.5, col="red")
abline(v=-1.5, col="red")
abline(h=1.3, col="green")
dev.off()

# DGE analysis of lumA(RD) vs lumA(pCR)
lumA <- lumA[order(lumA$pcr_rd),]   # sort groups in response to pcr / rd

rd_lumA <- lumA%>%filter(pcr_rd == 'RD') #261
pcr_lumA <- lumA%>%filter(pcr_rd  == 'pCR') #37
lumA_response <- rbind(rd_lumA, pcr_lumA)
other <- sdata[!sdata$sample.no.%in%lumA_response$sample.no.,]   #210 samples
other$subtypes <- 'other'

lumA_response <- rbind(rd_lumA, pcr_lumA,other) #508 samples
rawdata_lumA_response  <- read.celfiles(lumA_response$celfiles) 
normdata_lumA_response<-rma(rawdata_lumA_response)

mypch2_lumA_response <- c(rep(1,261),rep(2,37),rep(3,210))	 
design_lumA_response<- model.matrix(~-1+factor(mypch2_lumA_response))
colnames(design_lumA_response)<- c('rd', 'pcr','o')  #  group labels
top0_lumA_response<- eBayes(contrasts.fit(lmFit(normdata_lumA_response, design_lumA_response), makeContrasts(rd - pcr, levels=design_lumA_response)))
top1_lumA_response<- topTable(top0_lumA_response, coef='rd - pcr',adjust='BH',number=22283,sort.by='P') 

#upregulated probes in lumA(RD) #no significant probe
up.reg1_lumA_response <- top1_lumA_response%>%filter(top1_lumA_response$adj.P.Val <0.05)
up.reg1_lumA_response <- up.reg1_lumA_response%>%filter(up.reg1_lumA_response$logFC > 1.5)
names_upreg_lumA_rd <- rownames(up.reg1_lumA_response)
#upregulated probes in lumA(pCR)
up.reg2_lumA_response <- top1_lumA_response%>%filter(top1_lumA_response$adj.P.Val<0.05)
up.reg2_lumA_response <- up.reg2_lumA_response%>%filter(up.reg2_lumA_response$logFC < -1.5)
names_upreg_lumA_pcr <- rownames(up.reg2_lumA_response)

png(file="volcano-lumA(RD)_vs_lumA(pCR).png", units="px")
plot(top1_lumA_response$logFC,-log10(top1_lumA_response$adj.P.Val), main= 'Differential Gene Expression \n  between LumA(RD) and LumA(pCR)',xlab='log2FC(LumA(RD)/LumA(pCR))', ylab='-log10(BH P value)', cex=0.5, ,cex.main=1.2, cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.5,6,label='Upregulated in \n LumA(RD)',cex=1)
text(-1.5,6,label='Upregulated in \n LumA(pCR)',cex=1)
abline(v=1.5, col="red")
abline(v=-1.5, col="red")
abline(h=1.3, col="green")
dev.off()

#------------------------------------------------ Fold change testing ------------------------------------------------
#lum A
#lum A(RD vs pCR): change fold change threshold to 1.2
#upregulated probes in lumA(RD) 
up.reg1_lumA_response <- top1_lumA_response%>%filter(top1_lumA_response$adj.P.Val<0.05)
up.reg1_lumA_response <- up.reg1_lumA_response%>%filter(up.reg1_lumA_response$logFC > 1.2)
names_upreg_lumA_rd <- rownames(up.reg1_lumA_response) #2 significant probe
#upregulated probes in lumA(pCR)
up.reg2_lumA_response <- top1_lumA_response%>%filter(top1_lumA_response$adj.P.Val<0.05)
up.reg2_lumA_response <- up.reg2_lumA_response%>%filter(up.reg2_lumA_response$logFC < -1.2)
names_upreg_lumA_pcr <- rownames(up.reg2_lumA_response) #3 significant probes

#lum A(RD vs pCR): change fold change threshold to 1.1
#upregulated probes in lumA(RD) 
up.reg1_lumA_response <- top1_lumA_response%>%filter(top1_lumA_response$adj.P.Val<0.05)
up.reg1_lumA_response <- up.reg1_lumA_response%>%filter(up.reg1_lumA_response$logFC > 1.1)
length(names_upreg_lumA_rd <- rownames(up.reg1_lumA_response)) #4 significant probes
#upregulated probes in lumA(pCR)
up.reg2_lumA_response <- top1_lumA_response%>%filter(top1_lumA_response$adj.P.Val<0.05)
up.reg2_lumA_response <- up.reg2_lumA_response%>%filter(up.reg2_lumA_response$logFC < -1.1)
length(names_upreg_lumA_pcr <- rownames(up.reg2_lumA_response)) #5 significant probes

#lum A(RD vs pCR): change fold change threshold to 1.0
#upregulated probes in lumA(RD) 
up.reg1_lumA_response <- top1_lumA_response%>%filter(top1_lumA_response$adj.P.Val<0.05)
up.reg1_lumA_response <- up.reg1_lumA_response%>%filter(up.reg1_lumA_response$logFC > 1)
length(names_upreg_lumA_rd <- rownames(up.reg1_lumA_response)) #9 significant probes
#upregulated probes in lumA(pCR)
up.reg2_lumA_response <- top1_lumA_response%>%filter(top1_lumA_response$adj.P.Val<0.05)
up.reg2_lumA_response <- up.reg2_lumA_response%>%filter(up.reg2_lumA_response$logFC < -1)
length(names_upreg_lumA_pcr <- rownames(up.reg2_lumA_response)) #7 significant probes

# TNBC(RD vs pCR)
# TNBC(RD vs pCR): change fold change threshold to 1.2
#upregulated probes in TNBC(RD)
up.reg1_tnbc_respond <- top1_tnbc_response%>%filter(top1_tnbc_response$ adj.P.Val<0.05)
up.reg1_tnbc_respond <- up.reg1_tnbc_respond%>%filter(up.reg1_tnbc_respond$logFC > 1.2)
names_upreg_tnbc_rd <- rownames(up.reg1_tnbc_respond)
length(names_upreg_tnbc_rd) #0 significant probe
#upregulated probes in TNBC(pCR)
up.reg2_tnbc_respond <- top1_tnbc_response%>%filter(top1_tnbc_response$ adj.P.Val<0.05)
up.reg2_tnbc_respond <- up.reg2_tnbc_respond%>%filter(up.reg2_tnbc_respond$logFC < -1.2)
names_upreg_tnbc_pcr <- rownames(up.reg2_tnbc_respond)
length(names_upreg_tnbc_pcr) #0 significant probe

# TNBC(RD vs pCR): change fold change threshold to 1.1
#upregulated probes in TNBC(RD)
up.reg1_tnbc_respond <- top1_tnbc_response%>%filter(top1_tnbc_response$ adj.P.Val<0.05)
up.reg1_tnbc_respond <- up.reg1_tnbc_respond%>%filter(up.reg1_tnbc_respond$logFC > 1.1)
names_upreg_tnbc_rd <- rownames(up.reg1_tnbc_respond)
length(names_upreg_tnbc_rd)  #0 significant probe
#upregulated probes in TNBC(pCR)
up.reg2_tnbc_respond <- top1_tnbc_response%>%filter(top1_tnbc_response$ adj.P.Val<0.05)
up.reg2_tnbc_respond <- up.reg2_tnbc_respond%>%filter(up.reg2_tnbc_respond$logFC < -1.1)
names_upreg_tnbc_pcr <- rownames(up.reg2_tnbc_respond)
length(names_upreg_tnbc_pcr)  #0 significant probe

# TNBC(RD vs pCR): change fold change threshold to 1.0
#upregulated probes in TNBC(RD)
up.reg1_tnbc_respond <- top1_tnbc_response%>%filter(top1_tnbc_response$ adj.P.Val<0.05)
up.reg1_tnbc_respond <- up.reg1_tnbc_respond%>%filter(up.reg1_tnbc_respond$logFC > 1.0)
names_upreg_tnbc_rd <- rownames(up.reg1_tnbc_respond)
length(names_upreg_tnbc_rd)  #0 significant probe
#upregulated probes in TNBC(pCR)
up.reg2_tnbc_respond <- top1_tnbc_response%>%filter(top1_tnbc_response$ adj.P.Val<0.05)
up.reg2_tnbc_respond <- up.reg2_tnbc_respond%>%filter(up.reg2_tnbc_respond$logFC < -1.0)
names_upreg_tnbc_pcr <- rownames(up.reg2_tnbc_respond)
length(names_upreg_tnbc_pcr)  #0 significant probe

#nonTNBC(RD) vs nonTNBC(pCR)
#nonTNBC(RD) vs nonTNBC(pCR): change fold change threshold to 1.2
#upregulated probe_id in non-TNBC(RD)
up.reg1_nontnbc_respond  <- top1_nontnbc_response%>%filter(top1_nontnbc_response $logFC < -1.2)
up.reg1_nontnbc_respond  <-  up.reg1_nontnbc_respond%>%filter(up.reg1_nontnbc_respond$adj.P.Val<0.05)
names_upreg_rd_nonTNBC <- rownames(up.reg1_nontnbc_respond )
length(names_upreg_rd_nonTNBC) #5
#upregulated probe_id in non-TNBC(pCR)
up.reg2_nontnbc_respond  <- top1_nontnbc_response %>%filter(top1_nontnbc_response $logFC > 1.2)
up.reg2_nontnbc_respond  <- up.reg2_nontnbc_respond %>%filter(up.reg2_nontnbc_respond$adj.P.Val<0.05)
names_upreg_pcr_nonTNBC <- rownames(up.reg2_nontnbc_respond )
length(names_upreg_pcr_nonTNBC) #4

#nonTNBC(RD) vs nonTNBC(pCR): change fold change threshold to 1.1
#upregulated probe_id in non-TNBC(RD)
up.reg1_nontnbc_respond  <- top1_nontnbc_response %>%filter(top1_nontnbc_response $logFC < -1.1)
up.reg1_nontnbc_respond  <-  up.reg1_nontnbc_respond%>%filter(up.reg1_nontnbc_respond$adj.P.Val<0.05)
names_upreg_rd_nonTNBC <- rownames(up.reg1_nontnbc_respond )
length(names_upreg_rd_nonTNBC) #9
#upregulated probe_id in non-TNBC(pCR)
up.reg2_nontnbc_respond  <- top1_nontnbc_response %>%filter(top1_nontnbc_response $logFC > 1.1)
up.reg2_nontnbc_respond  <- up.reg2_nontnbc_respond %>%filter(up.reg2_nontnbc_respond$adj.P.Val<0.05)
names_upreg_pcr_nonTNBC <- rownames(up.reg2_nontnbc_respond )
length(names_upreg_pcr_nonTNBC) #10

#nonTNBC(RD) vs nonTNBC(pCR): change fold change threshold to 1.0
#upregulated probe_id in non-TNBC(RD)
up.reg1_nontnbc_respond  <- top1_nontnbc_response %>%filter(top1_nontnbc_response $logFC < -1.0)
up.reg1_nontnbc_respond  <-  up.reg1_nontnbc_respond%>%filter(up.reg1_nontnbc_respond$adj.P.Val<0.05)
names_upreg_rd_nonTNBC <- rownames(up.reg1_nontnbc_respond )
length(names_upreg_rd_nonTNBC) #15
#upregulated probe_id in non-TNBC(pCR)
up.reg2_nontnbc_respond  <- top1_nontnbc_response %>%filter(top1_nontnbc_response $logFC > 1.0)
up.reg2_nontnbc_respond  <- up.reg2_nontnbc_respond %>%filter(up.reg2_nontnbc_respond$adj.P.Val<0.05)
names_upreg_pcr_nonTNBC <- rownames(up.reg2_nontnbc_respond )
length(names_upreg_pcr_nonTNBC) #15

#TNBC(RD) vs nonTNBC(RD) 
#TNBC(RD) vs nonTNBC(RD): change fold change threshold to 1.2
#upregulated probe_id in non-TNBC(RD)
up.reg1_rd_tnbc_nontnbc  <- top1_rd_tnbc_nontnbc%>%filter(top1_rd_tnbc_nontnbc$logFC < -1.2)
up.reg1_rd_tnbc_nontnbc <- up.reg1_rd_tnbc_nontnbc%>%filter(up.reg1_rd_tnbc_nontnbc$ adj.P.Val<0.05)
names_upreg_rd_nontnbc<- rownames(up.reg1_rd_tnbc_nontnbc)
length(names_upreg_rd_nontnbc) #69
#upregulated probe_id in TNBC(RD)
up.reg2_rd_tnbc_nontnbc  <- top1_rd_tnbc_nontnbc %>%filter(top1_rd_tnbc_nontnbc$logFC > 1.2)
up.reg2_rd_tnbc_nontnbc  <- up.reg2_rd_tnbc_nontnbc%>%filter(up.reg2_rd_tnbc_nontnbc$ adj.P.Val<0.05) > 1)
names_upreg_rd_tnbc <- rownames(up.reg2_rd_tnbc_nontnbc)
length(names_upreg_rd_tnbc)  #25

#TNBC(RD) vs nonTNBC(RD): change fold change threshold to 1.1
#upregulated probe_id in non-TNBC(RD)
up.reg1_rd_tnbc_nontnbc  <- top1_rd_tnbc_nontnbc%>%filter(top1_rd_tnbc_nontnbc$logFC < -1.1)
up.reg1_rd_tnbc_nontnbc <- up.reg1_rd_tnbc_nontnbc%>%filter(up.reg1_rd_tnbc_nontnbc$ adj.P.Val<0.05) > 1)
names_upreg_rd_nontnbc<- rownames(up.reg1_rd_tnbc_nontnbc)
length(names_upreg_rd_nontnbc) #94
#upregulated probe_id in TNBC(RD)
up.reg2_rd_tnbc_nontnbc  <- top1_rd_tnbc_nontnbc %>%filter(top1_rd_tnbc_nontnbc$logFC > 1.1)
up.reg2_rd_tnbc_nontnbc  <- up.reg2_rd_tnbc_nontnbc%>%filter(up.reg2_rd_tnbc_nontnbc$ adj.P.Val<0.05) > 1)
names_upreg_rd_tnbc <- rownames(up.reg2_rd_tnbc_nontnbc)
length(names_upreg_rd_tnbc)  #36

#TNBC(RD) vs nonTNBC(RD): change fold change threshold to 1.0
#upregulated probe_id in non-TNBC(RD)
up.reg1_rd_tnbc_nontnbc  <- top1_rd_tnbc_nontnbc%>%filter(top1_rd_tnbc_nontnbc$logFC < -1.0)
up.reg1_rd_tnbc_nontnbc <- up.reg1_rd_tnbc_nontnbc%>%filter(up.reg1_rd_tnbc_nontnbc$ adj.P.Val<0.05) > 1)
names_upreg_rd_nontnbc<- rownames(up.reg1_rd_tnbc_nontnbc)
length(names_upreg_rd_nontnbc) #112
#upregulated probe_id in TNBC(RD)
up.reg2_rd_tnbc_nontnbc  <- top1_rd_tnbc_nontnbc %>%filter(top1_rd_tnbc_nontnbc$logFC > 1.0)
up.reg2_rd_tnbc_nontnbc  <- up.reg2_rd_tnbc_nontnbc%>%filter(up.reg2_rd_tnbc_nontnbc$ adj.P.Val<0.05) > 1)
names_upreg_rd_tnbc <- rownames(up.reg2_rd_tnbc_nontnbc)
length(names_upreg_rd_tnbc)  #55

#TNBC(pCR) vs nonTNBC(pCR) 
#TNBC(pCR) vs nonTNBC(pCR): change fold change threshold to 1.2
#upregulated probe_id in non-TNBC(pCR)
up.reg1_pcr_tnbc_nontnbc <- top1_pcr_tnbc_nontnbc%>%filter(top1_pcr_tnbc_nontnbc$logFC < -1.2)
up.reg1_pcr_tnbc_nontnbc <- up.reg1_pcr_tnbc_nontnbc%>%filter(up.reg1_pcr_tnbc_nontnbc$ adj.P.Val<0.05)
names_pcr_nontnbc  <- rownames(up.reg1_pcr_tnbc_nontnbc)
length(names_pcr_nontnbc) #22
#upregulated probe_id inTNBC(pCR)
up.reg2_pcr_tnbc_nontnbc <- top1_pcr_tnbc_nontnbc%>%filter(top1_pcr_tnbc_nontnbc$logFC > 1.2)
up.reg2_pcr_tnbc_nontnbc <- up.reg2_pcr_tnbc_nontnbc%>%filter(up.reg2_pcr_tnbc_nontnbc$ adj.P.Val<0.05)
names_pcr_tnbc <- rownames(up.reg2_pcr_tnbc_nontnbc)
length(names_pcr_tnbc) #15

#TNBC(pCR) vs nonTNBC(pCR): change fold change threshold to 1.1
#upregulated probe_id in non-TNBC(pCR)
up.reg1_pcr_tnbc_nontnbc <- top1_pcr_tnbc_nontnbc%>%filter(top1_pcr_tnbc_nontnbc$logFC < -1.1)
up.reg1_pcr_tnbc_nontnbc <- up.reg1_pcr_tnbc_nontnbc%>%filter(up.reg1_pcr_tnbc_nontnbc$ adj.P.Val<0.05)
names_pcr_nontnbc  <- rownames(up.reg1_pcr_tnbc_nontnbc)
length(names_pcr_nontnbc) #28
#upregulated probe_id inTNBC(pCR)
up.reg2_pcr_tnbc_nontnbc <- top1_pcr_tnbc_nontnbc%>%filter(top1_pcr_tnbc_nontnbc$logFC > 1.1)
up.reg2_pcr_tnbc_nontnbc <- up.reg2_pcr_tnbc_nontnbc%>%filter(up.reg2_pcr_tnbc_nontnbc$ adj.P.Val<0.05)
names_pcr_tnbc <- rownames(up.reg2_pcr_tnbc_nontnbc)
length(names_pcr_tnbc) #16

#TNBC(pCR) vs nonTNBC(pCR): change fold change threshold to 1.0
#upregulated probe_id in non-TNBC(pCR)
up.reg1_pcr_tnbc_nontnbc <- top1_pcr_tnbc_nontnbc%>%filter(top1_pcr_tnbc_nontnbc$logFC < -1.0)
up.reg1_pcr_tnbc_nontnbc <- up.reg1_pcr_tnbc_nontnbc%>%filter(up.reg1_pcr_tnbc_nontnbc$ adj.P.Val<0.05)
names_pcr_nontnbc  <- rownames(up.reg1_pcr_tnbc_nontnbc)
length(names_pcr_nontnbc) #36
#upregulated probe_id inTNBC(pCR)
up.reg2_pcr_tnbc_nontnbc <- top1_pcr_tnbc_nontnbc%>%filter(top1_pcr_tnbc_nontnbc$logFC > 1.0)
up.reg2_pcr_tnbc_nontnbc <- up.reg2_pcr_tnbc_nontnbc%>%filter(up.reg2_pcr_tnbc_nontnbc$ adj.P.Val<0.05)
names_pcr_tnbc <- rownames(up.reg2_pcr_tnbc_nontnbc)
length(names_pcr_tnbc) #25

# TNBC  vs nonTNBC
# TNBC  vs nonTNBC: change fold change threshold to 1.2
#upregulated probes in non-TNBC
up.reg1_tnbc_nontnbc  <- top1_tnbc_nontnbc%>%filter(top1_tnbc_nontnbc$logFC < -1.2) up.reg1_tnbc_nontnbc  <- up.reg1_tnbc_nontnbc %>%filter(up.reg1_tnbc_nontnbc$ adj.P.Val <0.05)
names_upreg_nontnbc <- rownames(up.reg1_tnbc_nontnbc)
length(names_upreg_nontnbc) #31
#upregulated probes inTNBC
up.reg2_tnbc_nontnbc <- top1_tnbc_nontnbc%>%filter(top1_tnbc_nontnbc$logFC > 1.2)
up.reg2_tnbc_nontnbc <- up.reg2_tnbc_nontnbc%>%filter(up.reg2_tnbc_nontnbc$ adj.P.Val <0.05)
names_upreg_tnbc <- rownames(up.reg2_tnbc_nontnbc)
length(names_upreg_tnbc) #26

# TNBC  vs nonTNBC: change fold change threshold to 1.1
#upregulated probes in non-TNBC
up.reg1_tnbc_nontnbc  <- top1_tnbc_nontnbc%>%filter(top1_tnbc_nontnbc$logFC < -1.1)
up.reg1_tnbc_nontnbc  <- up.reg1_tnbc_nontnbc %>%filter(up.reg1_tnbc_nontnbc$ adj.P.Val <0.05)
names_upreg_nontnbc <- rownames(up.reg1_tnbc_nontnbc)
length(names_upreg_nontnbc) #83
#upregulated probes inTNBC
up.reg2_tnbc_nontnbc <- top1_tnbc_nontnbc%>%filter(top1_tnbc_nontnbc$logFC > 1.1)
up.reg2_tnbc_nontnbc <- up.reg2_tnbc_nontnbc%>%filter(up.reg2_tnbc_nontnbc$ adj.P.Val <0.05)
names_upreg_tnbc <- rownames(up.reg2_tnbc_nontnbc)
length(names_upreg_tnbc) #37

# TNBC  vs nonTNBC: change fold change threshold to 1.0
#upregulated probes in non-TNBC
up.reg1_tnbc_nontnbc  <- top1_tnbc_nontnbc%>%filter(top1_tnbc_nontnbc$logFC < -1)
up.reg1_tnbc_nontnbc  <- up.reg1_tnbc_nontnbc %>%filter(up.reg1_tnbc_nontnbc$ adj.P.Val <0.05)
names_upreg_nontnbc <- rownames(up.reg1_tnbc_nontnbc)
length(names_upreg_nontnbc) #104
#upregulated probes inTNBC
up.reg2_tnbc_nontnbc <- top1_tnbc_nontnbc%>%filter(top1_tnbc_nontnbc$logFC > 1)
up.reg2_tnbc_nontnbc <- up.reg2_tnbc_nontnbc%>%filter(up.reg2_tnbc_nontnbc$ adj.P.Val <0.05)
names_upreg_tnbc <- rownames(up.reg2_tnbc_nontnbc)
length(names_upreg_tnbc) #56

#-------------------------------------------------- Probes annotation -----------------------------------------------
# TNBC vs non-TNBC (40)
# extracting significant probe ids
dd_tnbc_nontnbc <- subset(top1_tnbc_nontnbc, adj.P.Val<0.05 & abs(logFC)>1.5) # subset extreme values
str(dd_tnbc_nontnbc) # for TNBC vs non-TNBC 40 obs. of  6 variables 
out_tnbc_nontnbc <- data.frame(xxhuex[rownames(dd_tnbc_nontnbc)] ) # search for rownames example

gene_labs_tnbc_nontnbc <- NULL # gene ids
for (i in 1:length(rownames(dd_tnbc_nontnbc))){ gene_labs_tnbc_nontnbc <- c(gene_labs_tnbc_nontnbc, out_tnbc_nontnbc[, paste0("X", rownames(dd_tnbc_nontnbc)[i])]) }
dd_tnbc_nontnbc <- cbind(gene_labs_tnbc_nontnbc, dd_tnbc_nontnbc) # get probe and gene IDs
dd_tbl_tnbc_nontnbc <- as.data.table(dd_tnbc_nontnbc) # print as table 
df_tnbc_nontnbc <- data.frame(dd_tnbc_nontnbc)
df_tnbc_nontnbc$DGE <- 'TNBC vs nonTNBC'
table_tnbc_nontnbc_1.5 <- subset(df_tnbc_nontnbc, select=c('gene_labs_tnbc_nontnbc','logFC','adj.P.Val'))
info <- select(hgu133a.db, rownames(table_tnbc_nontnbc_1.5), c("SYMBOL","GENENAME"))
info <- subset(info,select=c("SYMBOL","GENENAME"))
table_tnbc_nontnbc_1.5 <- cbind(table_tnbc_nontnbc_1.5,info)
write.csv(table_tnbc_nontnbc_1.5, 'table_tnbc_nontnbc_probe_ann.csv')

# TNBC(RD) vs nonTNBC(RD) (43)
# extracting significant probe ids
dd_rd_tnbc_nontnbc <- subset(top1_rd_tnbc_nontnbc, adj.P.Val<0.05 & abs(logFC)>1.5) # subset extreme values
str(dd_rd_tnbc_nontnbc) # for (RD) TNBC vs non-TNBC 43 obs. of  6 variables
out_rd_tnbc_nontnbc <- data.frame(xxhuex[rownames(dd_rd_tnbc_nontnbc)] ) # search for rownames example

gene_labs_rd_tnbc_nontnbc <- NULL # gene ids
for (i in 1:length(rownames(dd_rd_tnbc_nontnbc))){ gene_labs_rd_tnbc_nontnbc <- c(gene_labs_rd_tnbc_nontnbc, out_rd_tnbc_nontnbc[, paste0("X", rownames(dd_rd_tnbc_nontnbc)[i])]) }
dd_rd_tnbc_nontnbc <- cbind(gene_labs_rd_tnbc_nontnbc, dd_rd_tnbc_nontnbc) # get probe and gene IDs
dd_tbl_rd_tnbc_nontnbc <- as.data.table(dd_rd_tnbc_nontnbc) # print as table 
df_rd_tnbc_nontnbc <- data.frame(dd_rd_tnbc_nontnbc)
df_rd_tnbc_nontnbc$DGE <- 'TNBC-RD vs nonTNBC-RD'
table_rd_tnbc_nontnbc_1.5 <- subset(df_rd_tnbc_nontnbc, select=c('gene_labs_rd_tnbc_nontnbc','logFC','adj.P.Val'))
info <- select(hgu133a.db, rownames(table_rd_tnbc_nontnbc_1.5), c("SYMBOL","GENENAME"))
table_rd_tnbc_nontnbc_1.5 <- cbind(table_rd_tnbc_nontnbc_1.5,info)
write.csv(table_rd_tnbc_nontnbc_1.5, 'table_rd_tnbc_nontnbc_probe_ann .csv')

# TNBC(pCR) vs nonTNBC(pCR) (15)
# extracting significant probe ids
dd_pcr_tnbc_nontnbc <- subset(top1_pcr_tnbc_nontnbc, adj.P.Val<0.05 & abs(logFC)>1.5) # subset extreme values
str(dd_pcr_tnbc_nontnbc) # for (pCR) TNBC vs non-TNBC 15 obs. of  6 variables
out_pcr_tnbc_nontnbc <- data.frame(xxhuex[rownames(dd_pcr_tnbc_nontnbc)] ) # search for rownames example

gene_labs_pcr_tnbc_nontnbc <- NULL # gene ids
for (i in 1:length(rownames(dd_pcr_tnbc_nontnbc))){ gene_labs_pcr_tnbc_nontnbc <- c(gene_labs_pcr_tnbc_nontnbc, out_pcr_tnbc_nontnbc[, paste0("X", rownames(dd_pcr_tnbc_nontnbc)[i])]) }
dd_pcr_tnbc_nontnbc <- cbind(gene_labs_pcr_tnbc_nontnbc, dd_pcr_tnbc_nontnbc) # get probe and gene IDs
dd_tbl_pcr_tnbc_nontnbc <- as.data.table(dd_pcr_tnbc_nontnbc) # print as table 
df_pcr_tnbc_nontnbc <- data.frame(dd_pcr_tnbc_nontnbc)
df_pcr_tnbc_nontnbc$DGE <- 'TNBC-pCR vs nonTNBC-pCR'
table_pcr_tnbc_nontnbc_1.5 <- subset(df_pcr_tnbc_nontnbc, select=c('gene_labs_pcr_tnbc_nontnbc','logFC','adj.P.Val'))
info <- select(hgu133a.db, rownames(table_pcr_tnbc_nontnbc_1.5), c("SYMBOL","GENENAME"))
table_pcr_tnbc_nontnbc_1.5 <- cbind(table_pcr_tnbc_nontnbc_1.5,info)
write.csv(table_pcr_tnbc_nontnbc_1.5, 'table_pcr_tnbc_nontnbc_probe_ann .csv')

# non-TNBC(RD) vs non-TNBC(pCR) (1)
# extracting significant probe ids
dd_nontnbc_response <- subset(top1_nontnbc_response, adj.P.Val<0.05 & abs(logFC)>1.5) # subset extreme values 
str(dd_nontnbc_response) # for non-TNBC(RD) vs(pCR) 1 obs. of  6 variables
out_nontnbc_response <- data.frame(xxhuex[rownames(dd_nontnbc_response)] ) # search for rownames example

gene_labs_nontnbc_respond <- NULL # gene ids
for (i in 1:length(rownames(dd_nontnbc_response))){ gene_labs_nontnbc_respond <- c(gene_labs_nontnbc_respond, out_nontnbc_response[, paste0("X", rownames(dd_nontnbc_response)[i])]) }
dd_nontnbc_response <- cbind(gene_labs_nontnbc_respond, dd_nontnbc_response) # get probe and gene IDs
dd_tbl_nontnbc_respond <- as.data.table(dd_nontnbc_response) # print as table 
df_nontnbc_response <- data.frame(dd_nontnbc_response)
df_nontnbc_response$DGE <- 'nonTNBC-RD vs nonTNBC-pCR'
table_nontnbc_1.5 <- subset(df_nontnbc_response, select=c('gene_labs_nontnbc_respond','logFC','adj.P.Val'))
info <- select(hgu133a.db, rownames(table_nontnbc_1.5), c("SYMBOL","GENENAME"))
table_nontnbc_1.5 <- cbind(table_nontnbc_1.5,info)
write.csv(table_nontnbc_1.5, 'table_nontnbc_probe_ann.csv')

# lumA(RD) vs lumA(pCR) with logFC>1.2 (5)
# extracting significant probe ids
dd_lumA_response <- subset(top1_lumA_response, adj.P.Val<0.05 & abs(logFC)>1.2) # subset extreme values 
str(dd_lumA_response) # for lumA(RD) vs(pCR) 4 obs. of  6 variables
out_lumA_response <- data.frame(xxhuex[rownames(dd_lumA_response)] ) # search for rownames example

gene_labs_lumA_response <- NULL # gene ids
for (i in 1:length(rownames(dd_lumA_response))){ gene_labs_lumA_response <- c(gene_labs_lumA_response, out_lumA_response[, paste0("X", rownames(dd_lumA_response)[i])]) }
dd_lumA_response <- cbind(gene_labs_lumA_response, dd_lumA_response) # get probe and gene IDs
dd_tbl_lumA_response <- as.data.table(dd_lumA_response) # print as table 
df_lumA_response <- data.frame(dd_lumA_response)
df_lumA_response$DGE <- 'lumA-RD vs lumA-pCR'
table_lumA_1.2 <- subset(df_lumA_response, select=c('gene_labs_lumA_response','logFC','adj.P.Val'))
info <- select(hgu133a.db, rownames(table_lumA_1.2), c("SYMBOL","GENENAME"))
table_lumA_1.2 <- cbind(table_lumA_1.2,info)

#-------------------------------------------------- venn diagram -------------------------------------------------- 
#showing the overlapping of significant probes of each comparison
x <- list(
TNBC_vs_NonTNBC = rownames(dd_tnbc_nontnbc), 
TNBC.RD_vs_NonTNBC.RD = rownames(dd_rd_tnbc_nontnbc), 
TNBC.pCR_vs_NonTNBC.pCR = rownames(dd_pcr_tnbc_nontnbc),  
NonTNBC.RD_vs_pCR = rownames(dd_nontnbc_response))
png('Venn.diagram_probes.png',height= 1000, width= 1900, units='px') 
ggvenn(x, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"), stroke_size = 1, set_name = FALSE,show_percentage = FALSE, text_size = 8)
dev.off()

#-------------------------------------------------- annotate all probes -------------------------------------------------- 
# link all probes to corresponding gene names (to check Ca2+ gene expression & identify the top differentially expressed probes)
# from DGE TNBC(RD) vs TNBC(pCR)
dd_tnbc_respond <- top1_tnbc_response
str(dd_tnbc_respond) 
acc.no <-  data.frame(unlist(xxhuex[rownames(dd_tnbc_respond)]))

gene_labs_tnbc_respond <- acc.no[,1] # gene ids
dd_tnbc_respond <- cbind(gene_labs_tnbc_respond, dd_tnbc_respond) # get probe and gene IDs

table_tnbc_respond <- subset(dd_tnbc_respond, select=c('gene_labs_tnbc_respond','logFC','adj.P.Val'))
info <- select(hgu133a.db, rownames(table_tnbc_respond), c("SYMBOL", "GENENAME")) # 24464 rows

info_temp <- data.frame(matrix( rep(0, length(info$PROBEID)), nrow=length(info$PROBEID), ncol=2))
colnames(info_temp) <- c("SYMBOL", "GENENAME")
str(info_temp)
for (i in 1:length( rownames(table_tnbc_respond) ) ) { # all 22283 probes
	rownames(info_temp)[i] <- rownames(table_tnbc_respond)[i] 
	info_temp[i,1] <-  select(hgu133a.db, rownames(table_tnbc_respond)[i], c("SYMBOL"))[2]
	info_temp[i,2] <-  select(hgu133a.db, rownames(table_tnbc_respond)[i], c("GENENAME"))[2]
} # end for 
head(info_temp) # “1:1 mapping” message not important
str(info_temp) # 24464 rows for each unique gene

link probes to gene name
gene.name <- subset(info_temp, rownames(info_temp) == rownames(table_tnbc_respond), select=c('SYMBOL','GENENAME'))
table_tnbc_respond <- cbind(table_tnbc_respond,gene.name)
write.csv(table_tnbc_respond,'tnbc_probes_annotation.csv')

#Annotate all probes in TNBC vs Non-TNBC
#extracting probe ids
dd_tnbc_nontnbc <- top1_tnbc_nontnbc
str(dd_tnbc_nontnbc) 
acc.no_tnbc_nontnbc <-  data.frame(unlist(xxhuex[rownames(dd_tnbc_nontnbc)])) # search for rownames 
gene_labs_tnbc_nontnbc <- acc.no_tnbc_nontnbc[,1] # gene ids
dd_tnbc_nontnbc <- cbind(gene_labs_tnbc_nontnbc, dd_tnbc_nontnbc) # get probe and gene IDs
df_tnbc_nontnbc$DGE <- 'TNBC vs nonTNBC'
table_tnbc_nontnbc_1.5 <- subset(dd_tnbc_nontnbc, select=c('gene_labs_tnbc_nontnbc','logFC','adj.P.Val'))
info <- select(hgu133a.db, rownames(table_tnbc_nontnbc_1.5), c("SYMBOL", "GENENAME"))

info_temp <- data.frame(matrix( rep(0, length(info$PROBEID)), nrow=length(info$PROBEID), ncol=2))
colnames(info_temp) <- c("SYMBOL", "GENENAME")
str(info_temp)
for (i in 1:length( rownames(table_tnbc_nontnbc_1.5 ) ) ) { # all 22283 probes
	rownames(info_temp)[i] <- rownames(table_tnbc_nontnbc_1.5 )[i] 
	info_temp[i,1] <-  select(hgu133a.db, rownames(table_tnbc_nontnbc_1.5 )[i], c("SYMBOL"))[2]
	info_temp[i,2] <-  select(hgu133a.db, rownames(table_tnbc_nontnbc_1.5 )[i], c("GENENAME"))[2]
} # end for 
head(info_temp) # “1:1 mapping” message not important
str(info_temp) # 24464 rows for each unique gene

link probes to gene name
gene.name <- subset(info_temp, rownames(info_temp) == rownames(table_tnbc_nontnbc_1.5), select=c('SYMBOL','GENENAME'))
table_tnbc_nontnbc_1.5 <- cbind(table_tnbc_nontnbc_1.5,gene.name)
write.csv(table_tnbc_nontnbc_1.5,'tnbc_nontnbc_1.5_probes_annotation_all.csv')

#Annotate all probes in TNBC(RD) vs Non-TNBC(RD)
dd_rd_tnbc_nontnbc <-top1_rd_tnbc_nontnbc
str(dd_rd_tnbc_nontnbc)
acc.no_rd_tnbc_nontnbc <-  data.frame(unlist(xxhuex[rownames(dd_rd_tnbc_nontnbc)])) # search for rownames 
gene_labs_rd_tnbc_nontnbc <- acc.no_rd_tnbc_nontnbc[,1] # gene ids
dd_rd_tnbc_nontnbc <- cbind(gene_labs_rd_tnbc_nontnbc, dd_rd_tnbc_nontnbc) # get probe and gene IDs
#df_rd_tnbc_nontnbc$DGE <- 'TNBC vs nonTNBC'
table_rd_tnbc_nontnbc_1.5 <- subset(dd_rd_tnbc_nontnbc, select=c('gene_labs_rd_tnbc_nontnbc','logFC','adj.P.Val'))
info <- select(hgu133a.db, rownames(table_rd_tnbc_nontnbc_1.5), c("SYMBOL", "GENENAME"))

info_temp <- data.frame(matrix( rep(0, length(info$PROBEID)), nrow=length(info$PROBEID), ncol=2))
colnames(info_temp) <- c("SYMBOL", "GENENAME")
str(info_temp)
for (i in 1:length( rownames(table_rd_tnbc_nontnbc_1.5 ) ) ) { # all 22283 probes
	rownames(info_temp)[i] <- rownames(table_rd_tnbc_nontnbc_1.5 )[i] 
	info_temp[i,1] <-  select(hgu133a.db, rownames(table_rd_tnbc_nontnbc_1.5 )[i], c("SYMBOL"))[2]
	info_temp[i,2] <-  select(hgu133a.db, rownames(table_rd_tnbc_nontnbc_1.5 )[i], c("GENENAME"))[2]
} # end for 
head(info_temp) # “1:1 mapping” message not important
str(info_temp) # 24464 rows for each unique gene

link probes to gene name
gene.name <- subset(info_temp, rownames(info_temp) == rownames(table_rd_tnbc_nontnbc_1.5), select=c('SYMBOL','GENENAME'))
table_rd_tnbc_nontnbc_1.5 <- cbind(table_rd_tnbc_nontnbc_1.5,gene.name)
write.csv(table_rd_tnbc_nontnbc_1.5,'table_rd_tnbc_nontnbc_1.5_probes_annotation_all2.csv')

# Annotate all probes in TNBC(pCR) vs Non-TNBC(pCR)
dd_pcr_tnbc_nontnbc <-top1_pcr_tnbc_nontnbc
str(dd_pcr_tnbc_nontnbc)
acc.no_pcr_tnbc_nontnbc <-  data.frame(unlist(xxhuex[rownames(dd_pcr_tnbc_nontnbc)])) # search for rownames 
gene_labs_pcr_tnbc_nontnbc <- acc.no_pcr_tnbc_nontnbc[,1] # gene ids
dd_pcr_tnbc_nontnbc <- cbind(gene_labs_pcr_tnbc_nontnbc, dd_pcr_tnbc_nontnbc) # get probe and gene IDs
#df_pcr_tnbc_nontnbc$DGE <- 'TNBC vs nonTNBC'
table_pcr_tnbc_nontnbc_1.5 <- subset(dd_pcr_tnbc_nontnbc, select=c('gene_labs_pcr_tnbc_nontnbc','logFC','adj.P.Val'))
info <- select(hgu133a.db, rownames(table_pcr_tnbc_nontnbc_1.5), c("SYMBOL", "GENENAME"))

info_temp <- data.frame(matrix( rep(0, length(info$PROBEID)), nrow=length(info$PROBEID), ncol=2))
colnames(info_temp) <- c("SYMBOL", "GENENAME")
str(info_temp)
for (i in 1:length( rownames(table_pcr_tnbc_nontnbc_1.5 ) ) ) { # all 22283 probes
	rownames(info_temp)[i] <- rownames(table_pcr_tnbc_nontnbc_1.5 )[i] 
	info_temp[i,1] <-  select(hgu133a.db, rownames(table_pcr_tnbc_nontnbc_1.5 )[i], c("SYMBOL"))[2]
	info_temp[i,2] <-  select(hgu133a.db, rownames(table_pcr_tnbc_nontnbc_1.5 )[i], c("GENENAME"))[2]
} # end for 
head(info_temp) # “1:1 mapping” message not important
str(info_temp) # 24464 rows for each unique gene

link probes to gene name
gene.name <- subset(info_temp, rownames(info_temp) == rownames(table_pcr_tnbc_nontnbc_1.5), select=c('SYMBOL','GENENAME'))
table_pcr_tnbc_nontnbc_1.5 <- cbind(table_pcr_tnbc_nontnbc_1.5,gene.name)
write.csv(table_pcr_tnbc_nontnbc_1.5,'table_pcr_tnbc_nontnbc_1.5_probes_annotation_all.csv')


#----------------------------------------------- DGE analysis (based on chemosensitivity) ------------------------------------------
# sort by sensitivity prediction
sen_tnbc <-  tnbc%>%filter(chemosensitivity.prediction == 'Rx Sensitive') #43
insen_tnbc <-  tnbc%>%filter(chemosensitivity.prediction == 'Rx Insensitive') #135

sen_nontnbc <-  nontnbc%>%filter(chemosensitivity.prediction == 'Rx Sensitive') #126
insen_nontnbc <-  nontnbc%>%filter(chemosensitivity.prediction == 'Rx Insensitive') #204

sensitivity<- rbind(sen_tnbc,insen_tnbc,sen_nontnbc,insen_nontnbc)
rawdata_sensitivity  <- read.celfiles(sensitivity$celfiles) 
normdata_sensitivity<-rma(rawdata_sensitivity)

# TNBC
mypch2_sensitivity <- c(rep(1,43),rep(2,135),rep(3,126),rep(4,204))	
design_sensitivity  <- model.matrix(~-1+factor(mypch2_sensitivity))
colnames(design_sensitivity)<- c('s.t', 'i.t','s.n','i.n')  #  group labels

top0_sensitivity_tnbc <- eBayes(contrasts.fit(lmFit(normdata_sensitivity, design_sensitivity),makeContrasts(s.t - i.t, levels=design_sensitivity)))
top1_sensitivity_tnbc <- topTable(top0_sensitivity_tnbc , coef='s.t - i.t',adjust='BH',number=22283,sort.by='P') 

#upregulated probe_id in sensitive group
up.reg1_sensitivity_tnbc <- top1_sensitivity_tnbc%>%filter(top1_sensitivity_tnbc$logFC > 1.5)
up.reg1_sensitivity_tnbc  <-  up.reg1_sensitivity_tnbc%>%filter(up.reg1_sensitivity_tnbc$adj.P.Val<0.05)
names_sensitive_tnbc <- rownames(up.reg1_sensitivity_tnbc)
length(names_sensitive_tnbc) #0
#upregulated probe_id in insensitive group
up.reg2_sensitivity_tnbc  <- top1_sensitivity_tnbc%>%filter(top1_sensitivity_tnbc$logFC < -1.5)
up.reg2_sensitivity_tnbc <- up.reg2_sensitivity_tnbc %>%filter(up.reg2_sensitivity_tnbc$adj.P.Val<0.05)
names_insensitive_tnbc <- rownames(up.reg2_sensitivity_tnbc)
length(names_insensitive_tnbc) #0

png(file="volcano-chemosensitivity_tnbc_sensitive-insensitive.png",units="px")
plot(top1_sensitivity_tnbc$logFC,-log10(top1_sensitivity_tnbc$adj.P.Val), main='Differential Gene Expression \n  between TNBC Chemo-sensitive and Chemo-insensitive group', xlab='log2FC(Sensitive/Insensitive)',ylab='-log10(BH P value)', cex=0.5, cex.main=1.2,cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.5,2.5,label='Upregulated in \n Sensitive group',cex=1)
text(-1.5,2.5,label='Upregulated in \n Insensitive group',cex=1)
abline(v=-1.5, col="red")
abline(v=1.5, col="red")
abline(h=1.3, col="green")
dev.off()

# Non-TNBC
mypch2_sensitivity <- c(rep(1,43),rep(2,135),rep(3,126),rep(4,204))	
design_sensitivity  <- model.matrix(~-1+factor(mypch2_sensitivity))
colnames(design_sensitivity)<- c('s.t', 'i.t','s.n','i.n')  #  group labels

top0_sensitivity_nontnbc <- eBayes(contrasts.fit(lmFit(normdata_sensitivity, design_sensitivity), makeContrasts(s.n - i.n, levels=design_sensitivity)))
top1_sensitivity_nontnbc <- topTable(top0_sensitivity_nontnbc , coef='s.n - i.n',adjust='BH',number=22283,sort.by='P') 

#upregulated probe_id in sensitive group
up.reg1_sensitivity_nontnbc <- top1_sensitivity_nontnbc%>%filter(top1_sensitivity_nontnbc$logFC > 1.5)
up.reg1_sensitivity_nontnbc  <-  up.reg1_sensitivity_nontnbc%>%filter(up.reg1_sensitivity_nontnbc$adj.P.Val<0.05)
names_sensitive_nontnbc <- rownames(up.reg1_sensitivity_nontnbc)
length(names_sensitive_nontnbc) #0
#upregulated probe_id in insensitive group
up.reg2_sensitivity_nontnbc  <- top1_sensitivity_nontnbc%>%filter(top1_sensitivity_nontnbc$logFC < -1.5)
up.reg2_sensitivity_nontnbc <- up.reg2_sensitivity_nontnbc %>%filter(up.reg2_sensitivity_nontnbc$adj.P.Val<0.05)
names_insensitive_nontnbc <- rownames(up.reg2_sensitivity_nontnbc)
length(names_insensitive_nontnbc) #0

png(file="volcano-nonTNBC_chemosensitivity.png",units="px")
plot(top1_sensitivity_nontnbc$logFC,-log10(top1_sensitivity_nontnbc$adj.P.Val), main='Differential Gene Expression \n between Non-TNBC Chemo-sensitive and Chemo-insensitive group', xlab='log2FC(Sensitive/Insensitive)',ylab='-log10(BH P value)', cex=0.5, cex.main=1,cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.5,4.5,label='Upregulated in \n Sensitive group',cex=1)
text(-1.5,4.5,label='Upregulated in \n Insensitive group',cex=1)
abline(v=-1.5, col="red")
abline(v=1.5, col="red")
abline(h=1.3, col="green")
dev.off()

# TNBC insensitive vs Non-TNBC insensitive
mypch2_sensitivity <- c(rep(1,43),rep(2,135),rep(3,126),rep(4,204))	
design_sensitivity  <- model.matrix(~-1+factor(mypch2_sensitivity))
colnames(design_sensitivity)<- c('s.t', 'i.t','s.n','i.n')  #  group labels

top0_sensitivity <- eBayes(contrasts.fit(lmFit(normdata_sensitivity, design_sensitivity),makeContrasts(i.t - i.n, levels=design_sensitivity)))
top1_insensitive <- topTable(top0_sensitivity , coef='i.t - i.n',adjust='BH',number=22283,sort.by='P') 

#upregulated probe_id in TNBC insensitive group
up.reg1_tnbc <- top1_insensitive%>%filter(top1_insensitive$logFC > 1.5)
up.reg1_tnbc  <-  up.reg1_tnbc%>%filter(up.reg1_tnbc$adj.P.Val<0.05)
names_insensitive_tnbc <- rownames(up.reg1_tnbc)
length(names_insensitive_tnbc) #12
#upregulated probe_id in nonTNBC  insensitive group
up.reg2_nontnbc  <- top1_insensitive%>%filter(top1_insensitive$logFC < -1.5)
up.reg2_nontnbc <- up.reg2_nontnbc %>%filter(up.reg2_nontnbc$adj.P.Val<0.05)
names_insensitive_nontnbc <- rownames(up.reg2_nontnbc)
length(names_insensitive_nontnbc) #29

png(file="volcano-insensitive_tnbc_nontnbc.png",units="px")
plot(top1_insensitive$logFC,-log10(top1_insensitive$adj.P.Val), main='Differential Gene Expression \n between Insensitive TNBC and Insensitive Non-TNBC group', xlab='log2FC(Insensitive TNBC/Insensitive Non-TNBC)',ylab='-log10(BH P value)', cex=0.5, cex.main=1,cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.5,50,label='Upregulated in \n Insensitive TNBC group',cex=0.9)
text(-1.3,50,label='Upregulated in \n Insensitive Non-TNBC group',cex=0.9)
abline(v=-1.5, col="red")
abline(v=1.5, col="red")
abline(h=1.3, col="green")
dev.off()

#probe annotation
dd_insensitive<- subset(top1_insensitive, adj.P.Val<0.05 & abs(logFC)>1.5) # subset extreme values
str(dd_insensitive) 
out_insensitive <- data.frame(xxhuex[rownames(dd_insensitive)] ) # search for rownames example

acc_no<- NULL # gene ids
for (i in 1:length(rownames(dd_insensitive))){ acc_no <- c(acc_no, out_insensitive[, paste0("X", rownames(dd_insensitive)[i])]) }
dd_insensitive <- cbind(acc_no, dd_insensitive) # get probe and gene IDs
dd_tbl_insensitive  <- as.data.table(dd_insensitive) # print as table 
df_insensitive <- data.frame(dd_insensitive)
df_insensitive$DGE <- 'insensitive'
table_insensitive <- subset(df_insensitive, select=c('acc_no','logFC','adj.P.Val'))
info <- select(hgu133a.db, rownames(table_insensitive), c("SYMBOL","GENENAME"))
info <- subset(info,select=c("SYMBOL","GENENAME"))
table_insensitive <- cbind(table_insensitive,info)
write.csv(table_insensitive, 'table_insensitive_probe_ann.csv')

# TNBC sensitive vs Non-TNBC sensitive
mypch2_sensitivity <- c(rep(1,43),rep(2,135),rep(3,126),rep(4,204))	
design_sensitivity  <- model.matrix(~-1+factor(mypch2_sensitivity))
colnames(design_sensitivity)<- c('s.t', 'i.t','s.n','i.n')  #  group labels

top0_sensitivity <- eBayes(contrasts.fit(lmFit(normdata_sensitivity, design_sensitivity),makeContrasts(s.t - s.n, levels=design_sensitivity)))
top1_sensitive <- topTable(top0_sensitivity , coef='s.t - s.n',adjust='BH',number=22283,sort.by='P') 

#upregulated probe_id in TNBC sensitive group
up.reg1_tnbc <- top1_sensitive%>%filter(top1_sensitive$logFC > 1.5)
up.reg1_tnbc  <-  up.reg1_tnbc%>%filter(up.reg1_tnbc$adj.P.Val<0.05)
names_sensitive_tnbc <- rownames(up.reg1_tnbc)
length(names_sensitive_tnbc) #18
#upregulated probe_id in nonTNBC  sensitive group
up.reg2_nontnbc  <- top1_sensitive%>%filter(top1_sensitive$logFC < -1.5)
up.reg2_nontnbc <- up.reg2_nontnbc %>%filter(up.reg2_nontnbc$ adj.P.Val<0.05)
names_sensitive_nontnbc <- rownames(up.reg2_nontnbc)
length(names_sensitive_nontnbc) #53

png(file="volcano-sensitive_tnbc_nontnbc.png",units="px")
plot(top1_sensitive$logFC,-log10(top1_sensitive$adj.P.Val), main='Differential Gene Expression \n between Sensitive TNBC and Sensitive Non-TNBC group', xlab='log2FC(Sensitive TNBC/Sensitive Non-TNBC)',ylab='-log10(BH P value)', cex=0.5, cex.main=1,cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.5,25,label='Upregulated in \n Sensitive TNBC group',cex=0.9)
text(-1.4,25,label='Upregulated in \n Sensitive Non-TNBC group',cex=0.9)
abline(v=-1.5, col="red")
abline(v=1.5, col="red")
abline(h=1.3, col="green")
dev.off()

#probes annotation
dd_sensitive<- subset(top1_sensitive, adj.P.Val<0.05 & abs(logFC)>1.5) # subset extreme values
str(dd_sensitive) 
out_sensitive <- data.frame(xxhuex[rownames(dd_sensitive)] ) # search for rownames example

acc_no<- NULL # gene ids
for (i in 1:length(rownames(dd_sensitive))){ acc_no <- c(acc_no, out_sensitive[, paste0("X", rownames(dd_sensitive)[i])]) }
dd_sensitive <- cbind(acc_no, dd_sensitive) # get probe and gene IDs
dd_tbl_sensitive  <- as.data.table(dd_sensitive) # print as table 
df_sensitive <- data.frame(dd_sensitive)
df_sensitive$DGE <- 'sensitive'
table_sensitive <- subset(df_sensitive, select=c('acc_no','logFC','adj.P.Val'))
info <- select(hgu133a.db, rownames(table_sensitive), c("SYMBOL","GENENAME"))
table_sensitive <- cbind(table_sensitive,info)
write.csv(table_sensitive, 'table_sensitive_probe_ann.csv')

#----------------------------------------------- DGE analysis (based on AJCC) ------------------------------------------
tnbc <- tnbc[order(tnbc$ajcc),] 
nontnbc <- nontnbc[order(nontnbc$ajcc),]

#filter by stage number
stage1t <- tnbc %>% filter(ajcc == "I" | ajcc == "IA" | ajcc == "IB" | ajcc == "IC") #4 samples
stage2t <- tnbc %>% filter(ajcc == "II" | ajcc == "IIA" | ajcc == "IIB" | ajcc == "IIC") #85
stage3t <- tnbc %>% filter(ajcc == "III" | ajcc == "IIIA" | ajcc == "IIIB" | ajcc == "IIIC") #87
stage1n <- nontnbc %>% filter(ajcc == "I" | ajcc == "IA" | ajcc == "IB" | ajcc == "IC") #4 cases
stage2n <- nontnbc %>% filter(ajcc == "II" | ajcc == "IIA" | ajcc == "IIB" | ajcc == "IIC") #187
stage3n <- nontnbc %>% filter(ajcc == "III" | ajcc == "IIIA" | ajcc == "IIIB" | ajcc == "IIIC") #137

ajcc <- rbind(stage1t, stage2t, stage3t,stage1n, stage2n, stage3n)
other <- sdata[!sdata$sample.no.%in%ajcc$sample.no.,]   #4
other$subtypes = 'other'

ajcc <- rbind(ajcc,other)
rawdata_ajcc  <- read.celfiles(ajcc$celfiles)
normdata_ajcc <-rma(rawdata_ajcc) #normalise the rawdata
eset_ajcc <- exprs(normdata_ajcc)

# compare stage2 (TNBC) to stage 2 (non-TNBC)
# draw volcano plots for each pairwise comparison
mypch_ajcc <- c(rep(1,4),rep(2,85),rep(3,87),rep(4,4),rep(5,187),rep(6,137),rep(7,4)) 
design_ajcc <- model.matrix(~-1+factor(mypch_ajcc))
colnames(design_ajcc)<- c('Stage1t', 'Stage2t','Stage3t','Stage1n', 'Stage2n','Stage3n','other') #  group labels

top0_ajcc <-  eBayes(contrasts.fit(lmFit(normdata_ajcc, design_ajcc), makeContrasts(Stage2t - Stage2n, levels=design_ajcc))) #compare stage2 of tnbc to stage2 of nontnbc
top1_stage2 <- topTable(top0_ajcc, coef='Stage2t - Stage2n',adjust='BH',number=22283,sort.by='P')

#upregulated probes in stage 2t
up.reg_tnbc_ajcc <- top1_stage2 %>%filter(top1_stage2$logFC > 1.5) 
up.reg_tnbc_ajcc  <- up.reg_tnbc_ajcc %>%filter(up.reg_tnbc_ajcc$ adj.P.Val<0.05)
names_upreg_tnbc_ajcc <- rownames(up.reg_tnbc_ajcc) #9
#upregulated probes in stage 2n
up.reg_nontnbc_ajcc <- top1_stage2%>%filter(top1_stage2 $logFC < -1.5)
up.reg_nontnbc_ajcc  <- up.reg_nontnbc_ajcc %>%filter(up.reg_nontnbc_ajcc$ adj.P.Val<0.05)
names_upreg_nontnbc_ajcc <- rownames(up.reg_nontnbc_ajcc) #31

png(file="volcano-TNBC-ajcc-stage2t_vs_stage2n.png",units="px")
plot(top1_stage2$logFC,-log10(top1_stage2$adj.P.Val), main='Differential Gene Expression \n between TNBC Stage2 and Non-TNBC Stage2', xlab='log2FC(TNBC Stage2/Non-TNBC Stage2)',ylab='-log10(BH P value)', cex=0.5, cex.main=1.2,cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.5,40,label='Upregulated in \n TNBC Stage2',cex=1)
text(-1.5,40,label='Upregulated in \n Non-TNBC Stage2',cex=1)
abline(v=-1.5, col="red")
abline(v=1.5, col="red")
abline(h=1.3, col="green")
dev.off()

#probes annotation
dd_stage2 <- subset(top1_stage2, adj.P.Val<0.05 & abs(logFC)>1.5) # subset extreme values
str(dd_stage2) 
out_stage2 <- data.frame(xxhuex[rownames(dd_stage2)] ) # search for rownames example

acc_no_stage2 <- NULL # gene ids
for (i in 1:length(rownames(dd_stage2))){ acc_no_stage2 <- c(acc_no_stage2, out_stage2[, paste0("X", rownames(dd_stage2)[i])]) }
dd_stage2 <- cbind(acc_no_stage2, dd_stage2) # get probe and gene IDs
dd_tbl_stage2  <- as.data.table(dd_stage2) # print as table 
df_stage2 <- data.frame(dd_stage2)
df_stage2 $DGE <- 'AJCC-stage2'
table_tnbc_nontnbc_stage2 <- subset(df_stage2, select=c('acc_no_stage2','logFC','adj.P.Val'))
info <- select(hgu133a.db, rownames(table_tnbc_nontnbc_stage2), c("SYMBOL","GENENAME"))
table_tnbc_nontnbc_stage2 <- cbind(table_tnbc_nontnbc_stage2,info)
write.csv(table_tnbc_nontnbc_stage2, 'table_ajcc_stage2_probe_ann.csv')

# compare stage3 (TNBC)  to stage 3 (non-TNBC)
# draw volcano plots for each pairwise comparison
mypch_ajcc <- c(rep(1,4),rep(2,85),rep(3,87),rep(4,4),rep(5,187),rep(6,137),rep(7,4)) 
design_ajcc <- model.matrix(~-1+factor(mypch_ajcc))
colnames(design_ajcc)<- c('Stage1t', 'Stage2t','Stage3t','Stage1n', 'Stage2n','Stage3n','other') #  group labels

top0_ajcc <-  eBayes(contrasts.fit(lmFit(normdata_ajcc, design_ajcc), makeContrasts(Stage3t - Stage3n, levels=design_ajcc))) #compare stage3 of tnbc to stage3 of nontnbc
top1_stage3 <- topTable(top0_ajcc, coef='Stage3t - Stage3n',adjust='BH',number=22283,sort.by='P')

#upregulated probes in stage 3t
up.reg_tnbc_ajcc <- top1_stage3 %>%filter(top1_stage3$logFC > 1.5) 
up.reg_tnbc_ajcc  <- up.reg_tnbc_ajcc %>%filter(up.reg_tnbc_ajcc$ adj.P.Val<0.05)
names_upreg_tnbc_ajcc <- rownames(up.reg_tnbc_ajcc) #10
#upregulated probes in stage 3n
up.reg_nontnbc_ajcc <- top1_stage3%>%filter(top1_stage3 $logFC < -1.5)
up.reg_nontnbc_ajcc  <- up.reg_nontnbc_ajcc %>%filter(up.reg_nontnbc_ajcc$ adj.P.Val<0.05)
names_upreg_nontnbc_ajcc <- rownames(up.reg_nontnbc_ajcc) #33

png(file="volcano-TNBC-ajcc-stage3t_vs_stage3n.png",units="px")
plot(top1_stage3$logFC,-log10(top1_stage3$adj.P.Val), main='Differential Gene Expression \n between TNBC Stage3 and Non-TNBC Stage3', xlab='log2FC(TNBC Stage3/Non-TNBC Stage3)',ylab='-log10(BH P value)', cex=0.5, cex.main=1.2,cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.5,40,label='Upregulated in \n TNBC Stage3',cex=1)
text(-1.5,40,label='Upregulated in \n Non-TNBC Stage3',cex=1)
abline(v=-1.5, col="red")
abline(v=1.5, col="red")
abline(h=1.3, col="green")
dev.off()

# probes annotation
dd_stage3 <- subset(top1_stage3, adj.P.Val<0.05 & abs(logFC)>1.5) # subset extreme values
str(dd_stage3) 
out_stage3 <- data.frame(xxhuex[rownames(dd_stage3)] ) # search for rownames example

acc_no_stage3 <- NULL # gene ids
for (i in 1:length(rownames(dd_stage3))){ acc_no_stage3 <- c(acc_no_stage3, out_stage3[, paste0("X", rownames(dd_stage3)[i])]) }
dd_stage3 <- cbind(acc_no_stage3, dd_stage3) # get probe and gene IDs
dd_tbl_stage3  <- as.data.table(dd_stage3) # print as table 
df_stage3 <- data.frame(dd_stage3)
df_stage3$DGE <- 'AJCC-stage3'
table_tnbc_nontnbc_stage3 <- subset(df_stage3, select=c('acc_no_stage3','logFC','adj.P.Val'))
info <- select(hgu133a.db, rownames(table_tnbc_nontnbc_stage3), c("SYMBOL","GENENAME"))
table_tnbc_nontnbc_stage3 <- cbind(table_tnbc_nontnbc_stage3,info)
write.csv(table_tnbc_nontnbc_stage3, 'table_ajcc_stage3_probe_ann.csv')

#compare stage2 (TNBC)  to stage 3 (TNBC)
# do volcano plots for each pairwise comparison
mypch_ajcc <- c(rep(1,4),rep(2,85),rep(3,87),rep(4,4),rep(5,187),rep(6,137),rep(7,4)) 
design_ajcc <- model.matrix(~-1+factor(mypch_ajcc))
colnames(design_ajcc)<- c('Stage1t', 'Stage2t','Stage3t','Stage1n', 'Stage2n','Stage3n','other') #  group labels

top0_tnbc <-  eBayes(contrasts.fit(lmFit(normdata_ajcc, design_ajcc), makeContrasts(Stage2t - Stage3t, levels=design_ajcc))) #compare stage2 of tnbc to stage3 of tnbc
top1_tnbc <- topTable(top0_tnbc, coef='Stage2t - Stage3t',adjust='BH',number=22283,sort.by='P')

#upregulated probes in stage 2t
up.reg2_tnbc_ajcc <- top1_tnbc %>%filter(top1_tnbc$logFC > 1.5) 
up.reg2_tnbc_ajcc  <- up.reg2_tnbc_ajcc %>%filter(up.reg2_tnbc_ajcc$adj.P.Val<0.05)
names_upreg2_tnbc_ajcc <- rownames(up.reg2_tnbc_ajcc) #9
#upregulated probes in stage 3t
up.reg3_tnbc_ajcc <- top1_tnbc %>%filter(top1_tnbc$logFC > 1.5) 
up.reg3_tnbc_ajcc  <- up.reg3_tnbc_ajcc %>%filter(up.reg3_tnbc_ajcc$adj.P.Val<0.05)
names_upreg3_tnbc_ajcc <- rownames(up.reg3_tnbc_ajcc) #9

png(file="volcano-TNBC-ajcc-stage2t_vs_stage3t.png",units="px")
plot(top1_tnbc$logFC,-log10(top1_tnbc$adj.P.Val), main='Differential Gene Expression of TNBC Stage2 vs TNBC Stage3', xlab='log2FC(TNBC Stage2/TNBC Stage3)',ylab='-log10(BH P value)', cex=0.5, cex.main=1.2,cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.5,40,label='Upregulated in \n Stage2',cex=1)
text(-1.5,40,label='Upregulated in \n Stage3',cex=1)
abline(v=-1.5, col="red")
abline(v=1.5, col="red")
abline(h=1, col="green")
dev.off()

# compare stage2 (Non-TNBC)  to stage 3 (Non-TNBC)
# do volcano plots for each pairwise comparison
mypch_ajcc <- c(rep(1,4),rep(2,85),rep(3,87),rep(4,4),rep(5,187),rep(6,137),rep(7,4)) 
design_ajcc <- model.matrix(~-1+factor(mypch_ajcc))
colnames(design_ajcc)<- c('Stage1t', 'Stage2t','Stage3t','Stage1n', 'Stage2n','Stage3n','other') #  group labels

top0_nontnbc <-  eBayes(contrasts.fit(lmFit(normdata_ajcc, design_ajcc), makeContrasts(Stage2n - Stage3n, levels=design_ajcc))) #compare stage2 of nontnbc to stage3 of nontnbc
top1_nontnbc <- topTable(top0_nontnbc, coef='Stage2n - Stage3n',adjust='BH',number=22283,sort.by='P')

#upregulated probes in stage 2n
up.reg2_nontnbc_ajcc <- top1_nontnbc %>%filter(top1_nontnbc$logFC > 1.5) 
up.reg2_nontnbc_ajcc  <- up.reg2_nontnbc_ajcc %>%filter(up.reg2_nontnbc_ajcc$adj.P.Val<0.05)
names_upreg2_nontnbc_ajcc <- rownames(up.reg2_nontnbc_ajcc) #9
#upregulated probes in stage 3n
up.reg3_nontnbc_ajcc <- top1_nontnbc %>%filter(top1_nontnbc$logFC > 1.5) 
up.reg3_nontnbc_ajcc  <- up.reg3_nontnbc_ajcc %>%filter(up.reg3_nontnbc_ajcc$adj.P.Val<0.05)
names_upreg3_nontnbc_ajcc <- rownames(up.reg3_nontnbc_ajcc) #9

png(file="volcano-NonTNBC-ajcc-stage2n_vs_stage3n.png",units="px")
plot(top1_nontnbc$logFC,-log10(top1_nontnbc$adj.P.Val), main='Differential Gene Expression of \n Non-TNBC Stage2 vs Non-TNBC Stage3', xlab='log2FC(Non-TNBC Stage2/Non-TNBC Stage3)',ylab='-log10(BH P value)', cex=0.5, cex.main=1,cex.axis=1, cex.lab=1, pch=21,xlim=c(-2,2))
text(1.2,40,label='Upregulated in \n Stage2',cex=0.8)
text(-1.2,40,label='Upregulated in \n Stage3',cex=0.8)
abline(v=-1.5, col="red")
abline(v=1.5, col="red")
abline(h=1, col="green")
dev.off()

#----------------------------------------------- Ca2+ channel mRNA expression ------------------------------------------
#ORAI3 expression in TNBC and NonTNBC
orai3 <- subset(eset2,  rownames(eset2) == '221864_at')
orai3_tnbc <-subset(orai3,select=c(tnbc$celfiles))
orai3_tnbc <- t(orai3_tnbc)
tnbc <- cbind(tnbc,orai3_tnbc)
x <- tnbc$'221864_at'

orai3 <- subset(eset2,  rownames(eset2) == '221864_at')
orai3_nontnbc <-subset(orai3,select=c(nontnbc$celfiles))
orai3_nontnbc <- t(orai3_nontnbc)
nontnbc <- cbind(nontnbc,orai3_nontnbc)
y <-nontnbc$'221864_at'

x <- data.frame(tnbc$'221864_at')
x$groups = 'TNBC'
names(x)[1] <- "value"
y <- data.frame(nontnbc$'221864_at')
y$groups = 'Non-TNBC'
names(y)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

png('boxplot_ORAI3_tnbc_nontnbc.png')
boxplot(value ~ groups, data = all, main = "ORAI3 Expression Value in Breast Cancer Subtypes", xlab = "Subtypes", ylab = "ORAI3 Expression Value", col = c("#00AFBB", "#E7B800"))
points(1:2, means$value, col = "red") 
text(1:2, means$value + 0.06, labels = round(means$value,3),cex = 0.9)
dev.off()

# ORAI3 in TNBC rd and TNBC pcr
rd_tnbc <- tnbc%>%filter(pcr_rd == 'RD') #113
pcr_tnbc <- tnbc%>%filter(pcr_rd == 'pCR') #57
rd_nontnbc <- nontnbc%>%filter(pcr_rd == 'RD') #276
pcr_nontnbc <- nontnbc%>%filter(pcr_rd == 'pCR') #42

orai3 <- subset(eset2,  rownames(eset2) == '221864_at')
orai3_tnbc_rd <-subset(orai3,select=c(rd_tnbc$celfiles))
orai3_tnbc_rd <- t(orai3_tnbc_rd)
tnbc_rd <- cbind(rd_tnbc,orai3_tnbc_rd)
orai3_tnbc_rd  <- tnbc_rd$'221864_at'

orai3 <- subset(eset2,  rownames(eset2) == '221864_at')
orai3_tnbc_pcr <-subset(orai3,select=c(pcr_tnbc$celfiles))
orai3_tnbc_pcr <- t(orai3_tnbc_pcr)
tnbc_pcr <- cbind(pcr_tnbc,orai3_tnbc_pcr)
orai3_tnbc_pcr <- tnbc_pcr$'221864_at'

x <- data.frame(tnbc_rd$'221864_at')
x$groups = 'RD'
names(x)[1] <- "value"
y <- data.frame(tnbc_pcr$'221864_at')
y$groups = 'pCR'
names(y)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

png('boxplot_orai3_tnbc.png')
boxplot(value ~ groups, data = all, main = "ORAI3 Expression Value in TNBC", xlab = "TNBC", ylab = "ORAI3 Expression Value", col = c("#00AFBB", "#E7B800"))
points(1:2, means$value, col = "red") 
text(1:2, means$value + 0.06, labels = round(means$value,3),cex = 0.9)
dev.off()

# ORAI3 in Non-TNBC rd and Non-TNBC pcr
orai3_nontnbc_rd <-subset(orai3,select=c(rd_nontnbc$celfiles))
orai3_nontnbc_rd <- t(orai3_nontnbc_rd)
nontnbc_rd <- cbind(rd_nontnbc,orai3_nontnbc_rd)
orai3_nontnbc_rd  <- nontnbc_rd$'221864_at'

orai3_nontnbc_pcr <-subset(orai3,select=c(pcr_nontnbc$celfiles))
orai3_nontnbc_pcr <- t(orai3_nontnbc_pcr)
nontnbc_pcr <- cbind(pcr_nontnbc,orai3_nontnbc_pcr)
orai3_nontnbc_pcr <- nontnbc_pcr$'221864_at'

x <- data.frame(nontnbc_rd$'221864_at')
x$groups = 'RD'
names(x)[1] <- "value"
y <- data.frame(nontnbc_pcr$'221864_at')
y$groups = 'pCR'
names(y)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

png('boxplot_orai3_nontnbc.png')
boxplot(value ~ groups, data = all, main = "ORAI3 Expression Value in Non-TNBC", xlab = "Non-TNBC", ylab = "ORAI3 Expression Value", col = c("#00AFBB", "#E7B800"))
points(1:2, means$value, col = "red") 
text(1:2, means$value + 0.12, labels = round(means$value,3),cex = 0.9)
dev.off()

# ORAI3 in TNBC rd and Non-TNBC rd
orai3_tnbc_rd <-subset(orai3,select=c(rd_tnbc$celfiles))
orai3_tnbc_rd <- t(orai3_tnbc_rd)
tnbc_rd <- cbind(rd_tnbc,orai3_tnbc_rd)
orai3_tnbc_rd  <- tnbc_rd$'221864_at'

orai3_nontnbc_rd <-subset(orai3,select=c(rd_nontnbc$celfiles))
orai3_nontnbc_rd <- t(orai3_nontnbc_rd)
nontnbc_rd <- cbind(rd_nontnbc,orai3_nontnbc_rd)
orai3_nontnbc_rd  <- nontnbc_rd$'221864_at'

y <- data.frame(nontnbc_rd$'221864_at')
y$groups = 'Non-TNBC(RD)'
names(y)[1] <- "value"
x <- data.frame(tnbc_rd$'221864_at')
x$groups = 'TNBC(RD)'
names(x)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

png('boxplot_orai3_rd.png')
boxplot(value ~ groups, data = all, main = "ORAI3 Expression Value in RD groups", xlab = "RD group in different subtypes", ylab = "ORAI3 Expression Value", col = c("#00AFBB", "#E7B800"))
points(1:2, means$value, col = "red") 
text(1:2, means$value + 0.06, labels = round(means$value,3),cex = 0.9)
dev.off()

#--------------- TRPC1 in TNBC rd and TNBC pcr ----------------
TRPC1 <- subset(eset2,  rownames(eset2) == '205802_at')
TRPC1_tnbc_rd <-subset(TRPC1,select=c(rd_tnbc$celfiles))
TRPC1_tnbc_rd <- t(TRPC1_tnbc_rd)
tnbc_rd <- cbind(rd_tnbc,TRPC1_tnbc_rd)
#TRPC1_tnbc_rd  <- tnbc_rd$'205802_at'

TRPC1 <- subset(eset2,  rownames(eset2) == '205802_at')
TRPC1_tnbc_pcr <-subset(TRPC1,select=c(pcr_tnbc$celfiles))
TRPC1_tnbc_pcr <- t(TRPC1_tnbc_pcr)
tnbc_pcr <- cbind(pcr_tnbc,TRPC1_tnbc_pcr)
#TRPC1_tnbc_pcr <- tnbc_pcr$'205802_at'

x <- data.frame(tnbc_rd$'205802_at')
x$groups = 'RD'
names(x)[1] <- "value"
y <- data.frame(tnbc_pcr$'205802_at')
y$groups = 'pCR'
names(y)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

png('boxplot_TRPC1_tnbc.png')
boxplot(value ~ groups, data = all, main = "TRPC1 Expression Value in TNBC",  xlab = "TNBC", ylab = "TRPC1 Expression Value", col = c("#00AFBB", "#E7B800"))
points(1:2, means$value, col = "red") 
text(1:2, means$value + 0.06, labels = round(means$value,3),cex = 0.9)
dev.off()

#--------------- TRPC5 in TNBC rd and TNBC pcr ----------------
TRPC5 <- subset(eset2,  rownames(eset2) == '220552_at')
TRPC5_tnbc_rd <-subset(TRPC5,select=c(rd_tnbc$celfiles))
TRPC5_tnbc_rd <- t(TRPC5_tnbc_rd)
tnbc_rd <- cbind(rd_tnbc,TRPC5_tnbc_rd)
#TRPC5_tnbc_rd  <- tnbc_rd$'220552_at'

TRPC5 <- subset(eset2,  rownames(eset2) == '220552_at')
TRPC5_tnbc_pcr <-subset(TRPC5,select=c(pcr_tnbc$celfiles))
TRPC5_tnbc_pcr <- t(TRPC5_tnbc_pcr)
tnbc_pcr <- cbind(pcr_tnbc,TRPC5_tnbc_pcr)
#TRPC5_tnbc_pcr <- tnbc_pcr$'220552_at'

x <- data.frame(tnbc_rd$'220552_at')
x$groups = 'RD'
names(x)[1] <- "value"
y <- data.frame(tnbc_pcr$'220552_at')
y$groups = 'pCR'
names(y)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

png('boxplot_TRPC5_tnbc.png')
boxplot(value ~ groups, data = all, main = "TRPC5 Expression Value in TNBC",  xlab = "TNBC", ylab = "TRPC5 Expression Value", col = c("#00AFBB", "#E7B800"))
points(1:2, means$value, col = "red") 
text(1:2, means$value + 0.04, labels = round(means$value,3),cex = 0.9)
dev.off()

#--------------- TRPC6 in TNBC rd and TNBC pcr ----------------
TRPC6 <- subset(eset2,  rownames(eset2) == '217287_s_at')
TRPC6_tnbc_rd <-subset(TRPC6,select=c(rd_tnbc$celfiles))
TRPC6_tnbc_rd <- t(TRPC6_tnbc_rd)
tnbc_rd <- cbind(rd_tnbc,TRPC6_tnbc_rd)
#TRPC6_tnbc_rd  <- tnbc_rd$'220552_at'

TRPC6 <- subset(eset2,  rownames(eset2) == '217287_s_at')
TRPC6_tnbc_pcr <-subset(TRPC6,select=c(pcr_tnbc$celfiles))
TRPC6_tnbc_pcr <- t(TRPC6_tnbc_pcr)
tnbc_pcr <- cbind(pcr_tnbc,TRPC6_tnbc_pcr)
#TRPC6_tnbc_pcr <- tnbc_pcr$'217287_s_at'

x <- data.frame(tnbc_rd$'217287_s_at')
x$groups = 'RD'
names(x)[1] <- "value"
y <- data.frame(tnbc_pcr$'217287_s_at')
y$groups = 'pCR'
names(y)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

png('boxplot_TRPC6_tnbc.png')
boxplot(value ~ groups, data = all, main = "TRPC6 Expression Value in TNBC",  xlab = "TNBC", ylab = "TRPC6 Expression Value", col = c("#00AFBB", "#E7B800"))
points(1:2, means$value, col = "red") 
text(1:2, means$value + 0.04, labels = round(means$value,3),cex = 0.9)
dev.off()

#--------------- ORAI2 in TNBC rd and TNBC pcr ----------------
ORAI2 <- subset(eset2,  rownames(eset2) == '217529_at')
ORAI2_tnbc_rd <-subset(ORAI2,select=c(rd_tnbc$celfiles))
ORAI2_tnbc_rd <- t(ORAI2_tnbc_rd)
tnbc_rd <- cbind(rd_tnbc,ORAI2_tnbc_rd)
#ORAI2_tnbc_rd  <- tnbc_rd$'220552_at'

ORAI2 <- subset(eset2,  rownames(eset2) == '217529_at')
ORAI2_tnbc_pcr <-subset(ORAI2,select=c(pcr_tnbc$celfiles))
ORAI2_tnbc_pcr <- t(ORAI2_tnbc_pcr)
tnbc_pcr <- cbind(pcr_tnbc,ORAI2_tnbc_pcr)
#ORAI2_tnbc_pcr <- tnbc_pcr$'220552_at'

x <- data.frame(tnbc_rd$'217529_at')
x$groups = 'RD'
names(x)[1] <- "value"
y <- data.frame(tnbc_pcr$'217529_at')
y$groups = 'pCR'
names(y)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

png('boxplot_ORAI2_tnbc.png')
boxplot(value ~ groups, data = all, main = "ORAI2 Expression Value in TNBC",  xlab = "TNBC", ylab = "ORAI2 Expression Value", col = c("#00AFBB", "#E7B800"))
points(1:2, means$value, col = "red") 
text(1:2, means$value + 0.06, labels = round(means$value,3),cex = 0.9)
dev.off()

#--------------- STIM1 in TNBC rd and TNBC pcr ----------------
STIM1 <- subset(eset2,  rownames(eset2) == '202764_at')
STIM1_tnbc_rd <-subset(STIM1,select=c(rd_tnbc$celfiles))
STIM1_tnbc_rd <- t(STIM1_tnbc_rd)
tnbc_rd <- cbind(rd_tnbc,STIM1_tnbc_rd)
#STIM1_tnbc_rd  <- tnbc_rd$'220552_at'

STIM1 <- subset(eset2,  rownames(eset2) == '202764_at')
STIM1_tnbc_pcr <-subset(STIM1,select=c(pcr_tnbc$celfiles))
STIM1_tnbc_pcr <- t(STIM1_tnbc_pcr)
tnbc_pcr <- cbind(pcr_tnbc,STIM1_tnbc_pcr)
#STIM1_tnbc_pcr <- tnbc_pcr$'220552_at'

x <- data.frame(tnbc_rd$'202764_at')
x$groups = 'RD'
names(x)[1] <- "value"
y <- data.frame(tnbc_pcr$'202764_at')
y$groups = 'pCR'
names(y)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

png('boxplot_STIM1_tnbc.png')
boxplot(value ~ groups, data = all, main = "STIM1 Expression Value in TNBC",  xlab = "TNBC", ylab = "STIM1 Expression Value", col = c("#00AFBB", "#E7B800"))
points(1:2, means$value, col = "red") 
text(1:2, means$value + 0.04, labels = round(means$value,3),cex = 0.9)
dev.off()

#--------------- TRPV6 in TNBC rd and TNBC pcr ----------------
TRPV6 <- subset(eset2,  rownames(eset2) == '206827_s_at')
TRPV6_tnbc_rd <-subset(TRPV6,select=c(rd_tnbc$celfiles))
TRPV6_tnbc_rd <- t(TRPV6_tnbc_rd)
tnbc_rd <- cbind(rd_tnbc,TRPV6_tnbc_rd)
#TRPV6_tnbc_rd  <- tnbc_rd$'206827_s_at'

TRPV6 <- subset(eset2,  rownames(eset2) == '206827_s_at')
TRPV6_tnbc_pcr <-subset(TRPV6,select=c(pcr_tnbc$celfiles))
TRPV6_tnbc_pcr <- t(TRPV6_tnbc_pcr)
tnbc_pcr <- cbind(pcr_tnbc,TRPV6_tnbc_pcr)
#TRPV6_tnbc_pcr <- tnbc_pcr$'206827_s_at'

x <- data.frame(tnbc_rd$'206827_s_at')
x$groups = 'RD'
names(x)[1] <- "value"
y <- data.frame(tnbc_pcr$'206827_s_at')
y$groups = 'pCR'
names(y)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

#--------------- CACNA1D in TNBC rd and TNBC pcr ----------------
CACNA1D <- subset(eset2,  rownames(eset2) == '210108_at')
CACNA1D_tnbc_rd <-subset(CACNA1D,select=c(rd_tnbc$celfiles))
CACNA1D_tnbc_rd <- t(CACNA1D_tnbc_rd)
tnbc_rd <- cbind(rd_tnbc,CACNA1D_tnbc_rd)
#CACNA1D_tnbc_rd  <- tnbc_rd$'210108_at'

CACNA1D <- subset(eset2,  rownames(eset2) == '210108_at')
CACNA1D_tnbc_pcr <-subset(CACNA1D,select=c(pcr_tnbc$celfiles))
CACNA1D_tnbc_pcr <- t(CACNA1D_tnbc_pcr)
tnbc_pcr <- cbind(pcr_tnbc,CACNA1D_tnbc_pcr)
#CACNA1D_tnbc_pcr <- tnbc_pcr$'210108_at'

x <- data.frame(tnbc_rd$'210108_at')
x$groups = 'RD'
names(x)[1] <- "value"
y <- data.frame(tnbc_pcr$'210108_at')
y$groups = 'pCR'
names(y)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

#--------------- CACNA1G in TNBC rd and TNBC pcr ----------------
CACNA1G <- subset(eset2,  rownames(eset2) == '211314_at')
CACNA1G_tnbc_rd <-subset(CACNA1G,select=c(rd_tnbc$celfiles))
CACNA1G_tnbc_rd <- t(CACNA1G_tnbc_rd)
tnbc_rd <- cbind(rd_tnbc,CACNA1G_tnbc_rd)
#CACNA1G_tnbc_rd  <- tnbc_rd$'211314_at'

CACNA1G <- subset(eset2,  rownames(eset2) == '211314_at')
CACNA1G_tnbc_pcr <-subset(CACNA1G,select=c(pcr_tnbc$celfiles))
CACNA1G_tnbc_pcr <- t(CACNA1G_tnbc_pcr)
tnbc_pcr <- cbind(pcr_tnbc,CACNA1G_tnbc_pcr)
#CACNA1G_tnbc_pcr <- tnbc_pcr$'211314_at'

x <- data.frame(tnbc_rd$'211314_at')
x$groups = 'RD'
names(x)[1] <- "value"
y <- data.frame(tnbc_pcr$'211314_at')
y$groups = 'pCR'
names(y)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

#--------------- CACNA1H in TNBC rd and TNBC pcr ----------------
CACNA1H <- subset(eset2,  rownames(eset2) == '205845_at')
CACNA1H_tnbc_rd <-subset(CACNA1H,select=c(rd_tnbc$celfiles))
CACNA1H_tnbc_rd <- t(CACNA1H_tnbc_rd)
tnbc_rd <- cbind(rd_tnbc,CACNA1H_tnbc_rd)
#CACNA1H_tnbc_rd  <- tnbc_rd$'205845_at'

CACNA1H <- subset(eset2,  rownames(eset2) == '205845_at')
CACNA1H_tnbc_pcr <-subset(CACNA1H,select=c(pcr_tnbc$celfiles))
CACNA1H_tnbc_pcr <- t(CACNA1H_tnbc_pcr)
tnbc_pcr <- cbind(pcr_tnbc,CACNA1H_tnbc_pcr)
#CACNA1H_tnbc_pcr <- tnbc_pcr$'205845_at'

x <- data.frame(tnbc_rd$'205845_at')
x$groups = 'RD'
names(x)[1] <- "value"
y <- data.frame(tnbc_pcr$'205845_at')
y$groups = 'pCR'
names(y)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

#--------------- CACNA1I in TNBC rd and TNBC pcr ----------------
CACNA1I <- subset(eset2,  rownames(eset2) == '208299_at')
CACNA1I_tnbc_rd <-subset(CACNA1I,select=c(rd_tnbc$celfiles))
CACNA1I_tnbc_rd <- t(CACNA1I_tnbc_rd)
tnbc_rd <- cbind(rd_tnbc,CACNA1I_tnbc_rd)
#CACNA1I_tnbc_rd  <- tnbc_rd$'208299_at'

CACNA1I <- subset(eset2,  rownames(eset2) == '208299_at')
CACNA1I_tnbc_pcr <-subset(CACNA1I,select=c(pcr_tnbc$celfiles))
CACNA1I_tnbc_pcr <- t(CACNA1I_tnbc_pcr)
tnbc_pcr <- cbind(pcr_tnbc,CACNA1I_tnbc_pcr)
#CACNA1I_tnbc_pcr <- tnbc_pcr$'208299_at'

x <- data.frame(tnbc_rd$'208299_at')
x$groups = 'RD'
names(x)[1] <- "value"
y <- data.frame(tnbc_pcr$'208299_at')
y$groups = 'pCR'
names(y)[1] <- "value"
all <- rbind(x,y)

# t test
group_by(all, groups) %>% summarise(count = n(),mean = mean(value, na.rm = TRUE),sd = sd(value, na.rm = TRUE))
t.test <- t.test(value ~ groups, data = all)
means <- aggregate(value ~ groups, all, mean)

#--------------------- correct the p-values ---------------------
#adjust p-values 
pvalues <- c(0.379,0.126,0.394,0.965,0.224,0.258,0.288,0.159,0.086,0.052,0.029)
p.adjust(pvalues,method="BH")
