#############################################################
### Example: Three BCR-ABL1 vairants: p230, p210, and p185
#############################################################
#
rm(list=ls())	# clear the worksapce
#
Home_Dir<-getwd()
library("FusionPathway")

## Input Data: Gene ID, protein domain of two fusion parental genes
Gene1="BCR"
GeneID1<-613
PFAM1=c("PF09036","PF00621","PF00169","PF00168","PF00620")
Gene2="ABL1"
GeneID2<-25
PFAM2=c("PF00018","PF00017","PF07714","PF08919")

## Result folder
Home_Dir<-file.path(Home_Dir,paste(Gene1,Gene2,sep="_"))
dir.create(Home_Dir)
setwd(Home_Dir)

##
for (Variant in 1:3) {
	# variants and their domain composition
	if(Variant==1)  {
		# p230
		PFAM1_lost<-"PF00620"	
		PFAM1_kept<-setdiff(PFAM1,PFAM1_lost)
		PFAM2_lost<-""
		PFAM2_kept<-setdiff(PFAM2,PFAM2_lost)
		Sub_dir="p230"
	} else if (Variant==2) {	
		#p210
		PFAM1_lost<-c("PF00620","PF00168")	
		PFAM1_kept<-setdiff(PFAM1,PFAM1_lost)
		PFAM2_lost<-""
		PFAM2_kept<-setdiff(PFAM2,PFAM2_lost)
		Sub_dir="p210"
	} else if (Variant==3) {	
		#p185
		PFAM1_lost<-c("PF00620","PF00168","PF00169","PF00621") 
		PFAM1_kept<-setdiff(PFAM1,PFAM1_lost)
		PFAM2_lost<-""
		PFAM2_kept<-setdiff(PFAM2,PFAM2_lost)
		Sub_dir="p185"
	} 
	
	## Folder for variants of the BCR-ABL1 variants
	print(paste("Fusion ",paste(Gene1,Gene2,sep="-"),": ",Sub_dir,sep=""))
	Dir2<-file.path(Home_Dir,Sub_dir)
	dir.create(Dir2)
	setwd(Dir2)
	
	## run FusionPathway
	DomainData<-list(pfam1_kept=PFAM1_kept,pfam1_lost=PFAM1_lost,pfam2_kept=PFAM2_kept,pfam2_lost=PFAM2_lost)
	GeneData<-data.frame(Gene1=Gene1,Gene2=Gene2,GeneID1=GeneID1,GeneID2=GeneID2)
	Result_List<-FusionPathway(GeneData,DomainData)
	save(file=file.path(Dir2,"Ranked_Result.RData"),Result_List)
		
	## Prediction evaluation using different benchmark gene sets
	Score<-as.numeric(Result_List$Ranked_Result$Score)
	GeneID_ranked<-as.numeric(Result_List$Ranked_Result$GeneID_ranked)
	load(file.path(system.file("data", package = "FusionPathway"), "BenchmarkSets_CML.RData"))
	for (i in 1:5) {
		##
		if (i==1) {
			BenchmarkName="CML_Genes"
			PubmedN<-BenchmarkSets_CML[["CML_Genes"]]$PubmedN
			GeneIDs<-BenchmarkSets_CML[["CML_Genes"]]$GeneIDs
			Ncitation_threshold=2
			#
			BenchID_pos<-GeneIDs[which(PubmedN>=Ncitation_threshold)]
			#
			ROC_M<-NULL
		} else if (i==2) {
			BenchmarkName="BCR_ABL1_Genes"
			PubmedN<-BenchmarkSets_CML[["BCR_ABL1_Genes"]]$PubmedN
			GeneIDs<-BenchmarkSets_CML[["BCR_ABL1_Genes"]]$GeneIDs
			Ncitation_threshold=2
			BenchID_pos<-GeneIDs[which(PubmedN>=Ncitation_threshold)]
		} else if (i==3) {
			BenchmarkName="WNT-CA-NFAT_Genes"
			BenchID_pos<-BenchmarkSets_CML[["WntCaNfat_Genes"]]$GeneIDs
		} else if (i==4) {
			BenchmarkName="CancerPathway_Genes"
			BenchID_pos<-BenchmarkSets_CML[["CancerPathway_Genes"]]$GeneIDs
		} else if (i==5) {
			BenchmarkName="Drug_TargetGenes"
			BenchID_pos<-BenchmarkSets_CML[["Drug_Targets"]]$GeneID
		}
		##
		BenchID_neg<-setdiff(GeneID_ranked,BenchID_pos)
		BenchID<-c(BenchID_pos,BenchID_neg)
	
		## Calculate AUC values for ROC curves
		loc<-match(GeneID_ranked,BenchID_pos)
		outcome=rep(0,length(GeneID_ranked))
		outcome[which(loc>0)]=1
		ROC<-ROCF(outcome,Score)
		ROC_M[[BenchmarkName]]$TPR<-ROC$TPR
		ROC_M[[BenchmarkName]]$FPR<-ROC$FPR
		ROC_M[[BenchmarkName]]$AUC<-round(ROC$AUC,3)
		print(paste("AUC using ",BenchmarkName," genes: ",round(ROC$AUC,3),sep=""))
		
		## ROC curve 
		outputfile=paste("AUC_Benchmark.",BenchmarkName,".tiff",sep="")
		tiff(outputfile,width = 510, height = 510, units = "px")
		pp<-plot(ROC$FPR,ROC$TPR,xlab="FPR",ylab="TPR",col = "red",lwd=3,"l",
			cex=1.45,cex.axis=1.45,cex.lab=1.45,font.lab=2.3)
		LabelS<-paste(BenchmarkName,"(AUC=",round(ROC$AUC,3),")",sep="")
		legend("bottomright",col="red",legend=LabelS,lty=rep(1,5),lwd=3,cex=1.4,text.font=1.5)
		dev.off()
	}
	
	## Save ROC
	output=file.path(Dir2,"ROC.RData")
	save(file=output,ROC_M)
	print ('')
	print ('')
}



###########################################################################################
### Figure 2: Evaluation of p210 BCR-ABL1 prediction 
###########################################################################################
###
### Combine p210 GSEA pathway pavlues of GSEA association & GSEA deregulation (for Figure 2a)
Work_Dir=file.path(Home_Dir,"p210")
setwd(Work_Dir)
#
## Gene expression data (GSE24493; Duy et al., 2011)
## differntially expression anlysis of three CML cell lines and their imatinib-treated CML cells
DataM1 <- read.delim(file.path(system.file("testData", package = "FusionPathway"), "ArrayData_GSE24493_Duy.2011.txt"), sep = "\t", header = TRUE, as.is = TRUE)
Gene.symbol<-as.character(DataM1$Gene.symbol)
# delete no symbols
pos<-which(Gene.symbol!="")
DataM1 <- DataM1[pos,]
P.Value<-as.numeric(as.character(DataM1$P.Value))
# sort genes based on p.value
AA<-sort(P.Value,decreasing =F, index.return = TRUE)
indx<-as.numeric(AA$ix)
DataM1 <- DataM1[indx,]
Gene.symbol<-as.character(DataM1$Gene.symbol)
P.Value<-as.numeric(as.character(DataM1$P.Value))
logFC<-as.numeric(as.character(DataM1$logFC))
#
## run fGSEA (GSEA deregulation)
library("fgsea")
load(file.path(system.file("data", package = "FusionPathway"), "GSEAPathway_Data.RData"))
Pathway_Names<-names(Pathway_M)
Score_ranked<-(length(Gene.symbol):1)
names(Score_ranked)<-Gene.symbol
fgseaRes <- fgsea(pathways = Pathway_M, stats = Score_ranked,
        minSize=15, maxSize=500, nperm=10000)
# pathways and write the results
FDR_threshold=1
Result_GSEA<-fgseaRes[fgseaRes[, padj < FDR_threshold],]
Result_GSEA = format(Result_GSEA)
fileoutput1="GSE24493_fGSEA_Result.xls"
write.table(Result_GSEA,fileoutput1, sep = "\t",append=FALSE,row.names = FALSE,col.names = TRUE,quote=FALSE)
#
## Pathway GSEA association & GSEA deregulation
filein1="BCR-ABL1_GSEA_Result.xls"
PathwayM1 <- read.table(filein1,sep="\t",header=T, quote="")
Pathways1<-as.character(PathwayM1$pathway)
filein2="GSE24493_fGSEA_Result.xls"
PathwayM2 <- read.table(filein2,sep="\t",header=T, quote="")
Pathways2<-as.character(PathwayM2$pathway)
loc<-match(Pathways1,Pathways2)
pos<-which(loc>0)
loc<-loc[pos]
PathwayM1<- PathwayM1[pos,]
PathwayM2<- PathwayM2[loc,]
#
## Combine pathway pvalues of GSEA association & GSEA deregulation
outfile="BCR-ABL1_GSEA_CombinedResult.xls"
ResultM<-Combine_Pathway_Pvalues(PathwayM1,PathwayM2,threshold_ind_pvalue=0.05,outfile)



###
### Figure 2a: Bar Chart of p-value of selected associated pathways for prediction of p210 BCR-ABL1
# output file: "Enrichment_AssocaitedPathways.tiff" in the p210 folder
#
## Selected Pathway
Pathway_Selected=c("REGULATION_OF_CELL_CYCLE","KEGG_JAK_STAT_SIGNALING_PATHWAY",
	"KEGG_P53_SIGNALING_PATHWAY","REGULATION_OF_APOPTOSIS",
	"BIOCARTA_GLEEVEC_PATHWAY",
	"KEGG_PATHWAYS_IN_CANCER","DNA_DAMAGE_CHECKPOINT",
	"KEGG_MAPK_SIGNALING_PATHWAY",
	"REACTOME_PI3K_AKT_ACTIVATION","BIOCARTA_MTOR_PATHWAY",
	"PID_HEDGEHOG_2PATHWAY","WNT_SIGNALING","PID_MYC_ACTIVPATHWAY","BIOCARTA_RAS_PATHWAY")	
#
## GSEA results (GSEA association, GSEA deregulation, and combined result) of the selected pathways
file_AsscoatedPathways<-"BCR-ABL1_GSEA_CombinedResult.xls"
DataM1 <- read.table(file_AsscoatedPathways,sep="\t",header=TRUE, quote="")
Pathways1<-as.character(DataM1[[1]])
loc<-match(Pathway_Selected,Pathways1)
DataM1 <- DataM1 [loc,]
Pathway_Selected<-as.character(DataM1[[1]])
pval1<-log10(as.numeric(DataM1$pval1))*(-1)	
pval2<-log10(as.numeric(DataM1$pval2))*(-1)		
combined_pval<-log10(as.numeric(DataM1$combined_pval))*(-1)	
combined_adj<-log10(as.numeric(DataM1$combined_adj))*(-1)	
#
##
ES_data<-combined_pval
Order<-(1:length(pval2))
DD<-sort(ES_data,decreasing =F,index.return=TRUE)
indx<-as.numeric(DD$ix)
pval1<-pval1[indx]
pval2<-pval2[indx]
combined_adj<-combined_adj[indx]
combined_pval<-combined_pval[indx]
Pathway_Selected<-Pathway_Selected[indx]
DataP<-data.frame(pathways=Pathway_Selected,pvalue1=pval1,pvalue2=pval2,combined_pval=combined_pval)
DataP$pathways= factor(DataP$pathways,Pathway_Selected) # change order
#
##
library(ggplot2)
library(gridExtra)
outputfile="Enrichment_AssocaitedPathways.tiff"
tiff(outputfile,width = 610, height = 300, units = "px")
g1<-ggplot(data = DataP, aes(x = pathways, y = pvalue1)) +
  geom_bar(stat = "identity", fill="#E69F00") + ggtitle("GSEA association")+
  labs(y="-log10(p value)", x = "", size=10)+ 
  theme(plot.title = element_text(face="bold"),
  axis.title.x =element_text(face="bold"), axis.text.x = element_text(size=10)
  , axis.title.y = element_blank(), axis.text.y =element_text(face="bold", size=10)
  , plot.margin = unit(c(1,0,1,0), "mm")) + coord_flip()
#  
g2<-ggplot(data = DataP, aes(x = pathways, y = pvalue2)) +
  geom_bar(stat = "identity", fill="#0072B2") + ggtitle("GSEA deregulation")+
  labs(y="-log10(p value)", x = "", size=10) +
  theme(plot.title = element_text(face="bold"), 
  axis.title.x =element_text(face="bold"), axis.text.x = element_text(size=10)
  , axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()
  , plot.margin = unit(c(1,0,1,0), "mm")) +coord_flip()
#
g3<-ggplot(data = DataP, aes(x = pathways, y = combined_pval)) +
  geom_bar(stat = "identity", fill="red") + ggtitle("Combined p values ")+
  labs(y="-log10(p value)", x = "", size=10)+ 
  theme(plot.title = element_text(face="bold"),
  axis.title.x =element_text(face="bold"), axis.text.x = element_text(size=10)
  , axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()
  , plot.margin = unit(c(1,0,1,0), "mm")) + coord_flip()
#
gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg3 <- ggplot_gtable(ggplot_build(g3))
grid.arrange(g1,g2,g3,ncol=3,widths=c(7.5/10,3.0/10,3.0/10))
dev.off()



###
### Figure 2c: Combine ROC plots for different benchmarks 
# output file: "AUC_Different_Benchmarks_v1.tiff" in the p210 folder
#
## Load ROC data
load("ROC.RData")
#
outputfile="AUC_Different_Benchmarks_v1.tiff"
tiff(outputfile,width = 510, height = 510, units = "px")
par(mar=c(4,4.2,0.5,0.5))
pp<-plot(ROC_M$CML_Genes$FPR,ROC_M$CML_Genes$TPR,xlab="FPR",ylab="TPR",col = "red",lwd=3,"l",
	cex=1.45,cex.axis=1.45,cex.lab=1.45,font.lab=2.3)
points(ROC_M$BCR_ABL1_Genes$FPR,ROC_M$BCR_ABL1_Genes$TPR,type="l",col="black",lwd=3)
points(ROC_M$CancerPathway_Genes$FPR,ROC_M$CancerPathway_Genes$TPR,type="l",col="blue",lwd=3)
points(ROC_M[["WNT-CA-NFAT_Genes"]]$FPR,ROC_M[["WNT-CA-NFAT_Genes"]]$TPR,type="l",col="purple",lwd=3)
points(ROC_M$Drug_TargetGenes$FPR,ROC_M$Drug_TargetGenes$TPR,type="l",col="green",lwd=3)
LabelS<-c(paste("CML Genes(AUC=",ROC_M$CML_Genes$AUC,")",sep=""),
	paste("BCR-ABL Genes(AUC=",ROC_M$BCR_ABL1_Genes$AUC,")",sep=""),
	paste("Cancer Pathway Genes(AUC=",ROC_M$CancerPathway_Genes$AUC,")",sep=""),
	paste("WNT_CA+_NFAT Genes(AUC=",ROC_M[["WNT-CA-NFAT_Genes"]]$AUC,")",sep=""),
	paste("Drug Targets(AUC=",ROC_M$Drug_TargetGenes$AUC,")",sep=""))
legend(0.145,0.27,col=c("red","black","blue","purple","green"),
	legend=LabelS,lty=rep(1,5),lwd=3,cex=1.4,text.font=2)
dev.off()




###########################################################################################
### Figure 3: Prediction evaluation of three known BCR-ABL1 variants (p185, p210, and p230) 
###########################################################################################
###
## Figure 3a: Compare number of PPI partners
# output file: "Number_PPIs.tiff" in the Comparison folder
#
setwd(Home_Dir)
Result_Dir=file.path(Home_Dir,"Comparison")
dir.create(Result_Dir)
#
##
Sub_dir=c("p230","p210","p185")
N1_kept<-rep(0,length(Sub_dir))
N2_kept<-rep(0,length(Sub_dir))
for (i in 1:length(Sub_dir)) {
	filein=file.path(Home_Dir,Sub_dir[i],"Ranked_Result.RData")
	load(filein)
	N1_kept[i]<-length(Result_List$Interaction_Data$GeneID1_kept)
	N2_kept[i]<-length(Result_List$Interaction_Data$GeneID2_kept)
}
PPI_data<-data.frame(N1_kept, N2_kept)  
#
## plot stacked histogram
outputfile=file.path(Result_Dir,"Number_PPIs.tiff")
tiff(outputfile,width = 510, height = 510, units = "px")
par(xpd=T, mar=par()$mar+c(0,0,0,4))
barplot(t(PPI_data),  ylab="Number of PPIs", 
   col=c(2,4), space=0.1, las=1,ylim=c(0,310), names.arg=Sub_dir
   ,cex=1.4,cex.axis=1.4,cex.lab=1.4,font.lab=2.3)	# font.axis=2  
legend(2.2, 351, legend = c("inherited from BCR","inherited from ABL1"), names(PPI_data), cex=1.4, fill=c(2,4),text.font=2.3)
dev.off()



### Figure 3b: correlation of pathway prediction of different BCR-ABL1 variants (fgsea;color)
# output file: "CorrealtionPathway_Variants_v1.tiff" in the Comparison folder
#
##
Sub_dir=c("p230","p210","p185")
for (i in 1:length(Sub_dir)) {
	## load ranked gene list
	fileoutput2="BCR-ABL1_GSEA_Result.xls"
	fileoutput2=file.path(Home_Dir,Sub_dir[i],fileoutput2)
	Pathway_Result<- read.table(fileoutput2,sep="\t",header=TRUE, quote="") # not select column 
	
	##
	Pathways<-as.character(Pathway_Result$pathway)
	NES<-as.numeric(Pathway_Result$NES)
	ES<-as.numeric(Pathway_Result$ES)
	qval<-as.numeric(Pathway_Result$padj)
	pval<-as.numeric(Pathway_Result$pval)
	#score2<-NES
	score2<-pval
	if (i==1)
	{
		Pathways_All<-Pathways
		pval_M<-matrix(0,nrow=length(Pathways_All), ncol=length(Sub_dir))
		colnames(pval_M)<-Sub_dir
		#rownames(pval_M)<-Pathways_All
	}
	
	## assign p-value
	Score<-rep(0,length(Pathways_All))
	loc<-match(Pathways_All,Pathways)
	pos<-which(loc>0)
	loc<-loc[pos]
	Score[pos]<-score2[loc]
	pval_M[,i]<-Score
}
#
## Correaltion plot
fileoutput3=file.path(Result_Dir,"CorrealtionPathway_Variants_v1.tiff")
tiff(filename = fileoutput3, ,width = 410, height = 410, units = "px");
#par(mar=c(0.2,0.2,0.2,0.2))
library(corrplot)
M <- cor(pval_M, method="spearman")
corrplot.mixed(M,number.cex=1.1,cl.cex = 1.1)
dev.off()


### Figure 3c & Figure 2S: Compare different ROC for different fusion variants
# output file (Figure 3c): "AUC_Comparison_CancerPathway_v1.tiff" in the Comparison folder
# output files (Figure 2S): "AUC_Comparison_BCR-ABL1_v1.tiff", "AUC_Comparison_CML_v1.tiff"
#		, "AUC_Comparison_Drug-Targets_v1.tiff", "AUC_Comparison_NFAT_v1.tiff"
##
Sub_dir=c("p230","p210","p185")
Col<-c(1,2,3)
# Names of behcnmark gene sets
Ben_Types=c("Bechmark: CML Associated Genes",
	"Bechmark: Genes Co-Cited with BCR-ABL1","Bechmark: Wnt/Ca+/NFAT Pathway Genes",
	"Bechmark: Cancer Pathway Genes",
	"Bechmark: Target Genes of Drugs")
B_Types=c("CML","BCR-ABL1","NFAT","CancerPathway","Drug-Targets")
#
##
for (tt in 1:length(B_Types)) {
	##
	LEG<-rep("",length(Sub_dir))
	COL<-rep(0,length(Sub_dir))
	outputfile=file.path(Result_Dir,paste("AUC_Comparison_",B_Types[tt],"_v1.tiff",sep=""))
	tiff(outputfile,width = 510, height = 510, units = "px")
	par(mar=c(4,4.2,2.0,0.5))
	for (i in 1:length(Sub_dir)) {
		filein=file.path(Home_Dir,Sub_dir[i],"ROC.RData")
		load(filein)
		TPR<-ROC_M[[tt]]$TPR
		FPR<-ROC_M[[tt]]$FPR
		auc<-ROC_M[[tt]]$AUC
		#
		## legend
		LEG[i]<-paste(Sub_dir[i],"(AUC=",auc,")")
		#
		##
		if (i>1) {par(new = TRUE)}
		plot(FPR,TPR,xlab="FPR",ylab="TPR",col = i,lwd=3,"l",cex=1.45,cex.axis=1.45,cex.lab=1.45,font.lab=2.3)
	}
	title(main = Ben_Types[tt],cex.main=1.5)
	legend(0.49,0.19, legend = LEG, lty = 1,lwd=3,col=Col,cex=1.4,text.font=2)	# for 4 groups
	#title(main = list(Ben_Types[tt]),cex.main = 1.3, font.main= 2)
	dev.off()
}



###########################################################################################
### Figure S1: Prediction evaluation of p210 BCR-ABL1 using data-driven gene signatures 
###########################################################################################
# outputfile: "AUC_Different_GeneSignatures_v1.tiff" in the p210 folder
#
##
setwd(Home_Dir)
#
## Load result of p210
Sub_dir=c("p210")
filein=file.path(Home_Dir,Sub_dir,"Ranked_Result.RData")
load(filein)
association_socre<-as.numeric(Result_List$Ranked_Result$Score)
Gene_Ranked<-Result_List$Ranked_Result$GeneSymbol_ranked
ranks<-(1:length(Gene_Ranked))
#
## Klein_2006 data
load(file.path(system.file("data", package = "FusionPathway"), "GSEA_ChemicalPurtubation_Data.RData"))
PathwayNames<-names(Pathway_M)
pos<-which(PathwayNames=="KLEIN_TARGETS_OF_BCR_ABL1_FUSION")
PathwayNames<-PathwayNames[pos]
Bench_Genes<-unique(Pathway_M[[pos[1]]])
loc<-match(Gene_Ranked,Bench_Genes)
outcome=rep(0,length(Gene_Ranked))
outcome[which(loc>0)]=1
ROC<-ROCF(outcome,association_socre)
auc<-round(ROC$AUC,3)
print(paste("AUC using BCR-ABL1 Signature",auc,sep="="))
ROC_M<-NULL
ROC_M$Klein_2006$TPR<-ROC$TPR
ROC_M$Klein_2006$FPR<-ROC$FPR
ROC_M$Klein_2006$AUC<-auc
#
## Diaz-Blanco_2007 data
load(file.path(system.file("data", package = "FusionPathway"), "GSEA_ChemicalPurtubation_Data.RData"))
PathwayNames<-names(Pathway_M)
pos<-grep("DIAZ_CHRONIC_MEYLOGENOUS_LEUKEMIA_",PathwayNames)
PathwayNames<-PathwayNames[pos]
Bench_Genes<-unique(c(Pathway_M[[pos[1]]],Pathway_M[[pos[2]]]))
loc<-match(Gene_Ranked,Bench_Genes)
outcome=rep(0,length(Gene_Ranked))
outcome[which(loc>0)]=1
ROC<-ROCF(outcome,association_socre)
auc<-round(ROC$AUC,3)
print(paste("AUC using CML Signature",auc,sep="="))
ROC_M$DiazBlanco_2007$TPR<-ROC$TPR
ROC_M$DiazBlanco_2007$FPR<-ROC$FPR
ROC_M$DiazBlanco_2007$AUC<-auc
#
## Ray_2004 data
load(file.path(system.file("data", package = "FusionPathway"), "GSEA_ChemicalPurtubation_Data.RData"))
PathwayNames<-names(Pathway_M)
pos<-grep("RAY_TARGETS_OF_P210_BCR_ABL_FUSION_",PathwayNames)
PathwayNames<-PathwayNames[pos]
Bench_Genes<-unique(c(Pathway_M[[pos[1]]],Pathway_M[[pos[2]]]))
loc<-match(Gene_Ranked,Bench_Genes)
outcome=rep(0,length(Gene_Ranked))
outcome[which(loc>0)]=1
ROC<-ROCF(outcome,association_socre)
auc<-round(ROC$AUC,3)
print(paste("AUC using CML Signature",auc,sep="="))
ROC_M$Ray_2004$TPR<-ROC$TPR
ROC_M$Ray_2004$FPR<-ROC$FPR
ROC_M$Ray_2004$AUC<-auc
#
## Save ROC
output="./p210/ROC_GeneSignatures.RData"
save(file=output,ROC_M)
#
## ROC curves
outputfile<-file.path(Home_Dir,Sub_dir,"AUC_Different_GeneSignatures_v1.tiff")
tiff(outputfile,width = 510, height = 510, units = "px")
par(mar=c(4,4.2,0.5,0.5))
pp<-plot(ROC_M$Klein_2006$FPR,ROC_M$Klein_2006$TPR,xlab="FPR",ylab="TPR",col = "red",lwd=3,"l",
	cex=1.4,cex.axis=1.4,cex.lab=1.4,font.lab=2)
points(ROC_M$DiazBlanco_2007$FPR,ROC_M$DiazBlanco_2007$TPR,type="l",col="blue",lwd=3)
points(ROC_M$Ray_2004$FPR,ROC_M$Ray_2004$TPR,type="l",col="black",lwd=3)
LabelS<-c(paste("Klein_2006 (AUC=",ROC_M$Klein_2006$AUC,")",sep=""),
	paste("Diaz-Blanco_2007 (AUC=",ROC_M$DiazBlanco_2007$AUC,")",sep=""),
	paste("Ray_2004 (AUC=",ROC_M$Ray_2004$AUC,")",sep="")
	)
legend(0.32,0.165,col=c("red","blue","black"),
	legend=LabelS,lty=rep(1,5),lwd=3,cex=1.3)
dev.off()
 


#########################################################################
### Table S3: Combine array data and rank genes to select top genes 
#########################################################################
# output file: "ArrayData_GSE24493_SelectedGenes.xls" in the p210 folder
#
Work_Dir=file.path(Home_Dir,"p210")
setwd(Work_Dir)
#
## gene expression GSEA (GSE24493)
DataM1 <- read.delim(file.path(system.file("testData", package = "FusionPathway"), "ArrayData_GSE24493_Duy.2011.txt"), sep = "\t", header = TRUE, as.is = TRUE)
Gene.symbol<-as.character(DataM1$Gene.symbol)
# delete no symbols
pos<-which(Gene.symbol!="")
DataM1 <- DataM1[pos,]
# select diff expressed genes
P.Value<-as.numeric(as.character(DataM1$P.Value))
logFC<-as.numeric(as.character(DataM1$logFC))
pos1<-which(P.Value<=0.05)
pos2<-which(logFC>=0.5)
pos<-intersect(pos1,pos2)
DataM1 <- DataM1[pos,]
Gene.symbol<-as.character(DataM1$Gene.symbol)
P.Value<-as.numeric(as.character(DataM1$P.Value))
logFC<-as.numeric(as.character(DataM1$logFC))
#
## Load ranked genes
filein="BCR.ABL1_RankedGenes.rnk"
DataM2 <- read.table(filein,sep="\t",header=F, quote="")
Ranked_Genes<-as.character(DataM2[[1]])
L<-round(length(Ranked_Genes)*0.1)
Selected_Ranked_Genes<-Ranked_Genes[1:L]
#
## combine
loc<-match(toupper(Gene.symbol),toupper(Selected_Ranked_Genes))
pos<-which(loc>0)
Ranks<-loc[pos]
DataM1b <- DataM1[pos,]
DataM1b$ranks<-Ranks
outputfile1="ArrayData_GSE24493_SelectedGenes.xls"
write.table(DataM1b, outputfile1, sep = "\t",row.names = FALSE, quote=FALSE)
#
## Map to drugs
Selected_Genes<-unique(Gene.symbol[pos])
Ranks<-match(Selected_Genes,Selected_Ranked_Genes)
load(file.path(system.file("data", package = "FusionPathway"), "DrugDatabase_DGIdb.RData"))
#data("DrugDatabase_DGIdb")
DataM<-DrugDatabase_DGIdb
Drug_GeneNames<-as.character(DataM$entrez_gene_symbol)
#	
## Mapping to drug-target information
loc<-match(Drug_GeneNames,Selected_Genes)
pos<-which(loc>0)
loc<-loc[pos]
DataM<-cbind(DataM[pos,],Ranks=Ranks[loc],Rank_Percentage=Ranks[loc]*100/length(Ranked_Genes))
# Write
fileoutput1="ArrayData_GSE24493__MappedDrugs_DGIdb.xls"
write.table(DataM, fileoutput1, sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)















