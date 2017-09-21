###########################################################################################
### Example: EWS-FLI1 Fusion
###########################################################################################
#
rm(list=ls())	# clear the worksapce
#
Home_Dir<-getwd()
library("FusionPathway")


## Input Data: Gene ID, protein domain of two fusion parental genes
Gene1="EWSR1"
GeneID1=2130
PFAM1=c("PF00641","PF00076")
Gene2="FLI1"
GeneID2=2313
PFAM2=c("PF00178","PF02198")
## Domain composition of the fusion
PFAM1_lost<-c("PF00641","PF00076")
PFAM1_kept<-setdiff(PFAM1,PFAM1_lost)
PFAM2_lost<-"PF02198"
PFAM2_kept<-setdiff(PFAM2,PFAM2_lost)

## Result folder
Home_Dir<-file.path(Home_Dir,paste(Gene1,Gene2,sep="_"))
dir.create(Home_Dir)
setwd(Home_Dir)

## run FusionPathway
DomainData<-list(pfam1_kept=PFAM1_kept,pfam1_lost=PFAM1_lost,pfam2_kept=PFAM2_kept,pfam2_lost=PFAM2_lost)
GeneData<-data.frame(Gene1=Gene1,Gene2=Gene2,GeneID1=GeneID1,GeneID2=GeneID2)
print(paste("Fusion ",paste(Gene1,Gene2,sep="-"),sep=""))
Result_List<-FusionPathway(GeneData,DomainData)
save(Result_List,file="Ranked_Result.RData")

## Evaluation using different benchmark gene sets
Score<-as.numeric(Result_List$Ranked_Result$Score)
GeneID_ranked<-as.numeric(Result_List$Ranked_Result$GeneID_ranked)
load(file.path(system.file("data", package = "FusionPathway"), "BenchmarkSets_Ewing.RData"))
for (i in 1:5) {
	##
	if (i==1) {
		BenchmarkName="Ewing_Genes"
		PubmedN<-BenchmarkSets_Ewing[["Ewing_Genes"]]$PubmedN
		GeneIDs<-BenchmarkSets_Ewing[["Ewing_Genes"]]$GeneIDs
		Ncitation_threshold=2
		#
		BenchID_pos<-GeneIDs[which(PubmedN>=Ncitation_threshold)]
		#
		ROC_M<-NULL
	} else if (i==2) {
		BenchmarkName="EWS_FLI1_Genes"
		PubmedN<-BenchmarkSets_Ewing[["EWS_FLI1_Genes"]]$PubmedN
		GeneIDs<-BenchmarkSets_Ewing[["EWS_FLI1_Genes"]]$GeneIDs
		Ncitation_threshold=2
		BenchID_pos<-GeneIDs[which(PubmedN>=Ncitation_threshold)]
	} else if (i==3) {
		BenchmarkName="CancerPathway_Genes"
		BenchID_pos<-BenchmarkSets_Ewing[["CancerPathway_Genes"]]$GeneIDs
	} else if (i==4) {
		BenchmarkName="Drug_TargetGenes"
		BenchID_pos<-BenchmarkSets_Ewing[["Drug_Targets"]]$GeneIDs
	} else if (i==5) {
		BenchmarkName="Drug_Screening"
		BenchID_pos<-BenchmarkSets_Ewing[["Drug_Screening"]]$GeneIDs
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
	
	## ROC curves 
	outputfile=paste("AUC_Benchmark.",BenchmarkName,".tiff",sep="")
	tiff(outputfile,width = 510, height = 510, units = "px")
	pp<-plot(ROC$FPR,ROC$TPR,xlab="FPR",ylab="TPR",col = "red",lwd=3,"l",
		cex=1.45,cex.axis=1.45,cex.lab=1.45,font.lab=2.3)
	LabelS<-paste(BenchmarkName,"(AUC=",round(ROC$AUC,3),")",sep="")
	legend("bottomright",col="red",legend=LabelS,lty=rep(1,5),lwd=3,cex=1.4,text.font=1.5)
	dev.off()
}
## Save ROC
output="ROC.RData"
save(file=output,ROC_M)



###########################################################################################
### Figure 4a: Combine ROC plots for different benchmarks (color setting)
###########################################################################################
### output file: "AUC_Different_Benchmarks_v1.tiff" 
###
setwd(Home_Dir)
#
### ROC curve
load("ROC.RData")
outputfile="AUC_Different_Benchmarks_v1.tiff"
tiff(outputfile,width = 510, height = 510, units = "px")
#outputfile="AUC_Different_Benchmarks.pdf"
#pdf(outputfile,7.2,7.2)
par(mar=c(4,4.2,0.5,0.5))
pp<-plot(ROC_M$Ewing_Genes$FPR,ROC_M$Ewing_Genes$TPR,xlab="FPR",ylab="TPR",col = "red",lwd=3,"l",
	cex=1.45,cex.axis=1.45,cex.lab=1.45,font.lab=2.3)
points(ROC_M$EWS_FLI1_Genes$FPR,ROC_M$EWS_FLI1_Genes$TPR,type="l",col="black",lwd=3)
points(ROC_M$CancerPathway_Genes$FPR,ROC_M$CancerPathway_Genes$TPR,type="l",col="blue",lwd=3)
points(ROC_M$Drug_TargetGenes$FPR,ROC_M$Drug_TargetGenes$TPR,type="l",col="green",lwd=3)
points(ROC_M$Drug_Screening$FPR,ROC_M$Drug_Screening$TPR,type="l",col="purple",lwd=3)
LabelS<-c(paste("Ewing Genes(AUC=",ROC_M$Ewing_Genes$AUC,")",sep=""),
	paste("EWS-FLI1 Genes(AUC=",ROC_M$EWS_FLI1_Genes$AUC,")",sep=""),
	paste("Cancer Pathway Genes(AUC=",ROC_M$CancerPathway_Genes$AUC,")",sep=""),
	paste("Drug Targets(AUC=",ROC_M$Drug_TargetGenes$AUC,")",sep=""),
	paste("Drug Screening(AUC=",ROC_M$Drug_Screening$AUC,")",sep=""))
legend(0.16,0.28,col=c("red","black","blue","green","purple"),
	legend=LabelS,lty=rep(1,5),lwd=3,cex=1.4,text.font=2)
dev.off()



###########################################################################################
### Figure S3: Bar plot of some pathways that are associated with EWS-FLI1 in the prediction
###########################################################################################
#
### Combine GSEA pathway pavlues of GSEA association & GSEA deregulation
setwd(Home_Dir)
#
##  Gene expression data (GSE27524; Bilke et al., 2013)
Array_name<-"GSE27524"
DataM1 <- read.delim(file.path(system.file("testData", package = "FusionPathway"), "ArrayData_GSE27524_siEWS_Bilke.2013.txt"), sep = "\t", header = TRUE, as.is = TRUE)
Gene.symbol<-as.character(DataM1$Gene.symbol)
# delete no symbols
pos<-which(Gene.symbol!="")
DataM1 <- DataM1[pos,]
Gene.symbol<-as.character(DataM1$Gene.symbol)
P.Value<-as.numeric(as.character(DataM1$P.Value))
logFC<-as.numeric(as.character(DataM1$logFC))
# sort genes
AA<-sort(P.Value,decreasing =F, index.return = TRUE)
indx<-as.numeric(AA$ix)
P.Value<-P.Value[indx]
Gene.symbol<-Gene.symbol[indx]		
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
fileoutput1=paste(Array_name,"_fGSEA_Result.xls",sep="")
write.table(Result_GSEA,fileoutput1, sep = "\t",append=FALSE,row.names = FALSE,col.names = TRUE,quote=FALSE)
#
## Pathway GSEA association & GSEA deregulation
filein1="EWSR1-FLI1_GSEA_Result.xls"
PathwayM1 <- read.table(filein1,sep="\t",header=T, quote="")
Pathways1<-as.character(PathwayM1$pathway)
filein2=paste(Array_name,"_fGSEA_Result.xls",sep="")
PathwayM2 <- read.table(filein2,sep="\t",header=T, quote="")
Pathways2<-as.character(PathwayM2$pathway)
loc<-match(Pathways1,Pathways2)
pos<-which(loc>0)
loc<-loc[pos]
PathwayM1<- PathwayM1[pos,]
PathwayM2<- PathwayM2[loc,]
#
## Combine pvalue
outfile="EWSR1-FLI1_GSEA_CombinedResult.xls"
ResultM<-Combine_Pathway_Pvalues(PathwayM1,PathwayM2,threshold_ind_pvalue=0.05,outfile)


### Figure S3: Bar Chart of p-value of selected associated pathways for prediction of EWS-FLI1
# output file: "Enrichment_AssocaitedPathways.tiff"
#
## Selected Pathway
Pathway_Selected=c("KEGG_CELL_CYCLE","REACTOME_SIGNALING_BY_WNT",
	"BIOCARTA_IGF1R_PATHWAY",
	"CHROMATIN_MODIFICATION","CHROMATIN_REMODELING",
	"KEGG_P53_SIGNALING_PATHWAY","APOPTOSIS_GO",
	"KEGG_SPLICEOSOME","REACTOME_SIGNALING_BY_PDGF",
	"KEGG_PATHWAYS_IN_CANCER","REACTOME_DNA_REPAIR",
	"BIOCARTA_MAPK_PATHWAY",
	"REACTOME_PI3K_AKT_ACTIVATION")	
#
## GSEA results (GSEA association, GSEA deregulation, and combined result) of the selected pathways
file_AsscoatedPathways<-"EWSR1-FLI1_GSEA_CombinedResult.xls"
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
#outputfile="Enrichment_AssocaitedPathways.pdf"
#pdf(outputfile,10,5)
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



###########################################################################################
### Figure S4: Prediction evaluation of EWS-FLI1 using gene signatures that are associated with EWS-FLI1 or Ewing’s sarcoma
###########################################################################################
# output file: "AUC_Different_GeneSignatures_v1.tiff"
#
setwd(Home_Dir)
#
##
filein="Ranked_Result.RData"
load(filein)
association_socre<-as.numeric(Result_List$Ranked_Result$Score)
Gene_Ranked<-Result_List$Ranked_Result$GeneSymbol_ranked
ranks<-(1:length(Gene_Ranked))
#
## Kinsey_2006 
load(file.path(system.file("data", package = "FusionPathway"), "GSEA_ChemicalPurtubation_Data.RData"))
PathwayNames<-names(Pathway_M)
pos<-grep("KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_",PathwayNames)
PathwayNames<-PathwayNames[pos]
Bench_Genes<-unique(c(Pathway_M[[pos[1]]],Pathway_M[[pos[2]]]))
loc<-match(Gene_Ranked,Bench_Genes)
outcome=rep(0,length(Gene_Ranked))
outcome[which(loc>0)]=1
ROC<-ROCF(outcome,association_socre)
auc<-round(ROC$AUC,3)
print(paste("AUC using Ewing Gene Signature",auc,sep="="))
ROC_M$Kinsey_2006$TPR<-ROC$TPR
ROC_M$Kinsey_2006$FPR<-ROC$FPR
ROC_M$Kinsey_2006$AUC<-auc
#
## Rorie_2004 Signature
load(file.path(system.file("data", package = "FusionPathway"), "GSEA_ChemicalPurtubation_Data.RData"))
PathwayNames<-names(Pathway_M)
pos<-grep("RORIE_TARGETS_OF_EWSR1_FLI1_FUSION_",PathwayNames)
PathwayNames<-PathwayNames[pos]
Bench_Genes<-unique(c(Pathway_M[[pos[1]]],Pathway_M[[pos[2]]]))
loc<-match(Gene_Ranked,Bench_Genes)
outcome=rep(0,length(Gene_Ranked))
outcome[which(loc>0)]=1
ROC<-ROCF(outcome,association_socre)
auc<-round(ROC$AUC,3)
print(paste("AUC using Ewing Gene Signature",auc,sep="="))
ROC_M$Rorie_2004$TPR<-ROC$TPR
ROC_M$Rorie_2004$FPR<-ROC$FPR
ROC_M$Rorie_2004$AUC<-auc
#
## Siligan_2005 Signature
load(file.path(system.file("data", package = "FusionPathway"), "GSEA_ChemicalPurtubation_Data.RData"))
PathwayNames<-names(Pathway_M)
pos<-grep("SILIGAN_TARGETS_OF_EWS_FLI1_FUSION_",PathwayNames)
PathwayNames<-PathwayNames[pos]
Bench_Genes<-unique(c(Pathway_M[[pos[1]]],Pathway_M[[pos[2]]]))
loc<-match(Gene_Ranked,Bench_Genes)
outcome=rep(0,length(Gene_Ranked))
outcome[which(loc>0)]=1
ROC<-ROCF(outcome,association_socre)
auc<-round(ROC$AUC,3)
print(paste("AUC using Ewing Gene Signature",auc,sep="="))
ROC_M$Siligan_2005$TPR<-ROC$TPR
ROC_M$Siligan_2005$FPR<-ROC$FPR
ROC_M$Siligan_2005$AUC<-auc
#
## Hu-Lieskovan_2005
load(file.path(system.file("data", package = "FusionPathway"), "GSEA_ChemicalPurtubation_Data.RData"))
PathwayNames<-names(Pathway_M)
pos<-which(PathwayNames=="ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION")
PathwayNames<-PathwayNames[pos]
Bench_Genes<-unique(Pathway_M[[pos[1]]])
loc<-match(Gene_Ranked,Bench_Genes)
outcome=rep(0,length(Gene_Ranked))
outcome[which(loc>0)]=1
ROC<-ROCF(outcome,association_socre)
auc<-round(ROC$AUC,3)
print(paste("AUC using Trabectedin Ewing Signature",auc,sep="="))
ROC_M$HuLieskovan_2005$TPR<-ROC$TPR
ROC_M$HuLieskovan_2005$FPR<-ROC$FPR
ROC_M$HuLieskovan_2005$AUC<-auc
#
## Save ROC
output="ROC_GeneSignatures.RData"
save(file=output,ROC_M)
#
##
outputfile<-"AUC_Different_GeneSignatures_v1.tiff"
tiff(outputfile,width = 610, height = 610, units = "px")
#outputfile<-"AUC_Different_GeneSignatures.pdf"
#pdf(outputfile,7.2,7.2)
pp<-plot(ROC_M$Siligan_2005$FPR,ROC_M$Siligan_2005$TPR,xlab="FPR",ylab="TPR",col = "red",lwd=3,"l",
	cex=1.4,cex.axis=1.4,cex.lab=1.4,font.lab=2)
points(ROC_M$Rorie_2004$FPR,ROC_M$Rorie_2004$TPR,type="l",col="purple",lwd=3)
points(ROC_M$HuLieskovan_2005$FPR,ROC_M$HuLieskovan_2005$TPR,type="l",col="cyan",lwd=3)
points(ROC_M$Kinsey_2006$FPR,ROC_M$Kinsey_2006$TPR,type="l",col="darkorange",lwd=3)
LabelS<-c(paste("Siligan 2005 (AUC=",ROC_M$Siligan_2005$AUC,")",sep=""),
	paste("Rorie_2004 (AUC=",ROC_M$Rorie_2004$AUC,")",sep=""),
	paste("Hu-Lieskovan_2005 (AUC=",ROC_M$HuLieskovan_2005$AUC,")",sep=""),
	paste("Kinsey 2006 (AUC=",ROC_M$Kinsey_2006$AUC,")",sep="")
	)
legend(0.41,0.19,col=c("red","purple","cyan","darkorange","blue","black","brown","green"),
	legend=LabelS,lty=rep(1,5),lwd=3,cex=1.3)
dev.off()
 
 
 
##########################################################################################
# Compare predicitons with Drug Screening Results
###########################################################################################
setwd(Home_Dir)

## Prediction of gene ranks
filein="EWSR1.FLI1_RankedGenes.rnk"
DataM<- read.table(filein,header = F,sep="\t") # not select column 1
Gene_Ranked<-as.character(DataM[[1]])
association_socre<-as.character(DataM[[2]])
All_GeneRanks<-(1:length(Gene_Ranked))
#
## Sensitive compounds
DataM1 <- read.delim(file.path(system.file("testData", package = "FusionPathway"), "Ewing_SensitiveDrugs.xls"), sep = "\t", header = TRUE, as.is = TRUE)
Sensitive_Drugs<-unique(as.character(DataM1[[1]]))
pos<-which(as.character(DataM1[[3]])!="-")
DataM1b<-DataM1[pos,]
Sensitive_Drugs2<-unique(as.character(DataM1b[[1]]))
#
## Map gene ranks to targets of sensitive compounds
loc<-match(as.character(DataM1[[3]]),Gene_Ranked)
pos<-which(loc>0)
loc<-loc[pos]
Ranks<-All_GeneRanks[loc]
#fileoutput="Ewing_PassDrugScreen_Mapped.xls"
#write.table(cbind(DataM1,Ranks),fileoutput, sep = "\t",append=FALSE,row.names = FALSE,col.names = T,quote=FALSE)
#
## Select drugs with top 10% ranked targets
rank_threshold=0.10
pos<-which(Ranks<=(length(All_GeneRanks)*rank_threshold))
TopRanked_Drugs<-unique(as.character(DataM1[pos,]$Drug))
TopRanked_Target<-unique(as.character(DataM1[pos,]$Target))
#
## Print
print(paste("Number of sensitive drugs",length(Sensitive_Drugs),sep="="))	
print(paste("Number of sensitive drugs with known targets",length(Sensitive_Drugs2),sep="="))	
print(paste("Number of sensitive drugs with top 10% targets",length(TopRanked_Drugs),sep="="))	
#print(paste("Top ranked target genes",TopRanked_Target,sep="="))
Percentage<-round(100*length(TopRanked_Drugs)/length(Sensitive_Drugs2),2)	
print(paste("% of sensitive drugs with top 10% targets = ",Percentage,"%",sep=""))	



#########################################################################
### Combine array data and rank genes to select top genes 
#########################################################################
# output file: "ArrayData_GSE24493_SelectedGenes.xls" in the p210 folder
#
setwd(Home_Dir)
#
## gene expression GSEA (GSE24493)
DataM1 <- read.delim(file.path(system.file("testData", package = "FusionPathway"), "ArrayData_GSE27524_siEWS_Bilke.2013.txt"), sep = "\t", header = TRUE, as.is = TRUE)
Gene.symbol<-as.character(DataM1$Gene.symbol)
# delete no symbols
pos<-which(Gene.symbol!="")
DataM1 <- DataM1[pos,]
# select diff expressed genes
P.Value<-as.numeric(as.character(DataM1$P.Value))
logFC<-as.numeric(as.character(DataM1$logFC))
pos1<-which(P.Value<=0.05)
pos2<-which(logFC<=(-1)*0.5)
pos<-intersect(pos1,pos2)
DataM1 <- DataM1[pos,]
Gene.symbol<-as.character(DataM1$Gene.symbol)
P.Value<-as.numeric(as.character(DataM1$P.Value))
logFC<-as.numeric(as.character(DataM1$logFC))
#
## Load ranked genes
filein="EWSR1.FLI1_RankedGenes.rnk"
DataM2 <- read.table(filein,sep="\t",header=F, quote="")
Ranked_Genes<-as.character(DataM2[[1]])
L<-round(length(Ranked_Genes)*0.1)	# select top 10% genes
Selected_Ranked_Genes<-Ranked_Genes[1:L]
#
## combine
loc<-match(toupper(Gene.symbol),toupper(Selected_Ranked_Genes))
pos<-which(loc>0)
Ranks<-loc[pos]
DataM1b <- DataM1[pos,]
DataM1b$ranks<-Ranks
outputfile1="ArrayData_GSE27524_SelectedGenes.xls"
write.table(DataM1b, outputfile1, sep = "\t",row.names = FALSE, quote=FALSE)
#
## Map to drugs
Selected_Genes<-unique(Gene.symbol[pos])
Ranks<-match(Selected_Genes,Selected_Ranked_Genes)
load(file.path(system.file("data", package = "FusionPathway"), "DrugDatabase_DGIdb.RData"))
DataM<-DrugDatabase_DGIdb
Drug_GeneNames<-as.character(DataM$entrez_gene_symbol)
#	
## Mapping to drug-target information
loc<-match(Drug_GeneNames,Selected_Genes)
pos<-which(loc>0)
loc<-loc[pos]
DataM<-cbind(DataM[pos,],Ranks=Ranks[loc],Rank_Percentage=Ranks[loc]*100/length(Ranked_Genes))
# Write
fileoutput1="ArrayData_GSE27524_MappedDrugs_DGIdb.xls"
write.table(DataM, fileoutput1, sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
















