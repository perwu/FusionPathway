###########################################################################################
### Example: FUS-DDIT3 Fusion
###########################################################################################
#
rm(list=ls())	# clear the worksapce
#
Home_Dir<-getwd()
library("FusionPathway")

## Input Data: Gene ID, protein domain of two fusion parental genes
Gene1="FUS"
GeneID1=2521
PFAM1=c("PF00076","PF00641")
Gene2="DDIT3"
GeneID2=1649
PFAM2=c("PF07716")
## Domain composition of the fusion
PFAM1_lost<-""	
PFAM1_kept<-setdiff(PFAM1,PFAM1_lost)
PFAM2_lost<-""
PFAM2_kept<-setdiff(PFAM2,PFAM2_lost)

## Result folder
Home_Dir<-file.path(Home_Dir,paste(Gene1,Gene2,sep="_"))
dir.create(Home_Dir)
setwd(Home_Dir)

## run FusionPathway
DomainData<-list(pfam1_kept=PFAM1_kept,pfam1_lost=PFAM1_lost,pfam2_kept=PFAM2_kept,pfam2_lost=PFAM2_lost)
GeneData<-list(Gene1=Gene1,Gene2=Gene2,GeneID1=GeneID1,GeneID2=GeneID2)
print(paste("Fusion ",paste(Gene1,Gene2,sep="-"),sep=""))
Result_List<-FusionPathway(GeneData,DomainData)
save(Result_List,file="Ranked_Result.RData")

## Evaluation using different benchmark gene sets
Score<-as.numeric(Result_List$Ranked_Result$Score)
GeneID_ranked<-as.numeric(Result_List$Ranked_Result$GeneID_ranked)
load(file.path(system.file("data", package = "FusionPathway"), "BenchmarkSets_Myxoid.RData"))
for (i in 1:4) {
	##
	if (i==1) {
		BenchmarkName="Myxoid_Genes"
		PubmedN<-BenchmarkSets_Myxoid[["Myxoid_Genes"]]$PubmedN
		GeneIDs<-BenchmarkSets_Myxoid[["Myxoid_Genes"]]$GeneIDs
		Ncitation_threshold=1
		#
		BenchID_pos<-GeneIDs[which(PubmedN>=Ncitation_threshold)]
		#
		ROC_M<-NULL
	} else if (i==2) {
		BenchmarkName="FUS_DDIT3_Genes"
		PubmedN<-BenchmarkSets_Myxoid[["FUS_DDIT3_Genes"]]$PubmedN
		GeneIDs<-BenchmarkSets_Myxoid[["FUS_DDIT3_Genes"]]$GeneIDs
		Ncitation_threshold=1
		BenchID_pos<-GeneIDs[which(PubmedN>=Ncitation_threshold)]
	} else if (i==3) {
		BenchmarkName="CancerPathway_Genes"
		BenchID_pos<-BenchmarkSets_Myxoid[["CancerPathway_Genes"]]$GeneIDs
	} else if (i==4) {
		BenchmarkName="Drug_Screening"
		BenchID_pos<-BenchmarkSets_Myxoid[["Drug_Screening"]]$GeneIDs
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
### Figure 4b: Combine ROC plots for different benchmarks (color setting)
###########################################################################################
### output file: "AUC_Different_Benchmarks_v1.tiff" 
##
setwd(Home_Dir)
#
## ROC curve
load("ROC.RData")
outputfile="AUC_Different_Benchmarks_v1.tiff"
tiff(outputfile,width = 510, height = 510, units = "px")
#outputfile="AUC_Different_Benchmarks.pdf"
#pdf(outputfile,7.2,7.2)
par(mar=c(4,4.2,0.5,0.5))
pp<-plot(ROC_M$Myxoid_Genes$FPR,ROC_M$Myxoid_Genes$TPR,xlab="FPR",ylab="TPR",col = "red",lwd=3,"l",
	cex=1.45,cex.axis=1.45,cex.lab=1.45,font.lab=2.3)
points(ROC_M$FUS_DDIT3_Genes$FPR,ROC_M$FUS_DDIT3_Genes$TPR,type="l",col="black",lwd=3)
points(ROC_M$CancerPathway_Genes$FPR,ROC_M$CancerPathway_Genes$TPR,type="l",col="blue",lwd=3)
points(ROC_M$Drug_Screening$FPR,ROC_M$Drug_Screening$TPR,type="l",col="purple",lwd=3)
LabelS<-c(paste("Myxoid Genes(AUC=",ROC_M$Myxoid_Genes$AUC,")",sep=""),
	paste("FUS-DDIT3 Genes(AUC=",ROC_M$FUS_DDIT3_Genes$AUC,")",sep=""),
	paste("Cancer Pathway Genes(AUC=",ROC_M$CancerPathway_Genes$AUC,")",sep=""),
	paste("Drug Screening(AUC=",ROC_M$Drug_Screening$AUC,")",sep=""))
legend(0.16,0.235,col=c("red","black","blue","purple"),
	legend=LabelS,lty=rep(1,5),lwd=3,cex=1.40,text.font=2)
dev.off()



###########################################################################################
### Figure S5: Bar plot of some pathways that are associated with FUS-DDIT3 in the prediction
###########################################################################################
#
### Combine GSEA pathway pavlues of GSEA association & GSEA deregulation
setwd(Home_Dir)
#
##  Gene expression data (GSE33616; Oikawa  et al., 2012)
Array_name<-"GSE48030"
DataM <- read.delim(file.path(system.file("testData", package = "FusionPathway"), "GSE48030_platform_annotation.xls"), sep = "\t", header = TRUE, as.is = TRUE)
GENE_SYMBOL<-as.character(DataM$GENE_SYMBOL)
pos<-which(GENE_SYMBOL!="")
DataM <- DataM[pos,]
GENE_SYMBOL<-as.character(DataM$GENE_SYMBOL)
ID<-as.character(DataM$ID)
#	
DataM1 <- read.delim(file.path(system.file("testData", package = "FusionPathway"), "ArrayData_GSE48030_FUS-CHOP_Rodriguez.2013.txt"), sep = "\t", header = TRUE, as.is = TRUE)
ID2<-as.character(DataM1$ID)
loc<-match(ID2,ID)
pos<-which(loc>0)
loc<-loc[pos]
Gene.symbol<-GENE_SYMBOL[loc]
P.Value<-as.numeric(as.character(DataM1$P.Value))
logFC<-as.numeric(as.character(DataM1$logFC))
## sort genes
AA<-sort(P.Value,decreasing =F, index.return = TRUE)
indx<-as.numeric(AA$ix)
P.Value<-P.Value[indx]
Gene.symbol<-Gene.symbol[indx]		
#
## run fGSEA
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
filein1="FUS-DDIT3_GSEA_Result.xls"
PathwayM1 <- read.table(filein1,sep="\t",header=T, quote="")
Pathways1<-as.character(PathwayM1$pathway)
filein2=paste(Array_name,"_fGSEA_Result.xls",sep="")
PathwayM2 <- read.table(filein2,sep="\t",header=T, quote="")
Pathways2<-as.character(PathwayM2$pathway)
#
loc<-match(Pathways1,Pathways2)
pos<-which(loc>0)
loc<-loc[pos]
PathwayM1<- PathwayM1[pos,]
PathwayM2<- PathwayM2[loc,]
#
## Combine pvalue
outfile="FUS-DDIT3_GSEA_CombinedResult.xls"
ResultM<-Combine_Pathway_Pvalues(PathwayM1,PathwayM2,threshold_ind_pvalue=0.05,outfile)


### Figure S5: Bar Chart of p-value of selected associated pathways for prediction of EWS-FLI1
# output file: "Enrichment_AssocaitedPathways.tiff"
#
## Selected Pathway	
Pathway_Selected=c("KEGG_CELL_CYCLE","REACTOME_SIGNALING_BY_WNT","KEGG_WNT_SIGNALING_PATHWAY",
	"REACTOME_DNA_REPAIR",
	"KEGG_ADIPOCYTOKINE_SIGNALING_PATHWAY",
	"KEGG_P53_SIGNALING_PATHWAY","APOPTOSIS_GO",
	"KEGG_PATHWAYS_IN_CANCER","RESPONSE_TO_DNA_DAMAGE_STIMULUS",
	"KEGG_MAPK_SIGNALING_PATHWAY",
	"PID_PI3KCIPATHWAY",
	"I_KAPPAB_KINASE_NF_KAPPAB_CASCADE",
	"REGULATION_OF_CELL_DIFFERENTIATION","REGULATION_OF_DEVELOPMENTAL_PROCESS",
	"REACTOME_SIGNALING_BY_FGFR_IN_DISEASE")	
## GSEA results (GSEA association, GSEA deregulation, and combined result) of the selected pathways
file_AsscoatedPathways<-"FUS-DDIT3_GSEA_CombinedResult.xls"
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
# pathway name is too long
pos<-which(Pathway_Selected=="REACTOME_TRANSCRIPTIONAL_REGULATION_OF_WHITE_ADIPOCYTE_DIFFERENTIATION")
Pathway_Selected[pos]<-"REACTOME_WHITE_ADIPOCYTE_DIFFERENTIATION"
DataP<-data.frame(pathways=Pathway_Selected,pvalue1=pval1,pvalue2=pval2,combined_pval=combined_pval)
DataP$pathways= factor(DataP$pathways,Pathway_Selected) # change order
##
library(ggplot2)
library(gridExtra)
outputfile="Enrichment_AssocaitedPathways.tiff"
tiff(outputfile,width = 750, height = 300, units = "px")
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
grid.arrange(g1,g2,g3,ncol=3,widths=c(8.5/10,3.0/10,3.0/10))
dev.off()



##########################################################################################
### Figure S6: Evaluation of FUS-DDIT3 prediction using gene signatures that are associated with Trabectedin
###########################################################################################
# output file: "AUC_Different_Trabectedin-Signatures_v1.tiff"
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
## Martinez_2005
load(file.path(system.file("data", package = "FusionPathway"), "GSEA_ChemicalPurtubation_Data.RData"))
PathwayNames<-names(Pathway_M)
pos<-grep("MARTINEZ_RESPONSE_TO_TRABECTEDIN_",PathwayNames)
PathwayNames<-PathwayNames[pos]
Bench_Genes<-unique(c(Pathway_M[[pos[1]]],Pathway_M[[pos[2]]]))
loc<-match(Gene_Ranked,Bench_Genes)
outcome=rep(0,length(Gene_Ranked))
outcome[which(loc>0)]=1
ROC<-ROCF(outcome,association_socre)
auc<-round(ROC$AUC,3)
print(paste("AUC using Trabectedin Gene Signature",auc,sep="="))
ROC_M$Martinez_2005$TPR<-ROC$TPR
ROC_M$Martinez_2005$FPR<-ROC$FPR
ROC_M$Martinez_2005$AUC<-auc
#
## Marchini_2005
load(file.path(system.file("data", package = "FusionPathway"), "GSEA_ChemicalPurtubation_Data.RData"))
PathwayNames<-names(Pathway_M)
pos<-grep("MARCHINI_TRABECTEDIN_RESISTANCE",PathwayNames)
PathwayNames<-PathwayNames[pos]
Bench_Genes<-unique(c(Pathway_M[[pos[1]]],Pathway_M[[pos[2]]]))
loc<-match(Gene_Ranked,Bench_Genes)
outcome=rep(0,length(Gene_Ranked))
outcome[which(loc>0)]=1
ROC<-ROCF(outcome,association_socre)
auc<-round(ROC$AUC,3)
print(paste("AUC using Trabectedin Gene Signature",auc,sep="="))
ROC_M$Marchini_2005$TPR<-ROC$TPR
ROC_M$Marchini_2005$FPR<-ROC$FPR
ROC_M$Marchini_2005$AUC<-auc
#
## Martinez_2001
load(file.path(system.file("data", package = "FusionPathway"), "GSEA_ChemicalPurtubation_Data.RData"))
PathwayNames<-names(Pathway_M)
pos<-which(PathwayNames=="MARTINEZ_RESPONSE_TO_TRABECTEDIN")
PathwayNames<-PathwayNames[pos]
Bench_Genes<-unique(Pathway_M[[pos[1]]])
loc<-match(Gene_Ranked,Bench_Genes)
outcome=rep(0,length(Gene_Ranked))
outcome[which(loc>0)]=1
ROC<-ROCF(outcome,association_socre)
auc<-round(ROC$AUC,3)
print(paste("AUC using Trabectedin Gene Signature",auc,sep="="))
ROC_M$Martinez_2001$TPR<-ROC$TPR
ROC_M$Martinez_2001$FPR<-ROC$FPR
ROC_M$Martinez_2001$AUC<-auc
#
## Gajate_2002 Signature
load(file.path(system.file("data", package = "FusionPathway"), "GSEA_ChemicalPurtubation_Data.RData"))
PathwayNames<-names(Pathway_M)
pos<-grep("GAJATE_RESPONSE_TO_TRABECTEDIN",PathwayNames)
PathwayNames<-PathwayNames[pos]
Bench_Genes<-unique(c(Pathway_M[[pos[1]]],Pathway_M[[pos[2]]]))
loc<-match(Gene_Ranked,Bench_Genes)
outcome=rep(0,length(Gene_Ranked))
outcome[which(loc>0)]=1
ROC<-ROCF(outcome,association_socre)
auc<-round(ROC$AUC,3)
print(paste("AUC using Trabectedin Gene Signature",auc,sep="="))
ROC_M$Gajate_2002$TPR<-ROC$TPR
ROC_M$Gajate_2002$FPR<-ROC$FPR
ROC_M$Gajate_2002$AUC<-auc
#
## Save ROC
output="ROC_Trabectedin-GeneSignatures.RData"
save(file=output,ROC_M)
#
##
outputfile<-"AUC_Different_Trabectedin-Signatures_v1.tiff"
tiff(outputfile,width = 610, height = 610, units = "px")
#outputfile<-"AUC_Different_GeneSignatures.pdf"
#pdf(outputfile,7.2,7.2)
pp<-plot(ROC_M$Gajate_2002$FPR,ROC_M$Gajate_2002$TPR,xlab="FPR",ylab="TPR",col = "red",lwd=3,"l",
	cex=1.4,cex.axis=1.4,cex.lab=1.4,font.lab=2)
points(ROC_M$Martinez_2001$FPR,ROC_M$Martinez_2001$TPR,type="l",col="blue",lwd=3)
points(ROC_M$Marchini_2005$FPR,ROC_M$Marchini_2005$TPR,type="l",col="black",lwd=3)
points(ROC_M$Martinez_2005$FPR,ROC_M$Martinez_2005$TPR,type="l",col="purple",lwd=3)
LabelS<-c(paste("Gajate_2002(AUC=",ROC_M$Gajate_2002$AUC,")",sep=""),
	paste("Martinez_2001(AUC=",ROC_M$Martinez_2001$AUC,")",sep=""),
	paste("Marchini_2005(AUC=",ROC_M$Marchini_2005$AUC,")",sep=""),
	paste("Martinez_2005(AUC=",ROC_M$Martinez_2005$AUC,")",sep=""))
legend(0.49,0.19,col=c("red","blue","black","purple","brown"),
	legend=LabelS,lty=rep(1,5),lwd=3,cex=1.3)
dev.off()
 

 
##########################################################################################
# Compare predicitons with Drug Screening Results
###########################################################################################
setwd(Home_Dir)

### Prediction of gene ranks
filein="FUS.DDIT3_RankedGenes.rnk"
DataM<- read.table(filein,header = F,sep="\t") # not select column 1
Gene_Ranked<-as.character(DataM[[1]])
association_socre<-as.character(DataM[[2]])
All_GeneRanks<-(1:length(Gene_Ranked))
# 
## Sensitive compounds
DataM1 <- read.delim(file.path(system.file("testData", package = "FusionPathway"), "Myxoid_SensitiveDrugs.xls"), sep = "\t", header = TRUE, as.is = TRUE)
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
##  gene expression GSEA (GSE48030)
Array_name<-"GSE48030"
DataM <- read.delim(file.path(system.file("testData", package = "FusionPathway"), "GSE48030_platform_annotation.xls"), sep = "\t", header = TRUE, as.is = TRUE)
GENE_SYMBOL<-as.character(DataM$GENE_SYMBOL)
pos<-which(GENE_SYMBOL!="")
DataM <- DataM[pos,]
GENE_SYMBOL<-as.character(DataM$GENE_SYMBOL)
ID<-as.character(DataM$ID)
#	
DataM1 <- read.delim(file.path(system.file("testData", package = "FusionPathway"), "ArrayData_GSE48030_FUS-CHOP_Rodriguez.2013.txt"), sep = "\t", header = TRUE, as.is = TRUE)
ID2<-as.character(DataM1$ID)
loc<-match(ID2,ID)
pos<-which(loc>0)
loc<-loc[pos]
Gene.symbol<-rep("",length(ID2))
Gene.symbol[pos]<-GENE_SYMBOL[loc]
DataM1<-cbind(DataM1,Gene.symbol)
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
filein="FUS.DDIT3_RankedGenes.rnk"
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
outputfile1="ArrayData_GSE48030_SelectedGenes.xls"
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
fileoutput1="ArrayData_GSE48030_MappedDrugs_DGIdb.xls"
write.table(DataM, fileoutput1, sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)




 
 

 



