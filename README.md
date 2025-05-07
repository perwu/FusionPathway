# FusionPathway

Understanding the molecular effects of a fusion is important in determining its biological function in oncogenesis and in developing therapeutic strategies for cancer patients with the fusion, especially if the fusion is undruggable. However, the molecular mechanisms have been elucidated for only a few fusions, in part because of the labor-intensive nature of the required functional assays. Therefore, we developed a domain-based network approach to infer the pathways that are functionally associated with a fusion and identify potential therapeutic targets in these pathways. The FusionPathway package provides the functionalities for implementing the developed domain-based network approach. 

Our proposed domain-based network approach is illustrated in the figure below. First, a fusion protein consists of the domains of its parental proteins, therefore we hypothesized that the fusion protein retains some of the molecular interactions of its parental proteins. Given the fusion’s domain composition, we can predict the protein-protein or protein-DNA interactions of the fusion (panel a in the figure). Second, we reasoned that genes that are functionally associated with these predicted molecular interactions are also functionally associated with the fusion (panel b in the figure). Therefore, we used the network “guilt-by-association” method, which has been used to prioritize novel disease genes and annotate gene function, to calculate the functional association of each gene with the fusion. Third, we considered pathways that enrich genes with high fusion-association scores are also functionally associated with the fusion (panel c in the figure). We ranked all the genes on the basis of their fusion-association scores. Gene set enrichment analysis (GSEA), which can evaluate the genes of a pathway for their distribution in the ordered gene list, was then used to identify pathways functionally associated with a fusion (termed GSEA association analysis). Forth, by leveraging gene expression data from samples that have or lack a fusion, we can identify pathways that are deregulated by the fusion (panel d in the figure). The deregulation level of an associated pathway can be also determined using GSEA (termed GSEA deregulation analysis). The truncated product method (Zaykin et al., 2002) was used to combine the p values of each pathway generated from the GSEA association and GSEA deregulation analyses to identify pathways that are both highly associated with the fusion and significantly deregulated by it in a specific cellular condition.

![figure1](https://user-images.githubusercontent.com/14062661/30707054-ded42b02-9ebf-11e7-872e-6d50b13302c7.jpg)



# Installation

```R
install.packages("devtools") 
library("devtools") 
devtools::install_github("perwu/FusionPathway") 
```


# Scripts for the three examples of fusions
This packages also contain scripts for the three examples mentioned in our paper, BCR-ABL1, EWSR1-FLI1, and FUS-DDIT3. The scripts can be accessed in the subfolder "scripts" under the package folder. 

```R
file.edit(file.path(system.file("scripts", package = "FusionPathway"), "Example_Fusion_BCR-ABL1.r")) 
file.edit(file.path(system.file("scripts", package = "FusionPathway"), "Example_Fusion_EWS_FLI1.r")) 
file.edit(file.path(system.file("scripts", package = "FusionPathway"), "Example_Fusion_FUS-DDIT3.r")) 
```


# Example: BCR-ABL1
We here applied our approach to the p210 BCR-ABL1. The ABL1 portion of p210 BCR-ABL1 contains tandem SRC homology (PF00017), the tyrosine kinase domains (PF07714), SH3 binding sites (PF00018), a DNA-binding domain, and an actin-binding domain (PF08919). The BCR portions of p210 BCR-ABL1 contains a coiled-coil oligomerisation domain (PF09036), Dbl homology domain (PF00621), and a pleckstrin homology domain (PF00169). 

First, provide gene IDs, gene symbols, and protein domains (Pfam IDs) of the two parental genes, BCR and ABL1.  

```R
library("FusionPathway")  
Gene1="BCR" 
GeneID1<-613 
PFAM1=c("PF09036","PF00621","PF00169","PF00168","PF00620") 
Gene2="ABL1" 
GeneID2<-25 
PFAM2=c("PF00018","PF00017","PF07714","PF08919") 
```

### Determine which proteins domains are retained in p210 BCR-ABL1: 

```R
 PFAM1_lost<-c("PF00620","PF00168") 
 PFAM1_kept<-setdiff(PFAM1,PFAM1_lost)
 PFAM2_lost<-"" 
 PFAM2_kept<-setdiff(PFAM2,PFAM2_lost) 
```

### Run FusionPathway with the data of parental genes and their domains: 

```R
DomainData<-list(pfam1_kept=PFAM1_kept,pfam1_lost=PFAM1_lost,pfam2_kept=PFAM2_kept,pfam2_lost=PFAM2_lost) 
GeneData<-data.frame(Gene1=Gene1,Gene2=Gene2,GeneID1=GeneID1,GeneID2=GeneID2) 
Result_List<-FusionPathway(GeneData,DomainData) 
save(file="Ranked_Result.RData",Result_List) 
```

The “Result_List” contains predicted the protein-protein interaction partners of p210 BCR-ABL1, an ordered gene list ranked by the fusion-association prediction, results of GSEA pathway association analysis. 

Prediction evaluation using 5 literature-based benchmark gene sets: 26 Wnt/Ca+/NFAT pathway genes 
(denoted WNT_CA+_NFAT genes), 1150 genes cited with BCR-ABL1 in the literature (denoted BCR-ABL genes), 
1240 genes cited with CML in the literature (denoted CML genes), 328 genes categorized in the 
KEGG cancer pathways (denoted cancer pathway genes), and 68 target genes of compounds that have been 
tested in clinical trials or used for the treatment of imatinib-resistant CML (denoted drug targets). 
Data file “BenchmarkSets_CML.RData” contains these benchmark gene sets.  

```R
Score<-as.numeric(Result_List$Ranked_Result$Score) 
GeneID_ranked<-as.numeric(Result_List$Ranked_Result$GeneID_ranked) 
BenchmarkNames<-c("CML_Genes","BCR_ABL1_Genes","WntCaNfat_Genes","CancerPathway_Genes","Drug_Targets") 
load(file.path(system.file("data", package = "FusionPathway"), "BenchmarkSets_CML.RData")) 
ROC_M<-NULL 
for (i in 1:5) { 
	if (i==1 | i==2) { 
		PubmedN<-BenchmarkSets_CML[[BenchmarkNames[i]]]$PubmedN 
		GeneIDs<-BenchmarkSets_CML[[BenchmarkNames[i]]]$GeneIDs 
		Ncitation_threshold=2 
		BenchID_pos<-GeneIDs[which(PubmedN>=Ncitation_threshold)] 
	} else { 
		BenchID_pos<-BenchmarkSets_CML[[BenchmarkNames[i]]]$GeneIDs	 
	}		 
	loc<-match(GeneID_ranked,BenchID_pos) 
	outcome=rep(0,length(GeneID_ranked)) 
	outcome[which(loc>0)]=1 
	ROC<-ROCF(outcome,Score) 
	ROC_M[[BenchmarkNames[i]]]$TPR<-ROC$TPR 
	ROC_M[[BenchmarkNames[i]]]$FPR<-ROC$FPR 
	ROC_M[[BenchmarkNames[i]]]$AUC<-round(ROC$AUC,3) 
	print(paste("AUC using ",BenchmarkNames[i]," genes: ",round(ROC$AUC,3),sep="")) 
} 
save(file="ROC.RData",ROC_M) 
```

### Plot ROC curves: 

```R
outputfile="AUC_Different_Benchmarks.jpeg" 
jpeg(outputfile,width = 510, height = 510, units = "px") 
par(mar=c(4,4.2,0.5,0.5)) 
pp<-plot(ROC_M$CML_Genes$FPR,ROC_M$CML_Genes$TPR,xlab="FPR",ylab="TPR", 
>	col = "red",lwd=3,"l",cex=1.45,cex.axis=1.45,cex.lab=1.45,font.lab=2.3) 
points(ROC_M$BCR_ABL1_Genes$FPR,ROC_M$BCR_ABL1_Genes$TPR,type="l",col="black",lwd=3) 
points(ROC_M$CancerPathway_Genes$FPR,ROC_M$CancerPathway_Genes$TPR,type="l",col="blue",lwd=3) 
points(ROC_M[["WntCaNfat_Genes"]]$FPR,ROC_M[["WntCaNfat_Genes"]]$TPR,type="l",col="purple",lwd=3) 
points(ROC_M$Drug_Targets$FPR,ROC_M$Drug_Targets$TPR,type="l",col="green",lwd=3) 
LabelS<-c(paste("CML Genes(AUC=",ROC_M$CML_Genes$AUC,")",sep=""), 
	paste("BCR-ABL Genes(AUC=",ROC_M$BCR_ABL1_Genes$AUC,")",sep=""), 
	paste("Cancer Pathway Genes(AUC=",ROC_M$CancerPathway_Genes$AUC,")",sep=""), 
	paste("WNT_CA_NFAT Genes(AUC=",ROC_M[["WntCaNfat_Genes"]]$AUC,")",sep=""), 
	paste("Drug Targets(AUC=",ROC_M$Drug_Targets$AUC,")",sep="")) 
legend(0.145,0.27,col=c("red","black","blue","purple","green"),  
>	legend=LabelS,lty=rep(1,5),lwd=3,cex=1.4,text.font=2) 
dev.off() 
```

![auc_different_benchmarks](https://user-images.githubusercontent.com/14062661/30716320-38c6b2aa-9edf-11e7-813e-e68d95b2f61e.jpeg)





