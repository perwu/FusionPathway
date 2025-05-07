#####
###################################################################
### FusionPathway: Run the domain-based network approach
###################################################################
#
FusionPathway<- function(GeneData,DomainData,Keep_AllInteractions="N"){
	##  Determine Interactions of the Fusion
	Interaction_Data<-Map_Interactions(GeneData,DomainData,Keep_AllInteractions=Keep_AllInteractions)

	## Rank all the genes based on association with the fusion
	Ranked_Result<-Association_Ranking_RWR(GeneData,Interaction_Data)
	
	## GSEA by fgsea
	Result_GSEA<-GSEA_fgsea(GeneData,Ranked_Result)
	
	## Mapping drugs
	Map_Drugs_DGIdb(Ranked_Result,GeneData) 
	
	## Result
	Result_List=list(Interaction_Data=Interaction_Data,Ranked_Result=Ranked_Result
		,Result_GSEA=Result_GSEA)
}



##########################################################################
## Mapping Molecular Interactions
##########################################################################
#
####### Protein and DNA interactions of The Fusion
Map_Interactions<-function (GeneData,DomainData,Keep_AllInteractions="N") 
{
	## Map Protein-Protein Interactions (PPI)
	PPI_Data<-Map_PPI(GeneData,DomainData)
	## Map Protein-DNA Interactions (PDI)
	PDI_Data<-Map_PDI(GeneData,DomainData)
	
	## kept and lost genes  (considering both PPI and PDI)
	if (Keep_AllInteractions=="N") {
		GeneID1_kept<-c(PPI_Data$GeneID1_kept,PDI_Data$GeneID1_kept)
		GeneID2_kept<-c(PPI_Data$GeneID2_kept,PDI_Data$GeneID2_kept)
		GeneID1_lost<-c(PPI_Data$GeneID1_lost,PDI_Data$GeneID1_lost)
		GeneID2_lost<-c(PPI_Data$GeneID2_lost,PDI_Data$GeneID2_lost)
		GeneSymbol1_kept<-c(PPI_Data$GeneSymbol1_kept,PDI_Data$GeneSymbol1_kept)
		GeneSymbol2_kept<-c(PPI_Data$GeneSymbol2_kept,PDI_Data$GeneSymbol2_kept)
		GeneSymbol1_lost<-c(PPI_Data$GeneSymbol1_lost,PDI_Data$GeneSymbol1_lost)
		GeneSymbol2_lost<-c(PPI_Data$GeneSymbol2_lost,PDI_Data$GeneSymbol2_lost)
	} else {
		GeneID1_kept<-c(PPI_Data$GeneID1_kept,PPI_Data$GeneID1_lost
			,PDI_Data$GeneID1_kept,PDI_Data$GeneID1_lost)
		GeneID2_kept<-c(PPI_Data$GeneID2_kept,PPI_Data$GeneID2_lost
			,PDI_Data$GeneID2_kept,PDI_Data$GeneID2_lost)
		GeneID1_lost<-NULL
		GeneID2_lost<-NULL
		GeneSymbol1_kept<-c(PPI_Data$GeneSymbol1_kept,PPI_Data$GeneSymbol1_lost
			,PDI_Data$GeneSymbol1_kept,PDI_Data$GeneSymbol1_lost)
		GeneSymbol2_kept<-c(PPI_Data$GeneSymbol2_kept,PPI_Data$GeneSymbol2_lost
			,PDI_Data$GeneSymbol2_kept,PDI_Data$GeneSymbol2_lost)
		GeneSymbol1_lost<-NULL
		GeneSymbol2_lost<-NULL
	}
	
	## Unique interacted genes
	GeneID1_kept_uniq<-unique(GeneID1_kept)
	GeneID2_kept_uniq<-unique(GeneID2_kept)
	GeneID1_lost_uniq<-unique(GeneID1_lost)
	GeneID2_lost_uniq<-unique(GeneID2_lost)
	#
	## Map Gene Symbol
	GeneSymbol1_kept_uniq<-GeneSymbol1_kept[match(GeneID1_kept_uniq,GeneID1_kept)]
	GeneSymbol2_kept_uniq<-GeneSymbol2_kept[match(GeneID2_kept_uniq,GeneID2_kept)]
	GeneSymbol1_lost_uniq<-GeneSymbol1_lost[match(GeneID1_lost_uniq,GeneID1_lost)]
	GeneSymbol2_lost_uniq<-GeneSymbol2_lost[match(GeneID2_lost_uniq,GeneID2_lost)]
	
	##
	Interaction_Data<-list(GeneID1_kept=GeneID1_kept_uniq,GeneID2_kept=GeneID2_kept_uniq,
		GeneID1_lost=GeneID1_lost_uniq,GeneID2_lost=GeneID2_lost_uniq,
		GeneSymbol1_kept=GeneSymbol1_kept_uniq,GeneSymbol2_kept=GeneSymbol2_kept_uniq,
		GeneSymbol1_lost=GeneSymbol1_lost_uniq,GeneSymbol2_lost=GeneSymbol2_lost_uniq)
	# Write the kept neighbor genes
	fileoutput1=paste(paste(GeneData$Gene1,GeneData$Gene2,sep="-")
		,"_KeptNeighbors.xls",sep="")
	write.table(union(Interaction_Data$GeneSymbol1_kept,Interaction_Data$GeneSymbol2_kept),
		fileoutput1, sep = "\t",append=FALSE,row.names = FALSE,col.names = FALSE,quote=FALSE)
	# Write the lost neighbor genes
	fileoutput2=paste(paste(GeneData$Gene1,GeneData$Gene2,sep="-")
		,"_LostNeighbors.xls",sep="")
	write.table(union(Interaction_Data$GeneSymbol1_lost,Interaction_Data$GeneSymbol2_lost),
		fileoutput2, sep = "\t",append=FALSE,row.names = FALSE,col.names = FALSE,quote=FALSE)
	
	##
	return (Interaction_Data)
}



####### Mapping Protein-DNA interactions
Map_PDI<-function (GeneData,DomainData) 
{
	## Extract PFAM Domain and Gene Information
	pfam1_kept=DomainData$pfam1_kept
	pfam2_kept=DomainData$pfam2_kept
	Gene1=GeneData$Gene1
	Gene2=GeneData$Gene2
	GeneID1=as.numeric(GeneData[["GeneID1"]])
	GeneID2=as.numeric(GeneData[["GeneID2"]])
	
	## Load DNA Binding Domain (PFAM database)
	load(file.path(system.file("data", package = "FusionPathway"), "DNAbindingDomain_pfam.RData"))
	PFAM_DNAbinding<-DNAbindingDomain_pfam[["PFAM_ID"]]
	
	## Check if DNA-bind Domains are kept]
	TFS1_kept="N"
	if(sum(pfam1_kept %in% PFAM_DNAbinding)>0) TFS1_kept="Y"
	TFS2_kept="N"
	if(sum(pfam2_kept %in% PFAM_DNAbinding)>0) TFS2_kept="Y"
	
	## transcription network
	load(file.path(system.file("data", package = "FusionPathway"), "Combined_TranscriptionNetwork.RData"))
	EntrezID_TFs<-as.numeric(Combined_TranscriptionNetwork$EntrezID_TFs)
	EntrezID_TG<-as.numeric(Combined_TranscriptionNetwork$EntrezID_TG)
	SYMBOL_TFs<-as.character(Combined_TranscriptionNetwork$SYMBOL_TFs)
	SYMBOL_TG<-as.character(Combined_TranscriptionNetwork$SYMBOL_TG)
	
	## Interaction kept and lost for first gene
	TFS1_pos<-which(EntrezID_TFs==GeneID1)
	if (length(TFS1_pos)>0 & TFS1_kept=="Y") {
		GeneID1_kept<-EntrezID_TG[TFS1_pos]
		GeneSymbol1_kept<-SYMBOL_TG[TFS1_pos]
		GeneID1_lost<-NULL
		GeneSymbol1_lost<-NULL
	} else if (length(TFS1_pos)>0 & TFS1_kept=="N") {
		GeneID1_kept<-NULL
		GeneSymbol1_kept<-NULL
		GeneID1_lost<-EntrezID_TG[TFS1_pos]
		GeneSymbol1_lost<-SYMBOL_TG[TFS1_pos]
	} else if (length(TFS1_pos)==0) {
		GeneID1_kept<-NULL
		GeneSymbol1_kept<-NULL
		GeneID1_lost<-NULL
		GeneSymbol1_lost<-NULL
	}
	## Interaction kept and lost for second gene
	TFS2_pos<-which(EntrezID_TFs==GeneID2)
	if (length(TFS2_pos)>0 & TFS2_kept=="Y") {
		GeneID2_kept<-EntrezID_TG[TFS2_pos]
		GeneSymbol2_kept<-SYMBOL_TG[TFS2_pos]
		GeneID2_lost<-NULL
		GeneSymbol2_lost<-NULL
	} else if (length(TFS2_pos)>0 & TFS2_kept=="N") {
		GeneID2_kept<-NULL
		GeneSymbol2_kept<-NULL
		GeneID2_lost<-EntrezID_TG[TFS2_pos]
		GeneSymbol2_lost<-SYMBOL_TG[TFS2_pos]
	} else if (length(TFS2_pos)==0) {
		GeneID2_kept<-NULL
		GeneSymbol2_kept<-NULL
		GeneID2_lost<-NULL
		GeneSymbol2_lost<-NULL
	}

	##
	Interaction_Data<-list(GeneID1_kept=GeneID1_kept,GeneID2_kept=GeneID2_kept,
		GeneID1_lost=GeneID1_lost,GeneID2_lost=GeneID2_lost,
		GeneSymbol1_kept=GeneSymbol1_kept,GeneSymbol2_kept=GeneSymbol2_kept,
		GeneSymbol1_lost=GeneSymbol1_lost,GeneSymbol2_lost=GeneSymbol2_lost)
	return (Interaction_Data)	
}


####### Mapping Protein-Protein interactions
Map_PPI<-function (GeneData,DomainData) 
{
	### Extract PFAM Domain and Gene Information
	pfam1_kept=DomainData$pfam1_kept
	pfam1_lost=DomainData$pfam1_lost
	pfam2_kept=DomainData$pfam2_kept
	pfam2_lost=DomainData$pfam2_lost	
	Gene1=GeneData$Gene1
	Gene2=GeneData$Gene2
	GeneID1=as.numeric(GeneData$GeneID1)
	GeneID2=as.numeric(GeneData$GeneID2)

	## Load Domain-Domain Interactions
	load(file.path(system.file("data", package = "FusionPathway"), "Combined_DomainInteraction.RData"))
	Inter_Dom1<-as.character(Combined_DomainInteraction$Inter_Dom1)
	Inter_Dom2<-as.character(Combined_DomainInteraction$Inter_Dom2)
	#
	## Domain interaction for the fusion partner gene 1
	DomB1_kept<-c(Inter_Dom2[which(Inter_Dom1 %in% pfam1_kept==TRUE)],
		Inter_Dom1[which(Inter_Dom2 %in% pfam1_kept==TRUE)])
	DomB1_lost<-c(Inter_Dom2[which(Inter_Dom1 %in% pfam1_lost==TRUE)],
		Inter_Dom1[which(Inter_Dom2 %in% pfam1_lost==TRUE)])
	## Domain interaction for the fusion partner gene 2
	DomB2_kept<-c(Inter_Dom2[which(Inter_Dom1 %in% pfam2_kept==TRUE)],
		Inter_Dom1[which(Inter_Dom2 %in% pfam2_kept==TRUE)])
	DomB2_lost<-c(Inter_Dom2[which(Inter_Dom1 %in% pfam2_lost==TRUE)],
		Inter_Dom1[which(Inter_Dom2 %in% pfam2_lost==TRUE)])

	## Load Protein Network
	load(file.path(system.file("data", package = "FusionPathway"), "Combined_PPI_Network.RData"))
	Network_EntrezID1=as.numeric(as.character(Combined_PPI_Network$Network_EntrezID1))
	Network_EntrezID2=as.numeric(as.character(Combined_PPI_Network$Network_EntrezID2))
	Network_GeneSymbol1=as.character(Combined_PPI_Network$Network_Symbol1)
	Network_GeneSymbol2=as.character(Combined_PPI_Network$Network_Symbol2)
	#
	# network partners of fusion partner gene 1
	GeneIDB_1<-c(Network_EntrezID2[which(Network_EntrezID1==GeneID1)],
		Network_EntrezID1[which(Network_EntrezID2==GeneID1)])
	GeneSymbolB_1<-c(Network_GeneSymbol2[which(Network_EntrezID1==GeneID1)],
		Network_GeneSymbol1[which(Network_EntrezID2==GeneID1)])
	# network partners of fusion partner gene 2
	GeneIDB_2<-c(Network_EntrezID2[which(Network_EntrezID1==GeneID2)],
		Network_EntrezID1[which(Network_EntrezID2==GeneID2)])
	GeneSymbolB_2<-c(Network_GeneSymbol2[which(Network_EntrezID1==GeneID2)],
		Network_GeneSymbol1[which(Network_EntrezID2==GeneID2)])
	#
	## Delete thos gene ID=NA
	GeneIDB_1<-GeneIDB_1[!is.na(GeneIDB_1)]
	GeneIDB_2<-GeneIDB_2[!is.na(GeneIDB_2)]
	GeneSymbolB_1<-GeneSymbolB_1[!is.na(GeneIDB_1)]
	GeneSymbolB_2<-GeneSymbolB_2[!is.na(GeneIDB_2)]
	# unique
	pos<-match(unique(GeneIDB_1),GeneIDB_1)
	GeneIDB_1<-GeneIDB_1[pos]
	GeneSymbolB_1<-GeneSymbolB_1[pos]
	pos<-match(unique(GeneIDB_2),GeneIDB_2)
	GeneIDB_2<-GeneIDB_2[pos]
	GeneSymbolB_2<-GeneSymbolB_2[pos]


	### Decide interaction partner of fusion partner genes are kept or lost in the new fusion based on domain-domain interaction
	data ("Biomart_Domain_AllGenes")
	PFAM<-Biomart_Genes_Domain$PFAM
	Homo_GeneID<-Biomart_Genes_Domain$Homo_GeneID
	## Partner of gene 1
	loc1<-match(Homo_GeneID,GeneIDB_1)
	pos1<-which(loc1>0)
	loc1<-loc1[pos1]
	pos2<-which(sapply(PFAM[pos1],function(x) sum(x %in% DomB1_kept))>0)
	pos_kept1<-unique(loc1[pos2])
	pos2<-which(sapply(PFAM[pos1],function(x) sum(x %in% DomB1_lost))>0)
	pos_lost1<-unique(loc1[pos2])
	#pos_part1<-setdiff(pos_kept1,pos_lost1)
	pos_part1<-pos_kept1
	## Partner of gene 2
	loc1<-match(Homo_GeneID,GeneIDB_2)
	pos1<-which(loc1>0)
	loc1<-loc1[pos1]
	pos2<-which(sapply(PFAM[pos1],function(x) sum(x %in% DomB2_kept))>0)
	pos_kept2<-unique(loc1[pos2])
	pos2<-which(sapply(PFAM[pos1],function(x) sum(x %in% DomB2_lost))>0)
	pos_lost2<-unique(loc1[pos2])
	#pos_part2<-setdiff(pos_kept2,pos_lost2)
	pos_part2<-pos_kept2

	### All kept interaction partners
	GeneID1_kept<-GeneIDB_1[pos_part1]	# interactions partners of gene 1
	GeneSymbol1_kept<-GeneSymbolB_1[pos_part1]
	GeneID2_kept<-GeneIDB_2[pos_part2]	# interactions partners of gene 2
	GeneSymbol2_kept<-GeneSymbolB_2[pos_part2]
	
	### All lost interaction partners 
	GeneID1_lost<-setdiff(GeneIDB_1,union(GeneID1_kept,GeneID2_kept))	# interactions partners of gene 1
	GeneSymbol1_lost<-GeneSymbolB_1[match(GeneID1_lost,GeneIDB_1)]
	GeneID2_lost<-setdiff(GeneIDB_2,union(GeneID1_kept,GeneID2_kept))	# interactions partners of gene 2
	GeneSymbol2_lost<-GeneSymbolB_2[match(GeneID2_lost,GeneIDB_2)]
	
	##
	Interaction_Data<-list(GeneID1_kept=GeneID1_kept,GeneID2_kept=GeneID2_kept,
		GeneID1_lost=GeneID1_lost,GeneID2_lost=GeneID2_lost,
		GeneSymbol1_kept=GeneSymbol1_kept,GeneSymbol2_kept=GeneSymbol2_kept,
		GeneSymbol1_lost=GeneSymbol1_lost,GeneSymbol2_lost=GeneSymbol2_lost)
	return (Interaction_Data)	
}




##########################################################################
### Rank Genes by network guilt-by-association 
##########################################################################
#
# random walk with restart (new version)
RandomWalk_Restart<- function(TranM, sourceGenes) {
    gamma <- 0.8
	N_Genes <- dim(TranM)[1]
	#
	if(sum(!sourceGenes %in% row.names(TranM))>0) {
		stop("sourceGenes contains genes not found in the network")
    }
    
    P_t0 <- rep(0,N_Genes)
    names(P_t0) <- row.names(TranM)
	#
	loc<-match(row.names(TranM),sourceGenes)
	selected_pos<-which(loc>0)

	# initial probability vector
    P_t0[selected_pos] <- 1
    P_t0 <- P_t0/sum(P_t0)
    
	## Random walk with restart (rwr)
	P_t0<-t(t(P_t0))
	P_ts1 <- P_t0
    delta <- 1
	delta2<-2
	Diff<-1
    while  ((delta > 1e-25) & (Diff> 1e-60)) {
	  #print(Diff)
      P_ts2 <- (((1-gamma)*t(TranM)) %*% P_ts1)+(gamma*P_t0)
      delta <- sum(abs(P_ts2 - P_ts1))
	  Diff<-delta2-delta
      P_ts1 <- P_ts2
	  delta2<-delta
    }
	#
	return(drop(P_ts1))
}


####### rank using random walk with restart
Association_Ranking_RWR<-function (GeneData,Interaction_Data) 
{
	## All Kept interactions
	GeneID_kept<-c(Interaction_Data$GeneID1_kept,Interaction_Data$GeneID2_kept)
	GeneSymbol_kept<-c(Interaction_Data$GeneSymbol1_kept,Interaction_Data$GeneSymbol2_kept)

	## Load Network Data (Transition Matrix)
	load(file.path(system.file("data", package = "FusionPathway"), "TranMatrix_undirected_sparse.RData"))
	row.names(TranMatrix)<-Network_GeneID

	## Rank gene based on random walk with restart
	if(sum(GeneID_kept %in% Network_GeneID)>0) {
		#random walk with restart
		Prob_V<-RandomWalk_Restart(TranMatrix,GeneID_kept)
		
		# Sort
		DD<-sort(Prob_V,decreasing=TRUE,index.return=TRUE)
		index<-as.numeric(DD$ix)
		GeneID_ranked<-Network_GeneID[index]
		Score<-Prob_V[index]
		GeneSymbol_ranked<-Network_GeneSymbol[index]
	} else {	# no neighbor genes in the network		
		Score<-seq(length(Network_GeneID),1,by=-1)
		GeneID_ranked<-Network_GeneID
		GeneSymbol_ranked=Network_GeneSymbol
	}
	
	##  rnk file for GSEA
	Ranked_Result<-list(GeneID_ranked=GeneID_ranked,GeneSymbol_ranked=GeneSymbol_ranked,Score=Score)
	fileoutput1=paste(paste(GeneData$Gene1,GeneData$Gene2,sep=".")
		,"_RankedGenes.rnk",sep="")
	write.table(cbind(Ranked_Result$GeneSymbol_ranked,seq(length(Ranked_Result$Score),1,by=-1))
		, fileoutput1, sep = "\t",row.names = FALSE,col.names = FALSE,quote=FALSE)
	
	##	
	return (Ranked_Result)
}



##########################################################################
### Pathway Analysis Based on Ranked Genes
##########################################################################
#
####### GSEA using "fgsea"
GSEA_fgsea <- function(GeneData,Ranked_Result,FDR_threshold=1){
	library("fgsea")

	## 
	load(file.path(system.file("data", package = "FusionPathway"), "GSEAPathway_Data.RData"))
	Pathway_Names<-names(Pathway_M)

	## association socre ranking
	Genes_ranked<-Ranked_Result$GeneSymbol_ranked
	Score_ranked<-as.numeric(Ranked_Result$Score)
	# uniuqe
	Genes_ranked_uniq<-unique(Genes_ranked)
	loc<-match(Genes_ranked_uniq,Genes_ranked)
	Genes_ranked<-Genes_ranked[loc]
	Score_ranked<-Score_ranked[loc]
	#
	Score_ranked[which(Score_ranked==Inf)]<-1000000
	names(Score_ranked)<-Genes_ranked
	# rank order
	Score_ranked2<-seq(length(Genes_ranked),1,by=-1)
	names(Score_ranked2)<-Genes_ranked
	
	## run fGSEA
    	fgseaRes <- fgsea(pathways = Pathway_M, stats = Score_ranked2,
                  minSize=15, maxSize=500, nperm=10000)
    
    ## Select top pathways and write the results
    #Result_GSEA<-fgseaRes[fgseaRes[, padj < FDR_threshold],]
    Result_GSEA<-fgseaRes
    Result_GSEA = format(Result_GSEA)
	fileoutput1=paste(paste(GeneData$Gene1,GeneData$Gene2,sep="-")
		,"_GSEA_Result",".xls",sep="")
	write.table(Result_GSEA,fileoutput1, sep = "\t",append=FALSE,row.names = FALSE,col.names = TRUE,quote=FALSE)
    return(Result_GSEA)
    
    ## list the top pathway
    #topPathways <- fgseaRes[ES > 0][head(order(pval), n=20), pathway]
	#plotGseaTable(Pathway_M[topPathways], Score_ranked2, fgseaRes, gseaParam = 0.5)
    
    # plot a pathway
    #library("ggplot2")
    #plotEnrichment(Pathway_M[["KEGG_JAK_STAT_SIGNALING_PATHWAY"]],
    #           Score_ranked2) + labs(title="JAK STAT Signaling")
}    


##########################################################################
## Combine pathway p-values
##########################################################################
#
#### Truncated product method (TPM)
TPM<-function(w,tau,L) {
	tpm_p=0
	for (k in 1:L) 
	{
		#if (k==1) {
		#	r1=2			
		#} else if (k==2) {
		#	r1=1
		#}
		DD=dim(combn(L, k))
		r1=DD[2]
		r2=(1-tau)^(L-k)
		#
		if (w<=(tau^k)) {
			innersum=0
			for (s in 1:k)
			{
				A=k*log(tau)-log(w)
				innersum=innersum+((A^(s-1))/gamma(s))
				#print(gamma(s))
			}
			innersum=w*innersum
		} else {
			innersum = tau^k
		}
		tpm_p=tpm_p+(r1*r2*innersum)
	}
	return (tpm_p)
}

###
Combine_Pathway_Pvalues<-function(PathwayM1,PathwayM2,threshold_ind_pvalue=0.05,outfile) {
	Pathways1<-PathwayM1[[1]]
	pvalue_1<-as.numeric(PathwayM1[[2]])
	Pathways2<-PathwayM2[[1]]
	pvalue_2<-as.numeric(PathwayM2[[2]])
	
	##
	Pathways<-intersect(Pathways1,Pathways2)
	Pathways_diff1<-setdiff(Pathways1,Pathways2)
	Pathways_diff2<-setdiff(Pathways2,Pathways1)
	#
	loc1<-match(Pathways,Pathways1)
	loc2<-match(Pathways,Pathways2)
	Pvalue_1<-pvalue_1[loc1]
	Pvalue_2<-pvalue_2[loc2]
	#
	loc3<-match(Pathways_diff1,Pathways1)
	p1a<-pvalue_1[loc3]
	p2a<-rep(1,length(p1a))
	#
	loc4<-match(Pathways_diff2,Pathways2)
	p1b<-pvalue_2[loc4]
	p2b<-rep(1,length(p1b))
	
	##
	Pathways<-c(Pathways,Pathways_diff1,Pathways_diff2)
	Pvalue_1<-c(Pvalue_1,p1a,p1b)
	Pvalue_2<-c(Pvalue_2,p2a,p2b)
	
	## combine two pavlues
	Combined_Pvalues<-rep(1,length(Pathways))
	for (i in 1:length(Pathways)) {
		p1=(Pvalue_1[i])^(ifelse(Pvalue_1[i]<=threshold_ind_pvalue,1,0)) 	
		p2=(Pvalue_2[i])^(ifelse(Pvalue_2[i]<=threshold_ind_pvalue,1,0)) 
		if (p1==0) {p1=0.000001}
		if (p2==0) {p2=0.000001}
		
		#
		if (p1<1 & p2<1 ) {				
			W=p1*p2
		} else{
			W=1
		}					
		if (W<1) {
			Combined_Pvalues[i]=TPM(W,threshold_ind_pvalue,2)
		} else {Combined_Pvalues[i]=1}
	}
	#
	Combined_Pvalues_adj<-p.adjust(Combined_Pvalues, method = "fdr", n = length(Combined_Pvalues))

	## 
	ResultM<-data.frame(Pathways=Pathways,pval1=Pvalue_1,pval2=Pvalue_2
		,combined_pval=Combined_Pvalues,combined_adj=Combined_Pvalues_adj)
	# sort
	DD<-sort(Combined_Pvalues_adj,decreasing=F,index.return=TRUE)
	index<-as.numeric(DD$ix)
	ResultM<-ResultM[index,]
	#
	write.table(ResultM,outfile, sep = "\t",append=FALSE,row.names = FALSE,col.names = TRUE,quote=FALSE)
}



##########################################################################
## Map Drugs 
##########################################################################
#
####### Map DGIdb
Map_Drugs_DGIdb<-function (Ranked_Result,GeneData) {	
	##
	load(file.path(system.file("data", package = "FusionPathway"), "DrugDatabase_DGIdb.RData"))
	DataM<-DrugDatabase_DGIdb
	Drug_GeneNames<-as.character(DataM$entrez_gene_symbol)
	
	## Ranked genes
	GeneNames_ranked<-Ranked_Result$GeneSymbol_ranked
	Score<-Ranked_Result$Score
	Ranks<-(1:length(GeneNames_ranked))
	Rank_Percentage<-Ranks*100/length(GeneNames_ranked)

	## Mapping to drug-target information
	loc<-match(Drug_GeneNames,GeneNames_ranked)
	pos<-which(loc>0)
	loc<-loc[pos]
	DataM<-cbind(DataM[pos,],Ranks=Ranks[loc],Rank_Percentage=Rank_Percentage[loc])

	
	### select the top 5% assocaited genes
	pos2<-which(DataM$Rank_Percentage<=100)
	DataM<-DataM[pos2,]
	# sort based on ranks
	DD<-sort(DataM$Rank_Percentage,index.return=TRUE)
	indx<-as.numeric(DD$ix)
	DataM<-DataM[indx,]
	
	## Write
	# Fusion partner genes
	fileoutput1=paste(paste(GeneData$Gene1,GeneData$Gene2,sep="-")
		,"_MappedDrugs_DGIdb.xls",sep="")
	write.table(DataM, fileoutput1, sep = "\t",row.names = FALSE,col.names = TRUE,quote=FALSE)
}

####### Calculate AUC
ROCF<-function(prediction,score)
{
	N<-length(score)
	N_pos<-sum(prediction)
	N_neg<-N-N_pos
	df<-data.frame(pred=prediction,score=score)
	df<-df[order(-df$score),]
	df$above=(1:N)-cumsum(df$pred)
	AUC<-1-sum(df$above*df$pred/(N_pos*(N-N_pos)))
	ranking<-seq(N,1,by=-1)
	TP<-cumsum(df$pred)
	FP =df$above
	threshold_indx = which(diff(ranking)!=0); 
	TPR = TP[threshold_indx]/N_pos;
	FPR = FP[threshold_indx]/N_neg;
	ROC=list(TPR=TPR,FPR=FPR,AUC=AUC)
	return(ROC)
}


