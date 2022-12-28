####################################################################################
# Analysis of epigenetic factors in TNBC Patient transcriptomic data to prioritize  
# TNBC Specific Epigenetic factor
# Patient data used : Metabric dataset ( downloaded from cbioportal )
####################################################################################


library(reshape)
library(pheatmap)
library(gplots)


################################################################################################
# Identified list of epigenetic factors from EpiFactors database
###################################################################################################

Epifactors_Genelist = read.csv('EpiFactors_Analysis/EpiGenes_main.csv')
Epifactors_Genelist = unique(Epifactors_Genelist$HGNC_symbol)


Epifactors_Metabric_genes = TNBC_METABRIC_Expression[which(as.matrix(TNBC_METABRIC_Expression$Hugo_Symbol) %in% as.matrix(Epifactors_Genelist)),]
colnames(Epifactors_Metabric_genes) = gsub("\\.","-",colnames(Epifactors_Metabric_genes))

TNBC_Samplenames = TNBC_METABRIC_Clinical[which(TNBC_METABRIC_Clinical$SAMPLE_ID %in% 
                                                  Results_From_Shira.selectTNBC$SAMPLEID),]
  
ER_Pos_Samplenames = TNBC_METABRIC_Clinical[which(TNBC_METABRIC_Clinical$ER_STATUS == "+" & TNBC_METABRIC_Clinical$HER2_STATUS == "-" & TNBC_METABRIC_Clinical$PR_STATUS == "-"),]
ER_Pos_Samplenames = ER_Pos_Samplenames[which(!ER_Pos_Samplenames$SAMPLE_ID %in% TNBC_Samplenames$SAMPLE_ID),]

HER_Pos_Samplenames = TNBC_METABRIC_Clinical[which(TNBC_METABRIC_Clinical$ER_STATUS == "-" & TNBC_METABRIC_Clinical$HER2_STATUS == "+" & TNBC_METABRIC_Clinical$PR_STATUS == "-"),]
HER_Pos_Samplenames = HER_Pos_Samplenames[which(!HER_Pos_Samplenames$SAMPLE_ID %in% TNBC_Samplenames$SAMPLE_ID),]

## Combine references ##
HER_Pos_Samplenames$MolecularSubtype = "HER2_Positive"
ER_Pos_Samplenames$MolecularSubtype = "ER_Positive"
TNBC_Samplenames$MolecularSubtype = "TNBC"

Metabric_Molecular_subtpyed_annotation = rbind(TNBC_Samplenames,HER_Pos_Samplenames,ER_Pos_Samplenames)
Epifactors_Metabric_genes.fil = Epifactors_Metabric_genes[,c(1,2,which(colnames(Epifactors_Metabric_genes) %in% as.matrix(Metabric_Molecular_subtpyed_annotation$SAMPLE_ID)))]

#### Compare T statistics across TNBC vs. HER2 ###
sel_comp = Metabric_Molecular_subtpyed_annotation[which(Metabric_Molecular_subtpyed_annotation$MolecularSubtype %in% c("TNBC","HER2_Positive")),]
sel.expr_all = Epifactors_Metabric_genes.fil[,c(1,which(colnames(Epifactors_Metabric_genes.fil) %in% sel_comp$SAMPLE_ID))]
sel.expr_all = sel.expr_all[which(rowSums(sel.expr_all[,-1] > 5) >=10),]
sel.expr = melt(sel.expr_all)
sel.expr = merge(sel.expr,sel_comp[,c("SAMPLE_ID","MolecularSubtype")],by.x="variable",by.y="SAMPLE_ID",all.x=TRUE)


################################################################################################
# Statistical comparison of epigenetic markers expressed between TNBC and HER2 Positive Samples
###################################################################################################

TNBC_vs_HER2_WilcoxRes =   data.frame(Gene=character(0),Comparison=character(0),GreaterIn=character(0),Wilcox_Pval=numeric(0))
Comparison = "TNBC vs. HER2+"
for (geneU in unique(sel.expr$Hugo_Symbol)){
  #geneU ="PRDM4"
  candidate.comp = sel.expr[which(sel.expr$Hugo_Symbol == geneU),]
  result_wilcox=wilcox.test(candidate.comp$value~candidate.comp$MolecularSubtype)
  wilcox_pval = result_wilcox$p.value
  Computed_mean = aggregate(candidate.comp, list(candidate.comp$MolecularSubtype), mean)
  Group_max=Computed_mean$Group.1[which(Computed_mean$value == max(Computed_mean$value))]
  TNBC_vs_HER2_WilcoxRes<-rbind(TNBC_vs_HER2_WilcoxRes, data.frame(Gene=geneU,Comparison=Comparison,GreaterIn=Group_max,Wilcox_Pval=wilcox_pval))

  
  
}

TNBC_vs_HER2_WilcoxRes = TNBC_vs_HER2_WilcoxRes[order(TNBC_vs_HER2_WilcoxRes$Wilcox_Pval),]
TNBC_vs_HER2_WilcoxRes$FDR = p.adjust(TNBC_vs_HER2_WilcoxRes$Wilcox_Pval,method = "BH")
Significant_TNBC_vs_HER2 = (TNBC_vs_HER2_WilcoxRes[which(TNBC_vs_HER2_WilcoxRes$FDR < 0.01),])
Significant_TNBC_vs_HER2 = Significant_TNBC_vs_HER2[which(Significant_TNBC_vs_HER2$GreaterIn== "TNBC"),]

################################################################################################
# Statistical comparison of epigenetic markers expressed between TNBC and  ER+  Samples
###################################################################################################

sel_comp = Metabric_Molecular_subtpyed_annotation[
  which(Metabric_Molecular_subtpyed_annotation$MolecularSubtype %in% c("TNBC","ER_Positive")),]
sel.expr_tnbc_er = Epifactors_Metabric_genes.fil[,c(1,which(colnames(Epifactors_Metabric_genes.fil) %in% as.matrix(sel_comp$SAMPLE_ID)))]

sel.expr = melt(sel.expr_tnbc_er)
sel.expr = merge(sel.expr,sel_comp[,c("SAMPLE_ID","MolecularSubtype")],by.x="variable",by.y="SAMPLE_ID",all.x=TRUE)

TNBC_vs_ER_WilcoxRes =   data.frame(Gene=character(0),Comparison=character(0),GreaterIn=character(0),Wilcox_Pval=numeric(0))
Comparison = "TNBC vs. ER+"
for (geneU in unique(sel.expr$Hugo_Symbol)){
  #geneU ="PRDM4" 
  candidate.comp = sel.expr[which(sel.expr$Hugo_Symbol == geneU),]
  result_wilcox=wilcox.test(candidate.comp$value~candidate.comp$MolecularSubtype)
  wilcox_pval = result_wilcox$p.value
  Computed_mean = aggregate(candidate.comp, list(candidate.comp$MolecularSubtype), mean)
  Group_max=Computed_mean$Group.1[which(Computed_mean$value == max(Computed_mean$value))]
  
  TNBC_vs_ER_WilcoxRes<-rbind(TNBC_vs_ER_WilcoxRes, data.frame(Gene=geneU,Comparison=Comparison,GreaterIn=Group_max,Wilcox_Pval=wilcox_pval))
  
}


TNBC_vs_ER_WilcoxRes = TNBC_vs_ER_WilcoxRes[order(TNBC_vs_ER_WilcoxRes$Wilcox_Pval),]
TNBC_vs_ER_WilcoxRes$FDR = p.adjust(TNBC_vs_ER_WilcoxRes$Wilcox_Pval,method = "BH")
Significant_TNBC_vs_ER = (TNBC_vs_ER_WilcoxRes[which(TNBC_vs_ER_WilcoxRes$FDR < 0.01),])
Significant_TNBC_vs_ER = Significant_TNBC_vs_ER[which(Significant_TNBC_vs_ER$GreaterIn =="TNBC"),]

#### Heatmap ####
anno_col= Metabric_Molecular_subtpyed_annotation
rownames(anno_col) = anno_col$SAMPLE_ID
anno_col=anno_col[,-1]


hm_input = sel.expr_all[which(sel.expr_all$Hugo_Symbol %in% Significant_TNBC_vs_HER2$Gene),]
rownames(hm_input) = hm_input$Hugo_Symbol
hm_input = hm_input[,-1]
anno_col.input = anno_col[which(anno_col$PATIENT_ID %in% colnames(hm_input)),]
anno_col.input = anno_col.input[,c(10),drop=FALSE]
Colors= colorRampPalette(c("green", "black", "red"))(n = 1000)

hm.obj_TNBC_HER2 = pheatmap(as.matrix(log2(hm_input+1)),
                  clustering_distance_rows = "correlation",
                  clustering_distance_cols = "euclidean",
                  cluster_cols = T,
                  color = Colors,
                  cluster_rows = T,
                  annotation_col = anno_col.input,
                  fontsize = 8,
                  scale = "row",
                  fontsize_row=2, 
                  fontsize_col = 6,
                  main = "195 Genes differentially expressed between TNBC and HER2 Positive (FDR <0.01)"
)

# Using hierarchical clustering to assign clustered TNBC Specific epifactors highly expressed compared to HER2+
TNBC_Specific_Factors_vs_HER2 = as.matrix(cutree(hm.obj_TNBC_HER2$tree_row,2)) 

# cluster 1 : genes up in TNBC vs HER2
TNBC_Specific_Factors_vs_HER2 = rownames(TNBC_Specific_Factors_vs_HER2)[which(TNBC_Specific_Factors_vs_HER2 == 1)]

######### Visualization of the subset of highly expressed epifactors in TNBC vs HER2 ##########

anno_col= Metabric_Molecular_subtpyed_annotation
rownames(anno_col) = anno_col$SAMPLE_ID
anno_col=anno_col[,-1]

hm_input = sel.expr_tnbc_er[which(sel.expr_tnbc_er$Hugo_Symbol %in% Significant_TNBC_vs_ER$Gene),]
rownames(hm_input) = hm_input$Hugo_Symbol
hm_input = hm_input[,-1]
anno_col.input = anno_col[which(anno_col$PATIENT_ID %in% as.matrix(colnames(hm_input))),]
anno_col.input = anno_col.input[,c(10),drop=FALSE]
Colors= colorRampPalette(c("green", "black", "red"))(n = 1000)

hm.obj_TNBC_ER_Pos = pheatmap(as.matrix(log2(hm_input+1)),
                  clustering_distance_rows = "correlation",
                  clustering_distance_cols = "euclidean",
                  cluster_cols = T,
                  color = Colors,
                  cluster_rows = T,
                  annotation_col = anno_col.input,
                  fontsize = 8,
                  #cutree_cols = 2,
                  # annotation_colors =  anno_colors,
                  scale = "row",
                  #breaks = myBreaks,
                  fontsize_row=0.4, 
                  fontsize_col = 2,
                  main = "377 Genes differentially expressed between TNBC and ER Positive (FDR <0.01)"
                  # treeheight_col = 0,border_color = 'NA'
)


# Using hierarchical clustering to assign clustered TNBC Specific epifactors highly expressed compared to ER+
TNBC_Specific_Factors_vs_ER= as.matrix(cutree(hm.obj_TNBC_ER_Pos$tree_row,3)) 
# cluster 2 : genes up in TNBC vs. ER+
TNBC_Specific_Factors_vs_ER = rownames(TNBC_Specific_Factors_vs_ER)[which(TNBC_Specific_Factors_vs_ER == 2)]




######### Combined enriched TNBC compared to ER+ and HER2 + ############ 
TO_PLot_Genes = unique(intersect(as.matrix(Significant_TNBC_vs_ER$Gene),as.matrix(Significant_TNBC_vs_HER2$Gene)))
TNBC_Specific_Epifactors = TO_PLot_Genes

anno_col= Metabric_Molecular_subtpyed_annotation
rownames(anno_col) = anno_col$SAMPLE_ID
anno_col=anno_col[,-1]

hm_input = Epifactors_Metabric_genes.fil[which(Epifactors_Metabric_genes.fil$Hugo_Symbol %in% TNBC_Specific_Epifactors),]
rownames(hm_input) = hm_input$Hugo_Symbol
hm_input = hm_input[,-1]
hm_input = hm_input[, which( colnames(hm_input) %in% as.matrix(anno_col$PATIENT_ID))]
anno_col.input = anno_col[which(anno_col$PATIENT_ID %in% as.matrix(colnames(hm_input))),]
anno_col.input = anno_col.input[,c(10),drop=FALSE]
paletteLength <- 7
breaks=seq(-4, 4, by=0.4)
breaks=append(breaks, 5)
breaks=append(breaks, -5,0)
mycol <- colorpanel(n=length(breaks)-1,low="green",mid="black",high="red")

Colors= colorRampPalette(c("green", "black", "red"))(n = 1000)
myBreaks <- c(seq(min(hm_input), 0, length.out=ceiling(paletteLength) + 1), 
              seq(max(hm_input)/paletteLength, max(hm_input), length.out=floor(paletteLength)))

hm.obj_all = pheatmap(as.matrix(log2(hm_input+1)),
                  clustering_distance_rows = "euclidean",
                  clustering_distance_cols = "euclidean",
                  cluster_cols = T,
                  color = mycol,
                  cluster_rows = T,
                  annotation_col = anno_col.input,
                  fontsize = 8,border_color = 'NA',
                  scale = "row",
                  breaks = breaks,
                  fontsize_row=5, 
                  fontsize_col = 2,
                  main = "84 Genes differentially expressed between TNBC and ER + & HER2 +(FDR <0.01)"
)

