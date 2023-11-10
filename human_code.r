###cellranger###
###
/data/biosoft/cellranger-4.0.0/cellranger count --id=HRR232020 \
--transcriptome=/data/reference/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=/data/ljs2/HRR232020 \
--sample=HRR232020 \
--expect-cells=3000 \
--localcores=20

K='HRR'
idx=232020 
#!/bin/bash 
start=232020 \
end=232047 \
for (( i=$start; i<=$end; i++ )) \
do \
id =HRR$id \
cellranger count --id=$id \
--transcriptome=/data/reference/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=/data/ljs2/HRR232020/$id \
--sample=$id \
--expect-cells=3000 \
--localcores=20 \
done

###
for(i in 0:27){
  address = paste('/data/ljs2/HRA000863/cellranger/HRR2320',i+20,'/outs/filtered_feature_bc_matrix',sep = '')
  data_conuts <- Read10X(address)
  projectname <- paste('HRR2320',i,sep = '')
  seruatobject <- CreateSeuratObject(data_conuts,project = projectname,assay = 'RNA',min.cells = 3, min.features = 200)
  rdsfile <- paste('/data/ljs2/HRA000863/HRA000863_data/',projectname,'.rds',sep = '')
  
  seruatobject[["percent.mt"]] <- PercentageFeatureSet(seruatobject, pattern = "^MT-")
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB.genes, rownames(seruatobject@assays$RNA)) 
  HB.genes <- rownames(seruatobject@assays$RNA)[HB_m] 
  HB.genes <- HB.genes[!is.na(HB.genes)] 
  seruatobject[["percent.HB"]]<-PercentageFeatureSet(seruatobject, features=HB.genes) 
  seruatobject <- subset(seruatobject, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20 & percent.HB < 20)
  seruatobject <- NormalizeData(seruatobject)
  seruatobject <- FindVariableFeatures(seruatobject, selection.method = "vst", nfeatures = 2000)
  seruatobject <- ScaleData(seruatobject, features = VariableFeatures(object = seruatobject))
  seruatobject <- RunPCA(seruatobject, verbose = FALSE)
  seruatobject <- RunUMAP(seruatobject, dims = 1:30)
  sweep.res.list <- paramSweep_v3(seruatobject, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- seruatobject@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)          
  nExp_poi <- round(0.075*nrow(seruatobject@meta.data)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seruatobject <- doubletFinder_v3(seruatobject, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  seruatobject <- doubletFinder_v3(seruatobject, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = colnames(seruatobject@meta.data)[6], sct = FALSE)
  colnames(seruatobject@meta.data)[ncol(seruatobject@meta.data)] = "doublet_info"
  seruatobject <- subset(seruatobject, subset = doublet_info == "Singlet")
  saveRDS(seruatobject,file = rdsfile)
}

###
HRA000863_icc <- merge(x = after_HRR232020, y = list (after_HRR232021,after_HRR232022,after_HRR232023,after_HRR232024,after_HRR232025,after_HRR232026,after_HRR232027,after_HRR232028,after_HRR232029,after_HRR232030,after_HRR232031
                                                      ,after_HRR232032,after_HRR232033,after_HRR232034,after_HRR232035,after_HRR232036,after_HRR232037,
                                                      after_HRR232038,after_HRR232039,after_HRR232040,after_HRR232041,after_HRR232042,after_HRR232043,after_HRR232044,after_HRR232045,after_HRR232046,after_HRR232047)) 
HRA000863_icc <- NormalizeData(HRA000863_icc, normalization.method = "LogNormalize", scale.factor = 10000)
HRA000863_icc <- FindVariableFeatures(HRA000863_icc, selection.method = "vst", nfeatures = 2000)
HRA000863_icc <-ScaleData(HRA000863_icc)%>%RunPCA(verbose=FALSE)
system.time({HRA000863_icc <- RunHarmony(HRA000863_icc, group.by.vars = "orig.ident")})
HRA000863_icc <-FindNeighbors(HRA000863_icc, reduction = "harmony",dims = 1:30)
HRA000863_icc <-FindClusters(HRA000863_icc, reduction = "harmony",resolution = 0.1)
HRA000863_icc <-RunUMAP(HRA000863_icc,dims = 1:30,reduction = "harmony")
HRA000863_icc <-RunTSNE(HRA000863_icc,dims = 1:30,reduction = "harmony")
DimPlot(HRA000863_icc, reduction = "umap", pt.size = 0.5,label = T)
HRA000863_icc$keep <- ifelse(HRA000863_icc$seurat_clusters %in% c(15),"not","keep")
HRA000863_icc<- subset(HRA000863_icc, subset = keep == "keep")
HRA000863_icc@meta.data$new.lables<-"T.cells"
HRA000863_icc@meta.data$new.lables[which(HRA000863_icc@meta.data$seurat_clusters %in% c(2,5,9,10,12,15))]<-"Epithelial.cells"#
HRA000863_icc@meta.data$new.lables[which(HRA000863_icc@meta.data$seurat_clusters %in% c(8))]<-"Endothelial.cells"#
HRA000863_icc@meta.data$new.lables[which(HRA000863_icc@meta.data$seurat_clusters %in% c(3))]<-"Neutrophils"
HRA000863_icc@meta.data$new.lables[which(HRA000863_icc@meta.data$seurat_clusters %in% c(4,11))]<-"Myeloid.cells"
HRA000863_icc@meta.data$new.lables[which(HRA000863_icc@meta.data$seurat_clusters %in% c(7))]<-"B.cells"
Idents(HRA000863_icc) <-'new.lables'
DimPlot(HRA000863_icc,reduction = "umap", pt.size = 0.5,label = T)

#######################
Tcells_H<-subset(bca.integrated,new.lables=='T.cells')

library("harmony")
Tcells_H <- NormalizeData(Tcells_H, normalization.method = "LogNormalize", scale.factor = 10000)
Tcells_H <- FindVariableFeatures(Tcells_H, selection.method = "vst", nfeatures = 2000)
Tcells_H<-ScaleData(Tcells_H)%>%RunPCA(verbose=FALSE)
system.time({Tcells_H <- RunHarmony(Tcells_H, group.by.vars = "orig.ident")})

Tcells_H<-FindNeighbors(Tcells_H, reduction = "harmony",dims = 1:30)
Tcells_H<-FindClusters(Tcells_H, reduction = "harmony",resolution = 0.3)
Tcells_H<-RunUMAP(Tcells_H,dims = 1:30,reduction = "harmony")
#Tcells_H<-RunTSNE(Tcells_H,dims = 1:30,reduction = "harmony")

DimPlot(Tcells_H, reduction = "umap", pt.size = 0.5,label = T)



Tcells_H$keep <- ifelse(Tcells_H$seurat_clusters %in% c(4),"not","keep")
Tcells_H<- subset(Tcells_H, subset = keep == "keep")



# DefaultAssay(CLUSTER4) <- "RNA"
# Idents(CLUSTER4)<-'seurat_clusters'
# CLUSTER4.markers <- FindAllMarkers(CLUSTER4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# write.csv(CLUSTER4.markers,"CLUSTER4.markers0.5.csv")




FeaturePlot(Tcells_H, features = c('CD3D','CD3E','CD4','CD8A'),reduction = "umap",label = T)  #CD8T
FeaturePlot(Tcells_H, features = c('PTPRC','CD8A','CD8B'),reduction = "umap",label = T)  #CD8T
FeaturePlot(Tcells_H, features = c('CD4','CD40LG'),reduction = "umap",label = T)  #cd4T
FeaturePlot(Tcells_H, features = c('IL2RA','FOXP3'),reduction = "umap",label = T)  #Treg
FeaturePlot(Tcells_H, features = c('TRDV2','TRGV9'),reduction = "umap",label = T)  #gdT
FeaturePlot(Tcells_H, features = c('NKG7', 'GNLY','KLRD1','KLRF1','NCAM1'),reduction = "umap",label = T)  #NK
FeaturePlot(Tcells_H, features = c('CD4', 'CD8A'),reduction = "umap",label = T)  #NK


Tcells_H@meta.data$new.lables<-"CD8.T.cells"
#Tcells_H@meta.data$new.lables[which(Tcells_H@meta.data$seurat_clusters %in% c(14))]<-"Stromal cells"#
Tcells_H@meta.data$new.lables[which(Tcells_H@meta.data$seurat_clusters %in% c(0,8,7))]<-"CD4.T.cells"#
Tcells_H@meta.data$new.lables[which(Tcells_H@meta.data$seurat_clusters %in% c(4))]<-"Tregs"#
Tcells_H@meta.data$new.lables[which(Tcells_H@meta.data$seurat_clusters %in% c(3,6))]<-"NKT.cells"#
Tcells_H@meta.data$new.lables[which(Tcells_H@meta.data$seurat_clusters %in% c(5,10))]<-"γδT.cells"#

table(Tcells_H$new.lables)
table(Tcells$new.lables)

#################################################Myeloid
MoMfDC_H<-subset(bca.integrated,new.lables=='Myeloid.cells')

library("harmony")
MoMfDC_H <- NormalizeData(MoMfDC_H, normalization.method = "LogNormalize", scale.factor = 10000)
MoMfDC_H <- FindVariableFeatures(MoMfDC_H, selection.method = "vst", nfeatures = 2000)
MoMfDC_H<-ScaleData(MoMfDC_H)%>%RunPCA(verbose=FALSE)
system.time({MoMfDC_H <- RunHarmony(MoMfDC_H, group.by.vars = "orig.ident")})
MoMfDC_H<-FindNeighbors(MoMfDC_H, reduction = "harmony",dims = 1:30)
MoMfDC_H<-FindClusters(MoMfDC_H, reduction = "harmony",resolution = 0.3)
MoMfDC_H<-RunUMAP(MoMfDC_H,dims = 1:30,reduction = "harmony")
DimPlot(MoMfDC_H, reduction = "umap", pt.size = 0.5,label = T)
MoMfDC_H$keep <- ifelse(MoMfDC_H$seurat_clusters %in% c(6),"not","keep")
MoMfDC_H<- subset(MoMfDC_H, subset = keep == "keep")
MoMfDC_H@meta.data$new.lables<-"DCs"
MoMfDC_H@meta.data$new.lables[which(MoMfDC_H@meta.data$seurat_clusters %in% c(0,4,8,2))]<-"Macrophages"#
MoMfDC_H@meta.data$new.lables[which(MoMfDC_H@meta.data$seurat_clusters %in% c(1,7))]<-"Monocytes"#
table(MoMfDC_H$new.lables)
table(Momfdc_newlabels_0709$new.lables)
table(MoMfDC_H$new.lables)
T_cell<-merge(MoMfDC,MoMfDC_H)
table(T_cell$orig.ident)
table(T_cell$new.lables)
table(MoMfDC$new.lables)
MoMfDC$new.lables<-ifelse(MoMfDC$new.lables%in%c('DCs'),'M_DCs',
                          ifelse(MoMfDC$new.lables%in%c('Macrophages'),'M_Macrophages','M_Monocytes'))
table(MoMfDC$new.lables)
table(MoMfDC_H$new.lables)
MoMfDC_H$new.lables<-ifelse(MoMfDC_H$new.lables%in%c('DCs'),'H_DCs',
                            ifelse(MoMfDC_H$new.lables%in%c('Macrophages'),'H_Macrophages','H_Monocytes'))
table(MoMfDC_H$new.lables)
Myeloids<-merge(MoMfDC,MoMfDC_H)
table(Tcells_H$new.lables)
DefaultAssay(MoMfDC) <- "RNA"
Idents(MoMfDC)<-'new.lables'
MoMfDC.markers <- FindAllMarkers(MoMfDC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(MoMfDC.markers,"MoMfDC.markers.csv")
top5_MoMfDC_H <- MoMfDC.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
DefaultAssay(MoMfDC_H) <- "RNA"
Idents(MoMfDC_H)<-'new.lables'
MoMfDC_H.markers <- FindAllMarkers(MoMfDC_H, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(MoMfDC_H.markers,"MoMfDC_H.markers.csv")
top5_MoMfDC_H <- MoMfDC_H.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(Tcells.markers,'Tcells.markers.csv')
write.csv(Tcells_H.markers,'Tcells_H.markers.csv')
table(MoMfDC_H$new.lables)
Idents(Tcells)<-Tcells$new.lables
mean_exp_matrix_MM<-AverageExpression(MoMfDC,return.seurat = T)
mean_exp_matrix_HS<-AverageExpression(MoMfDC_H,return.seurat = T)
df<-as.data.frame(mean_exp_matrix_MM@assays$RNA@counts)
df<-as.data.frame(mean_exp_matrix_HS@assays$RNA@counts)
MM_list<-read.csv('鼠髓系list.csv')
HS_list<-read.csv('人髓系list.csv')
df<-df[MM_list$gene,]
df<-df[HS_list$gene,]
df<-df[top5_MoMfDC_H$gene,]
colnames(df)<-factor(colnames(df),levels = c("CD8.T.cells","γδT.cells","NKT.cells","CD4.T.cells","Tregs","Cycling.T.cells"))
write.csv(df,'df.csv')
df<-read.csv('df.csv')
rownames(df)<-df[,1]
df<-df[,-1]
mat_MM_MM = scale(df, center = TRUE, scale = TRUE)
mat_HS_MM = scale(df, center = TRUE, scale = TRUE)
pheatmap(mat_MM_MM,
         show_rownames = F,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows=F,
         filename='test.pdf',#输出文件的名称
         fontsize_row=30, 
         fontsize_col =30,#行字体的大小
         height=15,  #输出图片的高度
         scale = "row",
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         clustering_distance_rows = 'euclidean', 
         clustering_method = 'single',
         #gaps_row = c(1:27),
         #gaps_col = c(1:8),
         #border_color = '#222222'
)
DefaultAssay(Tcells) <- "RNA"
Idents(Tcells)<-'new.lables'
Tcells.markers <- FindAllMarkers(Tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tcells.markers,"Tcells.markers.csv")
top5_Tcells <- Tcells.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
DefaultAssay(Tcells_H) <- "RNA"
Idents(Tcells_H)<-'new.lables'
Tcells_H.markers <- FindAllMarkers(Tcells_H, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Tcells_H.markers,"Tcells_H.markers.csv")
top5_Tcells_H <- Tcells_H.markers %>% group_by(cluster) %>% top_n(n=30, wt = avg_log2FC)
write.csv(top5_Tcells_H,'top5_Tcells_H.csv')
write.csv(Tcells.markers,'Tcells.markers.csv')
write.csv(Tcells_H.markers,'Tcells_H.markers.csv')
table(Tcells$new.lables)
Idents(Tcells)<-Tcells$new.lables
mean_exp_matrix_MM<-AverageExpression(Tcells,return.seurat = T)
mean_exp_matrix_HS<-AverageExpression(Tcells_H,return.seurat = T)
df<-as.data.frame(mean_exp_matrix_MM@assays$RNA@counts)
df<-as.data.frame(mean_exp_matrix_HS@assays$RNA@counts)
MM_list<-read.csv('鼠Tlist.csv')
HS_list<-read.csv('人Tlist.csv')
df<-df[MM_list$gene,]
df<-df[HS_list$gene,]
df<-df[top5_Tcells_H$gene,]
colnames(df)<-factor(colnames(df),levels = c("CD8.T.cells","γδT.cells","NKT.cells","CD4.T.cells","Tregs","Cycling.T.cells"))
write.csv(df,'df.csv')
df<-read.csv('df.csv')
rownames(df)<-df[,1]
df<-df[,-1]
mat_MM = scale(df, center = TRUE, scale = TRUE)
mat_HS = scale(df, center = TRUE, scale = TRUE)
write.csv(mat_HS,'mat_HS.csv')
write.csv(mat_MM,'mat_MM.csv')
mat_all<-read.csv('mat_HS.csv')
rownames(mat_all)<-mat_all[,1]
mat_all<-mat_all[,-1]
pheatmap(mat_HS,
         show_rownames = F,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows=F,
         filename='test.pdf',#输出文件的名称
         fontsize_row=30, 
         fontsize_col =30,#行字体的大小
         height=15,  #输出图片的高度
         scale = "row",
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         clustering_distance_rows = 'euclidean', 
         clustering_method = 'single',
         #gaps_row = c(1:27),
         #gaps_col = c(1:8),
         #border_color = '#222222'
)

table(Myeloids$orig.ident)
table(Myeloids$new.lables)
table(MoMfDC$new.lables)
MoMfDC$new.lables<-ifelse(MoMfDC$new.lables%in%c('CD4.T.cells'),'M_CD4.T.cells',
                          ifelse(MoMfDC$new.lables%in%c('CD8.T.cells'),'M_CD8.T.cells',
                                 ifelse(MoMfDC$new.lables%in%c('Cycling.T.cells'),'M_Cycling.T.cells',
                                        ifelse(MoMfDC$new.lables%in%c('NKT.cells'),'M_NKT.cells',
                                               ifelse(MoMfDC$new.lables%in%c('Tregs'),'M_Tregs','M_M_γδT.cells')))))
table(MoMfDC$new.lables)


table(MoMfDC_H$new.lables)
MoMfDC_H$new.lables<-ifelse(MoMfDC_H$new.lables%in%c('CD4.T.cells'),'H_CD4.T.cells',
                            ifelse(MoMfDC_H$new.lables%in%c('CD8.T.cells'),'H_CD8.T.cells',
                                   ifelse(MoMfDC_H$new.lables%in%c('Tregs'),'H_Tregs',
                                          ifelse(MoMfDC_H$new.lables%in%c('NKT.cells'),'H_NKT.cells','H_γδT.cells'))))
table(MoMfDC_H$new.lables)
Myeloids<-merge(MoMfDC,MoMfDC_H)
table(Myeloids$new.lables)

DefaultAssay(Myeloids) <- "RNA"
Idents(Myeloids)<-'new.lables'
Myeloids.markers <- FindAllMarkers(Myeloids, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- Myeloids.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
Myeloids@assays$RNA@scale.data <- scale(Myeloids@assays$RNA@data, scale = TRUE)
DoHeatmap(Myeloids, features = top5$gene, label=F,group.by = "new.lables") #+ NoLegend()
table(Neutrophils$new.lables)
Neutrophils$new.lables<-ifelse(Neutrophils$new.lables%in%c('Neutrophils'),'M_Neutrophils','NA')
table(Neutrophils$new.lables)
table(Neutrophils_H$new.lables)
Neutrophils_H$new.lables<-ifelse(Neutrophils_H$new.lables%in%c('Neutrophils'),'H_Neutrophils','NA')
table(Neutrophils_H$new.lables)
Neutrophils_all<-merge(Neutrophils_H,Neutrophils)
table(MoMfDC_H$new.lables)
MoMfDC_H$new.lables<-ifelse(MoMfDC_H$new.lables%in%c('CD4.T.cells'),'H_CD4.T.cells',
                            ifelse(MoMfDC_H$new.lables%in%c('CD8.T.cells'),'H_CD8.T.cells',
                                   ifelse(MoMfDC_H$new.lables%in%c('Tregs'),'H_Tregs',
                                          ifelse(MoMfDC_H$new.lables%in%c('NKT.cells'),'H_NKT.cells','H_γδT.cells'))))
table(MoMfDC_H$new.lables)
Myeloids<-merge(MoMfDC,MoMfDC_H)
table(Myeloids$new.lables)
DefaultAssay(Myeloids) <- "RNA"
Idents(Myeloids)<-'new.lables'
Myeloids.markers <- FindAllMarkers(Myeloids, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- Myeloids.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
Myeloids@assays$RNA@scale.data <- scale(Myeloids@assays$RNA@data, scale = TRUE)
DoHeatmap(Myeloids, features = top5$gene, label=F,group.by = "new.lables") #+ NoLegend()
table(Epi$new.lables)
table(Tcells$new.lables)
table(MoMfDC$new.lables)
table(B_cells$new.lables)
table(Neutrophils$new.lables)
table(Endothelial_cells$new.lables)
table(Epi_H$new.lables)
table(Tcells_H$new.lables)
table(MoMfDC_H$new.lables)
table(B_cells_H$new.lables)
table(Neutrophils_H$new.lables)
table(Endothelial_cells_H$new.lables)
MoMfDC$new.lables<-ifelse(MoMfDC$new.lables%in%c('DC1','DC2','Modc','MoMfDC'),'DCs',
                          ifelse(MoMfDC$new.lables%in%c('Monocyte'),'Monocytes','Macrophages'))
Idents(MoMfDC)<-MoMfDC$new.lables
DimPlot(MoMfDC,label = T)
Epi_H<-subset(bca.integrated,new.lables=='Epithelial.cells')
B_cells_H<-subset(bca.integrated,new.lables=='B.cells')
Neutrophils_H<-subset(bca.integrated,new.lables=='Neutrophils')
Endothelial_cells_H<-subset(bca.integrated,new.lables=='Endothelial.cells')
ICC_MM<-merge(x = Epi, y = list (Tcells,MoMfDC,B_cells,Neutrophils,Endothelial_cells))
table(ICC_MM$new.lables)
freq <- xtable(table(ICC_MM@meta.data$new.lables))
freq$terms<-rownames(freq)
colnames(freq)<-c("counts","terms")
freq<-as.data.frame(freq)
write.csv(freq,'ICC_MM.csv')
ICC_HS<-merge(x = Epi_H, y = list (Tcells_H,MoMfDC_H,B_cells_H,Neutrophils_H,Endothelial_cells_H))            
table(ICC_HS$new.lables)
table(bca.integrated$new.lables)
table(icc$new.lables)
freq <- xtable(table(ICC_HS@meta.data$new.lables))
freq$terms<-rownames(freq)
colnames(freq)<-c("counts","terms")
freq<-as.data.frame(freq)              
write.csv(freq,'ICC_HS.csv')             
install.packages('d3Network')
install.packages('xlsx')
library(d3Network)   
library(xlsx)
Sankey<-read.csv('ICC_ALL.csv')
Sankeylinks<-Sankey
Sankeynodes<-data.frame(name=unique(c(Sankeylinks$Source,Sankeylinks$Target)),stringsAsFactors=FALSE)  
Sankeynodes$index<-0:(nrow(Sankeynodes) - 1)
Sankeylinks<-merge(Sankeylinks,Sankeynodes,by.x="Source",by.y="name")
Sankeylinks<-merge(Sankeylinks,Sankeynodes,by.x="Target",by.y="name")
Sankeydata<-Sankeylinks[,c(4,5,3)];names(Sankeydata)<-c("Source","Target","Value")
Sankeyname<-Sankeynodes[,1,drop=FALSE]
d3Sankey(Links=Sankeydata,Nodes=Sankeyname,Source="Source",Target="Target",Value="Value",NodeID="name",  
         fontsize=12,nodeWidth=30,file="TestSankey.html")