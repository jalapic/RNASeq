#RNAseq analysis with DESeq 2 package 
# group analysis with both brain regions. 

# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering

#libraries 
library(tidyverse)
library(DESeq2)
library(pheatmap)
#devtools::install_github("stephenturner/annotables")
library(annotables)
grcm38 # mouse genes from annotables package.


#Getting metadata ready 

coldata #see other tutorial.


# First AMY data 

#Creating DESeq object for each brain region & getting rid of low count data 
# AMY
dag <- DESeqDataSetFromMatrix(countData = countData,
                              colData = coldata,
                              design = ~ status)

nrow(dag) #13913


# Normalized counts - but just for plotting.
amy<- estimateSizeFactors(dag)
sizeFactors(amy)  #base means

#Hierarchical heatmap
#Transform the normalized counts 
vsd_amy <-vst(amy, blind = T)

#Extract the matrix of transformed counts 
vsd_mat_amy <- assay(vsd_amy)

#Compute the correlation values between samples
vsd_cor_amy <- cor(vsd_mat_amy)

#Plot the heatmap based off group and condition
df1 <- as.data.frame(coldata[,c("status")], row.names=rownames(coldata))
rownames(df1) <- colnames(vsd_mat_amy)

#amy by group
pheatmap(vsd_cor_amy,
         cluster_rows = FALSE,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         annotation=df1)


#Principal component analysis 
#need transform normalized counts which we saved as vsd_pfc and vsd_amy
plotPCA(vsd_amy, intgroup = "status")


##################-----------------------------------------####################



#Differential expression analysis 

#DESeq model for amy 
des_amy <- DESeq(dag)

plotDispEsts(des_amy)#This looks terrible lol  
# good fit dispersion decreases with increasing mean and 
# the raw dispersion seem to cluster around the maximum likelihood line.

#DESeq model for pfc
des_pfc <- DESeq(dpg)
plotDispEsts(des_pfc) #seriously kill me

# get results and map to gene library
# not really sure if I should be a using lfcthreshold idk, video said it is not preferred but not all genes could be relevant?

#AMY
amy1 <- results(des_amy, contrast=c("group", "reorganized", "control"), alpha = 0.05)
summary(amy1)

res1 <- results(des_amy, contrast=c("group", "reorganized", "control"), alpha = 0.05)%>% 
  as.data.frame()%>%
  rownames_to_column('ensgene') %>% 
  as_tibble() %>% 
  select(ensgene,log2FoldChange,pvalue,padj) %>% 
  mutate(contrast = "Reorganized - Control") %>% 
  mutate(region = "AMY") 
head(res1)

#PFC
pfc1 <- results(des_pfc, contrast=c("group", "reorganized", "control"), alpha = 0.05)
summary(pfc1)

res2 <-results(des_pfc, contrast=c("group", "reorganized", "control"), alpha = 0.05)%>% 
  as.data.frame()%>%
  rownames_to_column('ensgene') %>% 
  as_tibble() %>% 
  select(ensgene,log2FoldChange,pvalue,padj) %>% 
  mutate(contrast = "Reorganized - Control") %>% 
  mutate(region = "PFC")
head(res2)

#Getting dataframe with all information 
group <- res1 %>% rbind(res2) %>% arrange(desc(abs(log2FoldChange)))
head(group)

group2 <- left_join(x=group, y =grcm38[,c("ensgene", "symbol", "description")],
          by = "ensgene")

group1 <- group2 %>% filter(pvalue <0.05)
groupadj <- group2 %>%  filter(padj <= 0.05)

#just getting them to csv for later use 
write.csv(group2,"manuscript/brain/genes/GenesbyGroup.csv")
write.csv(groupadj,"manuscript/brain/genes/GenesbyGroupAdj.csv")

#more plots
xx <- lfcShrink(des_amy, 
                      coef="group_reorganized_vs_control", 
                      type="apeglm")
# MA plot
plotMA(xx) #lol

xx1 <- lfcShrink(des_pfc, 
                coef="group_reorganized_vs_control", 
                type="apeglm")
# MA plot
plotMA(xx1) #lol

#count plot with normalized data 
a <- plotCounts(amy, gene=which.min(res1$padj), intgroup="group",returnData=TRUE)

ggplot(a, aes(x=group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ylab("Normalized Count")+
  xlab("Group")+
  theme_minimal()


p <- plotCounts(pfc, gene=which.min(res2$padj), intgroup="group",returnData=TRUE)

ggplot(a, aes(x=group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10()+
  ylab("Normalized Count")+
  xlab("Group")+
  theme_minimal()




##tick plot for group and reorganized I am going to use unadjusted pvalues for now. 
head(group1)

# code from Hannah and Won
# alright explored enough now make figure
group1 <- group1 %>% filter(pvalue <0.5)

group1%>% arrange(log2FoldChange) %>%
  split(.$region) %>%
  map(~mutate(.,real_LFC = log2FoldChange)) %>%
  map(~mutate(.,LFC_with_NAs = log2FoldChange)) %>%
  map(~select(.,ensgene, log2FoldChange,contrast)) %>%
  map(~pivot_wider(.,names_from = contrast,values_from=log2FoldChange)) %>%
  map(~arrange(.,`Reorganized - Control`))%>%
  map(~mutate(.,x_var = row_number())) %>%
  map(~gather(.,contrast,log2FoldChange,2)) %>%
  map2_df(.,names(.), ~mutate(.x, region = .y))-> DEG_result_for_plot


table(group1$region) %>%
  as.data.frame() %>%
  dplyr::mutate(facet_it = glue::glue('{Var1} ({Freq})')) %>%
  select(region = Var1, facet_it) -> freq1


table(group1$region, group1$contrast) %>%
  as.data.frame() %>%
  dplyr::rename(region = Var1) %>%
  mutate(each_region = glue::glue('({Freq})')) %>%
  left_join(freq1) %>%
  select(facet_it, contrast= Var2, each_region) -> freq2

DEG_result_for_plot %>%
  left_join(freq1) %>%
  ggplot(aes(x_var, contrast, fill = log2FoldChange))+
  geom_tile(color = NA, size = 10)+
  facet_grid(facet_it ~., switch = "y")+
  labs(x = "",
       y = "",
       fill = "Log2 Fold Change")+
  ggtitle("Reorganized - Control")+
  scale_fill_gradient2(low="purple4", mid="white", high="orange2", #colors in the scale
                       midpoint=0,    #same midpoint for plots (mean of the range)
                       limits=c(-2, 2),#same limits for plots
                       breaks=seq(-2,2,1),   #breaks in the scale bar
                       na.value = 'white')+
  theme(legend.position = "bottom",
        axis.text.x= element_blank(),
        axis.text.y= element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank(),
        # plot.background =  element_rect(color = "black"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        strip.text.y.left = element_text(angle = 0, size = 12))+
        theme(plot.title = element_text(hjust = 0.5)) -> RC_0.5


# ggsave(filename = "manuscript/brain/img/RC_0.5.png",
#        RC_0.5,
#        height = 6,
#        width = 8,
#        dpi = 100)
# 

##Volcano Plots 
# Generate logical column
res_all <- group2 %>% 
  mutate(threshold = padj <0.05)

# Create the volcano plot
ggplot(res_all) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  ggtitle("Reorganized - Control")+
  facet_wrap(~region) +
  theme_minimal()+
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25))) 

#labeled

library(ggrepel)
res_all$Significant <- ifelse(res_all$padj < 0.05, "FDR < 0.05", "Not Sig")
VolGroup <- ggplot(res_all, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("red", "grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  facet_wrap(~region) +
  geom_text_repel(
    data = subset(res_all, padj < 0.05),
    aes(label = symbol),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.63, "lines"))
# 
# ggsave(filename = "manuscript/brain/img/volcano_group.png",
#        VolGroup,
#        height = 6,
#        width = 6,
#        dpi = 100)


head(res_all)
#count data for top 11 genes 
head(group2)

a <- group2 %>% filter(region == "AMY" & padj <0.05)
p <- group2 %>% filter(region == "PFC" & padj <0.05)


amy_normalized <- counts(des_amy, normalized =T)
res_a <- res1 %>%  filter(padj<0.05)
counts1 <- amy_normalized[res_a$ensgene,]
meta_a <-amy_data1 %>%  rownames_to_column('samplename')

atop <- as.data.frame(counts1)%>% 
    rownames_to_column('ensgene') %>% 
    as_tibble(.) %>% 
    gather(.,key= "samplename", value="normalized_counts", 2:68) %>% 
    full_join(meta_a) %>% 
    full_join(.,a)

unique(atop$samplename)#correct. 


pfc_normalized <- counts(des_pfc, normalized =T)
res_p <- res2 %>%  filter(padj<0.05)
meta_p <-pfc_data1 %>%  rownames_to_column('samplename')

counts <- pfc_normalized[res_p$ensgene,] %>%
  as.data.frame(.) %>% rownames_to_column("samplename") %>% 
  as.tibble(.) %>% 
  select(samplename, normalized_counts = '.') %>% 
  mutate(ensgene = "ENSMUSG00000021250") %>% 
  select(ensgene,samplename,normalized_counts)

ptop <- counts %>% 
  full_join(meta_p, by ="samplename") %>% 
  full_join(.,p) 

unique(ptop$samplename) #correct


topgene_count <- atop %>% rbind(ptop) %>% select(normalized_counts,region,group,symbol)
head(topgene_count)
table(topgene_count$group,topgene_count$symbol)

FOS <- topgene_count %>% filter(symbol == "Fos")
AMY <- topgene_count %>% filter(symbol != "Fos")



diff_Fos <- ggplot(FOS, aes(x=region,y=normalized_counts, color=group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.8,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "",
       x = "",
       y = "Normalized Counts",
       color = "Group") +
  scale_y_log10()+
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  newggtheme +
  annotate(geom="text", x=.6, y=300, label="Fos",size = 8,
           color="Black")
# ggsave(filename = "manuscript/brain/img/FOS_group.png",
#        diff_Fos,
#        height = 6,
#        width = 8,
#        dpi = 100)


AMY_diff <- ggplot(AMY, aes(x=group,y=normalized_counts, color=group))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21,
                 alpha = 0.8,
                 jitter.height = 0.02, jitter.width = 0.030, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  labs(title = "",
       x = "",
       y = "Normalized Counts") +
  scale_y_log10()+
  scale_color_manual(values = viridis::viridis(4)) +
  scale_fill_manual(values = viridis::viridis(4))+
  facet_wrap(~symbol)+
  newggtheme+
  theme(legend.position = "none")

# ggsave(filename = "manuscript/brain/img/AMY_group.png",
#        AMY_diff,
#        height = 9,
#        width = 12,
#        dpi = 100)


########## pathway analysis
library(goseq)
library(clusterProfiler)
library(enrichplot)
library(biomaRt)
organism = 'org.Mm.eg.db'#obivously mouse
library(organism, character.only = TRUE)
library(DOSE)

a <- group1 %>% filter(region == "AMY")
p <- group1 %>% filter(region != "AMY")


#Go enrichment 
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

# AMY 1st
ent_gene_a <- getBM(attributes=c("entrezgene_id"),
                  filters=c('ensembl_gene_id'),
                  values= a$ensgene,
                  mart =ensembl )

ent_gene_a <-as.character(ent_gene_a$entrezgene_id)

#PFC
ent_gene_p <- getBM(attributes=c("entrezgene_id"),
                  filters=c('ensembl_gene_id'),
                  values= p$ensgene,
                  mart =ensembl )

ent_gene_p <-as.character(ent_gene_p$entrezgene_id)


#universe
ent_uni <- getBM(attributes=c("entrezgene_id"),
                 filters=c('ensembl_gene_id'),
                 values= group2$ensgene,
                 mart =ensembl )
ent_uni<-as.character(ent_uni$entrezgene_id)
#the universe might be wrong because I have pfc and amy genes in it but I still think it is better than using the whole mouse genome
AMY_ego_BP<- enrichGO(
  gene = ent_gene_a,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  universe = ent_uni)
View(summary(AMY_ego_BP)) 

PFC_ego_BP<- enrichGO(
  gene = ent_gene_p,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  universe = ent_uni)

View(summary(PFC_ego_BP))
dotplot(PFC_ego_BP, showCategory=20) + 
  ggtitle("Reorganized - Control (PFC)") +
  theme(plot.title = element_text(hjust = 0.5))
cnetplot(ego_BP)
goplot(ego_BP)


ego_CC<- enrichGO(
  gene = ent_gene_a,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  universe = ent_uni)

View(summary(ego_CC)) #nothing 


PFC_ego_CC<- enrichGO(
  gene = ent_gene_p,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  universe = ent_uni)

View(summary(PFC_ego_CC))
dotplot(PFC_ego_CC, showCategory=20)+ 
  ggtitle("Reorganized - Control (PFC)") +
  theme(plot.title = element_text(hjust = 0.5))
cnetplot(PFC_ego_CC)
goplot(PFC_ego_CC)


ego_MF<- enrichGO(
  gene = ent_gene_a,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  universe = ent_uni)
View(summary(ego_MF)) #nothing 

PFC_ego_MF<- enrichGO(
  gene = ent_gene_p,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  universe = ent_uni)

View(summary(PFC_ego_MF))
dotplot(PFC_ego_MF, showCategory=20)
cnetplot(PFC_ego_MF)
goplot(PFC_ego_MF)

###KEGG
ekga <- enrichKEGG(gene = ent_gene_p,
                    organism = 	'mmu',
                    universe = ent_uni)
View(summary(ekga))
dotplot(ekga, showCategory=20)
cnetplot(ekga)


