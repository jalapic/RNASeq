i = 4
filter_counts = 50

filter_counts = 10 # dds2

filter_counts = 20 # dds3

filter_counts = 100

library(DESeq2)

for (i in c(4,6)){
  


regions[[i]] -> my_region

df <- rawcount_list[[i]] %>% 
  dplyr::select(-ensgene)

coldata = behav %>% 
  filter(sampleID %in% colnames(df)) %>% 
  dplyr::select(sampleID, subjectID, status) %>% 
  # filter(status == "Alpha") %>% # # despotism - alpha only 
  left_join(alldata)

# # despotism - alpha only 
# df <- df %>% 
#   select_if(colnames(.) %in% coldata$sampleID)

countData <- as.matrix(df) 
rownames(countData) <- rawcount_list[[i]]$ensgene 


countData <- countData[!is.na(rowSums(countData)),]
countData %>% is.na() %>% sum

dim(countData)

countData <- countData[rowSums(countData > filter_counts) > round((length(df))*0.9), ]

dim(countData)




dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = coldata,
                              design = ~status)



# # despotism - alpha only 
# dds <- DESeqDataSetFromMatrix(countData = countData,
#                               colData = coldata,
#                               design = ~despotism)


dds <- DESeq(dds)

# saveRDS(dds, glue("results_RNAseqRDS/dds2_{my_region}.RDS"))


vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup = c("status"), returnData = T) -> d
attributes(d)
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
  geom_point(size = 3, alpha = 0.5) +
  # geom_text(aes(label = subjectID), vjust = -0.5)+
  scale_color_manual(values = c("purple4", "#21908CFF","orange")) +
  # scale_fill_manual(values = c("purple4", "#21908CFF","orange")) +
  labs(x = paste0("PC1: ",
                  round(attributes(d)$percentVar[1] * 100),
                  "% variance"),
       y = paste0("PC2: ",
                  round(attributes(d)$percentVar[2] * 100),
                  "% variance"),
       color = "Social status",
       fill = "Social status") +
  theme(legend.position = "right",
    # legend.position = c(0.11,0.87),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))+
  ggtitle(glue('{my_region} PCA plot')) -> p



theme_set(theme_bw(base_size = 10))

png(filename = glue("results_figures/PCA_{my_region}_filtercounts{filter_counts}.png"),
    width = 10.5, height = 7.3, units = "cm", res = 600)
print(p)
invisible(dev.off())


  # CORT 

# vsd <- vst(dds, blind=FALSE)
# plotPCA(vsd, intgroup = c("status"), returnData = T) %>%
#   rename(sampleID = name) %>% 
#   left_join(coldata)-> d
# 
# ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "cort_post")) +
#   geom_point(size = 3, alpha = 0.7) +
#   scale_color_distiller(palette = "Spectral") +
#   # scale_fill_distiller(values = c("purple4", "#21908CFF","orange")) +
#   labs(x = paste0("PC1: ",
#                   round(attributes(d)$percentVar[1] * 100),
#                   "% variance"),
#        y = paste0("PC2: ",
#                   round(attributes(d)$percentVar[2] * 100),
#                   "% variance"),
#        color = "Corticosterone (ng/ml)")+
#   theme(legend.position = "right",
#     # legend.position = c(0.11,0.87),
#         legend.key.size = unit(0.4, 'cm'),
#         legend.text = element_text(size=8),
#         legend.title = element_text(size=8))+
#   ggtitle(glue('{my_region} PCA plot')) -> p2
# 
# 
# 
# theme_set(theme_bw(base_size = 10))
# 
# png(filename = glue("results_figures/PCA_CORT_{my_region}.png"),
#     width = 12, height = 9, units = "cm", res = 600)
# print(p2)
# invisible(dev.off())

# # alpha only despotism
# vsd <- vst(dds, blind=FALSE)
# plotPCA(vsd, intgroup = c("despotism"), returnData = T) %>%
#   left_join(coldata %>%
#               rename(name = subjectID))-> d
# 
# ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "despotism")) +
#   geom_point(size = 3, alpha = 0.5) +
#   # geom_text(aes(label = subjectID), vjust = -0.5)+
#   scale_color_distiller(palette = "Spectral") +
#   labs(x = paste0("PC1: ",
#                   round(attributes(d)$percentVar[1] * 100),
#                   "% variance"),
#        y = paste0("PC2: ",
#                   round(attributes(d)$percentVar[2] * 100),
#                   "% variance"),
#        color = "Despotism",
#        fill = "Social status") +
#   theme(legend.position = "right",
#         # legend.position = c(0.11,0.87),
#         legend.key.size = unit(0.4, 'cm'),
#         legend.text = element_text(size=8),
#         legend.title = element_text(size=8))+
#   ggtitle(glue('{my_region} PCA plot (alpha only)')) -> p
# 
# 
# 
# theme_set(theme_bw(base_size = 10))
# 
# png(filename = glue("results_figures/PCA_despotism_{my_region}.png"),
#     width = 10.5, height = 9, units = "cm", res = 600)
# print(p)
# invisible(dev.off())


}







# for all region - sanity check ============================================


df <- counts 
filter_counts = 10

colnames(df)
coldata = data.frame(sample_id = as.character(colnames(df)[2:length(df)])) %>% 
  left_join (.,behav %>% 
               select(sample_id,region, status) %>% 
               mutate_if(is.character,factor), by = 'sample_id') %>% 
  column_to_rownames(var = 'sample_id')


# tidy count data format  
d.raw <- as.matrix(df[2:length(df)])
rownames(d.raw) <- df$ensgene

which(complete.cases(d.raw) == F) -> x

d.noNA <- d.raw[-x,]
countData <- d.noNA[rowSums(d.noNA > filter_counts) > 2, ]
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = coldata,
                              design = ~status+region)
dds_deseq <- DESeq(dds)

vsd <- vst(dds_deseq, blind=FALSE)

colnames(dds)
colnames(vsd)


plotPCA(vsd,intgroup = c("region"), returnData = T) %>% 
  rownames_to_column('sample_id') %>% 
  left_join(coldata %>% 
              rownames_to_column('sample_id'))-> d

ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
  geom_point(size = 3, alpha = 0.35) + 
  # geom_text(aes(label = subjectID), vjust = -0.5)+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  labs(x = paste0("PC1: ", 
                  round(attributes(d)$percentVar[1] * 100), 
                  "% variance"),
       y = paste0("PC2: ", 
                  round(attributes(d)$percentVar[2] * 100), 
                  "% variance"),
       color = "",
       fill = "") +
  theme(legend.position = "right",
        # legend.position = c(0.11,0.87),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))+
  ggtitle(glue('All regions PCA plot')) -> p



theme_set(theme_bw(base_size = 10))

png(filename = glue("results_figures/PCA_allregions.png"),
    width = 11.5, height = 9, units = "cm", res = 600)
print(p)
invisible(dev.off())



