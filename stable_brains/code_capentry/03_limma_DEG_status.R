i = 2

# How many random sampling
R = 5000


for (i in 3:length(regions)){

regions[[i]] -> my_region

dlNorm <- rawcount_list[[i]] %>% 
  column_to_rownames('ensgene')

var_info <- behav %>% 
  filter(sampleID %in% colnames(dlNorm)) %>% 
  select(sampleID, subjectID, status)

colnames(dlNorm)

var_info$status %>%
  factor(.,levels = c("Alpha","Subdominant","Subordinate")) -> group.dl

dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)

d0= DGEList(d, group = group.dl)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)


design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)

contrast.matrix <- makeContrasts(group.dlAlpha-group.dlSubdominant,
                                 group.dlAlpha-group.dlSubordinate, 
                                 group.dlSubdominant-group.dlSubordinate,
                                 levels=design.dl)

vfit.dl2 <- contrasts.fit(vfit.dl, contrast.matrix)

efit.dl2 = eBayes(vfit.dl2)

p.dl.limma2 = efit.dl2[["p.value"]]
head(p.dl.limma2)

# saveRDS(v.dl, glue("results_RNAseqRDS/limma_vdl_{my_region}"))

p.dl.rand = vector('list',length = R)

for(g in 1 : R){
  print(paste("Starting on Permutation", g))

  # Randomize the traits

  group.dl.rand = sample(group.dl)

  # Model
  design.dl.rand = model.matrix(~0 + group.dl.rand)
  colnames(design.dl.rand) <- mycolnames

  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  
  vfit.dl.rand2 <- contrasts.fit(vfit.dl.rand, contrast.matrix)
  
  efit.dl.rand2 = eBayes(vfit.dl.rand2)
  
  p.dl.rand[[g]] = efit.dl.rand2[["p.value"]]
  head(p.dl.rand[[g]])
}

q.dl = matrix(0, nrow = nrow(p.dl.limma2), ncol = ncol(p.dl.limma2))

for(h in 1 : R){
  print(paste("Calculating Permutation", h))

  temp = p.dl.rand[[h]]

  for(c in 1 : 3){
    for(r in 1 : nrow(p.dl.limma2)){
      if(temp[r, c] <= p.dl.limma2[r, c]){
        q.dl[r, c] = q.dl[r, c] + 1
      }
    }
  }
}

q.dl = q.dl / R
colnames(q.dl) <- mycolnames
q.dl = as.data.frame(q.dl)
row.names(q.dl) <- rownames(dge.dl)

saveRDS(q.dl,glue("results_RNAseqRDS/limma_eFDR_{my_region}_cutoff5_R{R}_contrast.RDS"))

q.dl <- readRDS(glue("results_RNAseqRDS/limma_eFDR_{my_region}_cutoff5_R{R}_contrast.RDS"))


png(filename = glue("results_figures/eFDR_hist_{my_region}R{R}_contrast.png"),
    width = 18, height = 17, units = "cm", res = 600)
hist(q.dl-p.dl.limma2, main = glue("{my_region}"))
invisible(dev.off())

efit.dl2[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl2$coefficients)))

tmp1 <- contrasts.fit(efit.dl2, coef = 1) # 
tmp2 <- contrasts.fit(efit.dl2, coef = 2) # 
tmp3 <- contrasts.fit(efit.dl2, coef = 3) # 

limma_list <- list()

topTable(tmp1, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$alphasubdom


topTable(tmp2, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$alphasub


topTable(tmp3, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$subdomsub


saveRDS(limma_list,glue("results_RNAseqRDS/limma_{my_region}.RDS"))

}



# ========================================================================
# replace p value with q value 


limma_list %>% map(head)
# ===========================================================================



i = 2
for (i in 1:length(regions)){
  
  regions[[i]] -> my_region
  
  dlNorm <- rawcount_list[[i]] %>% 
    column_to_rownames('ensgene')
  
  var_info <- behav %>% 
    filter(sample_id %in% colnames(dlNorm)) %>% 
    select(sample_id, subjectID, status)
  
  colnames(dlNorm)
  
  var_info$status %>%
    factor(.,levels = c("Alpha","Subdominant","Subordinate")) -> group.dl
  
  dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]
  
  d = apply(dlNorm, 2, as.numeric)
  dim(d)
  
  d0= DGEList(d, group = group.dl)
  dim(d0)
  rownames(d0) <- rownames(dlNorm)
  d0 <- calcNormFactors(d0)
  
  cutoff <- 5
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  dge.dl <- d0[-drop,]
  dim(dge.dl)
  
  
  design.dl <- model.matrix(~ 0 + group.dl)
  colnames(design.dl) -> mycolnames
  
  
  v.dl = voom(dge.dl, design.dl, plot = F)
  vfit.dl = lmFit(v.dl, design.dl)
  
  contrast.matrix <- makeContrasts(group.dlAlpha-group.dlSubdominant,
                                   group.dlAlpha-group.dlSubordinate, 
                                   group.dlSubdominant-group.dlSubordinate,
                                   levels=design.dl)
  
  vfit.dl2 <- contrasts.fit(vfit.dl, contrast.matrix)
  
  efit.dl = eBayes(vfit.dl)
  efit.dl2 = eBayes(vfit.dl2)
  p.dl.limma = efit.dl[["p.value"]]
  p.dl.limma2 = efit.dl2[["p.value"]]
  str(p.dl.limma2)
  glimpse(p.dl.limma2)
  topTable(efit.dl2, coef=1, adjust="BH")
  topTable(efit.dl2, coef=2, adjust="BH")
  topTable(efit.dl2, coef=3, adjust="BH")
  # saveRDS(v.dl, glue("results_RNAseqRDS/limma_vdl_{my_region}"))
  
  # 
  # How many random sampling
  R = 5000
  
  p.dl.rand = vector('list',length = R)
  
  for(i in 1 : R){
    print(paste("Starting on Permutation", i))
    
    # Randomize the traits
    
    group.dl.rand = sample(group.dl)
    
    # Model
    design.dl.rand = model.matrix(~group.dl.rand)
    colnames(design.dl.rand) <- mycolnames
    
    # Calculate p-values based on randomized traits
    v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
    vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
    efit.dl.rand = eBayes(vfit.dl.rand)
    p.dl.rand[[i]] = efit.dl.rand[["p.value"]]
  }
  
  q.dl = matrix(0, nrow = nrow(p.dl.limma), ncol = ncol(p.dl.limma))
  
  for(i in 1 : R){
    print(paste("Calculating Permutation", i))
    
    temp = p.dl.rand[[i]]
    
    for(c in 1 : 3){
      for(r in 1 : nrow(p.dl.limma)){
        if(temp[r, c] <= p.dl.limma[r, c]){
          q.dl[r, c] = q.dl[r, c] + 1
        }
      }
    }
  }
  
  q.dl = q.dl / R
  colnames(q.dl) <- mycolnames
  q.dl = as.data.frame(q.dl)
  row.names(q.dl) <- rownames(dge.dl)
  
  saveRDS(q.dl,glue("results_RNAseqRDS/limma_eFDR_{my_region}_cutoff5_R{R}_contrrast.RDS"))
  
  q.dl <- readRDS(glue("results_RNAseqRDS/limma_eFDR_{my_region}_cutoff5_R{R}.RDS"))
  
  
  png(filename = glue("results_figures/eFDR_hist_{my_region}.png"),
      width = 18, height = 17, units = "cm", res = 600)
  hist(q.dl-p.dl.limma, main = glue("{my_region}"))
  invisible(dev.off())
  
}

