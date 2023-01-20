
# tag-seq data cleaned from rna_seq github repo 

rawcount_list <- list()
temp = list.files(path = "data_clean",pattern="*_counts.csv")
temp0 <- paste0("data_clean/",temp)
rawcount_list = lapply(temp0, read_csv)
regions <- gsub("_counts.csv","",temp)
# behavior data cleaned from rna_seq github repo 
names(rawcount_list) <-regions

behavior <- readRDS('data_clean/all_behavior.RDS') %>% 
  mutate(cohort = as.numeric(as.factor(cohort))+104) %>% 
  mutate(subjectID = paste(cohort, mouseID, sep = '-'))

sample_id <- read_csv('data_raw/sample_id.csv') %>% 
  rename(sampleID = sample_id)

alldata <- readRDS("data_clean/alldata.RDS") %>% 
  mutate(cohort = as.numeric(factor(cohort))+104) %>% 
  mutate(subjectID = paste(cohort, mouseID, sep ="-")) %>% 
  mutate(bw_change_rate = (GD14_bw-GD01_BW)/GD01_BW) %>% 
  select(subjectID, cort_post,despotism, GD14_bw,bw_change_rate)



# tissue harvest time =====================================================
time_out <- read_csv("data_raw/end_bodyweight.csv") %>% 
  mutate(subjectID = glue("{cohort}-{mouseID}")) %>% 
  select(subjectID,time_out) %>% 
  filter(subjectID %in% unique(sample_id$subjectID)) %>% 
  mutate(timeout_minute = lubridate::minute(time_out)) 



# activity-regulated genes (ARGs) from Tyssowski et al. 2018 Neuron ========
# https://doi.org/10.1016/j.neuron.2018.04.001

readxl::read_excel("data_clean/Tyssowski_et_al_2018/mmc3.xlsx",
                   sheet = 2) %>% 
  rename(symbol = `gene name`,
         gene_class = `Gene Class`) %>% 
  select(symbol, gene_class) -> ARG_df


readxl::read_excel("data_clean/Tyssowski_et_al_2018/mmc6.xlsx",
                   sheet = 1) 


"data_clean/Tyssowski_et_al_2018/mmc6.xlsx" -> xl_data
tab_names <- readxl::excel_sheets(path = xl_data)

ARG_list <- lapply(tab_names, function(x) 
  readxl::read_excel(xl_data, sheet = x)) 
names(ARG_list) <- gsub(" ","_",tab_names)

ARG_df2 <- ARG_list %>% 
  map2_df(.,names(.), ~mutate(.x, gene_class = .y)) %>% 
  rename(symbol = `Gene ID`) %>% 
  select(symbol, gene_class)
  

# imprinted gene set from Higgs et al 2020 BioRxiv  =========================
# https://doi.org/10.1101/2020.07.27.222893

higgs_url = "https://raw.githubusercontent.com/MJHiggs/IG-Single-Cell-Enrichment/master/Imprinted_Gene_List.csv"

imprinted_genes  <- read_csv(url(higgs_url)) %>% 
  # found error, ensgene does not necessarily match with gene symbol. Going with gene symbol. 
  select(-Ensmbl) %>% 
  rename(symbol = Gene) %>% 
  left_join(grcm38 %>% select(symbol, ensgene, description), by = "symbol") %>% 
  mutate(Chrom = ifelse(symbol == "Tgfb1i1",7,Chrom)) # dou %>% ble check the literature


# read_csv(url(higgs_url)) %>% 
#   mutate(symbol = Gene) %>% 
#   left_join(grcm38 %>% select(ensgene, symbol)) %>% 
#   mutate(match_or_not = ifelse(ensgene == Ensmbl, "","DIFFERENT")) %>% 
#   arrange(Gene) %>% 
#   rename(grcm38_Ensmbl_ID = ensgene) %>% 
#   mutate(ensgene = Ensmbl) %>% 
#   select(-symbol) %>% 
#   left_join(grcm38 %>% select(ensgene, symbol)) %>% 
#   rename(grcm38_symbol = symbol) %>% 
#   mutate(grcm38_symbol = ifelse(grcm38_symbol == Gene, "",grcm38_symbol)) %>% 
#   write.csv("IGs_list_comparing_with_grcm38.csv", row.names = F)
