### DO NOT NEED THIS -  GENE SETS.


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
