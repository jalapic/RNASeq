
# tag-seq data cleaned from rna_seq github repo 

## JAMES EDIT - JUST ONE AT A TIME.

rawcount_list <- list()
temp = list.files(path = "stable_brains/raw_data",pattern="*_counts.csv")
temp0 <- paste0("stable_brains/raw_data/",temp)
rawcount_list = lapply(temp0, read_csv)
regions <- gsub("_counts.csv","",temp)
# behavior data cleaned from rna_seq github repo 
names(rawcount_list) <-regions


##

behavior <- readRDS('stable_brains/raw_data/all_behavior.RDS') %>% 
  mutate(cohort = as.numeric(as.factor(cohort))+104) %>% 
  mutate(subjectID = paste(cohort, mouseID, sep = '-'))

sample_id <- read_csv('stable_brains/raw_data/sample_id.csv') %>% 
  rename(sampleID = sample_id)

alldata <- readRDS("stable_brains/raw_data/alldata.RDS") %>% 
  mutate(cohort = as.numeric(factor(cohort))+104) %>% 
  mutate(subjectID = paste(cohort, mouseID, sep ="-")) %>% 
  mutate(bw_change_rate = (GD14_bw-GD01_BW)/GD01_BW) %>% 
  select(subjectID, cort_post,despotism, GD14_bw,bw_change_rate)



# tissue harvest time =====================================================
time_out <- read_csv("stable_brains/raw_data/end_bodyweight.csv") %>% 
  mutate(subjectID = paste(cohort,mouseID,sep="-")) %>% 
  select(subjectID,time_out) %>% 
  filter(subjectID %in% unique(sample_id$subjectID)) %>% 
  mutate(timeout_minute = lubridate::minute(time_out)) 

