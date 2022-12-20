library(tidyverse)
library(here)

# Read in data from CJS-ind traits study, which has more data than needed for Age of Return study
# data through 2022 for spring/summer/unknown runs of Chinook Salmon
DATA <- read.csv(here("DSYNTH2022", "dsynth_2022.spsu.csv"))

head (DATA)
str(DATA)

# Subset the data to include only wild Chinook detected at Lower Granite Dam as juveniles, 
# bypassed / run-of-river fish, and LGR-transported fish, and
# adult returns to BON and other mainstem dams upstream of BON
data_adults <- DATA %>% filter(adu_bon_first != "" | adu_mcn_first != "" | adu_ich_first != "" | adu_lgr_first != "" | adu_bon_last != "" | adu_mcn_last != "" | adu_ich_last != "" | adu_lgr_last != "") %>% 
  filter(t_rear_type == "W") %>% filter(juv_lgr_first != "") %>% 
  filter(last_trans_proj !="LMN") %>% filter(last_trans_proj !="LGS") %>% filter(last_trans_proj !="MCN" )
       
data_adults <- data_adults[,c("tag_id", "file_id", "t_species", "t_run", "t_rear_type", "pass_type", 
                            "tag_time", "rel_site", "rel_huc", "subbasin", "rel_rkm", "rel_time",
                            "length", "juv_lgr_first", "juv_lgr_last", "adu_bon_first", "adu_bon_last",
                            "adu_mcn_first", "adu_mcn_last", "adu_ich_first", "adu_ich_last",
                            "adu_lgr_first", "adu_lgr_last", "last_trans_proj", "last_trans_time", 
                            "totalrkm", "n_recaps")]

# Consolidate the transportation filter codes to two (T vs R) categories
# https://www.cbr.washington.edu/dart/metadata/pit#transport
unique(data_adults$pass_type)

dum <- str_replace(data_adults$pass_type, "LWG-T", "T")
dum <- str_replace(dum, "LWG-S", "ROR")
dum <- str_replace(dum, "LWG-SBD", "ROR")
dum <- str_replace(dum, "LWG-TB", "ROR")
dum <- str_replace(dum, "LWG-SD", "ROR")
dum <- str_replace(dum, "LWG-TD", "ROR")
# dum <- str_replace(dum, "RORBD", "ROR")
# dum <- str_replace(dum, "RORD", "ROR")
# dum <- str_replace(dum, "TB", "ROR")
# dum <- str_replace(dum, "TD", "ROR")

data_adults$pass_type_T_R <- dum


write.table(data_adults, file=here("data", "data_adults.csv"), col.names = TRUE, row.names = FALSE, sep = ",")




