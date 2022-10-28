
library(tidyverse)

### filling out table with RECIST info

clin <- readxl::read_xlsx("/Volumes/spakowicz/mitox/data/clinical/raw/Patient_demographics_baseline - MH 9-9-22.xlsx") %>%
  rename(record_id = "Patient Id")

samps <- read_csv("data/mitox_sample_matching.csv")

table(samps$Category)
# dis + val = Metastatic
# ipi_nivo_met = Metastatic - Ipi/Nivo

met_samps <- samps %>%
  filter(Category == "dis" | Category == "val") %>%
  select(record_id) %>%
  mutate(record_id = ifelse(record_id == "M100", "M0100", record_id))

met_nivo_samps <- samps %>%
  filter(Category == "ipi_nivo_met") %>%
  select(record_id)

adj <-  samps %>%
  filter(Category == "adjuvant") %>%
  select(record_id)

met_samps <- merge(met_samps, clin)
table(met_samps$...21)

met_nivo_samps <- merge(met_nivo_samps, clin)
table(met_nivo_samps$...21)

adj_samps <- merge(adj, clin)
table(adj_samps$...21)
table(adj_samps$RECIST)

# get response data for all
response <- samps %>%
  select(record_id)

response <- merge(response, clin) %>%
  mutate(response = ifelse(RECIST == "CR", RECIST, ...21)) %>%
  select(`MDSC #`, response) %>%
  drop_na() %>%
  distinct() %>%
  filter(`MDSC #` != 'N/A')

# save response variable
write.csv(response, "./data/MDSCresponse.csv", row.names = F)


