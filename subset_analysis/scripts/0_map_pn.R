library(tidyverse)





read_csv(here::here("data", "processed", "population_norm_mapped.csv")) %>% 
  left_join(list_studies) %>% 
  write_csv(here::here("data", "processed", "pn_mapped.csv"))
list_studies
