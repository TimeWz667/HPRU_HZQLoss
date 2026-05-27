library(targets)
library(tarchetypes)
library(rstan)
library(tidyverse)


tar_option_set(
  packages = c("tibble", "tidyverse", "tidybayes", "rstan", "readxl", "glue", "lme4")
)
options(dplyr.summarise.inform = FALSE)
options(mc.cores = 8)
rstan_options(auto_write = TRUE)
theme_set(theme_bw())


tar_source()


## Root settings

root_path <- here::here("docs", "gen")
fig_path <- here::here("docs", "gen", "figs")
fig_ext <- ".png"
tab_path <- here::here("docs", "gen", "tabs")
post_path <- here::here("posteriors", "gen")


dir.create(root_path, showWarnings = F)
dir.create(fig_path, showWarnings = F)
dir.create(tab_path, showWarnings = F)
dir.create(post_path, showWarnings = F)



## Meta parameters
n_mc <- 2000
n_seeds <- 10
meta_seed <- 11667



list_studies <- read_csv(here::here("data", "list_studies.csv")) %>% 
  select(study, SID, SIK)
list_studies


root_src <- here::here("data", "inputs_for_gen")







