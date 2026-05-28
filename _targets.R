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
dir.create(root_path, showWarnings = F)

fig_path <- here::here("docs", "gen", "figs")
fig_ext <- glue::as_glue(".png")
dir.create(fig_path, showWarnings = F)

tab_path <- here::here("docs", "gen", "tabs")
dir.create(tab_path, showWarnings = F)

post_path <- here::here("posteriors", "gen")
dir.create(post_path, showWarnings = F)


input_path <- here::here("data", "inputs_for_gen")




## Meta parameters
n_mc <- 2000
n_seeds <- 20

n_iter <- 2000
n_warmup <- 1500

tar_option_set(seed = 11667)


list_studies <- read_csv(here::here("data", "list_studies.csv")) %>% 
  select(study, SID, SIK)
list_studies



root_src <- here::here("data", "inputs_for_gen")

list(
  # Loading Data
  
  tar_target(file_norm, here::here("data", "processed", "pn_mapped.csv"), format = "file"),
  tar_target(data_norm, read_csv(file_norm)),
  
  tar_target(file_pars_qol, here::here(input_path, "pars_qol_km.csv"), format = "file"),
  tar_target(pars_qol, read_csv(file_pars_qol)),
  
  tar_target(file_pars_tte, here::here(input_path, "pars_tte_2gp.csv"), format = "file"),
  tar_target(pars_tte, read_csv(file_pars_tte)),
  
  tar_target(file_pars_tte_slopes, here::here(input_path, "pars_tte_slopes_2gp.csv"), format = "file"),
  tar_target(pars_tte_slopes, read_csv(file_pars_tte_slopes)),
  
  tar_target(time_steps, c(seq(7, 28, 7), 60, 90, 180, 270) / 365.25),
  
  tar_target(file_prof_a, here::here(input_path, "prof_a.csv"), format = "file"),
  tar_target(prof_a, read_csv(file_prof_a)),
  
  # Generating Data
  tar_target(seed, runif(n_seeds)),
  tar_target(data_gen, gen_patients(prof_a, time_steps, pars_qol, pars_tte, pars_tte_slopes, list_studies, seed = seed), pattern = map(seed)),
  
  # Model fitting
  tar_target(post_qol, fit_qol(data_gen, n_iter, n_warmup, n_attempts = 10), pattern = map(data_gen)),
  tar_target(post_tte, fit_tte_exact(data_gen, n_iter, n_warmup, n_attempts = 10), pattern = map(data_gen)),
  
  # Projection
  tar_target(summ_qloss, summarise_qls(post_qol, post_tte, batch = seed), pattern = map(post_qol, post_tte, seed)),
  
  tar_target(tab_time, summ_qloss$tab_time, pattern = map(summ_qloss)),
  tar_target(tab_age, summ_qloss$tab_age, pattern = map(summ_qloss)),
  
  tar_target(gs_qloss, vis_qls(tab_age)),

  tar_target(file_tab_time_batch, {
    f <- here::here(tab_path, "tab_time_by_batch.csv"); write_csv(tab_time, f); f
  }, format = "file"),
  tar_target(file_tab_time, {
    f <- here::here(tab_path, "tab_time.csv")
    tab_time %>% summarise(across(everything(), mean), .by = Time) %>% select(-Batch) %>% write_csv(f)
    f
  }, format = "file"),
  tar_target(file_tab_age_batch, {
    f <- here::here(tab_path, "tab_age_by_batch.csv"); write_csv(tab_age, f); f
  }, format = "file"),
  tar_target(file_tab_age, {
    f <- here::here(tab_path, "tab_age.csv")
    tab_age %>% summarise(across(everything(), mean), .by = Age) %>% select(-Batch) %>% write_csv(f)
    f
  }, format = "file"),
  tar_target(file_gs_age, {
    for (key in names(gs_qloss)) {
      ggsave(gs_qloss[[key]], filename = here::here(fig_path, paste0("g_", key, fig_ext)), width = 6, height = 5)
    }
    fig_path
  }, format = "file")
  ## extraction
)





